#include "internal/command/ResultDumpCommand.hpp"

#include <rhbm_gem/data/io/FileIO.hpp>
#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/GroupPotentialEntry.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/PotentialEntryQuery.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/ChimeraXHelper.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>

#include <algorithm>
#include <array>
#include <fstream>
#include <stdexcept>

namespace {
constexpr std::string_view kMapKey{ "map" };
constexpr std::string_view kPrinterOption{ "--printer" };
constexpr std::string_view kModelKeyListOption{ "--model-keylist" };
constexpr std::string_view kMapOption{ "--map" };

std::shared_ptr<rhbm_gem::ModelObject> LoadModelFromDatabase(
    rhbm_gem::DataObjectManager & data_manager,
    std::string_view key_tag,
    std::string_view label)
{
    try
    {
        data_manager.LoadDataObject(std::string(key_tag));
        return data_manager.GetTypedDataObject<rhbm_gem::ModelObject>(std::string(key_tag));
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to load " + std::string(label) + " with key tag '"
            + std::string(key_tag) + "' as ModelObject: " + ex.what());
    }
}

std::shared_ptr<rhbm_gem::MapObject> LoadMapFromFile(
    rhbm_gem::DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label)
{
    try
    {
        data_manager.ProcessFile(path, std::string(key_tag));
        return data_manager.GetTypedDataObject<rhbm_gem::MapObject>(std::string(key_tag));
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to process " + std::string(label) + " from '" + path.string()
            + "' as MapObject: " + ex.what());
    }
}
}

namespace rhbm_gem {

ResultDumpCommand::ResultDumpCommand() :
    CommandWithRequest<ResultDumpRequest>{},
    m_map_key_tag{ kMapKey }, m_map_object{ nullptr }
{
}

void ResultDumpCommand::NormalizeRequest()
{
    auto & request{ MutableRequest() };
    CoerceEnum(
        request.printer_choice,
        kPrinterOption,
        PrinterType::GAUS_ESTIMATES,
        "Printer choice");
    ValidateOptionalPath(request.map_file_path, kMapOption, "Map file");
    RequireNonEmptyList(request.model_key_tag_list, kModelKeyListOption, "Model key list");
}

void ResultDumpCommand::ValidateOptions()
{
    const auto & request{ RequestOptions() };
    RequireCondition(
        request.printer_choice != PrinterType::MAP_VALUE || !request.map_file_path.empty(),
        kMapOption,
        "A map file is required when '--printer map' is selected.");
}

void ResultDumpCommand::ResetRuntimeState()
{
    m_selected_atom_list_map.clear();
    m_model_object_list.clear();
    m_map_object.reset();
}

bool ResultDumpCommand::BuildDataObjectList()
{
    const auto & request{ RequestOptions() };
    ScopeTimer timer("ResultDumpCommand::BuildDataObjectList");
    try
    {
        m_data_manager.SetDatabaseManager(request.database_path);
        if (request.map_file_path.empty())
        {
            m_map_object.reset();
        }
        else
        {
            m_map_object = LoadMapFromFile(
                m_data_manager,
                request.map_file_path,
                m_map_key_tag,
                "map file");
        }
        m_selected_atom_list_map.clear();
        for (const auto & key : request.model_key_tag_list)
        {
            auto model_object{ LoadModelFromDatabase(
                m_data_manager,
                key,
                "model object") };
            m_model_object_list.emplace_back(model_object);
            const auto & selected_atom_list{ model_object->GetSelectedAtomList() };
            m_selected_atom_list_map[key] = {
                selected_atom_list.begin(),
                selected_atom_list.end()
            };
            Logger::Log(
                LogLevel::Info,
                "Selected atoms for key tag [" + key + "]: "
                    + std::to_string(m_selected_atom_list_map[key].size()));
        }
    }
    catch (const std::exception & e)
    {
        Logger::Log(
            LogLevel::Error,
            "ResultDumpCommand::BuildDataObjectList : " + std::string(e.what()));
        return false;
    }
    return true;
}

void ResultDumpCommand::RunAtomOutlierDumping()
{
    for (const auto & model_object : m_model_object_list)
    {
        const auto key_tag{ model_object->GetKeyTag() };
        std::unordered_map<std::string, std::vector<std::array<float, 3>>> outlier_position_map;
        const std::string file_name{ "atom_outlier_list_" + model_object->GetPdbID() };
        const auto output_path{ BuildOutputPath(file_name, ".csv") };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(
                LogLevel::Error,
                "Could not open file " + output_path.string() + " for writing.\n");
            return;
        }
        outfile << "SerialID,ClassType,Residue,Element,Spot\n";
        for (auto & atom : m_selected_atom_list_map.at(key_tag))
        {
            for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
            {
                const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
                if (!atom->GetLocalPotentialEntry()->GetOutlierTag(class_key)) continue;
                outfile << atom->GetSerialID() << ',' << class_key << ','
                        << ChemicalDataHelper::GetLabel(atom->GetResidue()) << ','
                        << ChemicalDataHelper::GetLabel(atom->GetElement()) << ','
                        << atom->GetAtomID() << '\n';
                if (atom->IsMainChainAtom())
                {
                    if (atom->GetElement() == Element::OXYGEN) continue;
                    outlier_position_map[class_key].emplace_back(atom->GetPosition());
                }
            }
        }
        outfile.close();
        Logger::Log(LogLevel::Info, "Output file: " + output_path.string());

        for (const auto & [class_key, position_list] : outlier_position_map)
        {
            const std::string class_file_name{
                "atom_outlier_position_" + class_key + "_" + model_object->GetPdbID()
            };
            const auto class_output_path{ BuildOutputPath(class_file_name, ".cmm") };
            const std::string marker_set_name{ class_key + "_" + model_object->GetPdbID() };
            if (!ChimeraXHelper::WriteCMMPoints(
                    position_list,
                    class_output_path.string(),
                    1.0f,
                    {},
                    marker_set_name))
            {
                Logger::Log(
                    LogLevel::Error,
                    "Could not open file " + class_output_path.string() + " for writing.\n");
                continue;
            }
            Logger::Log(LogLevel::Info, "Output file: " + class_output_path.string());
        }
    }
}

void ResultDumpCommand::RunAtomPositionDumping()
{
    for (const auto & model_object : m_model_object_list)
    {
        const auto key_tag{ model_object->GetKeyTag() };
        const std::string file_name{ "atom_position_list_" + model_object->GetPdbID() };
        const auto output_path{ BuildOutputPath(file_name, ".csv") };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(
                LogLevel::Error,
                "Could not open file " + output_path.string() + " for writing.\n");
            return;
        }
        outfile << "SerialID,X,Y,Z,Residue,Element,Spot\n";
        for (auto & atom : m_selected_atom_list_map.at(key_tag))
        {
            outfile << atom->GetSerialID() << ',' << atom->GetPosition().at(0) << ','
                    << atom->GetPosition().at(1) << ',' << atom->GetPosition().at(2) << ','
                    << ChemicalDataHelper::GetLabel(atom->GetResidue()) << ','
                    << ChemicalDataHelper::GetLabel(atom->GetElement()) << ','
                    << atom->GetAtomID() << '\n';
        }
        outfile.close();
        Logger::Log(LogLevel::Info, "Output file: " + output_path.string());
    }
}

void ResultDumpCommand::RunMapValueDumping()
{
    if (m_map_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Please give the path of map file via -m option.");
        return;
    }

    const auto margin{ 3.0f };

    for (const auto & model_object : m_model_object_list)
    {
        const auto key_tag{ model_object->GetKeyTag() };
        const auto atom_size{ m_selected_atom_list_map.at(key_tag).size() };
        if (atom_size == 0)
        {
            Logger::Log(
                LogLevel::Warning,
                "No selected atoms found for key tag [" + key_tag + "]. Skipping map-value dump.");
            continue;
        }
        std::array<float, 3> atom_range_min, atom_range_max;
        std::vector<float> x_list, y_list, z_list;
        x_list.reserve(atom_size);
        y_list.reserve(atom_size);
        z_list.reserve(atom_size);
        for (auto & atom : m_selected_atom_list_map.at(key_tag))
        {
            x_list.emplace_back(atom->GetPosition().at(0));
            y_list.emplace_back(atom->GetPosition().at(1));
            z_list.emplace_back(atom->GetPosition().at(2));
        }
        atom_range_min.at(0) = ArrayStats<float>::ComputeMin(x_list.data(), atom_size) - margin;
        atom_range_min.at(1) = ArrayStats<float>::ComputeMin(y_list.data(), atom_size) - margin;
        atom_range_min.at(2) = ArrayStats<float>::ComputeMin(z_list.data(), atom_size) - margin;
        atom_range_max.at(0) = ArrayStats<float>::ComputeMax(x_list.data(), atom_size) + margin;
        atom_range_max.at(1) = ArrayStats<float>::ComputeMax(y_list.data(), atom_size) + margin;
        atom_range_max.at(2) = ArrayStats<float>::ComputeMax(z_list.data(), atom_size) + margin;

        const std::string file_name{
            "map_value_list_" + model_object->GetEmdID() + "_" + model_object->GetPdbID()
        };
        const auto output_path{ BuildOutputPath(file_name, ".csv") };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(
                LogLevel::Error,
                "Could not open file " + output_path.string() + " for writing.\n");
            return;
        }
        outfile << "GridID,X,Y,Z,MapValue\n";
        auto count{ 0 };
        for (size_t i = 0; i < m_map_object->GetMapValueArraySize(); i++)
        {
            const auto grid_position{ m_map_object->GetGridPosition(i) };
            if (grid_position.at(0) < atom_range_min.at(0)) continue;
            if (grid_position.at(0) > atom_range_max.at(0)) continue;
            if (grid_position.at(1) < atom_range_min.at(1)) continue;
            if (grid_position.at(1) > atom_range_max.at(1)) continue;
            if (grid_position.at(2) < atom_range_min.at(2)) continue;
            if (grid_position.at(2) > atom_range_max.at(2)) continue;
            const auto map_value{ m_map_object->GetMapValue(i) };
            if (map_value <= 0.0f) continue;
            outfile << i << ',' << grid_position.at(0) << ',' << grid_position.at(1) << ','
                    << grid_position.at(2) << ',' << map_value << '\n';
            count++;
        }
        Logger::Log(
            LogLevel::Info,
            " Selected map grid size = " + std::to_string(count) + " / "
                + std::to_string(m_map_object->GetMapValueArraySize()));
    }
}

void ResultDumpCommand::RunGausEstimatesDumping()
{
    for (const auto & model_object : m_model_object_list)
    {
        const auto key_tag{ model_object->GetKeyTag() };
        const std::string csv_file_name{ "local_gaus_list_" + model_object->GetPdbID() };
        const auto output_csv_file{ BuildOutputPath(csv_file_name, ".csv") };
        std::ofstream outfile(output_csv_file);
        if (!outfile.is_open())
        {
            Logger::Log(
                LogLevel::Error,
                "Could not open file " + output_csv_file.string() + " for writing.\n");
            return;
        }

        outfile << "SerialID,Amplitude,Width,X,Y,Z,Residue,Element,Spot\n";
        for (auto & atom : m_selected_atom_list_map.at(key_tag))
        {
            auto entry{ atom->GetLocalPotentialEntry() };
            outfile << atom->GetSerialID() << ',' << entry->GetAmplitudeEstimateMDPDE() << ','
                    << entry->GetWidthEstimateMDPDE() << ',' << atom->GetPosition().at(0) << ','
                    << atom->GetPosition().at(1) << ',' << atom->GetPosition().at(2) << ','
                    << ChemicalDataHelper::GetLabel(atom->GetResidue()) << ','
                    << ChemicalDataHelper::GetLabel(atom->GetElement()) << ','
                    << atom->GetAtomID() << '\n';
        }
        outfile.close();
        Logger::Log(LogLevel::Info, "Output file: " + output_csv_file.string());
/*
        const std::string amplitude_file_name{ model_object->GetPdbID() + "_gaus_amplitude" };
        const auto output_amplitude_cif_file{ BuildOutputPath(amplitude_file_name, ".cif") };
        WriteModel(output_amplitude_cif_file, *model_object, 0);
        Logger::Log(LogLevel::Info, "Output file: " + output_amplitude_cif_file.string());

        const std::string width_file_name{ model_object->GetPdbID() + "_gaus_width" };
        const auto output_width_cif_file{ BuildOutputPath(width_file_name, ".cif") };
        WriteModel(output_width_cif_file, *model_object, 1);
        Logger::Log(LogLevel::Info, "Output file: " + output_width_cif_file.string());

        const std::string intensity_file_name{ model_object->GetPdbID() + "_gaus_intensity" };
        const auto output_intensity_cif_file{ BuildOutputPath(intensity_file_name, ".cif") };
        WriteModel(output_intensity_cif_file, *model_object, 2);
        Logger::Log(LogLevel::Info, "Output file: " + output_intensity_cif_file.string());*/
    }
}

void ResultDumpCommand::RunGroupGausEstimatesDumping()
{
    const auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };

    for (const auto & model_object : m_model_object_list)
    {
        auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object.get()) };

        const std::string file_name{ "group_gaus_list_" + model_object->GetPdbID() };
        const auto output_path{ BuildOutputPath(file_name, ".csv") };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(
                LogLevel::Error,
                "Could not open file " + output_path.string() + " for writing.\n");
            return;
        }
        outfile << "Residue,Spot,Amplitude,Width\n";
        for (auto & residue : ChemicalDataHelper::GetStandardAminoAcidList())
        {
            const auto residue_name{ ChemicalDataHelper::GetLabel(residue) };
            const auto component_key{ static_cast<ComponentKey>(residue) };
            for (auto & spot : ComponentHelper::GetSpotList(residue))
            {
                const auto atom_key{ static_cast<uint16_t>(spot) };
                const auto group_key{ AtomClassifier::GetGroupKeyInClass(component_key, atom_key) };
                const auto atom_id{ model_object->GetAtomKeySystemPtr()->GetAtomId(atom_key) };
                if (!entry_iter->IsAvailableAtomGroupKey(group_key, class_key)) continue;
                outfile << residue_name << ',' << atom_id << ','
                        << entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 0) << ','
                        << entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 1) << '\n';
            }
        }

        outfile.close();
        Logger::Log(LogLevel::Info, "Output file: " + output_path.string());
    }
}

void ResultDumpCommand::RunDumpWorkflow()
{
    ScopeTimer timer("ResultDumpCommand::RunResultDump");
    Logger::Log(
        LogLevel::Info,
        "Total number of model object sets to be dump: "
            + std::to_string(RequestOptions().model_key_tag_list.size()));

    const bool has_selected_atoms{
        std::any_of(
            m_selected_atom_list_map.begin(),
            m_selected_atom_list_map.end(),
            [](const auto & entry) { return !entry.second.empty(); })
    };
    if (!has_selected_atoms)
    {
        Logger::Log(
            LogLevel::Warning,
            "No selected atoms with local potential entries were found. Skipping dump.");
        return;
    }

    switch (RequestOptions().printer_choice)
    {
    case PrinterType::ATOM_POSITION:
        RunAtomPositionDumping();
        break;
    case PrinterType::MAP_VALUE:
        RunMapValueDumping();
        break;
    case PrinterType::GAUS_ESTIMATES:
        RunGroupGausEstimatesDumping();
        RunGausEstimatesDumping();
        break;
    case PrinterType::ATOM_OUTLIER:
        RunAtomOutlierDumping();
        break;
    default:
        Logger::Log(
            LogLevel::Warning,
            "Invalid printer choice input : ["
                + std::to_string(static_cast<int>(RequestOptions().printer_choice)) + "]");
        Logger::Log(
            LogLevel::Warning,
            "Available Printer Choices:\n"
            "  [0] AtomPositionDumping\n"
            "  [1] MapValueDumping\n"
            "  [2] GausEstimatesDumping\n"
            "  [3] AtomOutlierDumping");
        break;
    }
}

bool ResultDumpCommand::ExecuteImpl()
{
    if (!BuildDataObjectList())
    {
        return false;
    }
    RunDumpWorkflow();
    return true;
}

} // namespace rhbm_gem
