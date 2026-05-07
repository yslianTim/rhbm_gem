#include "detail/CommandBase.hpp"

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/ChimeraXHelper.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/KeyPacker.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace rhbm_gem {

class ResultDumpCommand final : public CommandBase<ResultDumpRequest>
{
public:
    ResultDumpCommand();

private:
    void NormalizeAndValidateRequest(ResultDumpRequest & request) override;
    void ValidatePreparedRequest(const ResultDumpRequest & request) override;
    bool ExecuteImpl(const ResultDumpRequest & request) override;
};

namespace {

struct ResultDumpInputs
{
    std::vector<std::unique_ptr<ModelObject>> model_objects;
    std::unique_ptr<MapObject> map_object;
};

std::filesystem::path BuildResultDumpOutputPath(
    const std::filesystem::path & output_dir,
    std::string_view stem,
    std::string_view extension)
{
    std::string normalized_extension{ extension };
    if (!normalized_extension.empty() && normalized_extension.front() != '.')
    {
        normalized_extension.insert(normalized_extension.begin(), '.');
    }
    return output_dir / (std::string(stem) + normalized_extension);
}

std::optional<ResultDumpInputs> LoadResultDumpInputs(const ResultDumpRequest & request)
{
    ScopeTimer timer("ResultDumpCommand::BuildDataObjectList");
    try
    {
        DataRepository repository{ request.database_path };
        ResultDumpInputs inputs;
        if (request.map_file_path.empty())
        {
            inputs.map_object.reset();
        }
        else
        {
            auto map_object{ ReadMap(request.map_file_path) };
            map_object->SetKeyTag("map");
            inputs.map_object = std::move(map_object);
        }
        inputs.model_objects.reserve(request.model_key_tag_list.size());
        for (const auto & key : request.model_key_tag_list)
        {
            auto model_object{ repository.LoadModel(key) };
            Logger::Log(LogLevel::Info,
                "Selected atoms for key tag [" + key + "]: "
                    + std::to_string(model_object->GetSelectedAtoms().size()));
            inputs.model_objects.emplace_back(std::move(model_object));
        }
        return inputs;
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "ResultDumpCommand::BuildDataObjectList : Failed to load dump inputs: "
                + std::string(e.what()));
        return std::nullopt;
    }
}

void RunAtomOutlierDumping(
    const std::vector<std::unique_ptr<ModelObject>> & model_object_list,
    const std::filesystem::path & output_dir)
{
    for (const auto & model_object : model_object_list)
    {
        std::unordered_map<std::string, std::vector<std::array<float, 3>>> outlier_position_map;
        const std::string file_name{ "atom_outlier_list_" + model_object->GetPdbID() };
        const auto output_path{ BuildResultDumpOutputPath(output_dir, file_name, ".csv") };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
                "Could not open file " + output_path.string() + " for writing.\n");
            return;
        }
        outfile << "SerialID,ClassType,Residue,Element,Spot\n";
        for (auto * atom : model_object->GetSelectedAtoms())
        {
            for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
            {
                const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
                const auto annotation{ LocalPotentialView::RequireFor(*atom).FindAnnotation(class_key) };
                if (!annotation.has_value() || !annotation->is_outlier) continue;
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
            const auto class_output_path{ BuildResultDumpOutputPath(output_dir, class_file_name, ".cmm") };
            const std::string marker_set_name{ class_key + "_" + model_object->GetPdbID() };
            if (!ChimeraXHelper::WriteCMMPoints(
                    position_list,
                    class_output_path.string(),
                    1.0f,
                    {},
                    marker_set_name))
            {
                Logger::Log(LogLevel::Error,
                    "Could not open file " + class_output_path.string() + " for writing.\n");
                continue;
            }
            Logger::Log(LogLevel::Info, "Output file: " + class_output_path.string());
        }
    }
}

void RunAtomPositionDumping(
    const std::vector<std::unique_ptr<ModelObject>> & model_object_list,
    const std::filesystem::path & output_dir)
{
    for (const auto & model_object : model_object_list)
    {
        const std::string file_name{ "atom_position_list_" + model_object->GetPdbID() };
        const auto output_path{ BuildResultDumpOutputPath(output_dir, file_name, ".csv") };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
                "Could not open file " + output_path.string() + " for writing.\n");
            return;
        }
        outfile << "SerialID,X,Y,Z,Residue,Element,Spot\n";
        for (auto * atom : model_object->GetSelectedAtoms())
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

void RunMapValueDumping(
    const std::vector<std::unique_ptr<ModelObject>> & model_object_list,
    const MapObject * map_object,
    const std::filesystem::path & output_dir)
{
    if (map_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Please give the path of map file via -m option.");
        return;
    }

    const auto margin{ 3.0f };

    for (const auto & model_object : model_object_list)
    {
        const auto key_tag{ model_object->GetKeyTag() };
        const auto & selected_atoms{ model_object->GetSelectedAtoms() };
        const auto atom_size{ selected_atoms.size() };
        if (atom_size == 0)
        {
            Logger::Log(LogLevel::Warning,
                "No selected atoms found for key tag [" + key_tag + "]. Skipping map-value dump.");
            continue;
        }
        std::array<float, 3> atom_range_min, atom_range_max;
        std::vector<float> x_list, y_list, z_list;
        x_list.reserve(atom_size);
        y_list.reserve(atom_size);
        z_list.reserve(atom_size);
        for (auto * atom : selected_atoms)
        {
            x_list.emplace_back(atom->GetPosition().at(0));
            y_list.emplace_back(atom->GetPosition().at(1));
            z_list.emplace_back(atom->GetPosition().at(2));
        }
        atom_range_min.at(0) = array_helper::ComputeMin(x_list.data(), atom_size) - margin;
        atom_range_min.at(1) = array_helper::ComputeMin(y_list.data(), atom_size) - margin;
        atom_range_min.at(2) = array_helper::ComputeMin(z_list.data(), atom_size) - margin;
        atom_range_max.at(0) = array_helper::ComputeMax(x_list.data(), atom_size) + margin;
        atom_range_max.at(1) = array_helper::ComputeMax(y_list.data(), atom_size) + margin;
        atom_range_max.at(2) = array_helper::ComputeMax(z_list.data(), atom_size) + margin;

        const std::string file_name{
            "map_value_list_" + model_object->GetEmdID() + "_" + model_object->GetPdbID()
        };
        const auto output_path{ BuildResultDumpOutputPath(output_dir, file_name, ".csv") };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
                "Could not open file " + output_path.string() + " for writing.\n");
            return;
        }
        outfile << "GridID,X,Y,Z,MapValue\n";
        auto count{ 0 };
        for (size_t i = 0; i < map_object->GetMapValueArraySize(); i++)
        {
            const auto grid_position{ map_object->GetGridPosition(i) };
            if (grid_position.at(0) < atom_range_min.at(0)) continue;
            if (grid_position.at(0) > atom_range_max.at(0)) continue;
            if (grid_position.at(1) < atom_range_min.at(1)) continue;
            if (grid_position.at(1) > atom_range_max.at(1)) continue;
            if (grid_position.at(2) < atom_range_min.at(2)) continue;
            if (grid_position.at(2) > atom_range_max.at(2)) continue;
            const auto map_value{ map_object->GetMapValue(i) };
            if (map_value <= 0.0f) continue;
            outfile << i << ',' << grid_position.at(0) << ',' << grid_position.at(1) << ','
                    << grid_position.at(2) << ',' << map_value << '\n';
            count++;
        }
        Logger::Log(LogLevel::Info,
            " Selected map grid size = " + std::to_string(count) + " / "
                + std::to_string(map_object->GetMapValueArraySize()));
    }
}

void RunGausEstimatesDumping(
    const std::vector<std::unique_ptr<ModelObject>> & model_object_list,
    const std::filesystem::path & output_dir)
{
    for (const auto & model_object : model_object_list)
    {
        const std::string csv_file_name{ "local_gaus_list_" + model_object->GetPdbID() };
        const auto output_csv_file{ BuildResultDumpOutputPath(output_dir, csv_file_name, ".csv") };
        std::ofstream outfile(output_csv_file);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
                "Could not open file " + output_csv_file.string() + " for writing.\n");
            return;
        }

        outfile << "SerialID,Amplitude,Width,X,Y,Z,Residue,Element,Spot\n";
        for (auto * atom : model_object->GetSelectedAtoms())
        {
            const auto entry{ LocalPotentialView::RequireFor(*atom) };
            const auto & estimate{ entry.GetEstimateMDPDE() };
            outfile << atom->GetSerialID() << ',' << estimate.GetAmplitude() << ','
                    << estimate.GetWidth() << ',' << atom->GetPosition().at(0) << ','
                    << atom->GetPosition().at(1) << ',' << atom->GetPosition().at(2) << ','
                    << ChemicalDataHelper::GetLabel(atom->GetResidue()) << ','
                    << ChemicalDataHelper::GetLabel(atom->GetElement()) << ','
                    << atom->GetAtomID() << '\n';
        }
        outfile.close();
        Logger::Log(LogLevel::Info, "Output file: " + output_csv_file.string());
/*
        const std::string amplitude_file_name{ model_object->GetPdbID() + "_gaus_amplitude" };
        const auto output_amplitude_cif_file{ BuildResultDumpOutputPath(output_dir, amplitude_file_name, ".cif") };
        WriteModel(output_amplitude_cif_file, *model_object, 0);
        Logger::Log(LogLevel::Info, "Output file: " + output_amplitude_cif_file.string());

        const std::string width_file_name{ model_object->GetPdbID() + "_gaus_width" };
        const auto output_width_cif_file{ BuildResultDumpOutputPath(output_dir, width_file_name, ".cif") };
        WriteModel(output_width_cif_file, *model_object, 1);
        Logger::Log(LogLevel::Info, "Output file: " + output_width_cif_file.string());

        const std::string intensity_file_name{ model_object->GetPdbID() + "_gaus_intensity" };
        const auto output_intensity_cif_file{ BuildResultDumpOutputPath(output_dir, intensity_file_name, ".cif") };
        WriteModel(output_intensity_cif_file, *model_object, 2);
        Logger::Log(LogLevel::Info, "Output file: " + output_intensity_cif_file.string());*/
    }
}

void RunGroupGausEstimatesDumping(
    const std::vector<std::unique_ptr<ModelObject>> & model_object_list,
    const std::filesystem::path & output_dir)
{
    const auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    for (const auto & model_object : model_object_list)
    {
        const ModelAnalysisView entry_view{ *model_object };
        const std::string file_name{ "group_gaus_list_" + model_object->GetPdbID() };
        const auto output_path{ BuildResultDumpOutputPath(output_dir, file_name, ".csv") };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
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
                const auto group_key{ KeyPackerComponentAtomClass::Pack(component_key, atom_key) };
                const auto atom_id{ model_object->FindAtomID(atom_key) };
                if (!entry_view.HasAtomGroup(group_key, class_key)) continue;
                outfile << residue_name << ',' << atom_id << ','
                        << entry_view.GetAtomGausEstimatePrior(group_key, class_key, 0) << ','
                        << entry_view.GetAtomGausEstimatePrior(group_key, class_key, 1) << '\n';
            }
        }

        outfile.close();
        Logger::Log(LogLevel::Info, "Output file: " + output_path.string());
    }
}

} // namespace

ResultDumpCommand::ResultDumpCommand() : CommandBase<ResultDumpRequest>{}
{
}

void ResultDumpCommand::NormalizeAndValidateRequest(ResultDumpRequest & request)
{
    CoerceEnum(request.printer_choice, "--printer", PrinterType::GAUS_ESTIMATES, "Printer choice");
    ValidateOptionalPath(request.map_file_path, "--map", "Map file");
    RequireNonEmptyList(request.model_key_tag_list, "--model-keylist", "Model key list");
}

void ResultDumpCommand::ValidatePreparedRequest(const ResultDumpRequest & request)
{
    RequirePrepareCondition(
        request.printer_choice != PrinterType::MAP_VALUE || !request.map_file_path.empty(),
        "--map",
        "A map file is required when '--printer map' is selected.");
}

bool ResultDumpCommand::ExecuteImpl(const ResultDumpRequest & request)
{
    auto inputs{ LoadResultDumpInputs(request) };
    if (!inputs.has_value())
    {
        return false;
    }

    ScopeTimer timer("ResultDumpCommand::RunResultDump");
    Logger::Log(LogLevel::Info,
        "Total number of model object sets to be dump: "
            + std::to_string(request.model_key_tag_list.size()));

    const bool has_selected_atoms{
        std::any_of(
            inputs->model_objects.begin(),
            inputs->model_objects.end(),
            [](const auto & model_object) { return !model_object->GetSelectedAtoms().empty(); })
    };
    if (!has_selected_atoms)
    {
        Logger::Log(
            LogLevel::Warning,
            "No selected atoms with local potential entries were found. Skipping dump.");
        return true;
    }

    switch (request.printer_choice)
    {
    case PrinterType::ATOM_POSITION:
        RunAtomPositionDumping(inputs->model_objects, request.output_dir);
        break;
    case PrinterType::MAP_VALUE:
        RunMapValueDumping(inputs->model_objects, inputs->map_object.get(), request.output_dir);
        break;
    case PrinterType::GAUS_ESTIMATES:
        RunGroupGausEstimatesDumping(inputs->model_objects, request.output_dir);
        RunGausEstimatesDumping(inputs->model_objects, request.output_dir);
        break;
    case PrinterType::ATOM_OUTLIER:
        RunAtomOutlierDumping(inputs->model_objects, request.output_dir);
        break;
    default:
        Logger::Log(
            LogLevel::Warning,
            "Invalid printer choice input : ["
                + std::to_string(static_cast<int>(request.printer_choice)) + "]");
        Logger::Log(
            LogLevel::Warning,
            "Available Printer Choices:\n"
            "  [0] AtomPositionDumping\n"
            "  [1] MapValueDumping\n"
            "  [2] GausEstimatesDumping\n"
            "  [3] AtomOutlierDumping");
        break;
    }
    return true;
}

namespace command_internal {

CommandResult ExecuteResultDumpCommand(const ResultDumpRequest & request)
{
    ResultDumpCommand command;
    return command.ExecuteRequest(request);
}

} // namespace command_internal

} // namespace rhbm_gem
