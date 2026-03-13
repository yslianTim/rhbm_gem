#include "internal/workflow/ResultDumpWorkflow.hpp"

#include "internal/CommandDataLoader.hpp"
#include "workflow/DataObjectWorkflowOps.hpp"
#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/io/FileIO.hpp>
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
#include <memory>
#include <string>
#include <unordered_map>

namespace rhbm_gem::detail {
namespace {

std::filesystem::path BuildOutputPath(
    const ResultDumpWorkflowContext & context,
    std::string_view stem,
    std::string_view extension)
{
    std::string normalized_extension{ extension };
    if (!normalized_extension.empty() && normalized_extension.front() != '.')
    {
        normalized_extension.insert(normalized_extension.begin(), '.');
    }
    return context.options.folder_path / (std::string(stem) + normalized_extension);
}

bool BuildDataObjectList(ResultDumpWorkflowContext & context)
{
    ScopeTimer timer("ResultDumpCommand::BuildDataObjectList");
    try
    {
        context.data_manager.SetDatabaseManager(context.options.database_path);
        context.map_object = command_data_loader::OptionalProcessMapFile(
            context.data_manager,
            context.options.map_file_path,
            context.map_key_tag,
            "map file");
        context.selected_atom_list_map.clear();
        for (auto & key : context.options.model_key_tag_list)
        {
            auto model_object{ command_data_loader::LoadModelObject(
                context.data_manager,
                key,
                "model object") };
            context.model_object_list.emplace_back(model_object);
            ModelAtomCollectorOptions collector_options;
            collector_options.selected_only = false;
            collector_options.require_local_potential_entry = true;
            context.selected_atom_list_map[key] = CollectModelAtoms(*model_object, collector_options);
            Logger::Log(
                LogLevel::Info,
                "Selected atoms for key tag [" + key + "]: "
                    + std::to_string(context.selected_atom_list_map[key].size()));
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

void RunAtomOutlierDumping(ResultDumpWorkflowContext & context)
{
    for (const auto & model_object : context.model_object_list)
    {
        const auto key_tag{ model_object->GetKeyTag() };
        std::unordered_map<std::string, std::vector<std::array<float, 3>>> outlier_position_map;
        const std::string file_name{ "atom_outlier_list_" + model_object->GetPdbID() };
        const auto output_path{ BuildOutputPath(context, file_name, ".csv") };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(
                LogLevel::Error,
                "Could not open file " + output_path.string() + " for writing.\n");
            return;
        }
        outfile << "SerialID,ClassType,Residue,Element,Spot\n";
        for (auto & atom : context.selected_atom_list_map.at(key_tag))
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
            const auto class_output_path{ BuildOutputPath(context, class_file_name, ".cmm") };
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

void RunAtomPositionDumping(ResultDumpWorkflowContext & context)
{
    for (const auto & model_object : context.model_object_list)
    {
        const auto key_tag{ model_object->GetKeyTag() };
        const std::string file_name{ "atom_position_list_" + model_object->GetPdbID() };
        const auto output_path{ BuildOutputPath(context, file_name, ".csv") };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(
                LogLevel::Error,
                "Could not open file " + output_path.string() + " for writing.\n");
            return;
        }
        outfile << "SerialID,X,Y,Z,Residue,Element,Spot\n";
        for (auto & atom : context.selected_atom_list_map.at(key_tag))
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

void RunMapValueDumping(ResultDumpWorkflowContext & context)
{
    if (context.map_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Please give the path of map file via -m option.");
        return;
    }

    const auto margin{ 3.0f };

    for (const auto & model_object : context.model_object_list)
    {
        const auto key_tag{ model_object->GetKeyTag() };
        const auto atom_size{ context.selected_atom_list_map.at(key_tag).size() };
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
        for (auto & atom : context.selected_atom_list_map.at(key_tag))
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
        const auto output_path{ BuildOutputPath(context, file_name, ".csv") };
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
        for (size_t i = 0; i < context.map_object->GetMapValueArraySize(); i++)
        {
            const auto grid_position{ context.map_object->GetGridPosition(i) };
            if (grid_position.at(0) < atom_range_min.at(0)) continue;
            if (grid_position.at(0) > atom_range_max.at(0)) continue;
            if (grid_position.at(1) < atom_range_min.at(1)) continue;
            if (grid_position.at(1) > atom_range_max.at(1)) continue;
            if (grid_position.at(2) < atom_range_min.at(2)) continue;
            if (grid_position.at(2) > atom_range_max.at(2)) continue;
            const auto map_value{ context.map_object->GetMapValue(i) };
            if (map_value <= 0.0f) continue;
            outfile << i << ',' << grid_position.at(0) << ',' << grid_position.at(1) << ','
                    << grid_position.at(2) << ',' << map_value << '\n';
            count++;
        }
        Logger::Log(
            LogLevel::Info,
            " Selected map grid size = " + std::to_string(count) + " / "
                + std::to_string(context.map_object->GetMapValueArraySize()));
    }
}

void RunGausEstimatesDumping(ResultDumpWorkflowContext & context)
{
    for (const auto & model_object : context.model_object_list)
    {
        const auto key_tag{ model_object->GetKeyTag() };
        const std::string csv_file_name{ "atom_gaus_list_" + model_object->GetPdbID() };
        const auto output_csv_file{ BuildOutputPath(context, csv_file_name, ".csv") };
        std::ofstream outfile(output_csv_file);
        if (!outfile.is_open())
        {
            Logger::Log(
                LogLevel::Error,
                "Could not open file " + output_csv_file.string() + " for writing.\n");
            return;
        }

        outfile << "SerialID,Amplitude,Width,X,Y,Z,Residue,Element,Spot\n";
        for (auto & atom : context.selected_atom_list_map.at(key_tag))
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

        const std::string amplitude_file_name{ model_object->GetPdbID() + "_gaus_amplitude" };
        const auto output_amplitude_cif_file{ BuildOutputPath(context, amplitude_file_name, ".cif") };
        WriteModel(output_amplitude_cif_file, *model_object, 0);
        Logger::Log(LogLevel::Info, "Output file: " + output_amplitude_cif_file.string());

        const std::string width_file_name{ model_object->GetPdbID() + "_gaus_width" };
        const auto output_width_cif_file{ BuildOutputPath(context, width_file_name, ".cif") };
        WriteModel(output_width_cif_file, *model_object, 1);
        Logger::Log(LogLevel::Info, "Output file: " + output_width_cif_file.string());

        const std::string intensity_file_name{ model_object->GetPdbID() + "_gaus_intensity" };
        const auto output_intensity_cif_file{ BuildOutputPath(context, intensity_file_name, ".cif") };
        WriteModel(output_intensity_cif_file, *model_object, 2);
        Logger::Log(LogLevel::Info, "Output file: " + output_intensity_cif_file.string());
    }
}

void RunGroupGausEstimatesDumping(ResultDumpWorkflowContext & context)
{
    const auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };

    for (const auto & model_object : context.model_object_list)
    {
        auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object.get()) };

        const std::string file_name{ "group_gaus_list_" + model_object->GetPdbID() };
        const auto output_path{ BuildOutputPath(context, file_name, ".csv") };
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

void RunResultDump(ResultDumpWorkflowContext & context)
{
    ScopeTimer timer("ResultDumpCommand::RunResultDump");
    Logger::Log(
        LogLevel::Info,
        "Total number of model object sets to be dump: "
            + std::to_string(context.options.model_key_tag_list.size()));

    const bool has_selected_atoms{
        std::any_of(
            context.selected_atom_list_map.begin(),
            context.selected_atom_list_map.end(),
            [](const auto & entry) { return !entry.second.empty(); })
    };
    if (!has_selected_atoms)
    {
        Logger::Log(
            LogLevel::Warning,
            "No selected atoms with local potential entries were found. Skipping dump.");
        return;
    }

    switch (context.options.printer_choice)
    {
    case PrinterType::ATOM_POSITION:
        RunAtomPositionDumping(context);
        break;
    case PrinterType::MAP_VALUE:
        RunMapValueDumping(context);
        break;
    case PrinterType::GAUS_ESTIMATES:
        RunGroupGausEstimatesDumping(context);
        RunGausEstimatesDumping(context);
        break;
    case PrinterType::ATOM_OUTLIER:
        RunAtomOutlierDumping(context);
        break;
    default:
        Logger::Log(
            LogLevel::Warning,
            "Invalid printer choice input : ["
                + std::to_string(static_cast<int>(context.options.printer_choice)) + "]");
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

} // namespace

bool ExecuteResultDumpWorkflow(ResultDumpWorkflowContext & context)
{
    if (!BuildDataObjectList(context))
    {
        return false;
    }
    RunResultDump(context);
    return true;
}

} // namespace rhbm_gem::detail
