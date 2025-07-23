#include "ResultDumpCommand.hpp"
#include "DataObjectManager.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "FilePathHelper.hpp"
#include "ScopeTimer.hpp"
#include "ArrayStats.hpp"
#include "AtomicPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "PotentialEntryIterator.hpp"
#include "AtomicInfoHelper.hpp"
#include "KeyPacker.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"

#include <memory>
#include <fstream>

namespace {
CommandRegistrar<ResultDumpCommand> registrar_result_dump{
    "result_dump",
    "Run result dump"};
}

ResultDumpCommand::ResultDumpCommand(void) :
    m_map_key_tag{"map"}, m_map_object{ nullptr }
{
    Logger::Log(LogLevel::Debug, "ResultDumpCommand::ResultDumpCommand() called");
}

void ResultDumpCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    Logger::Log(LogLevel::Debug, "ResultDumpCommand::RegisterCLIOptionsExtend() called");
    std::map<std::string, PrinterType> printer_map
    {
        {"0", PrinterType::ATOM_POSITION},  {"atom", PrinterType::ATOM_POSITION},
        {"1", PrinterType::MAP_VALUE},      {"map",  PrinterType::MAP_VALUE},
        {"2", PrinterType::GAUS_ESTIMATES}, {"gaus", PrinterType::GAUS_ESTIMATES}
    };
    cmd->add_option("-p,--printer", m_options.printer_choice,
        "Printer choice")
        ->required()
        ->transform(CLI::CheckedTransformer(printer_map, CLI::ignore_case));
    cmd->add_option("-k,--model-keylist", m_options.model_key_tag_list,
        "List of model key tag to be display")->required()->delimiter(',');
    cmd->add_option("-m,--map", m_options.map_file_path,
        "Map file path")->default_val(m_options.map_file_path.string());
}

bool ResultDumpCommand::Execute(void)
{
    Logger::Log(LogLevel::Debug, "ResultDumpCommand::Execute() called");
    Logger::Log(LogLevel::Info, "Total number of model object sets to be print: "
                + std::to_string(m_options.model_key_tag_list.size()));

    auto data_manager{ std::make_unique<DataObjectManager>(m_options.database_path) };
    try
    {
        if (!m_options.map_file_path.empty())
        {
            data_manager->ProcessFile(m_options.map_file_path, "map");
        }
        for (auto & key : m_options.model_key_tag_list)
        {
            data_manager->LoadDataObject(key);
        }
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error, "ResultDumpCommand::Execute() : " + std::string(e.what()));
        return false;
    }

    for (auto & key : m_options.model_key_tag_list)
    {
        auto object{ data_manager->GetTypedDataObjectPtr<ModelObject>(key) };
        m_model_object_list.emplace_back(object);
    }
    if (!m_options.map_file_path.empty())
    {
        m_map_object = data_manager->GetTypedDataObjectPtr<MapObject>("map");
    }

    BuildSelectedAtomList();
    switch (m_options.printer_choice)
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
        default:
            Logger::Log(LogLevel::Warning,
                        "Invalid printer choice input : ["
                        + std::to_string(static_cast<int>(m_options.printer_choice)) + "]");
            Logger::Log(LogLevel::Warning,
                        "Available Printer Choices:\n"
                        "  [0] AtomPositionDumping\n"
                        "  [1] MapValueDumping\n"
                        "  [2] GausEstimatesDumping");
            break;
    }
    return true;
}

bool ResultDumpCommand::ValidateOptions(void) const
{
    Logger::Log(LogLevel::Debug, "ResultDumpCommand::ValidateOptions() called");
    if (m_options.model_key_tag_list.empty())
    {
        Logger::Log(LogLevel::Error, "Model key list cannot be empty");
        return false;
    }
    if (!m_options.map_file_path.empty() &&
        !FilePathHelper::EnsureFileExists(m_options.map_file_path, "Map file"))
    {
        return false;
    }
    return true;
}

void ResultDumpCommand::SetModelKeyTagList(const std::string & value)
{
    m_options.model_key_tag_list = StringHelper::ParseListOption<std::string>(value, ',');
}

void ResultDumpCommand::BuildSelectedAtomList(void)
{
    Logger::Log(LogLevel::Debug, "ResultDumpCommand::BuildSelectedAtomList() called");
    m_selected_atom_list_map.clear();
    for (auto * model_object : m_model_object_list)
    {
        auto key_tag{ model_object->GetKeyTag() };
        for (auto & atom : model_object->GetComponentsList())
        {
            if (atom->GetAtomicPotentialEntry() == nullptr) continue;
            m_selected_atom_list_map[key_tag].emplace_back(atom.get());
        }
        Logger::Log(LogLevel::Info,
                    "Selected atoms for key tag [" + key_tag + "]: "
                    + std::to_string(m_selected_atom_list_map[key_tag].size()));
    }
}

void ResultDumpCommand::RunAtomPositionDumping(void)
{
    Logger::Log(LogLevel::Debug, "ResultDumpCommand::RunAtomPositionDumping() called");
    ScopeTimer timer("ResultDumpVisitor::RunAtomPositionDumping");

    for (auto * model_object : m_model_object_list)
    {
        auto key_tag{ model_object->GetKeyTag() };
        std::string file_name{ "atom_position_list_"+ model_object->GetPdbID() +".csv" };
        std::string output_path{
            FilePathHelper::EnsureTrailingSlash(m_options.folder_path.string()) + file_name
        };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
                        "Could not open file " + output_path + " for writing.\n");
            return;
        }
        outfile << "SerialID,X,Y,Z\n";
        for (auto & atom : m_selected_atom_list_map.at(key_tag))
        {
            outfile << atom->GetSerialID() <<','
                    << atom->GetPosition().at(0) <<','
                    << atom->GetPosition().at(1) <<','
                    << atom->GetPosition().at(2) <<'\n';
        }
        outfile.close();
        Logger::Log(LogLevel::Info, "Output file: " + output_path);
    }
}

void ResultDumpCommand::RunMapValueDumping(void)
{
    Logger::Log(LogLevel::Debug, "ResultDumpCommand::RunMapValueDumping() called");
    ScopeTimer timer("ResultDumpVisitor::RunMapValueDumping");
    if (m_map_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Please give the path of map file via -m option.");
        return;
    }

    auto margin{ 3.0f };

    for (auto * model_object : m_model_object_list)
    {
        auto key_tag{ model_object->GetKeyTag() };
        auto atom_size{ m_selected_atom_list_map.at(key_tag).size() };
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

        std::string file_name{
            "map_value_list_"+ model_object->GetEmdID() +"_"+ model_object->GetPdbID() +".csv"
        };
        std::string output_path{
            FilePathHelper::EnsureTrailingSlash(m_options.folder_path.string()) + file_name
        };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
                        "Could not open file " + output_path + " for writing.\n");
            return;
        }
        outfile << "GridID,X,Y,Z,MapValue\n";
        auto count{ 0 };
        for (size_t i = 0; i < m_map_object->GetMapValueArraySize(); i++)
        {
            auto grid_position{ m_map_object->GetGridPosition(i) };
            if (grid_position.at(0) < atom_range_min.at(0)) continue;
            if (grid_position.at(0) > atom_range_max.at(0)) continue;
            if (grid_position.at(1) < atom_range_min.at(1)) continue;
            if (grid_position.at(1) > atom_range_max.at(1)) continue;
            if (grid_position.at(2) < atom_range_min.at(2)) continue;
            if (grid_position.at(2) > atom_range_max.at(2)) continue;
            auto map_value{ m_map_object->GetMapValue(i) };
            if (map_value <= 0.0f) continue;
            outfile << i <<','
                    << grid_position.at(0) <<','<< grid_position.at(1) <<','<< grid_position.at(2) <<','
                    << map_value <<'\n';
            count++;
        }
        Logger::Log(LogLevel::Info,
                    " Selected map grid size = " + std::to_string(count) +
                    " / " + std::to_string(m_map_object->GetMapValueArraySize()));
    }
}

void ResultDumpCommand::RunGausEstimatesDumping(void)
{
    Logger::Log(LogLevel::Debug, "ResultDumpCommand::RunGausEstimatesDumping() called");
    ScopeTimer timer("ResultDumpCommand::RunGausEstimatesDumping");

    for (auto * model_object : m_model_object_list)
    {
        auto key_tag{ model_object->GetKeyTag() };
        std::string file_name{ "atom_gaus_list_"+ model_object->GetPdbID() +".csv" };
        std::string output_path{ m_options.folder_path.string() + file_name };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
                        "Could not open file " + output_path + " for writing.\n");
            return;
        }
        outfile << "SerialID,Amplitude,Width\n";
        for (auto & atom : m_selected_atom_list_map.at(key_tag))
        {
            auto entry{ atom->GetAtomicPotentialEntry() };
            outfile << atom->GetSerialID() <<','
                    << entry->GetAmplitudeEstimateMDPDE() <<','
                    << entry->GetWidthEstimateMDPDE() <<'\n';
        }
        outfile.close();
        Logger::Log(LogLevel::Info, "Output file: " + output_path);
    }
}

void ResultDumpCommand::RunGroupGausEstimatesDumping(void)
{
    Logger::Log(LogLevel::Debug, "ResultDumpCommand::RunGroupGausEstimatesDumping() called");
    ScopeTimer timer("ResultDumpCommand::RunGroupGausEstimatesDumping");

    auto class_key{ AtomicInfoHelper::GetResidueClassKey() };

    for (auto * model_object : m_model_object_list)
    {
        auto key_tag{ model_object->GetKeyTag() };
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

        std::string file_name{ "group_gaus_list_"+ model_object->GetPdbID() +".csv" };
        std::string output_path{ m_options.folder_path.string() + file_name };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
                        "Could not open file " + output_path + " for writing.\n");
            return;
        }
        outfile << "Residue,Element,Remoteness,Branch,Amplitude,Width\n";
        for (auto & residue : AtomicInfoHelper::GetStandardResidueList())
        {
            auto residue_name{ AtomicInfoHelper::GetLabel(residue) };
            for (auto & element : AtomicInfoHelper::GetStandardElementList())
            {
                auto element_name{ AtomicInfoHelper::GetLabel(element) };
                for (auto & remoteness : AtomicInfoHelper::GetStandardRemotenessList())
                {
                    auto remoteness_name{ AtomicInfoHelper::GetLabel(remoteness) };
                    for (auto & branch : AtomicInfoHelper::GetStandardBranchList())
                    {
                        auto branch_name{ AtomicInfoHelper::GetLabel(branch) };
                        auto group_key{
                            KeyPackerResidueClass::Pack(residue, element, remoteness, branch, false)
                        };
                        if (entry_iter->IsAvailableGroupKey(group_key, class_key) == false) continue;
                        outfile << residue_name <<','
                                << element_name <<','
                                << remoteness_name <<','
                                << branch_name <<','
                                << entry_iter->GetGausEstimatePrior(group_key, class_key, 0) <<','
                                << entry_iter->GetGausEstimatePrior(group_key, class_key, 1) <<'\n';
                    }
                }
            }
        }

        outfile.close();
        Logger::Log(LogLevel::Info, "Output file: " + output_path);
    }
}
