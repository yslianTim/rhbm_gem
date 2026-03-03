#include "ResultDumpCommand.hpp"
#include "DataObjectManager.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "FilePathHelper.hpp"
#include "ScopeTimer.hpp"
#include "ArrayStats.hpp"
#include "LocalPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "PotentialEntryIterator.hpp"
#include "ChemicalDataHelper.hpp"
#include "ComponentHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "AtomClassifier.hpp"
#include "StringHelper.hpp"
#include "ModelFileWriter.hpp"
#include "ChimeraXHelper.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"
#include "AtomKeySystem.hpp"

#include <memory>
#include <fstream>

namespace {
rhbm_gem::CommandRegistrar<rhbm_gem::ResultDumpCommand> registrar_result_dump{
    "result_dump",
    "Run result dump"};
}

namespace rhbm_gem {

ResultDumpCommand::ResultDumpCommand(void) :
    CommandBase(),
    m_map_key_tag{"map"}, m_map_object{ nullptr }
{
}

ResultDumpCommand::~ResultDumpCommand()
{
}

void ResultDumpCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    std::map<std::string, PrinterType> printer_map
    {
        {"0", PrinterType::ATOM_POSITION},  {"atom_pos", PrinterType::ATOM_POSITION},
        {"1", PrinterType::MAP_VALUE},      {"map",      PrinterType::MAP_VALUE},
        {"2", PrinterType::GAUS_ESTIMATES}, {"gaus",     PrinterType::GAUS_ESTIMATES},
        {"3", PrinterType::ATOM_OUTLIER},   {"atom_out", PrinterType::ATOM_OUTLIER}
    };
    cmd->add_option_function<PrinterType>("-p,--printer",
        [&](PrinterType value) { SetPrinterChoice(value); },
        "Printer choice")->required()
        ->transform(CLI::CheckedTransformer(printer_map, CLI::ignore_case));
    cmd->add_option_function<std::string>("-k,--model-keylist",
        [&](const std::string & value) { SetModelKeyTagList(value); },
        "List of model key tag to be display")->required();
    cmd->add_option_function<std::string>("-m,--map",
        [&](const std::string & value) { SetMapFilePath(value); },
        "Map file path")->default_val(m_options.map_file_path.string());
}

bool ResultDumpCommand::Execute(void)
{
    if (BuildDataObjectList() == false) return false;
    RunResultDump();
    return true;
}

void ResultDumpCommand::SetPrinterChoice(PrinterType value)
{
    m_options.printer_choice = value;
}

void ResultDumpCommand::SetMapFilePath(const std::filesystem::path & path)
{
    m_options.map_file_path = path;
    if (!m_options.map_file_path.empty() &&
        !FilePathHelper::EnsureFileExists(m_options.map_file_path, "Map file"))
    {
        Logger::Log(LogLevel::Error,
            "Map file does not exist: " + m_options.map_file_path.string());
        m_valiate_options = false;
    }
}

void ResultDumpCommand::SetModelKeyTagList(const std::string & value)
{
    m_options.model_key_tag_list = StringHelper::ParseListOption<std::string>(value, ',');
    if (m_options.model_key_tag_list.empty())
    {
        Logger::Log(LogLevel::Error, "Model key list cannot be empty");
        m_valiate_options = false;
    }
}

bool ResultDumpCommand::BuildDataObjectList(void)
{
    ScopeTimer timer("ResultDumpCommand::BuildDataObjectList");
    auto data_manager{ GetDataManagerPtr() };
    data_manager->SetDatabaseManager(m_options.database_path);
    try
    {
        if (!m_options.map_file_path.empty())
        {
            data_manager->ProcessFile(m_options.map_file_path, "map");
            m_map_object = data_manager->GetTypedDataObject<MapObject>("map");
        }
        m_selected_atom_list_map.clear();
        for (auto & key : m_options.model_key_tag_list)
        {
            data_manager->LoadDataObject(key);
            auto model_object{ data_manager->GetTypedDataObject<ModelObject>(key) };
            m_model_object_list.emplace_back(model_object);
            for (auto & atom : model_object->GetAtomList())
            {
                if (atom->GetLocalPotentialEntry() == nullptr) continue;
                m_selected_atom_list_map[key].emplace_back(atom.get());
            }
            Logger::Log(LogLevel::Info,
                "Selected atoms for key tag [" + key + "]: "
                + std::to_string(m_selected_atom_list_map[key].size()));
        }
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "ResultDumpCommand::BuildDataObjectList : " + std::string(e.what()));
        return false;
    }
    return true;
}

void ResultDumpCommand::RunResultDump(void)
{
    ScopeTimer timer("ResultDumpCommand::RunResultDump");
    Logger::Log(LogLevel::Info,
        "Total number of model object sets to be dump: "
        + std::to_string(m_options.model_key_tag_list.size()));
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
        case PrinterType::ATOM_OUTLIER:
            RunAtomOutlierDumping();
            break;
        default:
            Logger::Log(LogLevel::Warning,
                        "Invalid printer choice input : ["
                        + std::to_string(static_cast<int>(m_options.printer_choice)) + "]");
            Logger::Log(LogLevel::Warning,
                        "Available Printer Choices:\n"
                        "  [0] AtomPositionDumping\n"
                        "  [1] MapValueDumping\n"
                        "  [2] GausEstimatesDumping\n"
                        "  [3] AtomOutlierDumping");
            break;
    }
}

void ResultDumpCommand::RunAtomOutlierDumping(void)
{

    for (const auto & model_object : m_model_object_list)
    {
        auto key_tag{ model_object->GetKeyTag() };
        std::unordered_map<std::string, std::vector<std::array<float,3>>> outlier_position_map;
        std::string file_name{ "atom_outlier_list_"+ model_object->GetPdbID() +".csv" };
        std::filesystem::path output_path{ m_options.folder_path / file_name };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
                        "Could not open file " + output_path.string() + " for writing.\n");
            return;
        }
        outfile << "SerialID,ClassType,Residue,Element,Spot\n";
        for (auto & atom : m_selected_atom_list_map.at(key_tag))
        {
            for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
            {
                const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
                if (atom->GetLocalPotentialEntry()->GetOutlierTag(class_key) == false) continue;
                outfile << atom->GetSerialID() <<','
                        << class_key <<','
                        << ChemicalDataHelper::GetLabel(atom->GetResidue()) <<','
                        << ChemicalDataHelper::GetLabel(atom->GetElement()) <<','
                        << atom->GetAtomID() <<'\n';
                if (atom->IsMainChainAtom()) // Only output main chain atom for outlier display
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
            std::string class_file_name{ "atom_outlier_position_" + class_key + "_" + model_object->GetPdbID() + ".cmm" };
            std::filesystem::path class_output_path{ m_options.folder_path / class_file_name };
            std::string marker_set_name{ class_key + "_" + model_object->GetPdbID() };
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

void ResultDumpCommand::RunAtomPositionDumping(void)
{

    for (const auto & model_object : m_model_object_list)
    {
        auto key_tag{ model_object->GetKeyTag() };
        std::string file_name{ "atom_position_list_"+ model_object->GetPdbID() +".csv" };
        std::filesystem::path output_path{ m_options.folder_path / file_name };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
                        "Could not open file " + output_path.string() + " for writing.\n");
            return;
        }
        outfile << "SerialID,X,Y,Z,Residue,Element,Spot\n";
        for (auto & atom : m_selected_atom_list_map.at(key_tag))
        {
            outfile << atom->GetSerialID() <<','
                    << atom->GetPosition().at(0) <<','
                    << atom->GetPosition().at(1) <<','
                    << atom->GetPosition().at(2) <<','
                    << ChemicalDataHelper::GetLabel(atom->GetResidue()) <<','
                    << ChemicalDataHelper::GetLabel(atom->GetElement()) <<','
                    << atom->GetAtomID() <<'\n';
        }
        outfile.close();
        Logger::Log(LogLevel::Info, "Output file: " + output_path.string());
    }
}

void ResultDumpCommand::RunMapValueDumping(void)
{
    if (m_map_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Please give the path of map file via -m option.");
        return;
    }

    auto margin{ 3.0f };

    for (const auto & model_object : m_model_object_list)
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
        std::filesystem::path output_path{ m_options.folder_path / file_name };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
                        "Could not open file " + output_path.string() + " for writing.\n");
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

    for (const auto & model_object : m_model_object_list)
    {
        auto key_tag{ model_object->GetKeyTag() };
        std::string csv_file_name{ "atom_gaus_list_"+ model_object->GetPdbID() +".csv" };
        std::filesystem::path output_csv_file{ m_options.folder_path / csv_file_name };
        std::ofstream outfile(output_csv_file);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
                        "Could not open file " + output_csv_file.string() + " for writing.\n");
            return;
        }

        outfile << "SerialID,Amplitude,Width,X,Y,Z,Residue,Element,Spot\n";
        for (auto & atom : m_selected_atom_list_map.at(key_tag))
        {
            auto entry{ atom->GetLocalPotentialEntry() };
            outfile << atom->GetSerialID() <<','
                    << entry->GetAmplitudeEstimateMDPDE() <<','
                    << entry->GetWidthEstimateMDPDE() <<','
                    << atom->GetPosition().at(0) <<','
                    << atom->GetPosition().at(1) <<','
                    << atom->GetPosition().at(2) <<','
                    << ChemicalDataHelper::GetLabel(atom->GetResidue()) <<','
                    << ChemicalDataHelper::GetLabel(atom->GetElement()) <<','
                    << atom->GetAtomID() <<'\n';
        }
        outfile.close();
        Logger::Log(LogLevel::Info, "Output file: " + output_csv_file.string());

        // Output result to mmCIF file
        std::string amplitude_file_name{ model_object->GetPdbID() +"_gaus_amplitude.cif" };
        std::filesystem::path output_amplitude_cif_file{ m_options.folder_path / amplitude_file_name };
        ModelFileWriter amplitude_writer{ output_amplitude_cif_file.string(), model_object.get(), 0 };
        amplitude_writer.Write();
        Logger::Log(LogLevel::Info, "Output file: " + output_amplitude_cif_file.string());

        std::string width_file_name{ model_object->GetPdbID() +"_gaus_width.cif" };
        std::filesystem::path output_width_cif_file{ m_options.folder_path / width_file_name };
        ModelFileWriter width_writer{ output_width_cif_file.string(), model_object.get(), 1 };
        width_writer.Write();
        Logger::Log(LogLevel::Info, "Output file: " + output_width_cif_file.string());

        std::string intensity_file_name{ model_object->GetPdbID() +"_gaus_intensity.cif" };
        std::filesystem::path output_intensity_cif_file{ m_options.folder_path / intensity_file_name };
        ModelFileWriter intensity_writer{ output_intensity_cif_file.string(), model_object.get(), 2 };
        intensity_writer.Write();
        Logger::Log(LogLevel::Info, "Output file: " + output_intensity_cif_file.string());
    }
}

void ResultDumpCommand::RunGroupGausEstimatesDumping(void)
{

    auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };

    for (const auto & model_object : m_model_object_list)
    {
        auto key_tag{ model_object->GetKeyTag() };
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object.get()) };

        std::string file_name{ "group_gaus_list_"+ model_object->GetPdbID() +".csv" };
        std::filesystem::path output_path{ m_options.folder_path / file_name };
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
            auto residue_name{ ChemicalDataHelper::GetLabel(residue) };
            auto component_key{ static_cast<ComponentKey>(residue) };
            for (auto & spot : ComponentHelper::GetSpotList(residue))
            {
                auto atom_key{ static_cast<uint16_t>(spot) };
                auto group_key{ AtomClassifier::GetGroupKeyInClass(component_key, atom_key) };
                auto atom_id{ model_object->GetAtomKeySystemPtr()->GetAtomId(atom_key) };
                if (entry_iter->IsAvailableAtomGroupKey(group_key, class_key) == false) continue;
                outfile << residue_name <<','
                        << atom_id <<','
                        << entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 0) <<','
                        << entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 1) <<'\n';
            }
        }

        outfile.close();
        Logger::Log(LogLevel::Info, "Output file: " + output_path.string());
    }
}

} // namespace rhbm_gem
