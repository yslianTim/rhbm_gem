#include "ResultDumpVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectBase.hpp"
#include "ScopeTimer.hpp"
#include "FilePathHelper.hpp"
#include "ArrayStats.hpp"
#include "AtomicPotentialEntry.hpp"
#include "GroupPotentialEntry.hpp"
#include "PotentialEntryIterator.hpp"
#include "AtomicInfoHelper.hpp"
#include "KeyPacker.hpp"
#include "Logger.hpp"

#include <memory>
#include <fstream>

ResultDumpVisitor::ResultDumpVisitor(const ResultDumpCommand::Options & options) :
    m_printer_choice{ options.printer_choice },
    m_folder_path{ FilePathHelper::EnsureTrailingSlash(options.folder_path) },
    m_map_key_tag{ "map" },
    m_model_key_tag_list{ options.model_key_tag_list }, m_map_object{ nullptr }
{

}

ResultDumpVisitor::~ResultDumpVisitor()
{
    Logger::Log(LogLevel::Debug, "ResultDumpVisitor::~ResultDumpVisitor() called.");
}

void ResultDumpVisitor::Finalize()
{
    BuildSelectedAtomList();
    switch (m_printer_choice)
    {
        case 0:
            RunAtomPositionDumping();
            break;
        case 1:
            RunMapValueDumping();
            break;
        case 2:
            RunGroupGausEstimatesDumping();
            RunGausEstimatesDumping();
            break;
        default:
            Logger::Log(LogLevel::Warning,
                        "Invalid printer choice input : [" + std::to_string(m_printer_choice) + "]");
            Logger::Log(LogLevel::Warning,
                        "Available Printer Choices:\n"
                        "  [0] AtomPositionDumping\n"
                        "  [1] MapValueDumping\n"
                        "  [2] GausEstimatesDumping");
            break;
    }
}

void ResultDumpVisitor::VisitModelObject(ModelObject * data_object)
{
    if (data_object == nullptr) return;
    for (auto & key : m_model_key_tag_list)
    {
        if (data_object->GetKeyTag() == key)
        {
            m_model_object_list.emplace_back(data_object);
            return;
        }
    }
}

void ResultDumpVisitor::VisitMapObject(MapObject * data_object)
{
    if (data_object == nullptr) return;
    if (data_object->GetKeyTag() == m_map_key_tag)
    {
        m_map_object = data_object;
    }
}

void ResultDumpVisitor::BuildSelectedAtomList(void)
{
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

void ResultDumpVisitor::RunAtomPositionDumping(void)
{
    ScopeTimer timer("ResultDumpVisitor::RunAtomPositionDumping");

    for (auto * model_object : m_model_object_list)
    {
        auto key_tag{ model_object->GetKeyTag() };
        std::string file_name{ "atom_position_list_"+ model_object->GetPdbID() +".csv" };
        std::string output_path{ m_folder_path + file_name };
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

void ResultDumpVisitor::RunMapValueDumping(void)
{
    ScopeTimer timer("ResultDumpVisitor::RunMapValueDumping");
    if (m_map_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Please give the path of map file via -m option.");
        return;
    }
    MapObject * map_object{ m_map_object };

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

        auto emd_id{ model_object->GetEmdID() };
        std::string file_name{ "map_value_list_"+ emd_id +"_"+ model_object->GetPdbID() +".csv" };
        std::string output_path{ m_folder_path + file_name };
        std::ofstream outfile(output_path);
        if (!outfile.is_open())
        {
            Logger::Log(LogLevel::Error,
                        "Could not open file " + output_path + " for writing.\n");
            return;
        }
        outfile << "GridID,X,Y,Z,MapValue\n";
        auto count{ 0 };
        for (size_t i = 0; i < map_object->GetMapValueArraySize(); i++)
        {
            auto grid_position{ map_object->GetGridPosition(i) };
            if (grid_position.at(0) < atom_range_min.at(0)) continue;
            if (grid_position.at(0) > atom_range_max.at(0)) continue;
            if (grid_position.at(1) < atom_range_min.at(1)) continue;
            if (grid_position.at(1) > atom_range_max.at(1)) continue;
            if (grid_position.at(2) < atom_range_min.at(2)) continue;
            if (grid_position.at(2) > atom_range_max.at(2)) continue;
            auto map_value{ map_object->GetMapValue(i) };
            if (map_value <= 0.0f) continue;
            outfile << i <<','
                    << grid_position.at(0) <<','<< grid_position.at(1) <<','<< grid_position.at(2) <<','
                    << map_value <<'\n';
            count++;
        }
        Logger::Log(LogLevel::Info,
                    " Selected map grid size = " + std::to_string(count) +
                    " / " + std::to_string(map_object->GetMapValueArraySize()));
    }
}

void ResultDumpVisitor::RunGausEstimatesDumping(void)
{
    ScopeTimer timer("ResultDumpVisitor::RunGausEstimatesDumping");

    for (auto * model_object : m_model_object_list)
    {
        auto key_tag{ model_object->GetKeyTag() };
        std::string file_name{ "atom_gaus_list_"+ model_object->GetPdbID() +".csv" };
        std::string output_path{ m_folder_path + file_name };
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

void ResultDumpVisitor::RunGroupGausEstimatesDumping(void)
{
    ScopeTimer timer("ResultDumpVisitor::RunGroupGausEstimatesDumping");

    auto class_key{ AtomicInfoHelper::GetResidueClassKey() };

    for (auto * model_object : m_model_object_list)
    {
        auto key_tag{ model_object->GetKeyTag() };
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

        std::string file_name{ "group_gaus_list_"+ model_object->GetPdbID() +".csv" };
        std::string output_path{ m_folder_path + file_name };
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
                        auto group_key{ KeyPackerResidueClass::Pack(residue, element, remoteness, branch, false) };
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
