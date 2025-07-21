#include "ResultDumpVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
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

ResultDumpVisitor::ResultDumpVisitor() :
    m_printer_choice{ 0 }
{

}

void ResultDumpVisitor::VisitAtomObject(AtomObject * data_object)
{
    if (data_object == nullptr) return;
}

void ResultDumpVisitor::VisitModelObject(ModelObject * data_object)
{
    if (data_object == nullptr) return;
}

void ResultDumpVisitor::VisitMapObject(MapObject * data_object)
{
    if (data_object == nullptr) return;
}

void ResultDumpVisitor::Analysis(DataObjectManager * data_manager)
{
    if (data_manager == nullptr)
    {
        Logger::Log(LogLevel::Warning, "Data manager is null, please check the input.");
        return;
    }
    BuildSelectedAtomList(data_manager);
    
    switch (m_printer_choice)
    {
        case 0:
            RunAtomPositionDumping(data_manager);
            break;
        case 1:
            RunMapValueDumping(data_manager);
            break;
        case 2:
            RunGroupGausEstimatesDumping(data_manager);
            RunGausEstimatesDumping(data_manager);
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

void ResultDumpVisitor::SetFolderPath(const std::string & value)
{
    m_folder_path = FilePathHelper::EnsureTrailingSlash(value);
}

void ResultDumpVisitor::BuildSelectedAtomList(DataObjectManager * data_manager)
{
    if (data_manager == nullptr) return;
    for (auto & key_tag : m_model_key_tag_list)
    {
        auto data_object{ data_manager->GetDataObjectRef(key_tag) };
        auto model_object{ dynamic_cast<ModelObject *>(data_object) };
        if (model_object == nullptr) continue;
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

void ResultDumpVisitor::RunAtomPositionDumping(DataObjectManager * data_manager)
{
    ScopeTimer timer("ResultDumpVisitor::RunAtomPositionDumping");
    if (data_manager == nullptr) return;

    for (auto & key_tag : m_model_key_tag_list)
    {
        auto data_object{ data_manager->GetDataObjectRef(key_tag) };
        auto model_object{ dynamic_cast<ModelObject *>(data_object) };
        if (model_object == nullptr) continue;
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

void ResultDumpVisitor::RunMapValueDumping(DataObjectManager * data_manager)
{
    ScopeTimer timer("ResultDumpVisitor::RunMapValueDumping");
    if (data_manager == nullptr) return;

    MapObject * map_object{ nullptr };
    try
    {
        auto data_object{ data_manager->GetDataObjectRef(m_map_key_tag) };
        map_object = dynamic_cast<MapObject *>(data_object);
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error, e.what());
        Logger::Log(LogLevel::Error, "Please give the path of map file via -m option.");
        return;
    }
    if (map_object == nullptr) return;

    auto margin{ 3.0f };

    for (auto & key_tag : m_model_key_tag_list)
    {
        auto data_object{ data_manager->GetDataObjectRef(key_tag) };
        auto model_object{ dynamic_cast<ModelObject *>(data_object) };
        if (model_object == nullptr) continue;
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

void ResultDumpVisitor::RunGausEstimatesDumping(DataObjectManager * data_manager)
{
    ScopeTimer timer("ResultDumpVisitor::RunGausEstimatesDumping");
    if (data_manager == nullptr) return;

    for (auto & key_tag : m_model_key_tag_list)
    {
        auto data_object{ data_manager->GetDataObjectRef(key_tag) };
        auto model_object{ dynamic_cast<ModelObject *>(data_object) };
        if (model_object == nullptr) continue;
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

void ResultDumpVisitor::RunGroupGausEstimatesDumping(DataObjectManager * data_manager)
{
    ScopeTimer timer("ResultDumpVisitor::RunGroupGausEstimatesDumping");
    if (data_manager == nullptr) return;

    auto class_key{ AtomicInfoHelper::GetResidueClassKey() };

    for (auto & key_tag : m_model_key_tag_list)
    {
        auto data_object{ data_manager->GetDataObjectRef(key_tag) };
        auto model_object{ dynamic_cast<ModelObject *>(data_object) };
        if (model_object == nullptr) continue;
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
