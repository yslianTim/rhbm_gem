#include "MapSimulationVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "ScopeTimer.hpp"
#include "FilePathHelper.hpp"
#include "ElectricPotential.hpp"
#include "KDTreeAlgorithm.hpp"
#include "AminoAcidInfoHelper.hpp"
#include "ArrayStats.hpp"
#include "StringHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "Logger.hpp"

#include <algorithm>
#include <limits>

MapSimulationVisitor::MapSimulationVisitor(void) :
    m_thread_size{ 1 }, m_kd_tree_root{ nullptr }
{

}

MapSimulationVisitor::~MapSimulationVisitor()
{
    m_kd_tree_root.reset();
}

void MapSimulationVisitor::VisitAtomObject(AtomObject * data_object)
{
    data_object->SetSelectedFlag(true);
}

void MapSimulationVisitor::VisitModelObject(ModelObject * data_object)
{
    ScopeTimer timer("MapSimulationVisitor::VisitModelObject");
    auto selected_atom_size{ static_cast<size_t>(data_object->GetNumberOfSelectedAtom()) };
    const auto & atom_list{ data_object->GetComponentsList() };
    m_selected_atom_list.clear();
    m_selected_atom_list.reserve(selected_atom_size);
    m_atom_charge_map.clear();
    m_atom_charge_map.reserve(selected_atom_size);
    for (auto & atom : atom_list)
    {
        if (atom->IsUnknownAtom() == true)
        {
            Logger::Log(LogLevel::Warning,
                        "Unknown atom found in the model object: Serial-ID = "
                        + std::to_string(atom->GetSerialID()) +
                        ", this atom will be ignored in the simulation.");
            continue;
        }
        m_selected_atom_list.emplace_back(atom.get());
        auto charge{ 0.0 };
        switch (m_partial_charge_choice)
        {
            case 0: // Neutral
                charge = 0.0;
                break;
            case 1: // Partial Charge
                charge = AminoAcidInfoHelper::GetPartialCharge(
                    atom->GetResidue(),
                    atom->GetElement(),
                    atom->GetRemoteness(),
                    atom->GetBranch(),
                    atom->GetStructure());
                break;
            case 2: // Partial Charge (AMBER table)
                charge = AminoAcidInfoHelper::GetPartialCharge(
                    atom->GetResidue(),
                    atom->GetElement(),
                    atom->GetRemoteness(),
                    atom->GetBranch(),
                    atom->GetStructure(), true);
                break;
            default:
                Logger::Log(LogLevel::Error,
                            "Invalid partial charge choice: "
                            + std::to_string(m_partial_charge_choice));
                break;
        }
        m_atom_charge_map.emplace(atom->GetSerialID(), charge);
    }
    m_kd_tree_root = KDTreeAlgorithm<AtomObject>::BuildKDTree(m_selected_atom_list, 0);
    m_pdb_id = data_object->GetPdbID();
    Logger::Log(LogLevel::Info,
                "Number of selected atoms to be simulated = "
                + std::to_string(m_selected_atom_list.size()) +" / "+
                std::to_string(atom_list.size()) + " atoms.");
}

void MapSimulationVisitor::VisitMapObject(MapObject * data_object)
{
    (void) data_object;
}

void MapSimulationVisitor::Analysis(DataObjectManager * data_manager)
{
    Logger::Log(LogLevel::Info, "- MapSimulationVisitor::Analysis");
    try
    {
        auto model_object{ data_manager->GetDataObjectRef(m_model_key_tag) };
        model_object->Accept(this);

        for (auto & blurring_width : m_blurring_width_list)
        {
            auto map_key_tag{
                m_pdb_id + "_bw" +
                StringHelper::ToStringWithPrecision<double>(blurring_width, 2)
            };
            data_manager->AddDataObject(map_key_tag, CreateSimulatedMapObject(blurring_width));

            auto extension{ std::string(".map") };
            auto file_name{ m_map_file_name +"_"+ map_key_tag + extension };
            auto output_file_name{ FilePathHelper::EnsureTrailingSlash(m_folder_path) + file_name };
            data_manager->ProduceFile(output_file_name, map_key_tag);
        }
    }
    catch(const std::exception & e)
    {
        Logger::Log(LogLevel::Error, e.what());
    }
}

std::unique_ptr<MapObject> MapSimulationVisitor::CreateSimulatedMapObject(double blurring_width)
{
    ScopeTimer timer("MapSimulationVisitor::CreateSimulatedMapObject");
    Logger::Log(LogLevel::Info,
                std::string("  - Create simulated map object with blurring width = ") +
                StringHelper::ToStringWithPrecision<double>(blurring_width, 2));

    auto electric_potential{ std::make_unique<ElectricPotential>() };
    electric_potential->SetBlurringWidth(blurring_width);
    electric_potential->SetModelChoice(m_potential_model_choice);
    
    std::array<float, 3> grid_spacing{ m_grid_spacing, m_grid_spacing, m_grid_spacing };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    std::array<int, 3> grid_size{ CalculateGridSize(grid_spacing, origin) };
    auto map_object{ std::make_unique<MapObject>(grid_size, grid_spacing, origin) };
    map_object->SetThreadSize(static_cast<int>(m_thread_size));

    auto voxel_size{ map_object->GetMapValueArraySize() };
    auto map_value_array{ std::make_unique<float[]>(voxel_size) };
    std::fill_n(map_value_array.get(), voxel_size, 0.0f);

    Logger::Log(LogLevel::Info, "  - Start map value array production ...");
    std::vector<AtomObject*> in_range_atom_list;
    #ifdef USE_OPENMP
    #pragma omp parallel for num_threads(m_thread_size) private(in_range_atom_list)
    #endif
    for (size_t i = 0; i < voxel_size; i++)
    {
        AtomObject query_atom;
        query_atom.SetPosition(map_object->GetGridPosition(i));
        KDTreeAlgorithm<AtomObject>::RangeSearch(
            m_kd_tree_root.get(), &query_atom, m_cutoff_distance, in_range_atom_list);

        for (const auto & atom : in_range_atom_list)
        {
            auto distance{
                ArrayStats<float>::ComputeNorm(query_atom.GetPosition(), atom->GetPosition())
            };
            auto charge{ 0.0 };
            auto iter{ m_atom_charge_map.find(atom->GetSerialID()) };
            if (iter != m_atom_charge_map.end()) charge = iter->second;
            map_value_array[i] += static_cast<float>(
                electric_potential->GetPotentialValue(atom->GetElement(), distance, charge)
            );
        }
    }
    Logger::Log(LogLevel::Info, "  - End map value array production.");
    map_object->SetMapValueArray(std::move(map_value_array));
    map_object->Display();

    return map_object;
}

std::array<int, 3> MapSimulationVisitor::CalculateGridSize(
    const std::array<float, 3> & grid_spacing, const std::array<float, 3> & origin)
{
    auto selected_atom_size{ m_selected_atom_list.size() };
    if (selected_atom_size == 0)
    {
        Logger::Log(LogLevel::Warning, "Warning: no atoms selected. Grid size is set to [1,1,1].");
        return std::array{ 1, 1, 1 };
    }
    std::array<float, 3> atom_range_max{
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest()
    };
    for (auto & atom : m_selected_atom_list)
    {
        const auto & pos{ atom->GetPositionRef() };
        atom_range_max.at(0) = std::max(atom_range_max.at(0), pos.at(0));
        atom_range_max.at(1) = std::max(atom_range_max.at(1), pos.at(1));
        atom_range_max.at(2) = std::max(atom_range_max.at(2), pos.at(2));
    }
    atom_range_max.at(0) += static_cast<float>(m_cutoff_distance);
    atom_range_max.at(1) += static_cast<float>(m_cutoff_distance);
    atom_range_max.at(2) += static_cast<float>(m_cutoff_distance);

    auto grid_size_x{ static_cast<int>(std::ceil((atom_range_max.at(0) - origin.at(0)) / grid_spacing.at(0))) };
    auto grid_size_y{ static_cast<int>(std::ceil((atom_range_max.at(1) - origin.at(1)) / grid_spacing.at(1))) };
    auto grid_size_z{ static_cast<int>(std::ceil((atom_range_max.at(2) - origin.at(2)) / grid_spacing.at(2))) };
    Logger::Log(LogLevel::Info,
                "Grid size = [" + std::to_string(grid_size_x) + "," +
                std::to_string(grid_size_y) + "," +
                std::to_string(grid_size_z) + "]");
    return std::array{ grid_size_x, grid_size_y, grid_size_z };
}
