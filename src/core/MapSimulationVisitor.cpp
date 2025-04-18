#include "MapSimulationVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "AtomSelector.hpp"
#include "ScopeTimer.hpp"
#include "FilePathHelper.hpp"
#include "ElectricPotential.hpp"
#include "KDTreeAlgorithm.hpp"
#include "AminoAcidInfoHelper.hpp"
#include "ArrayStats.hpp"
#include "StringHelper.hpp"

MapSimulationVisitor::MapSimulationVisitor(AtomSelector * atom_selector) :
    m_thread_size{ 1 },
    m_atom_selector{ atom_selector }
{

}

MapSimulationVisitor::~MapSimulationVisitor()
{

}

void MapSimulationVisitor::VisitAtomObject(AtomObject * data_object)
{
    bool selected_flag
    {
        m_atom_selector->GetSelectionFlag(
            data_object->GetChainID(),
            data_object->GetIndicator(),
            data_object->GetResidue(),
            data_object->GetElement(),
            data_object->GetRemoteness(),
            data_object->GetBranch()
        )
    };
    data_object->SetSelectedFlag(selected_flag);
}

void MapSimulationVisitor::VisitModelObject(ModelObject * data_object)
{
    ScopeTimer timer("MapSimulationVisitor::VisitModelObject");
    auto selected_atom_size{ static_cast<size_t>(data_object->GetNumberOfSelectedAtom()) };
    const auto & atom_list{ data_object->GetComponentsList() };
    m_selected_atom_list.clear();
    m_selected_atom_list.reserve(selected_atom_size);
    for (auto & atom : atom_list)
    {
        if (atom->GetSelectedFlag() == false) continue;
        if (atom->IsUnknownAtom() == true) continue;
        m_selected_atom_list.emplace_back(atom.get());
    }
    std::cout <<" Number of selected atoms for simulation = "<< selected_atom_size << std::endl;
}

void MapSimulationVisitor::VisitMapObject(MapObject * data_object)
{
    (void) data_object;
}

void MapSimulationVisitor::Analysis(DataObjectManager * data_manager)
{
    std::cout <<"- Analysis..." << std::endl;
    try
    {
        const auto & model_object{ data_manager->GetDataObjectRef(m_model_key_tag) };
        model_object->Accept(this);

        for (auto & blurring_width : m_blurring_width_list)
        {
            auto map_key_tag{
                m_model_key_tag + "_bw" +
                StringHelper::ToStringWithPrecision<double>(blurring_width, 2)
            };
            data_manager->AddDataObject(map_key_tag, CreateSimulatedMapObject(blurring_width));

            auto extension{ std::string(".map") };
            auto file_name{ std::string("sim_map_") + map_key_tag + extension };
            auto output_file_name{ FilePathHelper::EnsureTrailingSlash(m_folder_path) + file_name };
            data_manager->ProduceFile(output_file_name, map_key_tag);
        }
    }
    catch(const std::exception & e)
    {
        std::cerr << e.what() << std::endl;
    }
}

std::unique_ptr<MapObject> MapSimulationVisitor::CreateSimulatedMapObject(double blurring_width)
{
    ScopeTimer timer("MapSimulationVisitor::CreateSimulatedMapObject");
    std::cout <<"  - Create simulated map object with blurring width = "<< blurring_width << std::endl;

    auto electric_potential{ std::make_unique<ElectricPotential>() };
    electric_potential->SetBlurringWidth(blurring_width);
    electric_potential->SetModelChoice(m_potential_model_choice);

    auto kd_tree_root{ KDTreeAlgorithm<AtomObject>::BuildKDTree(m_selected_atom_list, 0) };
    
    std::array<float, 3> grid_spacing{ m_grid_spacing, m_grid_spacing, m_grid_spacing };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    std::array<int, 3> grid_size{ CalculateGridSize(grid_spacing, origin) };
    auto map_object{ std::make_unique<MapObject>(grid_size, grid_spacing, origin) };

    auto voxel_size{ map_object->GetMapValueArraySize() };
    auto map_value_array{ std::make_unique<float[]>(voxel_size) };

    std::cout <<"  - Start map value array production ..."<< std::endl;
    #ifdef USE_OPENMP
    #pragma omp parallel for num_threads(m_thread_size)
    #endif
    for (size_t i = 0; i < voxel_size; i++)
    {
        auto query_pseudo_atom{ std::make_unique<AtomObject>() };
        query_pseudo_atom->SetPosition(map_object->GetGridPosition(i));
        auto in_range_atom_list{
            KDTreeAlgorithm<AtomObject>::RangeSearch(
                kd_tree_root.get(), query_pseudo_atom.get(), m_cutoff_distance)
        };

        for (const auto & atom : in_range_atom_list)
        {
            auto distance{
                ArrayStats<float>::ComputeNorm(query_pseudo_atom->GetPosition(), atom->GetPosition())
            };
            auto charge{ (m_partial_charge_choice == 0) ? 0.0 :
                AminoAcidInfoHelper::GetPartialCharge(
                    atom->GetResidueType(),
                    atom->GetElementType(),
                    atom->GetRemotenessType(),
                    atom->GetBranchType())
            };
            map_value_array[i] += static_cast<float>(
                electric_potential->GetPotentialValue(atom->GetElementType(), distance, charge)
            );
        }
    }
    std::cout <<"  - End map value array production."<< std::endl;
    map_object->SetMapValueArray(std::move(map_value_array));
    map_object->Display();

    return map_object;
}

std::array<int, 3> MapSimulationVisitor::CalculateGridSize(
    const std::array<float, 3> & grid_spacing, const std::array<float, 3> & origin)
{
    auto selected_atom_size{ m_selected_atom_list.size() };
    std::array<float, 3> atom_range_max;
    std::vector<float> x_list, y_list, z_list;
    x_list.reserve(selected_atom_size);
    y_list.reserve(selected_atom_size);
    z_list.reserve(selected_atom_size);
    for (auto & atom : m_selected_atom_list)
    {
        x_list.emplace_back(atom->GetPosition().at(0));
        y_list.emplace_back(atom->GetPosition().at(1));
        z_list.emplace_back(atom->GetPosition().at(2));
    }
    atom_range_max.at(0) = ArrayStats<float>::ComputeMax(x_list.data(), selected_atom_size) + static_cast<float>(m_cutoff_distance);
    atom_range_max.at(1) = ArrayStats<float>::ComputeMax(y_list.data(), selected_atom_size) + static_cast<float>(m_cutoff_distance);
    atom_range_max.at(2) = ArrayStats<float>::ComputeMax(z_list.data(), selected_atom_size) + static_cast<float>(m_cutoff_distance);

    auto grid_size_x{ static_cast<int>(std::ceil((atom_range_max.at(0) - origin.at(0)) / grid_spacing.at(0))) };
    auto grid_size_y{ static_cast<int>(std::ceil((atom_range_max.at(1) - origin.at(1)) / grid_spacing.at(1))) };
    auto grid_size_z{ static_cast<int>(std::ceil((atom_range_max.at(2) - origin.at(2)) / grid_spacing.at(2))) };
    std::cout <<"Grid size = ["<< grid_size_x << "," << grid_size_y << "," << grid_size_z << "]" << std::endl;
    return std::array{ grid_size_x, grid_size_y, grid_size_z };
}