#include "MapSimulationVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "AtomSelector.hpp"
#include "ScopeTimer.hpp"
#include "FilePathHelper.hpp"

#include <array>

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
        m_selected_atom_list.emplace_back(atom.get());
    }
    std::cout <<" Number of selected atom = "<< selected_atom_size << std::endl;
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
            auto map_key_tag{ m_model_key_tag + "_bw" + std::to_string(blurring_width) };
            data_manager->AddDataObject(map_key_tag, CreateSimulatedMapObject(blurring_width));

            auto output_file_name{ FilePathHelper::EnsureTrailingSlash(m_folder_path) + map_key_tag + ".mrc" };
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
    std::cout << "  - Create simulated map object with blurring width = "<< blurring_width << std::endl;
    
    std::array<int, 3> grid_size;
    std::array<float, 3> grid_spacing;
    std::array<float, 3> origin;
    std::unique_ptr<float[]> map_value_array;

    size_t voxel_size{ 0 };

    #ifdef USE_OPENMP
    #pragma omp parallel for num_threads(m_thread_size)
    #endif
    for (size_t i = 0; i < voxel_size; i++)
    {

    }

    return std::make_unique<MapObject>(grid_size, grid_spacing, origin, std::move(map_value_array));
}