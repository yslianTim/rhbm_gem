#pragma once

#include <memory>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include "DataObjectVisitorBase.hpp"
#include "KDTreeAlgorithm.hpp"
#include "MapSimulationCommand.hpp"

class AtomSelector;

class MapSimulationVisitor : public DataObjectVisitorBase
{
    unsigned int m_thread_size;
    int m_potential_model_choice;
    int m_partial_charge_choice;
    double m_cutoff_distance;
    float m_grid_spacing;
    std::string m_model_key_tag, m_folder_path, m_pdb_id, m_map_file_name;
    std::vector<double> m_blurring_width_list;
    std::vector<AtomObject *> m_selected_atom_list;
    std::unordered_map<int, double> m_atom_charge_map;
    std::unique_ptr<KDNode<AtomObject>> m_kd_tree_root;

public:
    explicit MapSimulationVisitor(const MapSimulationCommand::Options & options);
    ~MapSimulationVisitor();
    void VisitAtomObject(AtomObject * data_object) override;
    void VisitModelObject(ModelObject * data_object) override;
    void VisitDataObjectManager(DataObjectManager * data_manager) override;

private:
    std::unique_ptr<MapObject> CreateSimulatedMapObject(double blurring_width);
    std::array<int, 3> CalculateGridSize(
        const std::array<float, 3> & grid_spacing, const std::array<float, 3> & origin);

};
