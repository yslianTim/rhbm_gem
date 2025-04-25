#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include "DataObjectVisitorBase.hpp"
#include "AtomicInfoHelper.hpp"

class AtomSelector;

class MapSimulationVisitor : public DataObjectVisitorBase
{
    unsigned int m_thread_size;
    int m_potential_model_choice;
    int m_partial_charge_choice;
    double m_cutoff_distance;
    float m_grid_spacing;
    std::string m_model_key_tag, m_folder_path, m_pdb_id;
    std::vector<double> m_blurring_width_list;
    AtomSelector * m_atom_selector;
    std::vector<AtomObject *> m_selected_atom_list;

public:
    MapSimulationVisitor(AtomSelector * atom_selector);
    ~MapSimulationVisitor();
    void VisitAtomObject(AtomObject * data_object) override;
    void VisitModelObject(ModelObject * data_object) override;
    void VisitMapObject(MapObject * data_object) override;
    void Analysis(DataObjectManager * data_manager) override;

    void SetFolderPath(const std::string & path) { m_folder_path = path; }
    void SetThreadSize(unsigned int thread_size) { m_thread_size = thread_size; }
    void SetPotentialModelChoice(int value) { m_potential_model_choice = value; }
    void SetPartialChargeChoice(int value) { m_partial_charge_choice = value; }
    void SetCutoffDistance(double value) { m_cutoff_distance = value; }
    void SetGridSpacing(double value) { m_grid_spacing = static_cast<float>(value); }
    void SetBlurringWidthList(const std::vector<double> & value) { m_blurring_width_list = value; }
    void SetModelObjectKeyTag(const std::string & value) { m_model_key_tag = value; }

private:
    std::unique_ptr<MapObject> CreateSimulatedMapObject(double blurring_width);
    std::array<int, 3> CalculateGridSize(
        const std::array<float, 3> & grid_spacing, const std::array<float, 3> & origin);

};