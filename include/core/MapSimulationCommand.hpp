#pragma once

#include <memory>
#include <string>
#include <vector>
#include "CommandBase.hpp"

class AtomSelector;

class MapSimulationCommand : public CommandBase
{
    int m_thread_size;
    int m_potential_model_choice;
    int m_partial_charge_choice;
    double m_cutoff_distance;
    double m_grid_spacing;
    std::string m_model_file_path, m_folder_path, m_map_file_name;
    std::vector<double> m_blurring_width_list;
    std::unique_ptr<AtomSelector> m_atom_selector;
    
public:
    MapSimulationCommand(void);
    ~MapSimulationCommand();
    void Execute(void) override;

    void SetThreadSize(int value) { m_thread_size = value; }
    void SetPotentialModelChoice(int value) { m_potential_model_choice = value; }
    void SetPartialChargeChoice(int value) { m_partial_charge_choice = value; }
    void SetCutoffDistance(double value) { m_cutoff_distance = value; }
    void SetModelFilePath(const std::string & value) { m_model_file_path = value; }
    void SetFolderPath(const std::string & path) { m_folder_path = path; }
    void SetMapFileName(const std::string & value) { m_map_file_name = value; }
    void SetGridSpacing(double value) { m_grid_spacing = value; }
    void SetBlurringWidthList(const std::string & value);
    void SetPickChainID(const std::string & value);
    void SetPickResidueType(const std::string & value);
    void SetPickElementType(const std::string & value);
    void SetPickRemotenessType(const std::string & value);
    void SetVetoChainID(const std::string & value);
    void SetVetoResidueType(const std::string & value);
    void SetVetoElementType(const std::string & value);
    void SetVetoRemotenessType(const std::string & value);

};