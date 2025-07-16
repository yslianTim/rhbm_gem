#pragma once

#include <memory>
#include <string>
#include <vector>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"
#include "GlobalOptions.hpp"

class MapSimulationCommand : public CommandBase
{
    int m_thread_size;
    int m_potential_model_choice;
    int m_partial_charge_choice;
    double m_cutoff_distance;
    double m_grid_spacing;
    std::string m_model_file_path, m_folder_path, m_map_file_name;
    std::vector<double> m_blurring_width_list;
    
public:
    struct Options
    {
        std::string model_file_path;
        std::string map_file_name{"sim_map"};
        int potential_model_choice{ 1 };
        int partial_charge_choice{ 1 };
        double cutoff_distance{ 5.0 };
        double grid_spacing{ 0.5 };
        std::string blurring_width_list{"0.50"};
    };
    MapSimulationCommand(void);
    MapSimulationCommand(const Options & options, const GlobalOptions & globals);
    ~MapSimulationCommand() = default;
    void Execute(void) override;

    static void RegisterCLIOptions(CLI::App * cmd, Options & options);
    void SetThreadSize(int value) { m_thread_size = value; }
    void SetPotentialModelChoice(int value) { m_potential_model_choice = value; }
    void SetPartialChargeChoice(int value) { m_partial_charge_choice = value; }
    void SetCutoffDistance(double value) { m_cutoff_distance = value; }
    void SetModelFilePath(const std::string & value) { m_model_file_path = value; }
    void SetFolderPath(const std::string & path) { m_folder_path = path; }
    void SetMapFileName(const std::string & value) { m_map_file_name = value; }
    void SetGridSpacing(double value) { m_grid_spacing = value; }
    void SetBlurringWidthList(const std::string & value);

};
