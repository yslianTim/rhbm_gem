#pragma once

#include <memory>
#include <string>
#include <vector>
#include <CLI/CLI.hpp>

#include "CommandBase.hpp"

class MapSimulationCommand : public CommandBase
{
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

private:
    Options m_options{};
    std::vector<double> m_blurring_width_list;

public:
    MapSimulationCommand(void);
    ~MapSimulationCommand() = default;
    void Execute(void) override;
    void RegisterCLIOptions(CLI::App * cmd) override;

    void SetThreadSize(int value) { m_globals.thread_size = value; }
    void SetPotentialModelChoice(int value) { m_options.potential_model_choice = value; }
    void SetPartialChargeChoice(int value) { m_options.partial_charge_choice = value; }
    void SetCutoffDistance(double value) { m_options.cutoff_distance = value; }
    void SetModelFilePath(const std::string & value) { m_options.model_file_path = value; }
    void SetFolderPath(const std::string & path) { m_globals.folder_path = path; }
    void SetMapFileName(const std::string & value) { m_options.map_file_name = value; }
    void SetGridSpacing(double value) { m_options.grid_spacing = value; }
    void SetBlurringWidthList(const std::string & value);

};
