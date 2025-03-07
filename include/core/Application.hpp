#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <CLI/CLI.hpp>

class CommandBase;

class Application
{
    CLI::App * m_cli_app;
    CLI::App * m_potential_analysis_cmd;
    CLI::App * m_test_cmd;
    std::string m_selected_command;

    struct PotentialAnalysisOptions
    {
        int thread_size, sampling_size;
        double sampling_range_min, sampling_range_max;
        double fit_range_min, fit_range_max;
        double alpha_r, alpha_g;
        std::string database_path;
        std::string model_file_path;
        std::string map_file_path;
        std::string pick_chain_id;
        std::string pick_indicator;
        std::string pick_residue;
        std::string pick_element;
        std::string pick_remoteness;
        std::string pick_branch;
        std::string veto_chain_id;
        std::string veto_indicator;
        std::string veto_residue;
        std::string veto_element;
        std::string veto_remoteness;
        std::string veto_branch;
    } m_potential_analysis_options;

    struct TestOptions
    {
        int thread_size;
    } m_test_options;

public:
    Application(CLI::App * app);
    ~Application() = default;
    void Run(void);

private:
    std::unique_ptr<CommandBase> CreateCommand(void);
    void RegisterCommands(void);
    void RegisterPotentialAnalysisCommand(void);
    void RegisterTestCommand(void);

};