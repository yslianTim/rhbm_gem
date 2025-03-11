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
    CLI::App * m_potential_display_cmd;
    CLI::App * m_test_cmd;
    std::string m_selected_command;

    struct AtomSelectorOptions
    {
        std::string pick_chain_id, veto_chain_id;
        std::string pick_indicator, veto_indicator;
        std::string pick_residue, veto_residue;
        std::string pick_element, veto_element;
        std::string pick_remoteness, veto_remoteness;
        std::string pick_branch, veto_branch;
    } m_atom_selector_options;

    struct SphereSamplerOptions
    {
        int sampling_size;
        double sampling_range_min, sampling_range_max;
    } m_sphere_sampler_options;

    struct PotentialAnalysisOptions
    {
        double fit_range_min, fit_range_max;
        double alpha_r, alpha_g;
        std::string model_file_path;
        std::string map_file_path;
        std::string saved_key_tag;
    } m_potential_analysis_options;

    struct PotentialDisplayOptions
    {
        std::string model_key_tag;
        std::string folder_path;
    } m_potential_display_options;

    struct GlobalOptions
    {
        int thread_size, verbose_level;
        std::string database_path;
    } m_global_options;

public:
    Application(CLI::App * app);
    ~Application() = default;
    void Run(void);

private:
    std::unique_ptr<CommandBase> CreateCommand(void);
    void RegisterCommands(void);
    void RegisterPotentialAnalysisCommand(void);
    void RegisterPotentialDisplayCommand(void);
    void RegisterTestCommand(void);

};