#pragma once

#include <memory>
#include <string>

namespace CLI
{
    class App;
}

class CommandBase;

class Application
{
    CLI::App * m_cli_app;
    CLI::App * m_potential_analysis_cmd;
    CLI::App * m_potential_display_cmd;
    CLI::App * m_result_dump_cmd;
    CLI::App * m_map_simulation_cmd;
    std::string m_selected_command;

    struct AtomSelectorOptions
    {
        std::string pick_chain_id, veto_chain_id;
        std::string pick_residue, veto_residue;
        std::string pick_element, veto_element;
        std::string pick_remoteness, veto_remoteness;
    } m_atom_selector_options;

    struct SphereSamplerOptions
    {
        int sampling_size;
        double sampling_range_min, sampling_range_max;
    } m_sphere_sampler_options;

    struct PotentialAnalysisOptions
    {
        bool is_asymmetry, is_simulation;
        double fit_range_min, fit_range_max;
        double alpha_r, alpha_g;
        double resolution_simulation;
        std::string model_file_path;
        std::string map_file_path;
        std::string saved_key_tag;
    } m_potential_analysis_options;

    struct PotentialDisplayOptions
    {
        int painter_choice;
        std::string model_key_tag_list;
        std::string ref_model_key_tag_list;
    } m_potential_display_options;

    struct ResultDumpOptions
    {
        int printer_choice;
        std::string model_key_tag_list;
        std::string map_file_path;
    } m_result_dump_options;

    struct MapSimulationOptions
    {
        std::string model_file_path, map_file_name;
        int potential_model_choice;
        int partial_charge_choice;
        double cutoff_distance;
        double grid_spacing;
        std::string blurring_width_list;
    } m_map_simulation_options;

    struct GlobalOptions
    {
        int thread_size, verbose_level;
        std::string database_path, folder_path;
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
    void RegisterResultDumpCommand(void);
    void RegisterMapSimulationCommand(void);
    void RegisterGlobalOptions(CLI::App * command);

};
