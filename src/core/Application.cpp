#include "Application.hpp"
#include "CommandBase.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "ScopeTimer.hpp"

#include <CLI/CLI.hpp>

Application::Application(CLI::App * app) :
    m_cli_app{ app }, m_potential_analysis_cmd{ nullptr },
    m_potential_display_cmd{ nullptr }, m_result_dump_cmd{ nullptr },
    m_map_simulation_cmd{ nullptr },
    m_selected_command{ "" }
{
    if (m_cli_app == nullptr)
    {
        throw std::runtime_error("Application::Application() failed: CLI::App instance is nullptr.");
    }
    RegisterCommands();
}

void Application::Run(void)
{
    ScopeTimer timer("Application::Run");
    Logger::SetLogLevel(m_global_options.verbose_level);
    std::unique_ptr<CommandBase> command{ CreateCommand() };
    if (command)
    {
        command->Execute();
    }
}

std::unique_ptr<CommandBase> Application::CreateCommand(void)
{
    if (m_cli_app->got_subcommand(m_potential_analysis_cmd))
    {
        auto command{ std::make_unique<PotentialAnalysisCommand>() };
        command->SetThreadSize(m_global_options.thread_size);
        command->SetAsymmetryFlag(m_potential_analysis_options.is_asymmetry);
        command->SetSimulationFlag(m_potential_analysis_options.is_simulation);
        command->SetSimulatedMapResolution(m_potential_analysis_options.resolution_simulation);
        command->SetSavedKeyTag(m_potential_analysis_options.saved_key_tag);
        command->SetSamplingSize(m_sphere_sampler_options.sampling_size);
        command->SetSamplingRangeMinimum(m_sphere_sampler_options.sampling_range_min);
        command->SetSamplingRangeMaximum(m_sphere_sampler_options.sampling_range_max);
        command->SetDatabasePath(m_global_options.database_path);
        command->SetModelFilePath(m_potential_analysis_options.model_file_path);
        command->SetMapFilePath(m_potential_analysis_options.map_file_path);
        command->SetFitRangeMinimum(m_potential_analysis_options.fit_range_min);
        command->SetFitRangeMaximum(m_potential_analysis_options.fit_range_max);
        command->SetAlphaR(m_potential_analysis_options.alpha_r);
        command->SetAlphaG(m_potential_analysis_options.alpha_g);
        return command;
    }
    else if (m_cli_app->got_subcommand(m_potential_display_cmd))
    {
        auto command{ std::make_unique<PotentialDisplayCommand>() };
        command->SetPainterChoice(m_potential_display_options.painter_choice);
        command->SetDatabasePath(m_global_options.database_path);
        command->SetFolderPath(m_global_options.folder_path);
        command->SetModelKeyTagList(m_potential_display_options.model_key_tag_list);
        command->SetRefModelKeyTagListMap(m_potential_display_options.ref_model_key_tag_list);
        command->SetPickChainID(m_atom_selector_options.pick_chain_id);
        command->SetPickResidueType(m_atom_selector_options.pick_residue);
        command->SetPickElementType(m_atom_selector_options.pick_element);
        command->SetPickRemotenessType(m_atom_selector_options.pick_remoteness);
        command->SetVetoChainID(m_atom_selector_options.veto_chain_id);
        command->SetVetoResidueType(m_atom_selector_options.veto_residue);
        command->SetVetoElementType(m_atom_selector_options.veto_element);
        command->SetVetoRemotenessType(m_atom_selector_options.veto_remoteness);
        return command;
    }
    else if (m_cli_app->got_subcommand(m_result_dump_cmd))
    {
        auto command{ std::make_unique<ResultDumpCommand>() };
        command->SetPrinterChoice(m_result_dump_options.printer_choice);
        command->SetModelKeyTagList(m_result_dump_options.model_key_tag_list);
        command->SetMapFilePath(m_result_dump_options.map_file_path);
        command->SetDatabasePath(m_global_options.database_path);
        command->SetFolderPath(m_global_options.folder_path);
        return command;
    }
    else if (m_cli_app->got_subcommand(m_map_simulation_cmd))
    {
        auto command{ std::make_unique<MapSimulationCommand>() };
        command->SetModelFilePath(m_map_simulation_options.model_file_path);
        command->SetMapFileName(m_map_simulation_options.map_file_name);
        command->SetFolderPath(m_global_options.folder_path);
        command->SetThreadSize(m_global_options.thread_size);
        command->SetPotentialModelChoice(m_map_simulation_options.potential_model_choice);
        command->SetPartialChargeChoice(m_map_simulation_options.partial_charge_choice);
        command->SetCutoffDistance(m_map_simulation_options.cutoff_distance);
        command->SetGridSpacing(m_map_simulation_options.grid_spacing);
        command->SetBlurringWidthList(m_map_simulation_options.blurring_width_list);
        return command;
    }
    else
    {
        std::cerr <<"The sub-command is not defined!"<< std::endl;
        return nullptr;
    }
}

void Application::RegisterCommands(void)
{
    RegisterPotentialAnalysisCommand();
    RegisterPotentialDisplayCommand();
    RegisterResultDumpCommand();
    RegisterMapSimulationCommand();
}

void Application::RegisterPotentialAnalysisCommand(void)
{
    m_potential_analysis_cmd = m_cli_app->add_subcommand("potential_analysis", "Run potential analysis");
    m_potential_analysis_cmd->add_option(
        "-a,--model", m_potential_analysis_options.model_file_path,
        "Model file path")->required();
    m_potential_analysis_cmd->add_option(
        "-m,--map", m_potential_analysis_options.map_file_path,
        "Map file path")->required();
    m_potential_analysis_cmd->add_option(
        "--simulation", m_potential_analysis_options.is_simulation,
        "Simulation flag")->default_val(false);
    m_potential_analysis_cmd->add_option(
        "-r,--sim-resolution", m_potential_analysis_options.resolution_simulation,
        "Set simulated map's resolution (blurring width)")->default_val(0.0);
    m_potential_analysis_cmd->add_option(
        "-k,--save-key", m_potential_analysis_options.saved_key_tag,
        "New key tag for saving ModelObject results into database")->default_val("");
    m_potential_analysis_cmd->add_option(
        "--asymmetry", m_potential_analysis_options.is_asymmetry,
        "Turn On/Off asymmetry flag")->default_val(false);
    m_potential_analysis_cmd->add_option(
        "-s,--sampling", m_sphere_sampler_options.sampling_size,
        "Number of sampling points per atom")->default_val(1500);
    m_potential_analysis_cmd->add_option(
        "--sampling-min", m_sphere_sampler_options.sampling_range_min,
        "Minimum sampling range")->default_val(0.0);
    m_potential_analysis_cmd->add_option(
        "--sampling-max", m_sphere_sampler_options.sampling_range_max,
        "Maximum sampling range")->default_val(1.5);
    m_potential_analysis_cmd->add_option(
        "--fit-min", m_potential_analysis_options.fit_range_min,
        "Minimum fitting range")->default_val(0.0);
    m_potential_analysis_cmd->add_option(
        "--fit-max", m_potential_analysis_options.fit_range_max,
        "Maximum fitting range")->default_val(1.0);
    m_potential_analysis_cmd->add_option(
        "--alpha-r", m_potential_analysis_options.alpha_r,
        "Alpha value for R")->default_val(0.1);
    m_potential_analysis_cmd->add_option(
        "--alpha-g", m_potential_analysis_options.alpha_g,
        "Alpha value for G")->default_val(0.2);
    RegisterGlobalOptions(m_potential_analysis_cmd);
    
    m_cli_app->callback([&]()
    {
        m_selected_command = "potential_analysis";
    });
}

void Application::RegisterPotentialDisplayCommand(void)
{
    m_potential_display_cmd = m_cli_app->add_subcommand("potential_display", "Run potential display");
    m_potential_display_cmd->add_option(
        "-p,--painter", m_potential_display_options.painter_choice,
        "Painter choice")->required();
    m_potential_display_cmd->add_option(
        "-k,--model-keylist", m_potential_display_options.model_key_tag_list,
        "List of model key tag to be display")->required();
    m_potential_display_cmd->add_option(
        "-r,--ref-model-keylist", m_potential_display_options.ref_model_key_tag_list,
        "List of reference model key tag to be display")->default_val("");
    m_potential_display_cmd->add_option(
        "--pick-chain", m_atom_selector_options.pick_chain_id,
        "Pick chain ID")->default_val("");
    m_potential_display_cmd->add_option(
        "--pick-residue", m_atom_selector_options.pick_residue,
        "Pick residue type")->default_val("");
    m_potential_display_cmd->add_option(
        "--pick-element", m_atom_selector_options.pick_element,
        "Pick element type")->default_val("");
    m_potential_display_cmd->add_option(
        "--pick-remoteness", m_atom_selector_options.pick_remoteness,
        "Pick remoteness type")->default_val("");
    m_potential_display_cmd->add_option(
        "--veto-chain", m_atom_selector_options.veto_chain_id,
        "Veto chain ID")->default_val("");
    m_potential_display_cmd->add_option(
        "--veto-residue", m_atom_selector_options.veto_residue,
        "Veto residue type")->default_val("");
    m_potential_display_cmd->add_option(
        "--veto-element", m_atom_selector_options.veto_element,
        "Veto element type")->default_val("");
    m_potential_display_cmd->add_option(
        "--veto-remoteness", m_atom_selector_options.veto_remoteness,
        "Veto remoteness type")->default_val("");
    RegisterGlobalOptions(m_potential_display_cmd);

    m_cli_app->callback([&]()
    {
        m_selected_command = "potential_display";
    });
}

void Application::RegisterResultDumpCommand(void)
{
    m_result_dump_cmd = m_cli_app->add_subcommand("result_dump", "Run result dump");
    m_result_dump_cmd->add_option(
        "-p,--printer", m_result_dump_options.printer_choice,
        "Printer choice")->required();
    m_result_dump_cmd->add_option(
        "-k,--model-keylist", m_result_dump_options.model_key_tag_list,
        "List of model key tag to be display")->required();
    m_result_dump_cmd->add_option(
        "-m,--map", m_result_dump_options.map_file_path,
        "Map file path")->default_val("");
    RegisterGlobalOptions(m_result_dump_cmd);

    m_cli_app->callback([&]()
    {
        m_selected_command = "result_dump";
    });
}

void Application::RegisterMapSimulationCommand(void)
{
    m_map_simulation_cmd = m_cli_app->add_subcommand("map_simulation", "Run map simulation command");
    m_map_simulation_cmd->add_option(
        "-a,--model", m_map_simulation_options.model_file_path,
        "Model file path")->required();
    m_map_simulation_cmd->add_option(
        "-n,--name", m_map_simulation_options.map_file_name,
        "File name for output map files")->default_val("sim_map");
    m_map_simulation_cmd->add_option(
        "--potential-model", m_map_simulation_options.potential_model_choice,
        "Atomic potential model option")->default_val(1);
    m_map_simulation_cmd->add_option(
        "--charge", m_map_simulation_options.partial_charge_choice,
        "Partial charge table option")->default_val(1);
    m_map_simulation_cmd->add_option(
        "-c,--cut-off", m_map_simulation_options.cutoff_distance,
        "Cutoff distance")->default_val(5.0);
    m_map_simulation_cmd->add_option(
        "-g,--grid-spacing", m_map_simulation_options.grid_spacing,
        "Grid spacing")->default_val(0.5);
    m_map_simulation_cmd->add_option(
        "--blurring-width", m_map_simulation_options.blurring_width_list,
        "Blurring width (list) setting")->default_val("0.50");
    RegisterGlobalOptions(m_map_simulation_cmd);

    m_map_simulation_cmd->callback([&]()
    {
        m_selected_command = "map_simulation";
    });
}

void Application::RegisterGlobalOptions(CLI::App * command)
{
    command->add_option(
        "-d,--database", m_global_options.database_path,
        "Database file path")->default_val("database.sqlite");
    command->add_option(
        "-o,--folder", m_global_options.folder_path,
        "folder path for output files")->default_val("");
    command->add_option(
        "-j,--jobs", m_global_options.thread_size,
        "Number of threads")->default_val(1);
    command->add_option(
        "-v,--verbose", m_global_options.verbose_level,
        "Verbose level")->default_val(2);
}
