#include "Application.hpp"
#include "CommandBase.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "TestCommand.hpp"
#include "ScopeTimer.hpp"

Application::Application(CLI::App * app) :
    m_cli_app{ app }, m_potential_analysis_cmd{ nullptr }, m_test_cmd{ nullptr },
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
        command->SetPickChainID(m_atom_selector_options.pick_chain_id);
        command->SetPickIndicator(m_atom_selector_options.pick_indicator);
        command->SetPickResidueType(m_atom_selector_options.pick_residue);
        command->SetPickElementType(m_atom_selector_options.pick_element);
        command->SetPickRemotenessType(m_atom_selector_options.pick_remoteness);
        command->SetPickBranchType(m_atom_selector_options.pick_branch);
        command->SetVetoChainID(m_atom_selector_options.veto_chain_id);
        command->SetVetoIndicator(m_atom_selector_options.veto_indicator);
        command->SetVetoResidueType(m_atom_selector_options.veto_residue);
        command->SetVetoElementType(m_atom_selector_options.veto_element);
        command->SetVetoRemotenessType(m_atom_selector_options.veto_remoteness);
        command->SetVetoBranchType(m_atom_selector_options.veto_branch);

        return command;
    }
    else if (m_cli_app->got_subcommand(m_test_cmd))
    {
        auto command{ std::make_unique<TestCommand>() };
        command->SetThreadSize(m_global_options.thread_size);
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
    RegisterTestCommand();
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
        "-d,--database", m_global_options.database_path,
        "Database file path")->default_str("database.sqlite");
    m_potential_analysis_cmd->add_option(
        "-k,--save-key", m_potential_analysis_options.saved_key_tag,
        "Number of sampling points per atom")->default_str("");
    m_potential_analysis_cmd->add_option(
        "-j,--jobs", m_global_options.thread_size,
        "Number of threads")->default_val(1);
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
    m_potential_analysis_cmd->add_option(
        "--pick-chain", m_atom_selector_options.pick_chain_id,
        "Pick chain ID")->default_str("");
    m_potential_analysis_cmd->add_option(
        "--pick-indicator", m_atom_selector_options.pick_indicator,
        "Pick indicator")->default_str(".,A");
    m_potential_analysis_cmd->add_option(
        "--pick-residue", m_atom_selector_options.pick_residue,
        "Pick residue type")->default_str("");
    m_potential_analysis_cmd->add_option(
        "--pick-element", m_atom_selector_options.pick_element,
        "Pick element type")->default_str("");
    m_potential_analysis_cmd->add_option(
        "--pick-remoteness", m_atom_selector_options.pick_remoteness,
        "Pick remoteness type")->default_str("");
    m_potential_analysis_cmd->add_option(
        "--pick-branch", m_atom_selector_options.pick_branch,
        "Pick branch type")->default_str("");
    m_potential_analysis_cmd->add_option(
        "--veto-chain", m_atom_selector_options.veto_chain_id,
        "Veto chain ID")->default_str("");
    m_potential_analysis_cmd->add_option(
        "--veto-indicator", m_atom_selector_options.veto_indicator,
        "Veto indicator")->default_str("");
    m_potential_analysis_cmd->add_option(
        "--veto-residue", m_atom_selector_options.veto_residue,
        "Veto residue type")->default_str("");
    m_potential_analysis_cmd->add_option(
        "--veto-element", m_atom_selector_options.veto_element,
        "Veto element type")->default_str("");
    m_potential_analysis_cmd->add_option(
        "--veto-remoteness", m_atom_selector_options.veto_remoteness,
        "Veto remoteness type")->default_str("");
    m_potential_analysis_cmd->add_option(
        "--veto-branch", m_atom_selector_options.veto_branch,
        "Veto branch type")->default_str("");
    
    m_cli_app->callback([&]()
    {
        m_selected_command = "potential_analysis";
    });
}

void Application::RegisterTestCommand(void)
{
    m_test_cmd = m_cli_app->add_subcommand("test", "Run test command");
    m_test_cmd->add_option(
        "-j,--jobs", m_global_options.thread_size,
        "Number of threads");
    m_cli_app->callback([&]()
    {
        m_selected_command = "test";
    });
}