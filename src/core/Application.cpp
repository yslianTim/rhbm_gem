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
        command->SetThreadSize(m_potential_analysis_options.thread_size);
        command->SetDatabasePath(m_potential_analysis_options.database_path);
        command->SetModelFilePath(m_potential_analysis_options.model_file_path);
        command->SetMapFilePath(m_potential_analysis_options.map_file_path);
        command->SetPickChainID(m_potential_analysis_options.pick_chain_id);
        command->SetPickIndicator(m_potential_analysis_options.pick_indicator);
        command->SetPickResidueType(m_potential_analysis_options.pick_residue);
        command->SetPickElementType(m_potential_analysis_options.pick_element);
        command->SetPickRemotenessType(m_potential_analysis_options.pick_remoteness);
        command->SetPickBranchType(m_potential_analysis_options.pick_branch);
        command->SetVetoChainID(m_potential_analysis_options.veto_chain_id);
        command->SetVetoIndicator(m_potential_analysis_options.veto_indicator);
        command->SetVetoResidueType(m_potential_analysis_options.veto_residue);
        command->SetVetoElementType(m_potential_analysis_options.veto_element);
        command->SetVetoRemotenessType(m_potential_analysis_options.veto_remoteness);
        command->SetVetoBranchType(m_potential_analysis_options.veto_branch);


        return command;
    }
    else if (m_cli_app->got_subcommand(m_test_cmd))
    {
        auto command{ std::make_unique<TestCommand>() };
        command->SetThreadSize(m_test_options.thread_size);
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
        "-d,--database", m_potential_analysis_options.database_path,
        "Database file path");
    m_potential_analysis_cmd->add_option(
        "-a,--model", m_potential_analysis_options.model_file_path,
        "Model file path")->required();
    m_potential_analysis_cmd->add_option(
        "-m,--map", m_potential_analysis_options.map_file_path,
        "Map file path")->required();
    m_potential_analysis_cmd->add_option(
        "-j,--jobs", m_potential_analysis_options.thread_size,
        "Number of threads");
    m_potential_analysis_cmd->add_option(
        "-s,--sampling", m_potential_analysis_options.sampling_size,
        "Number of sampling points per atom");
    m_potential_analysis_cmd->add_option(
        "--sampling-min", m_potential_analysis_options.sampling_range_min,
        "Minimum sampling range");
    m_potential_analysis_cmd->add_option(
        "--sampling-max", m_potential_analysis_options.sampling_range_max,
        "Maximum sampling range");
    m_potential_analysis_cmd->add_option(
        "--fit-min", m_potential_analysis_options.fit_range_min,
        "Minimum fitting range");
    m_potential_analysis_cmd->add_option(
        "--fit-max", m_potential_analysis_options.fit_range_max,
        "Maximum fitting range");
    m_potential_analysis_cmd->add_option(
        "--alpha-r", m_potential_analysis_options.alpha_r,
        "Alpha value for R");
    m_potential_analysis_cmd->add_option(
        "--alpha-g", m_potential_analysis_options.alpha_g,
        "Alpha value for G");
    m_potential_analysis_cmd->add_option(
        "--pick-chain", m_potential_analysis_options.pick_chain_id,
        "Pick chain ID");
    m_potential_analysis_cmd->add_option(
        "--pick-indicator", m_potential_analysis_options.pick_indicator,
        "Pick indicator");
    m_potential_analysis_cmd->add_option(
        "--pick-residue", m_potential_analysis_options.pick_residue,
        "Pick residue type");
    m_potential_analysis_cmd->add_option(
        "--pick-element", m_potential_analysis_options.pick_element,
        "Pick element type");
    m_potential_analysis_cmd->add_option(
        "--pick-remoteness", m_potential_analysis_options.pick_remoteness,
        "Pick remoteness type");
    m_potential_analysis_cmd->add_option(
        "--pick-branch", m_potential_analysis_options.pick_branch,
        "Pick branch type");
    m_potential_analysis_cmd->add_option(
        "--veto-chain", m_potential_analysis_options.veto_chain_id,
        "Veto chain ID");
    m_potential_analysis_cmd->add_option(
        "--veto-indicator", m_potential_analysis_options.veto_indicator,
        "Veto indicator");
    m_potential_analysis_cmd->add_option(
        "--veto-residue", m_potential_analysis_options.veto_residue,
        "Veto residue type");
    m_potential_analysis_cmd->add_option(
        "--veto-element", m_potential_analysis_options.veto_element,
        "Veto element type");
    m_potential_analysis_cmd->add_option(
        "--veto-remoteness", m_potential_analysis_options.veto_remoteness,
        "Veto remoteness type");
    m_potential_analysis_cmd->add_option(
        "--veto-branch", m_potential_analysis_options.veto_branch,
        "Veto branch type");
    
    m_cli_app->callback([&]()
    {
        m_selected_command = "potential_analysis";
    });
}

void Application::RegisterTestCommand(void)
{
    m_test_cmd = m_cli_app->add_subcommand("test", "Run test command");
    m_test_cmd->add_option(
        "-j,--jobs", m_test_options.thread_size,
        "Number of threads");
    m_cli_app->callback([&]()
    {
        m_selected_command = "test";
    });
}