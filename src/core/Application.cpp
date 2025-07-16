#include "Application.hpp"
#include "CommandBase.hpp"
#include "ScopeTimer.hpp"

#include <CLI/CLI.hpp>

Application::Application(CLI::App & app) :
    m_cli_app{ app }, m_potential_analysis_cmd{ nullptr },
    m_potential_display_cmd{ nullptr }, m_result_dump_cmd{ nullptr },
    m_map_simulation_cmd{ nullptr }
{
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
    if (m_cli_app.got_subcommand(m_potential_analysis_cmd))
    {
        auto command
        {
            std::make_unique<PotentialAnalysisCommand>(m_potential_analysis_options, m_global_options)
        };
        return command;
    }
    else if (m_cli_app.got_subcommand(m_potential_display_cmd))
    {
        auto command
        {
            std::make_unique<PotentialDisplayCommand>(m_potential_display_options, m_global_options)
        };
        return command;
    }
    else if (m_cli_app.got_subcommand(m_result_dump_cmd))
    {
        auto command
        {
            std::make_unique<ResultDumpCommand>(m_result_dump_options, m_global_options)
        };
        return command;
    }
    else if (m_cli_app.got_subcommand(m_map_simulation_cmd))
    {
        auto command
        {
            std::make_unique<MapSimulationCommand>(m_map_simulation_options, m_global_options)
        };
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
    m_potential_analysis_cmd = m_cli_app.add_subcommand("potential_analysis", "Run potential analysis");
    PotentialAnalysisCommand::RegisterCLIOptions(m_potential_analysis_cmd, m_potential_analysis_options);
    RegisterGlobalOptions(m_potential_analysis_cmd);
}

void Application::RegisterPotentialDisplayCommand(void)
{
    m_potential_display_cmd = m_cli_app.add_subcommand("potential_display", "Run potential display");
    PotentialDisplayCommand::RegisterCLIOptions(m_potential_display_cmd, m_potential_display_options);
    RegisterGlobalOptions(m_potential_display_cmd);
}

void Application::RegisterResultDumpCommand(void)
{
    m_result_dump_cmd = m_cli_app.add_subcommand("result_dump", "Run result dump");
    ResultDumpCommand::RegisterCLIOptions(m_result_dump_cmd, m_result_dump_options);
    RegisterGlobalOptions(m_result_dump_cmd);
}

void Application::RegisterMapSimulationCommand(void)
{
    m_map_simulation_cmd = m_cli_app.add_subcommand("map_simulation", "Run map simulation command");
    MapSimulationCommand::RegisterCLIOptions(m_map_simulation_cmd, m_map_simulation_options);
    RegisterGlobalOptions(m_map_simulation_cmd);
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
        "Verbose level")->default_val(3);
}
