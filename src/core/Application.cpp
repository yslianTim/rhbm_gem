#include "Application.hpp"
#include "CommandBase.hpp"
#include "ScopeTimer.hpp"
#include "Logger.hpp"

#include <CLI/CLI.hpp>

Application::Application(CLI::App & app) :
    m_cli_app{ app }
{
    m_cli_app.require_subcommand(1);
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
    for (const auto & [command, factory] : m_command_map)
    {
        if (m_cli_app.got_subcommand(command))
        {
            Logger::Log(LogLevel::Info, "Subcommand found: " + command->get_name());
            return factory();
        }
    }
    // If no subcommand was found, log an error and return nullptr
    Logger::Log(LogLevel::Error, "No valid subcommand provided.");
    return nullptr;
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
    CLI::App * command{ m_cli_app.add_subcommand("potential_analysis", "Run potential analysis") };
    PotentialAnalysisCommand::RegisterCLIOptions(command, m_potential_analysis_options);
    RegisterGlobalOptions(command);
    m_command_map.emplace(command, [this]() {
        return std::make_unique<PotentialAnalysisCommand>(m_potential_analysis_options, m_global_options);
    });
}

void Application::RegisterPotentialDisplayCommand(void)
{
    CLI::App * command{ m_cli_app.add_subcommand("potential_display", "Run potential display") };
    PotentialDisplayCommand::RegisterCLIOptions(command, m_potential_display_options);
    RegisterGlobalOptions(command);
    m_command_map.emplace(command, [this]() {
        return std::make_unique<PotentialDisplayCommand>(m_potential_display_options, m_global_options);
    });
}

void Application::RegisterResultDumpCommand(void)
{
    CLI::App * command{ m_cli_app.add_subcommand("result_dump", "Run result dump") };
    ResultDumpCommand::RegisterCLIOptions(command, m_result_dump_options);
    RegisterGlobalOptions(command);
    m_command_map.emplace(command, [this]() {
        return std::make_unique<ResultDumpCommand>(m_result_dump_options, m_global_options);
    });
}

void Application::RegisterMapSimulationCommand(void)
{
    CLI::App * command{ m_cli_app.add_subcommand("map_simulation", "Run map simulation command") };
    MapSimulationCommand::RegisterCLIOptions(command, m_map_simulation_options);
    RegisterGlobalOptions(command);
    m_command_map.emplace(command, [this]() {
        return std::make_unique<MapSimulationCommand>(m_map_simulation_options, m_global_options);
    });
}

void Application::RegisterGlobalOptions(CLI::App * command)
{
    command->add_option("-d,--database", m_global_options.database_path,
        "Database file path")->default_val("database.sqlite");
    command->add_option("-o,--folder", m_global_options.folder_path,
        "folder path for output files")->default_val("");
    command->add_option("-j,--jobs", m_global_options.thread_size,
        "Number of threads")->default_val(1);
    command->add_option("-v,--verbose", m_global_options.verbose_level,
        "Verbose level")->default_val(3);
}
