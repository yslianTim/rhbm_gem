#include "Application.hpp"
#include "CommandBase.hpp"
#include "ScopeTimer.hpp"
#include "Logger.hpp"

#include <CLI/CLI.hpp>

Application::Application(CLI::App & app) :
    m_cli_app{ app }
{
    m_cli_app.require_subcommand(1);
    RegisterAllCommands();
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
    auto selected_command{ m_cli_app.get_subcommands() };
    if (selected_command.size() != 1)
    {
        Logger::Log(LogLevel::Error, "No valid subcommand provided.");
        return nullptr;
    }

    auto command_name{ selected_command[0]->get_name() };
    auto it{ m_command_map.find(command_name) };
    if (it != m_command_map.end())
    {
        Logger::Log(LogLevel::Info, "Subcommand found: " + command_name);
        return it->second();
    }

    Logger::Log(LogLevel::Error, "No valid subcommand provided.");
    return nullptr;
}

void Application::RegisterAllCommands(void)
{
    RegisterCommand<PotentialAnalysisCommand>("potential_analysis",
        "Run potential analysis", m_potential_analysis_options);
    RegisterCommand<PotentialDisplayCommand>("potential_display",
        "Run potential display", m_potential_display_options);
    RegisterCommand<ResultDumpCommand>("result_dump",
        "Run result dump", m_result_dump_options);
    RegisterCommand<MapSimulationCommand>("map_simulation",
        "Run map simulation command", m_map_simulation_options);
}

template<typename Type>
void Application::RegisterCommand(
    const std::string & name,
    const std::string & description,
    typename Type::Options & options)
{
    CLI::App * command{ m_cli_app.add_subcommand(name, description) };
    Type::RegisterCLIOptions(command, options);
    RegisterGlobalOptions(command);
    m_command_map.emplace(name, [this, &options]() {
        return std::make_unique<Type>(options, m_global_options);
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
