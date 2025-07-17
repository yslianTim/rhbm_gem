#include "Application.hpp"
#include "CommandBase.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "ScopeTimer.hpp"
#include "Logger.hpp"

#include <CLI/CLI.hpp>

Application::Application(CLI::App & app) :
    m_cli_app{ app }
{
    m_cli_app.require_subcommand(1);
    RegisterAllCommands();
}

void Application::RegisterAllCommands(void)
{
    RegisterCommand<PotentialAnalysisCommand>("potential_analysis",
        "Run potential analysis");
    RegisterCommand<PotentialDisplayCommand>("potential_display",
        "Run potential display");
    RegisterCommand<ResultDumpCommand>("result_dump",
        "Run result dump");
    RegisterCommand<MapSimulationCommand>("map_simulation",
        "Run map simulation command");
}

template<typename Type>
void Application::RegisterCommand(const std::string & name, const std::string & description)
{
    auto options{ std::make_shared<typename Type::Options>() };
    CLI::App * command{ m_cli_app.add_subcommand(name, description) };
    Type::RegisterCLIOptions(command, *options);
    RegisterGlobalOptions(command);

    command->callback([this, options]() {
        ScopeTimer timer("Command in Application");
        Logger::SetLogLevel(m_global_options.verbose_level);
        Type command_object(*options, m_global_options);
        command_object.Execute();
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
