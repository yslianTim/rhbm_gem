#include "Application.hpp"
#include "CommandBase.hpp"
#include "PotentialAnalysisCommand.hpp"
#include "PotentialDisplayCommand.hpp"
#include "ResultDumpCommand.hpp"
#include "MapSimulationCommand.hpp"
#include "ScopeTimer.hpp"
#include "Logger.hpp"

#include <memory>
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
    auto command_object{ std::make_unique<Type>() };
    CLI::App * command{ m_cli_app.add_subcommand(name, description) };
    command_object->RegisterCLIOptions(command);

    auto shared_cmd{ std::shared_ptr<Type>(std::move(command_object)) };
    command->callback([cmd = std::move(shared_cmd)]() {
        ScopeTimer timer("Command in Application");
        Logger::SetLogLevel(cmd->GetOptions().verbose_level);
        if (!cmd->Execute())
        {
            Logger::Log(LogLevel::Error, "Command execution failed");
        }
    });
}
