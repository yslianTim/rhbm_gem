#include "Application.hpp"
#include "CommandBase.hpp"
#include "ScopeTimer.hpp"
#include "Logger.hpp"

#include <memory>
#include <CLI/CLI.hpp>

namespace rhbm_gem {

Application::Application(CLI::App & app) :
    m_cli_app{ app }
{
    m_cli_app.require_subcommand(1);
    RegisterAllCommands();
}

void Application::RegisterAllCommands()
{
    const auto & commands{ CommandRegistry::Instance().GetCommands() };
    for (const auto & [name, info] : commands)
    {
        RegisterCommand(info.name, info.description, info.factory);
    }
}

void Application::RegisterCommand(
    const std::string & name,
    const std::string & description,
    std::function<std::unique_ptr<CommandBase>()> factory)
{
    auto command_object{ factory() };
    CLI::App * command{ m_cli_app.add_subcommand(name, description) };
    command_object->RegisterCLIOptions(command);

    auto shared_cmd{ std::shared_ptr<CommandBase>(std::move(command_object)) };
    command->callback([cmd = std::move(shared_cmd)]() {
        ScopeTimer timer("Command in Application");
        const auto & options{ cmd->GetOptions() };
        Logger::SetLogLevel(options.verbose_level);
        if (cmd->IsValidateOptions() == false)
        {
            Logger::Log(LogLevel::Error,
                "Invalid command options detected. Aborting command execution.");
            return;
        }
        if (cmd->Execute() == false)
        {
            Logger::Log(LogLevel::Error,
                "Command execution failed. Aborting command execution.");
        }
    });
}

} // namespace rhbm_gem
