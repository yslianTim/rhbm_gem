#include "Application.hpp"
#include "CommandBase.hpp"
#include "ScopeTimer.hpp"
#include "Logger.hpp"

#include <memory>
#include <CLI/CLI.hpp>

Application::Application(CLI::App & app) :
    m_cli_app{ app }
{
    Logger::Log(LogLevel::Debug, "Application::Application() called");
    m_cli_app.require_subcommand(1);
    RegisterAllCommands();
}

void Application::RegisterAllCommands(void)
{
    Logger::Log(LogLevel::Debug, "Application::RegisterAllCommands() called");
    const auto & commands{ CommandRegistry::Instance().GetCommands() };
    for (const auto & [name, info] : commands)
    {
        RegisterCommand(info.name, info.description, info.factory);
    }
}

void Application::RegisterCommand(
    const std::string & name,
    const std::string & description,
    CommandRegistry::FactoryFunc factory)
{
    Logger::Log(LogLevel::Debug, "Application::RegisterCommand() called");
    auto command_object{ factory() };
    CLI::App * command{ m_cli_app.add_subcommand(name, description) };
    command_object->RegisterCLIOptions(command);

    auto shared_cmd{ std::shared_ptr<CommandBase>(std::move(command_object)) };
    command->callback([cmd = std::move(shared_cmd)]() {
        ScopeTimer timer("Command in Application");
        Logger::SetLogLevel(cmd->GetOptions().verbose_level);
        if (!cmd->Execute())
        {
            Logger::Log(LogLevel::Error, "Command execution failed");
        }
    });
}
