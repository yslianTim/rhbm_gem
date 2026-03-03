#include "Application.hpp"
#include "CommandBase.hpp"
#include "Logger.hpp"
#include "RegisterBuiltInCommands.hpp"
#include "ScopeTimer.hpp"

#include <memory>
#include <CLI/CLI.hpp>

namespace rhbm_gem {

Application::Application(CLI::App & app) :
    m_cli_app{ app }
{
    m_cli_app.require_subcommand(1);
    RegisterBuiltInCommands();
    RegisterAllCommands();
}

void Application::RegisterAllCommands()
{
    const auto & commands{ CommandRegistry::Instance().GetCommands() };
    for (const auto & info : commands)
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
        cmd->Execute();
    });
}

} // namespace rhbm_gem
