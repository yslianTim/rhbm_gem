#include "Application.hpp"
#include "BuiltInCommandCatalogInternal.hpp"
#include "CommandBase.hpp"
#include "Logger.hpp"
#include "ScopeTimer.hpp"

#include <memory>
#include <CLI/CLI.hpp>

namespace rhbm_gem {

namespace {

void RegisterBuiltInCommand(CLI::App & cli_app, const CommandDescriptor & descriptor)
{
    auto command_object{ descriptor.factory() };
    CLI::App * command{
        cli_app.add_subcommand(
            std::string(descriptor.name),
            std::string(descriptor.description))
    };
    command_object->RegisterCLIOptions(command);

    auto shared_cmd{ std::shared_ptr<CommandBase>(std::move(command_object)) };
    command->callback([cmd = std::move(shared_cmd)]() {
        ScopeTimer timer("Command in Application");
        if (!cmd->Execute())
        {
            throw CLI::RuntimeError(1);
        }
    });
}

} // namespace

Application::Application(CLI::App & app) :
    m_cli_app{ app }
{
    m_cli_app.require_subcommand(1);
    RegisterAllCommands();
}

void Application::RegisterAllCommands()
{
    for (const auto & descriptor : BuiltInCommandCatalog())
    {
        RegisterBuiltInCommand(m_cli_app, descriptor);
    }
}

} // namespace rhbm_gem
