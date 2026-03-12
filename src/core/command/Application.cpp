#include <rhbm_gem/core/command/Application.hpp>
#include "internal/CommandCatalogInternal.hpp"
#include <rhbm_gem/core/command/CommandBase.hpp>
#include <rhbm_gem/data/io/DataIoServices.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>

#include <memory>
#include <CLI/CLI.hpp>

namespace rhbm_gem {

Application::Application(CLI::App & app, const DataIoServices & data_io_services) :
    m_cli_app{ app },
    m_data_io_services{ data_io_services }
{
    m_cli_app.require_subcommand(1);
    RegisterAllCommands();
}

void Application::RegisterAllCommands()
{
    for (const auto & descriptor : CommandCatalog())
    {
        auto command_object{ descriptor.factory(m_data_io_services) };
        CLI::App * command{
            m_cli_app.add_subcommand(
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
}

} // namespace rhbm_gem
