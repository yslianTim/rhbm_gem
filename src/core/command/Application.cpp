#include <rhbm_gem/core/command/Application.hpp>
#include "internal/CommandCatalog.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>

#include <CLI/CLI.hpp>
#include <utility>

namespace rhbm_gem {

Application::Application(CLI::App & app) :
    m_cli_app{ app }
{
    m_cli_app.require_subcommand(1);
    RegisterAllCommands();
}

void Application::RegisterAllCommands()
{
    for (const auto & descriptor : CommandCatalog())
    {
        CLI::App * command{
            m_cli_app.add_subcommand(
                std::string(descriptor.name),
                std::string(descriptor.description))
        };
        auto run_command{ descriptor.bind_runtime(command) };

        command->callback([run = std::move(run_command)]() {
            ScopeTimer timer("Command in Application");
            const auto report{ run() };
            if (!report.executed)
            {
                throw CLI::RuntimeError(1);
            }
        });
    }
}

} // namespace rhbm_gem
