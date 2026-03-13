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
    RegisterCommandSubcommands(m_cli_app);
}

} // namespace rhbm_gem
