#pragma once

#include <string>

#include "CommandRegistry.hpp"

namespace CLI
{
    class App;
}

namespace rhbm_gem {

class Application
{
    ::CLI::App & m_cli_app;

public:
    Application(::CLI::App & app);
    ~Application() = default;

private:
    void RegisterAllCommands();
    void RegisterCommand(
        const std::string & name,
        const std::string & description,
        std::function<std::unique_ptr<CommandBase>()> factory);
};

} // namespace rhbm_gem
