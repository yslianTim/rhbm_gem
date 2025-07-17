#pragma once

#include <memory>
#include <string>
#include <vector>

#include "GlobalOptions.hpp"

namespace CLI
{
    class App;
}

class CommandBase;

class Application
{
    CLI::App & m_cli_app;
    GlobalOptions m_global_options{};
    std::vector<std::unique_ptr<CommandBase>> m_commands;

public:
    Application(CLI::App & app);
    ~Application();

private:
    void RegisterAllCommands(void);
    template<typename Type>
    void RegisterCommand(const std::string & name, const std::string & description);
    void RegisterGlobalOptions(CLI::App * command);

};
