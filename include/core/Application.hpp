#pragma once

#include <memory>
#include <string>

#include "GlobalOptions.hpp"

namespace CLI
{
    class App;
}

class Application
{
    CLI::App & m_cli_app;
    GlobalOptions m_global_options{};

public:
    Application(CLI::App & app);
    ~Application() = default;

private:
    void RegisterAllCommands(void);
    template<typename Type>
    void RegisterCommand(const std::string & name, const std::string & description);
    void RegisterGlobalOptions(CLI::App * command);

};
