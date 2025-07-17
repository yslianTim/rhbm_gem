#pragma once

#include <string>

namespace CLI
{
    class App;
}

class Application
{
    CLI::App & m_cli_app;

public:
    Application(CLI::App & app);
    ~Application() = default;

private:
    void RegisterAllCommands(void);
    template<typename Type>
    void RegisterCommand(const std::string & name, const std::string & description);

};
