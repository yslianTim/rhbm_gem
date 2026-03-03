#pragma once

namespace CLI
{
    class App;
}

namespace rhbm_gem {

struct CommandDescriptor;

class Application {
public:
    Application(::CLI::App & app);
    ~Application() = default;

private:
    ::CLI::App & m_cli_app;

    void RegisterAllCommands();
    void RegisterCommand(const CommandDescriptor & descriptor);
};

} // namespace rhbm_gem
