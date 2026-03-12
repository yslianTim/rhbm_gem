#pragma once

namespace CLI
{
    class App;
}

namespace rhbm_gem {

class Application {
public:
    explicit Application(::CLI::App & app);
    ~Application() = default;

private:
    ::CLI::App & m_cli_app;

    void RegisterAllCommands();
};

} // namespace rhbm_gem
