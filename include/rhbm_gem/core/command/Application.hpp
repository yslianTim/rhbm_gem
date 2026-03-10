#pragma once

#include <rhbm_gem/data/io/DataIoServices.hpp>

namespace CLI
{
    class App;
}

namespace rhbm_gem {

class Application {
public:
    Application(::CLI::App & app, const DataIoServices & data_io_services);
    ~Application() = default;

private:
    ::CLI::App & m_cli_app;
    DataIoServices m_data_io_services;

    void RegisterAllCommands();
};

} // namespace rhbm_gem
