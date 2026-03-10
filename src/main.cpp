#include <iostream>
#include <CLI/CLI.hpp>

#include <rhbm_gem/core/command/Application.hpp>
#include <rhbm_gem/data/io/DataIoServices.hpp>

int main(int argc, char * argv[])
{
    CLI::App app{"RHBM-GEM"};
    const auto data_io_services{ rhbm_gem::DataIoServices::BuildDefault() };

    rhbm_gem::Application app_controller(app, data_io_services);

    CLI11_PARSE(app, argc, argv);

    return 0;
}
