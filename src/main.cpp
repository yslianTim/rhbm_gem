#include <iostream>
#include <CLI/CLI.hpp>

#include <rhbm_gem/core/Application.hpp>

int main(int argc, char * argv[])
{
    CLI::App app{"RHBM-GEM"};

    rhbm_gem::Application app_controller(app);

    CLI11_PARSE(app, argc, argv);

    return 0;
}
