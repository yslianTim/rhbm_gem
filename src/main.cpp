#include <iostream>
#include <CLI/CLI.hpp>

#include "Application.hpp"

int main(int argc, char * argv[])
{
    CLI::App app{"RHBM-GEM"};

    Application app_controller(app);

    CLI11_PARSE(app, argc, argv);

    return 0;
}
