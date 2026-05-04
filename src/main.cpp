#include <iostream>
#include <CLI/CLI.hpp>

#include <rhbm_gem/core/command/CommandSystem.hpp>

int main(int argc, char * argv[])
{
    CLI::App app{"RHBM-GEM"};

    rhbm_gem::ConfigureCommandCLI(app);

    CLI11_PARSE(app, argc, argv);

    return 0;
}
