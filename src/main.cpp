#include <iostream>
#include <CLI/CLI.hpp>

#include "core/command/detail/CommandCli.hpp"

int main(int argc, char * argv[])
{
    CLI::App app{"RHBM-GEM"};

    rhbm_gem::ConfigureCommandCli(app);

    CLI11_PARSE(app, argc, argv);

    return 0;
}
