#include <iostream>
#include <CLI/CLI.hpp>

#include <rhbm_gem/core/command/CommandApi.hpp>

int main(int argc, char * argv[])
{
    CLI::App app{"RHBM-GEM"};

    rhbm_gem::ConfigureCommandCli(app);

    CLI11_PARSE(app, argc, argv);

    return 0;
}
