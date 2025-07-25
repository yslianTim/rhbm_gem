#include <iostream>
#include <CLI/CLI.hpp>

#include "Application.hpp"
#include "FileProcessFactoryRegistry.hpp"

int main(int argc, char* argv[])
{
    CLI::App app{"App with subcommands"};

    FileProcessFactoryRegistry::Instance().RegisterDefaultFactories();

    Application app_controller(app);

    CLI11_PARSE(app, argc, argv);

    return 0;
}
