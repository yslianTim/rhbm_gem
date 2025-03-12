#include <iostream>
#include <memory>
#include <CLI/CLI.hpp>
#include "Application.hpp"

using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
    CLI::App app{"App with subcommands"};

    Application app_controller(&app);

    CLI11_PARSE(app, argc, argv);

    app_controller.Run();

    return 0;
}