#pragma once

#include "CommandOptions.hpp"

namespace CLI
{
    class App;
}

class CommandBase
{
public:
    virtual ~CommandBase() = default;
    virtual bool Execute(void) = 0;
    virtual void RegisterCLIOptions(CLI::App * command) = 0;
    virtual CommandOptions & GetOptions(void) = 0;

protected:
    static void RegisterCommandOptions(CLI::App * command, CommandOptions & options);
};
