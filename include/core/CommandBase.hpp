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
    virtual bool ValidateOptions(void) const = 0;
    void RegisterCLIOptions(CLI::App * command);
    virtual void RegisterCLIOptionsExtend(CLI::App * command) = 0;
    virtual const CommandOptions & GetOptions(void) const = 0;
    virtual CommandOptions & GetOptions(void)
    {
        return const_cast<CommandOptions &>(static_cast<const CommandBase &>(*this).GetOptions());
    }

protected:
    void RegisterCLIOptionsBasic(CLI::App * command);
};
