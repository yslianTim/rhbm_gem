#pragma once

namespace CLI
{
    class App;
}
class GlobalOptions;

class CommandBase
{
public:
    virtual ~CommandBase() = default;
    virtual void Execute(void) = 0;
    virtual void RegisterCLIOptions(CLI::App * command) = 0;
    virtual void SetGlobalOptions(const GlobalOptions & globals) = 0;
};
