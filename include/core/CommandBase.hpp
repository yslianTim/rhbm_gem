#pragma once

#include "GlobalOptions.hpp"

namespace CLI
{
    class App;
}

class CommandBase
{
protected:
    GlobalOptions m_globals{};

public:
    virtual ~CommandBase() = default;
    virtual void Execute(void) = 0;
    virtual void RegisterCLIOptions(CLI::App * command) = 0;
    virtual void SetGlobalOptions(const GlobalOptions & globals) { m_globals = globals; }
};
