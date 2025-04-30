#pragma once

class CommandBase
{
public:
    virtual ~CommandBase() = default;
    virtual void Execute(void) = 0;
};