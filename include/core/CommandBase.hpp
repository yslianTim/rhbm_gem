#pragma once

#include <iostream>

class CommandBase
{
public:
    virtual ~CommandBase() = default;
    virtual void Execute(void) = 0;
};