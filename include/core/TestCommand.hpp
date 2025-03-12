#pragma once

#include <iostream>
#include "CommandBase.hpp"

class TestCommand : public CommandBase
{
    int m_thread_size;
    
public:
    TestCommand(void) = default;
    ~TestCommand() = default;
    void Execute(void) override;

    void SetThreadSize(int value) { m_thread_size = value; }
};