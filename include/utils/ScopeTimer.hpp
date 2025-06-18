#pragma once

#include "Timer.hpp"

class ScopeTimer
{
    Timer<> m_timer;
    std::string m_label;

public:
    explicit ScopeTimer(const std::string & label) : m_timer(), m_label(label) {}
    
    ~ScopeTimer()
    {
        m_timer.PrintFormatted(m_label);
    }
};
