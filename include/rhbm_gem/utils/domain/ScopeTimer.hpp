#pragma once

#include <rhbm_gem/utils/Timer.hpp>

class ScopeTimer
{
    Timer<> m_timer;
    std::string m_label;
    LogLevel m_level;

public:
    explicit ScopeTimer(const std::string & label, LogLevel level = LogLevel::Debug) :
        m_timer(), m_label(label), m_level(level) {}
    
    ~ScopeTimer()
    {
        m_timer.PrintFormatted(m_label, m_level);
    }
};
