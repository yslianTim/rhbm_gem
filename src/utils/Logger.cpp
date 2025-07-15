#include "Logger.hpp"

#include <iostream>

int Logger::m_current_level{ static_cast<int>(LogLevel::Info) };

void Logger::SetLogLevel(int level)
{
    m_current_level = level;
}

void Logger::Log(LogLevel level, const std::string & message)
{
    if (static_cast<int>(level) <= m_current_level)
    {
        if (level == LogLevel::Error || level == LogLevel::Warning)
        {
            std::cerr << message << std::endl;
        }
        else
        {
            std::cout << message << std::endl;
        }
    }
}
