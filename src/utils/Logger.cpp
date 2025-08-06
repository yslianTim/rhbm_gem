#include "Logger.hpp"

#include <iostream>

LogLevel Logger::m_current_level{ LogLevel::Info };

void Logger::SetLogLevel(int level)
{
    if (level < static_cast<int>(LogLevel::Error) ||
        level > static_cast<int>(LogLevel::Debug))
    {
        m_current_level = LogLevel::Info;
        return;
    }
    m_current_level = static_cast<LogLevel>(level);
}

void Logger::SetLogLevel(LogLevel level)
{
    if (level < LogLevel::Error || level > LogLevel::Debug)
    {
        m_current_level = LogLevel::Info;
        return;
    }
    m_current_level = level;
}

LogLevel Logger::GetLogLevel(void)
{
    return m_current_level;
}

void Logger::Log(LogLevel level, const std::string & message)
{
    if (level > m_current_level) return;

    switch (level)
    {
        case LogLevel::Error:
            std::cerr << "[Error] " << message << std::endl;
            break;
        case LogLevel::Warning:
            std::cerr << "[Warning] " << message << std::endl;
            break;
        case LogLevel::Notice:
            std::cout << "[Notice] " << message << std::endl;
            break;
        case LogLevel::Info:
            std::cout << message << std::endl;
            break;
        case LogLevel::Debug:
            std::cout << "[Debug] " << message << std::endl;
            break;
        default:
            std::cerr << "[Unknown] " << message << std::endl;
            break;
    }
}
