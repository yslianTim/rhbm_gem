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
        switch (level)
        {
            case LogLevel::Error:
                std::cerr << "[Error] " << message << std::endl;
                break;
            case LogLevel::Warning:
                std::cerr << "[Warning] " << message << std::endl;
                break;
            case LogLevel::Info:
                std::cout << message << std::endl;
                break;
            case LogLevel::Debug:
                std::cout << "[Debug] " << message << std::endl;
                break;
        }
    }
}
