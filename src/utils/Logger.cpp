#include "Logger.hpp"

#include <iostream>

std::atomic<LogLevel> Logger::m_current_level{ LogLevel::Info };
std::mutex Logger::m_stream_mutex{};

LogLevel Logger::NormalizeLevel(LogLevel level)
{
    if (level < LogLevel::Error || level > LogLevel::Debug)
    {
        return LogLevel::Info;
    }
    return level;
}

void Logger::SetLogLevel(int level)
{
    m_current_level = NormalizeLevel(static_cast<LogLevel>(level));
}

void Logger::SetLogLevel(LogLevel level)
{
    m_current_level = NormalizeLevel(level);
}

void Logger::Log(LogLevel level, const std::string & message)
{
    if (level < LogLevel::Error || level > LogLevel::Debug)
    {
        std::cerr << "[Unknown] " << message << std::endl;
        return;
    }

    std::lock_guard<std::mutex> lock(m_stream_mutex);
    if (level > m_current_level.load()) return;

    switch (level)
    {
        case LogLevel::Error:
            std::cerr << "[Error] " << message << '\n' << std::flush;
            break;
        case LogLevel::Warning:
            std::cerr << "[Warning] " << message << '\n' << std::flush;
            break;
        case LogLevel::Notice:
            std::cout << "[Notice] " << message << '\n';
            break;
        case LogLevel::Info:
            std::cout << message << '\n';
            break;
        case LogLevel::Debug:
            std::cout << "[Debug] " << message << '\n';
            break;
    }
}
