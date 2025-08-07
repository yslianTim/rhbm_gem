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
    const auto normalized_level{ NormalizeLevel(level) };
    std::lock_guard<std::mutex> lock(m_stream_mutex);

    if (normalized_level != level)
    {
        if (normalized_level > m_current_level.load()) return;
        std::cerr << "[Unknown] " << message << '\n' << std::flush;
        return;
    }

    if (normalized_level > m_current_level.load()) return;

    switch (normalized_level)
    {
        case LogLevel::Error:
            std::cerr << "[Error] " << message << '\n' << std::flush;
            break;
        case LogLevel::Warning:
            std::cerr << "[Warning] " << message << '\n' << std::flush;
            break;
        case LogLevel::Notice:
            std::cout << "[Notice] " << message << '\n' << std::flush;
            break;
        case LogLevel::Info:
            std::cout << message << '\n' << std::flush;
            break;
        case LogLevel::Debug:
            std::cout << "[Debug] " << message << '\n' << std::flush;
            break;
    }
}
