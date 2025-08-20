#include "Logger.hpp"

#include <iostream>
#include <sstream>

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
    const auto current_level{ m_current_level.load() };

    if (normalized_level != level)
    {
        std::lock_guard<std::mutex> lock(m_stream_mutex);
        std::cerr << "[Unknown] " << message << '\n' << std::flush;
        return;
    }

    if (normalized_level > current_level) return;

    std::lock_guard<std::mutex> lock(m_stream_mutex);
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

void Logger::Log(LogLevel level, std::string_view message)
{
    Log(level, std::string{ message });
}

void Logger::Log(LogLevel level, const char * message)
{
    Log(level, std::string_view{ message });
}

void Logger::ProgressBar(std::size_t current, std::size_t total, std::size_t bar_width)
{
    if (total == 0) return;
    auto ratio{ static_cast<double>(current) / static_cast<double>(total) };
    auto pos{ static_cast<std::size_t>(ratio * static_cast<double>(bar_width)) };
    std::ostringstream oss;
    oss << '[';
    for (std::size_t i = 0; i < bar_width; i++)
    {
        if (i < pos) oss << '=';
        else if (i == pos && current < total) oss << '>';
        else oss << '.';
    }
    auto percent{ static_cast<int>(ratio * 100.0) };
    oss << "] " << percent << "% (" << current << '/' << total << ')';

    std::lock_guard<std::mutex> lock(m_stream_mutex);
    std::cout << '\r' << oss.str() << std::flush;
    if (current >= total)
    {
        std::cout << std::endl;
    }
}
