#pragma once

#include <chrono>
#include <string>
#include <sstream>
#include <memory>
#include <type_traits>

#include "Logger.hpp"

template<typename DurationUnit = std::chrono::microseconds>
class Timer final
{
    std::chrono::time_point<std::chrono::steady_clock> start_time;

public:
    Timer(void) : start_time{ std::chrono::steady_clock::now() } {}

    void Reset(void) { start_time = std::chrono::steady_clock::now(); }

    void Print(const std::string & text, LogLevel level = LogLevel::Debug) const
    {
        std::ostringstream oss;
        oss << text << " took " << Elapsed() << " " << UnitString();
        Logger::Log(level, oss.str());
    }

    double Elapsed(void) const
    {
        auto end_time{ std::chrono::steady_clock::now() };
        auto duration{ std::chrono::duration_cast<DurationUnit>(end_time - start_time) };
        return static_cast<double>(duration.count());
    }

    std::string UnitString(void) const
    {
        if constexpr (std::is_same_v<DurationUnit, std::chrono::nanoseconds>) return "ns";
        else if constexpr (std::is_same_v<DurationUnit, std::chrono::microseconds>) return "us";
        else if constexpr (std::is_same_v<DurationUnit, std::chrono::milliseconds>) return "ms";
        else if constexpr (std::is_same_v<DurationUnit, std::chrono::seconds>) return "s";
        else if constexpr (std::is_same_v<DurationUnit, std::chrono::minutes>) return "min";
        else if constexpr (std::is_same_v<DurationUnit, std::chrono::hours>) return "h";
        else return "Unknown";
    }

    void PrintFormatted(const std::string & text, LogLevel level = LogLevel::Debug) const
    {
        auto now{ std::chrono::steady_clock::now() };
        auto elapsed_time{ now - start_time };
        auto minutes_part{ std::chrono::duration_cast<std::chrono::minutes>(elapsed_time) };
        elapsed_time -= minutes_part;
        auto seconds_part{ std::chrono::duration_cast<std::chrono::seconds>(elapsed_time) };
        elapsed_time -= seconds_part;
        auto milliseconds_part{ std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_time) };
        elapsed_time -= milliseconds_part;
        auto microseconds_part{ std::chrono::duration_cast<std::chrono::microseconds>(elapsed_time) };

        auto printed_any{ false };
        std::ostringstream oss;
        oss << text << " took ";
        if (minutes_part.count() != 0)
        {
            oss << minutes_part.count() << " min " << seconds_part.count() << " s";
            printed_any = true;
        }
        if (seconds_part.count() != 0 && printed_any == false)
        {
            oss << seconds_part.count() << '.' << milliseconds_part.count() << " s";
            printed_any = true;
        }
        if (milliseconds_part.count() != 0 && printed_any == false)
        {
            oss << milliseconds_part.count() << '.' << microseconds_part.count() << " ms";
            printed_any = true;
        }
        if (printed_any == false)
        {
            oss << microseconds_part.count() << " us";
        }
        Logger::Log(level, oss.str());
    }
};
