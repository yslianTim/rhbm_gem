#pragma once

#include <chrono>
#include <string>
#include <iostream>
#include <memory>
#include <type_traits>

using namespace std::chrono;

template<typename DurationUnit = microseconds>
class Timer final
{
    time_point<steady_clock> start_time;

public:
    Timer(void) : start_time{ steady_clock::now() } {}

    void Reset(void) { start_time = steady_clock::now(); }

    void Print(const std::string & text) const
    {
        std::cout << text << " took " << Elapsed() << " " << UnitString() << std::endl;
    }

    double Elapsed(void) const
    {
        auto end_time{ steady_clock::now() };
        auto duration{ duration_cast<DurationUnit>(end_time - start_time) };
        return static_cast<double>(duration.count());
    }

    std::string UnitString(void) const
    {
        if constexpr (std::is_same_v<DurationUnit, nanoseconds>) return "ns";
        else if constexpr (std::is_same_v<DurationUnit, microseconds>) return "us";
        else if constexpr (std::is_same_v<DurationUnit, milliseconds>) return "ms";
        else if constexpr (std::is_same_v<DurationUnit, seconds>) return "s";
        else if constexpr (std::is_same_v<DurationUnit, minutes>) return "min";
        else if constexpr (std::is_same_v<DurationUnit, hours>) return "h";
        else return "Unknown";
    }

    void PrintFormatted(const std::string & text) const
    {
        auto now{ steady_clock::now() };
        auto elapsed_time{ now - start_time };
        auto minutes_part{ duration_cast<minutes>(elapsed_time) };
        elapsed_time -= minutes_part;
        auto seconds_part{ duration_cast<seconds>(elapsed_time) };
        elapsed_time -= seconds_part;
        auto milliseconds_part{ duration_cast<milliseconds>(elapsed_time) };
        elapsed_time -= milliseconds_part;
        auto microseconds_part{ duration_cast<microseconds>(elapsed_time) };

        auto printed_any{ false };
        std::cout << text << " took ";
        if (minutes_part.count() != 0)
        {
            std::cout << minutes_part.count() << " min " << seconds_part.count() << " s";
            printed_any = true;
        }
        if (seconds_part.count() != 0 && printed_any == false)
        {
            std::cout << seconds_part.count() << "." << milliseconds_part.count() << " s";
            printed_any = true;
        }
        if (milliseconds_part.count() != 0 && printed_any == false)
        {
            std::cout << milliseconds_part.count() << "." << microseconds_part.count() << " ms";
            printed_any = true;
        }
        if (printed_any == false)
        {
            std::cout << microseconds_part.count() << " us";
        }
        std::cout << std::endl;
    }
};

class TimedComponent
{
protected:
    Timer<> m_timer;

public:
    void ResetTimer(void) { m_timer.Reset(); }
    void PrintTimer(const std::string & label) const { m_timer.PrintFormatted(label); }
};