#pragma once

#include <string>

enum class LogLevel : int
{
    Error   = 0,   // [Error]   Only error messages are shown.
    Warning = 1,   // [Warning] Warnings and errors are shown.
    Info    = 2,   // [Info]    Informational messages are shown.
    Debug   = 3    // [Debug]   Detailed debug output.
};

class Logger
{
    static int m_current_level;

public:
    static void SetLogLevel(int level);
    static void Log(LogLevel level, const std::string & message);
};
