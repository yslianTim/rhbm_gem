#pragma once

#include <string>

enum class LogLevel : int
{
    Error   = 0,   // [Error]   Only error messages are shown.
    Warning = 1,   // [Warning] Warnings and errors are shown.
    Notice  = 2,   // [Notice]  Notice messages are shown.
    Info    = 3,   // [Info]    Informational messages are shown.
    Debug   = 4    // [Debug]   Detailed debug output.
};

class Logger
{
    static LogLevel m_current_level;

public:
    static void SetLogLevel(int level);
    static void SetLogLevel(LogLevel level);
    static LogLevel GetLogLevel(void);
    static void Log(LogLevel level, const std::string & message);
};
