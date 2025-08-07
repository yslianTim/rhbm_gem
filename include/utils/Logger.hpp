#pragma once

#include <string>
#include <atomic>
#include <mutex>

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
    static std::atomic<LogLevel> m_current_level;
    static std::mutex m_stream_mutex;

public:
    static LogLevel GetLogLevel(void) { return m_current_level.load(); }
    static void SetLogLevel(int level);
    static void SetLogLevel(LogLevel level);
    static void Log(LogLevel level, const std::string & message);

private:
    static LogLevel NormalizeLevel(LogLevel level);
};
