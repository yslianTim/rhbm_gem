#include <gtest/gtest.h>

#include "Logger.hpp"

class LoggerLogLevelTest : public ::testing::TestWithParam<LogLevel>
{
protected:
    void TearDown(void) override
    {
        Logger::SetLogLevel(LogLevel::Info);
    }
};

INSTANTIATE_TEST_SUITE_P(AllLogLevels, LoggerLogLevelTest,
                         ::testing::Values(LogLevel::Error, LogLevel::Warning,
                                           LogLevel::Notice, LogLevel::Info,
                                           LogLevel::Debug));

TEST_P(LoggerLogLevelTest, SetAndGetEnum)
{
    const auto level{ GetParam() };
    Logger::SetLogLevel(level);
    EXPECT_EQ(level, Logger::GetLogLevel());
}

TEST_P(LoggerLogLevelTest, SetAndGetInt)
{
    const auto level{GetParam()};
    Logger::SetLogLevel(static_cast<int>(level));
    EXPECT_EQ(level, Logger::GetLogLevel());
}

TEST(LoggerTest, DefaultLogLevelIsInfo)
{
    EXPECT_EQ(LogLevel::Info, Logger::GetLogLevel());
}

TEST(LoggerTest, InfoMessageIsSuppressedWhenLevelIsWarning)
{
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    Logger::Log(LogLevel::Info, "msg");
    const std::string out{ testing::internal::GetCapturedStdout() };
    const std::string err{ testing::internal::GetCapturedStderr() };
    EXPECT_TRUE(out.empty());
    EXPECT_TRUE(err.empty());
}

TEST(LoggerTest, DebugMessageIsSuppressedWhenLevelIsWarning)
{
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    Logger::Log(LogLevel::Debug, "msg");
    const std::string out{ testing::internal::GetCapturedStdout() };
    const std::string err{ testing::internal::GetCapturedStderr() };
    EXPECT_TRUE(out.empty());
    EXPECT_TRUE(err.empty());
}

TEST(LoggerTest, NoticeMessageIsSuppressedWhenLevelIsWarning)
{
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    Logger::Log(LogLevel::Notice, "msg");
    const std::string out{ testing::internal::GetCapturedStdout() };
    const std::string err{ testing::internal::GetCapturedStderr() };
    EXPECT_TRUE(out.empty());
    EXPECT_TRUE(err.empty());
}

TEST(LoggerTest, WarningMessageIsSuppressedWhenLevelIsError)
{
    Logger::SetLogLevel(LogLevel::Error);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    Logger::Log(LogLevel::Warning, "msg");
    const std::string out{ testing::internal::GetCapturedStdout() };
    const std::string err{ testing::internal::GetCapturedStderr() };
    EXPECT_TRUE(out.empty());
    EXPECT_TRUE(err.empty());
}

TEST(LoggerTest, SetLogLevelIntBelowRangeDefaultsToInfo)
{
    Logger::SetLogLevel(LogLevel::Debug);
    Logger::SetLogLevel(-1);
    EXPECT_EQ(LogLevel::Info, Logger::GetLogLevel());
}

TEST(LoggerTest, SetLogLevelIntAboveRangeDefaultsToInfo)
{
    Logger::SetLogLevel(LogLevel::Error);
    Logger::SetLogLevel(5);
    EXPECT_EQ(LogLevel::Info, Logger::GetLogLevel());
}

TEST(LoggerTest, SetLogLevelIntWithinRangeSetsLevel)
{
    Logger::SetLogLevel(static_cast<int>(LogLevel::Debug));
    EXPECT_EQ(LogLevel::Debug, Logger::GetLogLevel());
    Logger::SetLogLevel(LogLevel::Info);
}

TEST(LoggerTest, LogsToCorrectStreamWithPrefix)
{
    struct TestCase
    {
        LogLevel level;
        const char * prefix;
        bool to_stderr;
    };

    const std::array<TestCase, 5> cases{ {
        { LogLevel::Error,   "[Error] ",   true  },
        { LogLevel::Warning, "[Warning] ", true  },
        { LogLevel::Notice,  "[Notice] ",  false },
        { LogLevel::Info,    "",           false },
        { LogLevel::Debug,   "[Debug] ",   false }
    } };

    const auto prev_level{ Logger::GetLogLevel() };

    for (const auto & tc : cases)
    {
        Logger::SetLogLevel(LogLevel::Debug);

        testing::internal::CaptureStdout();
        testing::internal::CaptureStderr();

        Logger::Log(tc.level, "msg");

        const std::string stdout_output{ testing::internal::GetCapturedStdout() };
        const std::string stderr_output{ testing::internal::GetCapturedStderr() };

        const std::string expected{ std::string(tc.prefix) + "msg\n" };

        if (tc.to_stderr)
        {
            EXPECT_EQ(expected, stderr_output);
            EXPECT_TRUE(stdout_output.empty());
        }
        else
        {
            EXPECT_EQ(expected, stdout_output);
            EXPECT_TRUE(stderr_output.empty());
        }
    }

    Logger::SetLogLevel(prev_level);
}
