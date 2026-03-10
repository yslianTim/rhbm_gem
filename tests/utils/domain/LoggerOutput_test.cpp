#include <gtest/gtest.h>

#include <rhbm_gem/utils/domain/Logger.hpp>

struct LoggerOutputTestCase
{
    LogLevel level;
    const char * prefix;
    bool to_stderr;
};

class LoggerOutputTest : public ::testing::TestWithParam<LoggerOutputTestCase>
{
protected:
    void TearDown() override
    {
        Logger::SetLogLevel(LogLevel::Info);
    }
};

INSTANTIATE_TEST_SUITE_P(
    AllLogLevelsOutput,
    LoggerOutputTest,
    ::testing::Values(
        LoggerOutputTestCase{ LogLevel::Error, "[Error] ", true },
        LoggerOutputTestCase{ LogLevel::Warning, "[Warning] ", true },
        LoggerOutputTestCase{ LogLevel::Notice, "[Notice] ", false },
        LoggerOutputTestCase{ LogLevel::Info, "", false },
        LoggerOutputTestCase{ LogLevel::Debug, "[Debug] ", false }));

TEST_P(LoggerOutputTest, LogsToCorrectStreamWithPrefix)
{
    const auto tc{ GetParam() };
    Logger::SetLogLevel(tc.level);
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

TEST_P(LoggerOutputTest, MessagesAppearOnSeparateLines)
{
    const auto tc{ GetParam() };
    Logger::SetLogLevel(tc.level);
    const std::string expected{
        std::string(tc.prefix) + "msg1\n" + std::string(tc.prefix) + "msg2\n" };
    if (tc.to_stderr)
    {
        testing::internal::CaptureStderr();
        Logger::Log(tc.level, "msg1");
        Logger::Log(tc.level, "msg2");
        const std::string stderr_output{ testing::internal::GetCapturedStderr() };
        EXPECT_EQ(expected, stderr_output);
    }
    else
    {
        testing::internal::CaptureStdout();
        Logger::Log(tc.level, "msg1");
        Logger::Log(tc.level, "msg2");
        const std::string stdout_output{ testing::internal::GetCapturedStdout() };
        EXPECT_EQ(expected, stdout_output);
    }
}

struct LoggerLevelPairTestCase
{
    LogLevel current_level;
    LogLevel message_level;
    const char * prefix;
    bool to_stderr;
};

class LoggerLevelPairTest : public ::testing::TestWithParam<LoggerLevelPairTestCase>
{
protected:
    void TearDown() override
    {
        Logger::SetLogLevel(LogLevel::Info);
    }
};

INSTANTIATE_TEST_SUITE_P(
    CurrentGreaterThanMessage,
    LoggerLevelPairTest,
    ::testing::Values(
        LoggerLevelPairTestCase{ LogLevel::Warning, LogLevel::Error, "[Error] ", true },
        LoggerLevelPairTestCase{ LogLevel::Notice, LogLevel::Error, "[Error] ", true },
        LoggerLevelPairTestCase{ LogLevel::Notice, LogLevel::Warning, "[Warning] ", true },
        LoggerLevelPairTestCase{ LogLevel::Info, LogLevel::Error, "[Error] ", true },
        LoggerLevelPairTestCase{ LogLevel::Info, LogLevel::Warning, "[Warning] ", true },
        LoggerLevelPairTestCase{ LogLevel::Info, LogLevel::Notice, "[Notice] ", false },
        LoggerLevelPairTestCase{ LogLevel::Debug, LogLevel::Error, "[Error] ", true },
        LoggerLevelPairTestCase{ LogLevel::Debug, LogLevel::Warning, "[Warning] ", true },
        LoggerLevelPairTestCase{ LogLevel::Debug, LogLevel::Notice, "[Notice] ", false },
        LoggerLevelPairTestCase{ LogLevel::Debug, LogLevel::Info, "", false }));

TEST_P(LoggerLevelPairTest, LogsMessageWithCorrectPrefixAndStream)
{
    const auto tc{ GetParam() };
    Logger::SetLogLevel(tc.current_level);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    Logger::Log(tc.message_level, "msg");
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
