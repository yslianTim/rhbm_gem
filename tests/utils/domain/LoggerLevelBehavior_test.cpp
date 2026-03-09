#include <gtest/gtest.h>

#include <algorithm>
#include <thread>
#include <vector>

#include <rhbm_gem/utils/domain/Logger.hpp>

class LoggerLogLevelTest : public ::testing::TestWithParam<LogLevel>
{
protected:
    void TearDown() override
    {
        Logger::SetLogLevel(LogLevel::Info);
    }
};

INSTANTIATE_TEST_SUITE_P(
    AllLogLevels,
    LoggerLogLevelTest,
    ::testing::Values(
        LogLevel::Error,
        LogLevel::Warning,
        LogLevel::Notice,
        LogLevel::Info,
        LogLevel::Debug));

TEST_P(LoggerLogLevelTest, SetAndGetEnum)
{
    const auto level{ GetParam() };
    Logger::SetLogLevel(level);
    EXPECT_EQ(level, Logger::GetLogLevel());
}

TEST_P(LoggerLogLevelTest, SetAndGetInt)
{
    const auto level{ GetParam() };
    Logger::SetLogLevel(static_cast<int>(level));
    EXPECT_EQ(level, Logger::GetLogLevel());
}

struct LoggerSuppressionTestCase
{
    LogLevel current_level;
    LogLevel message_level;
};

class LoggerSuppressionTest : public ::testing::TestWithParam<LoggerSuppressionTestCase>
{
protected:
    void TearDown() override
    {
        Logger::SetLogLevel(LogLevel::Info);
    }
};

INSTANTIATE_TEST_SUITE_P(
    SuppressedPairs,
    LoggerSuppressionTest,
    ::testing::Values(
        LoggerSuppressionTestCase{ LogLevel::Error, LogLevel::Warning },
        LoggerSuppressionTestCase{ LogLevel::Error, LogLevel::Notice },
        LoggerSuppressionTestCase{ LogLevel::Error, LogLevel::Info },
        LoggerSuppressionTestCase{ LogLevel::Error, LogLevel::Debug },
        LoggerSuppressionTestCase{ LogLevel::Warning, LogLevel::Notice },
        LoggerSuppressionTestCase{ LogLevel::Warning, LogLevel::Info },
        LoggerSuppressionTestCase{ LogLevel::Warning, LogLevel::Debug },
        LoggerSuppressionTestCase{ LogLevel::Notice, LogLevel::Info },
        LoggerSuppressionTestCase{ LogLevel::Notice, LogLevel::Debug },
        LoggerSuppressionTestCase{ LogLevel::Info, LogLevel::Debug }));

TEST_P(LoggerSuppressionTest, SuppressesOutput)
{
    const auto tc{ GetParam() };
    Logger::SetLogLevel(tc.current_level);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    Logger::Log(tc.message_level, "msg");
    const std::string out{ testing::internal::GetCapturedStdout() };
    const std::string err{ testing::internal::GetCapturedStderr() };
    EXPECT_TRUE(out.empty());
    EXPECT_TRUE(err.empty());
}

class LoggerUnknownLevelAboveDebugTest : public ::testing::TestWithParam<LogLevel>
{
protected:
    void TearDown() override
    {
        Logger::SetLogLevel(LogLevel::Info);
    }
};

INSTANTIATE_TEST_SUITE_P(
    AllLogLevels,
    LoggerUnknownLevelAboveDebugTest,
    ::testing::Values(
        LogLevel::Error,
        LogLevel::Warning,
        LogLevel::Notice,
        LogLevel::Info,
        LogLevel::Debug));

TEST_P(LoggerUnknownLevelAboveDebugTest, LogsUnknownLevelAboveDebug)
{
    const auto level{ GetParam() };
    Logger::SetLogLevel(level);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    const auto unknown_level{ static_cast<LogLevel>(static_cast<int>(LogLevel::Debug) + 1) };
    Logger::Log(unknown_level, "msg");
    const std::string out{ testing::internal::GetCapturedStdout() };
    const std::string err{ testing::internal::GetCapturedStderr() };
    EXPECT_TRUE(out.empty());
    EXPECT_EQ(std::string("[Unknown] msg\n"), err);
}

TEST(LoggerLevelBehaviorTest, DefaultLogLevelIsInfo)
{
    EXPECT_EQ(LogLevel::Info, Logger::GetLogLevel());
}

TEST(LoggerLevelBehaviorTest, ThreadedLoggingAndLevelChanges)
{
    const int num_threads{ 4 };
    const int messages_per_thread{ 50 };
    testing::internal::CaptureStdout();
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int i = 0; i < num_threads; i++)
    {
        threads.emplace_back(
            [&]()
            {
                for (int j = 0; j < messages_per_thread; j++)
                {
                    Logger::SetLogLevel(LogLevel::Debug);
                    Logger::Log(LogLevel::Info, "msg");
                }
                Logger::SetLogLevel(LogLevel::Info);
            });
    }
    for (auto & thread : threads)
    {
        thread.join();
    }

    const std::string out{ testing::internal::GetCapturedStdout() };
    const auto line_count{ std::count(out.begin(), out.end(), '\n') };
    EXPECT_EQ(static_cast<std::size_t>(num_threads * messages_per_thread), line_count);
    EXPECT_EQ(LogLevel::Info, Logger::GetLogLevel());
}

TEST(LoggerLevelBehaviorTest, ThreadedErrorLoggingAndLevelChanges)
{
    const int num_threads{ 4 };
    const int messages_per_thread{ 50 };
    testing::internal::CaptureStderr();
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int i = 0; i < num_threads; i++)
    {
        threads.emplace_back(
            [&]()
            {
                for (int j = 0; j < messages_per_thread; j++)
                {
                    Logger::SetLogLevel(LogLevel::Debug);
                    Logger::Log(LogLevel::Error, "msg");
                }
                Logger::SetLogLevel(LogLevel::Info);
            });
    }
    for (auto & thread : threads)
    {
        thread.join();
    }

    const std::string err{ testing::internal::GetCapturedStderr() };
    const auto line_count{ std::count(err.begin(), err.end(), '\n') };
    EXPECT_EQ(static_cast<std::size_t>(num_threads * messages_per_thread), line_count);
    EXPECT_EQ(LogLevel::Info, Logger::GetLogLevel());
}

TEST(LoggerLevelBehaviorTest, SetLogLevelIntBelowRangeDefaultsToInfo)
{
    Logger::SetLogLevel(LogLevel::Debug);
    Logger::SetLogLevel(-1);
    EXPECT_EQ(LogLevel::Info, Logger::GetLogLevel());
}

TEST(LoggerLevelBehaviorTest, SetLogLevelIntAboveRangeDefaultsToInfo)
{
    Logger::SetLogLevel(LogLevel::Error);
    Logger::SetLogLevel(5);
    EXPECT_EQ(LogLevel::Info, Logger::GetLogLevel());
}

TEST(LoggerLevelBehaviorTest, SetLogLevelIntWithinRangeSetsLevel)
{
    Logger::SetLogLevel(static_cast<int>(LogLevel::Debug));
    EXPECT_EQ(LogLevel::Debug, Logger::GetLogLevel());
    Logger::SetLogLevel(LogLevel::Info);
}

TEST(LoggerLevelBehaviorTest, SetLogLevelEnumOutOfRangeDefaultsToInfo)
{
    Logger::SetLogLevel(LogLevel::Debug);
    Logger::SetLogLevel(static_cast<LogLevel>(-1));
    EXPECT_EQ(LogLevel::Info, Logger::GetLogLevel());
}

TEST(LoggerLevelBehaviorTest, SetLogLevelEnumAboveRangeDefaultsToInfo)
{
    Logger::SetLogLevel(LogLevel::Error);
    Logger::SetLogLevel(static_cast<LogLevel>(5));
    EXPECT_EQ(LogLevel::Info, Logger::GetLogLevel());
}

TEST(LoggerLevelBehaviorTest, UnknownLogLevelDefaultsToDiagnostic)
{
    const auto bogus_level{ static_cast<LogLevel>(-1) };
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    Logger::Log(bogus_level, "msg");
    const std::string stdout_output{ testing::internal::GetCapturedStdout() };
    const std::string stderr_output{ testing::internal::GetCapturedStderr() };
    EXPECT_TRUE(stdout_output.empty());
    EXPECT_EQ(std::string("[Unknown] msg\n"), stderr_output);
}

TEST(LoggerLevelBehaviorTest, InfoLogsAppearOnSeparateLines)
{
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    Logger::Log(LogLevel::Info, "first");
    Logger::Log(LogLevel::Info, "second");
    const std::string stdout_output{ testing::internal::GetCapturedStdout() };
    EXPECT_EQ(std::string("first\nsecond\n"), stdout_output);
    Logger::SetLogLevel(LogLevel::Info);
}

TEST(LoggerLevelBehaviorTest, LogAcceptsStringView)
{
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    std::string_view msg{ "sv_msg" };
    Logger::Log(LogLevel::Info, msg);
    const std::string out{ testing::internal::GetCapturedStdout() };
    EXPECT_EQ(std::string("sv_msg\n"), out);
}

TEST(LoggerLevelBehaviorTest, LogAcceptsStdString)
{
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    std::string msg{ "std_msg" };
    Logger::Log(LogLevel::Info, msg);
    const std::string out{ testing::internal::GetCapturedStdout() };
    EXPECT_EQ(std::string("std_msg\n"), out);
}

TEST(LoggerLevelBehaviorTest, UnknownNegativeLevelOutputsUnknownAtAllLevels)
{
    for (auto level : { LogLevel::Error,
                        LogLevel::Warning,
                        LogLevel::Notice,
                        LogLevel::Info,
                        LogLevel::Debug })
    {
        Logger::SetLogLevel(level);
        testing::internal::CaptureStdout();
        testing::internal::CaptureStderr();
        Logger::Log(static_cast<LogLevel>(-1), "msg");
        const std::string stdout_output{ testing::internal::GetCapturedStdout() };
        const std::string stderr_output{ testing::internal::GetCapturedStderr() };
        EXPECT_TRUE(stdout_output.empty());
        EXPECT_EQ(std::string("[Unknown] msg\n"), stderr_output);
    }
    Logger::SetLogLevel(LogLevel::Info);
}

TEST(LoggerLevelBehaviorTest, UnknownLevelAlwaysLogs)
{
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    Logger::Log(static_cast<LogLevel>(-1), "msg");
    const std::string stdout_output{ testing::internal::GetCapturedStdout() };
    const std::string stderr_output{ testing::internal::GetCapturedStderr() };
    EXPECT_TRUE(stdout_output.empty());
    EXPECT_EQ(std::string("[Unknown] msg\n"), stderr_output);

    Logger::SetLogLevel(LogLevel::Error);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    Logger::Log(static_cast<LogLevel>(-1), "msg");
    const std::string stdout_output2{ testing::internal::GetCapturedStdout() };
    const std::string stderr_output2{ testing::internal::GetCapturedStderr() };
    EXPECT_TRUE(stdout_output2.empty());
    EXPECT_EQ(std::string("[Unknown] msg\n"), stderr_output2);
    Logger::SetLogLevel(LogLevel::Info);
}
