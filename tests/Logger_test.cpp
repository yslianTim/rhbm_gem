#include <gtest/gtest.h>
#include <algorithm>
#include <thread>
#include <vector>

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
    ::testing::Values(
        LogLevel::Error,
        LogLevel::Warning,
        LogLevel::Notice,
        LogLevel::Info,
        LogLevel::Debug
    )
);

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

TEST(LoggerTest, ThreadedLoggingAndLevelChanges)
{
    const int num_threads{4};
    const int messages_per_thread{50};
    testing::internal::CaptureStdout();
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int i = 0; i < num_threads; i++)
    {
        threads.emplace_back([messages_per_thread]() {
            for (int j = 0; j < messages_per_thread; j++)
            {
                Logger::SetLogLevel(LogLevel::Debug);
                Logger::Log(LogLevel::Info, "msg");
            }
            Logger::SetLogLevel(LogLevel::Info);
        });
    }
    for (auto & t : threads) t.join();
    const std::string out{ testing::internal::GetCapturedStdout() };
    const auto line_count{ std::count(out.begin(), out.end(), '\n') };
    EXPECT_EQ(static_cast<std::size_t>(num_threads * messages_per_thread), line_count);
    EXPECT_EQ(LogLevel::Info, Logger::GetLogLevel());
}

struct LoggerSuppressionTestCase
{
    LogLevel current_level;
    LogLevel message_level;
};

class LoggerSuppressionTest : public ::testing::TestWithParam<LoggerSuppressionTestCase>
{
protected:
    void TearDown(void) override
    {
        Logger::SetLogLevel(LogLevel::Info);
    }
};

INSTANTIATE_TEST_SUITE_P(SuppressedPairs, LoggerSuppressionTest,
    ::testing::Values(
        LoggerSuppressionTestCase{ LogLevel::Error,   LogLevel::Warning },
        LoggerSuppressionTestCase{ LogLevel::Error,   LogLevel::Notice  },
        LoggerSuppressionTestCase{ LogLevel::Error,   LogLevel::Info    },
        LoggerSuppressionTestCase{ LogLevel::Error,   LogLevel::Debug   },
        LoggerSuppressionTestCase{ LogLevel::Warning, LogLevel::Notice  },
        LoggerSuppressionTestCase{ LogLevel::Warning, LogLevel::Info    },
        LoggerSuppressionTestCase{ LogLevel::Warning, LogLevel::Debug   },
        LoggerSuppressionTestCase{ LogLevel::Notice,  LogLevel::Info    },
        LoggerSuppressionTestCase{ LogLevel::Notice,  LogLevel::Debug   },
        LoggerSuppressionTestCase{ LogLevel::Info,    LogLevel::Debug   }
    )
);

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

TEST(LoggerTest, SetLogLevelEnumOutOfRangeDefaultsToInfo)
{
    Logger::SetLogLevel(LogLevel::Debug);
    Logger::SetLogLevel(static_cast<LogLevel>(-1));
    EXPECT_EQ(LogLevel::Info, Logger::GetLogLevel());
}

TEST(LoggerTest, SetLogLevelEnumAboveRangeDefaultsToInfo)
{
    Logger::SetLogLevel(LogLevel::Error);
    Logger::SetLogLevel(static_cast<LogLevel>(5));
    EXPECT_EQ(LogLevel::Info, Logger::GetLogLevel());
}

TEST(LoggerTest, UnknownLogLevelDefaultsToDiagnostic)
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

TEST(LoggerTest, LogsUnknownLevelAboveDebug)
{
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    const auto level_above_debug{ static_cast<LogLevel>(static_cast<int>(LogLevel::Debug) + 1) };
    Logger::Log(level_above_debug, "msg");
    const std::string out{ testing::internal::GetCapturedStdout() };
    const std::string err{ testing::internal::GetCapturedStderr() };
    EXPECT_TRUE(out.empty());
    EXPECT_EQ(std::string("[Unknown] msg\n"), err);
}

TEST(LoggerTest, InfoLogsAppearOnSeparateLines)
{
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    Logger::Log(LogLevel::Info, "first");
    Logger::Log(LogLevel::Info, "second");
    const std::string stdout_output{ testing::internal::GetCapturedStdout() };
    EXPECT_EQ(std::string("first\nsecond\n"), stdout_output);
    Logger::SetLogLevel(LogLevel::Info);
}

TEST(LoggerTest, UnknownNegativeLevelLogsAtError)
{
    Logger::SetLogLevel(LogLevel::Error);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    Logger::Log(static_cast<LogLevel>(-1), "msg");
    const std::string stdout_output{ testing::internal::GetCapturedStdout() };
    const std::string stderr_output{ testing::internal::GetCapturedStderr() };
    EXPECT_TRUE(stdout_output.empty());
    EXPECT_EQ(std::string("[Unknown] msg\n"), stderr_output);
    Logger::SetLogLevel(LogLevel::Info);
}

struct LoggerOutputTestCase
{
    LogLevel level;
    const char * prefix;
    bool to_stderr;
};

class LoggerOutputTest : public ::testing::TestWithParam<LoggerOutputTestCase>
{
protected:
    void TearDown(void) override
    {
        Logger::SetLogLevel(LogLevel::Info);
    }
};

INSTANTIATE_TEST_SUITE_P(AllLogLevelsOutput, LoggerOutputTest,
    ::testing::Values(
        LoggerOutputTestCase{ LogLevel::Error,   "[Error] ",   true  },
        LoggerOutputTestCase{ LogLevel::Warning, "[Warning] ", true  },
        LoggerOutputTestCase{ LogLevel::Notice,  "[Notice] ",  false },
        LoggerOutputTestCase{ LogLevel::Info,    "",           false },
        LoggerOutputTestCase{ LogLevel::Debug,   "[Debug] ",   false }
    )
);

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
        std::string(tc.prefix) + "msg1\n" +
        std::string(tc.prefix) + "msg2\n"
    };
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
    void TearDown(void) override
    {
        Logger::SetLogLevel(LogLevel::Info);
    }
};

INSTANTIATE_TEST_SUITE_P(
    CurrentGreaterThanMessage, LoggerLevelPairTest,
    ::testing::Values(
        LoggerLevelPairTestCase{ LogLevel::Warning, LogLevel::Error,   "[Error] ",   true  },
        LoggerLevelPairTestCase{ LogLevel::Notice,  LogLevel::Error,   "[Error] ",   true  },
        LoggerLevelPairTestCase{ LogLevel::Notice,  LogLevel::Warning, "[Warning] ", true  },
        LoggerLevelPairTestCase{ LogLevel::Info,    LogLevel::Error,   "[Error] ",   true  },
        LoggerLevelPairTestCase{ LogLevel::Info,    LogLevel::Warning, "[Warning] ", true  },
        LoggerLevelPairTestCase{ LogLevel::Info,    LogLevel::Notice,  "[Notice] ",  false },
        LoggerLevelPairTestCase{ LogLevel::Debug,   LogLevel::Error,   "[Error] ",   true  },
        LoggerLevelPairTestCase{ LogLevel::Debug,   LogLevel::Warning, "[Warning] ", true  },
        LoggerLevelPairTestCase{ LogLevel::Debug,   LogLevel::Notice,  "[Notice] ",  false },
        LoggerLevelPairTestCase{ LogLevel::Debug,   LogLevel::Info,    "",           false }
    )
);

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
