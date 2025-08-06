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
                             LoggerSuppressionTestCase{ LogLevel::Info,    LogLevel::Debug   }));

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
                         ::testing::Values(LoggerOutputTestCase{ LogLevel::Error,   "[Error] ",   true  },
                                           LoggerOutputTestCase{ LogLevel::Warning, "[Warning] ", true  },
                                           LoggerOutputTestCase{ LogLevel::Notice,  "[Notice] ",  false },
                                           LoggerOutputTestCase{ LogLevel::Info,    "",           false },
                                           LoggerOutputTestCase{ LogLevel::Debug,   "[Debug] ",   false }));

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
    CurrentGreaterThanMessage,
    LoggerLevelPairTest,
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
        LoggerLevelPairTestCase{ LogLevel::Debug,   LogLevel::Info,    "",           false }));

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
