#include <gtest/gtest.h>

#include <algorithm>

#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/domain/Timer.hpp>

TEST(LoggerProgressTimerTest, ProgressBarUpdatesAndFinishes)
{
    testing::internal::CaptureStdout();
    size_t bar_width{ 10 };
    Logger::ProgressBar(1, 10, bar_width);
    Logger::ProgressBar(4, 10, bar_width);
    Logger::ProgressBar(10, 10, bar_width);
    const std::string out{ testing::internal::GetCapturedStdout() };
    EXPECT_NE(out.find("[=>........] 10% (1/10)"), std::string::npos);
    EXPECT_NE(out.find("[====>.....] 40% (4/10)"), std::string::npos);
    EXPECT_NE(out.find("[==========] 100% (10/10)"), std::string::npos);
    EXPECT_EQ('\n', out.back());
    const auto newline_count{ std::count(out.begin(), out.end(), '\n') };
    EXPECT_EQ(static_cast<std::size_t>(1), newline_count);
}

TEST(LoggerProgressTimerTest, ProgressPercentOnlyUpdatesOnChange)
{
    testing::internal::CaptureStdout();
    size_t bar_width{ 10 };
    Logger::ProgressPercent(1, 10, bar_width);
    Logger::ProgressPercent(1, 10, bar_width);
    Logger::ProgressPercent(2, 10, bar_width);
    Logger::ProgressPercent(10, 10, bar_width);
    const std::string out{ testing::internal::GetCapturedStdout() };
    EXPECT_EQ(static_cast<std::size_t>(3), std::count(out.begin(), out.end(), '\r'));
    EXPECT_NE(out.find("[=>........] 10%"), std::string::npos);
    EXPECT_NE(out.find("[==>.......] 20%"), std::string::npos);
    EXPECT_NE(out.find("[==========] 100%"), std::string::npos);
    EXPECT_EQ(std::string::npos, out.find('('));
    EXPECT_EQ('\n', out.back());
    const auto newline_count{ std::count(out.begin(), out.end(), '\n') };
    EXPECT_EQ(static_cast<std::size_t>(1), newline_count);
}

TEST(LoggerProgressTimerTest, ScopeTimerDefaultIsSuppressedAtInfo)
{
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    {
        ScopeTimer timer("suppressed timer");
    }
    const std::string out{ testing::internal::GetCapturedStdout() };
    EXPECT_TRUE(out.empty());
}

TEST(LoggerProgressTimerTest, ScopeTimerDefaultAppearsAtDebug)
{
    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    {
        ScopeTimer timer("visible timer");
    }
    const std::string out{ testing::internal::GetCapturedStdout() };
    EXPECT_NE(out.find("visible timer took"), std::string::npos);
}

TEST(LoggerProgressTimerTest, TimerPrintFormattedHonorsExplicitLogLevel)
{
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    Timer<> timer;
    timer.PrintFormatted("explicit timer", LogLevel::Info);
    const std::string out{ testing::internal::GetCapturedStdout() };
    EXPECT_NE(out.find("explicit timer took"), std::string::npos);
}
