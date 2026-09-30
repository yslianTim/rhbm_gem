#include <gtest/gtest.h>

#include "core/detail/second_stage/observation/PerformanceCounters.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <initializer_list>
#include <string>

namespace {

TEST(SecondStagePerformanceTest, RuntimeSummaryRespectsQuietMode)
{
    using rhbm_gem::core::detail::PerformanceCounters;

    const auto previous_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    for (const bool quiet : { true, false })
    {
        testing::internal::CaptureStdout();
        {
            PerformanceCounters counters{ quiet };
            counters.RecordBoundaryReconciliation(4.25);
            counters.RecordBoundaryJointCorrection(3.0);
            counters.RecordDependencyPolish(8.5);
        }
        const auto output{ testing::internal::GetCapturedStdout() };
        if (quiet)
        {
            EXPECT_TRUE(output.empty());
        }
        else
        {
            EXPECT_NE(output.find("boundary_reconciliation_ms = 4.250"), std::string::npos);
            EXPECT_NE(output.find("dependency_polish_ms = 8.500"), std::string::npos);
        }
    }
    Logger::SetLogLevel(previous_level);
}

} // namespace
