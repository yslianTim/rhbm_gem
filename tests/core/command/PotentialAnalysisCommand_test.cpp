#include <gtest/gtest.h>

#include <filesystem>

#include <Eigen/Dense>

#include "CommandTestHelpers.hpp"
#include "CommandValidationAssertions.hpp"
#include "internal/PotentialAnalysisTrainingSupport.hpp"
#include <rhbm_gem/core/command/PotentialAnalysisCommand.hpp>

namespace rg = rhbm_gem;

TEST(PotentialAnalysisCommandTest, NegativeAlphaValuesBecomeValidationErrors)
{
    rg::PotentialAnalysisCommand command{ command_test::BuildDataIoServices() };
    command.SetAlphaR(-0.1);
    command.SetAlphaG(0.0);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(command, "--alpha-r", rg::ValidationPhase::Parse),
        nullptr);
    ASSERT_NE(
        command_test::FindValidationIssue(command, "--alpha-g", rg::ValidationPhase::Parse),
        nullptr);
}

TEST(PotentialAnalysisCommandTest, SimulationRequiresPositiveResolutionAtPrepare)
{
    rg::PotentialAnalysisCommand command{ command_test::BuildDataIoServices() };
    command.SetSimulationFlag(true);
    command.SetSimulatedMapResolution(0.0);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--sim-resolution",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

TEST(PotentialAnalysisCommandTest, NonPositiveSamplingHeightBecomesValidationError)
{
    rg::PotentialAnalysisCommand command{ command_test::BuildDataIoServices() };
    command.SetSamplingHeight(0.0);

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--sampling-height",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(PotentialAnalysisCommandTest, EmptySavedKeyBecomesValidationError)
{
    rg::PotentialAnalysisCommand command{ command_test::BuildDataIoServices() };
    command.SetSavedKeyTag("");

    EXPECT_FALSE(command.PrepareForExecution());
    ASSERT_NE(
        command_test::FindValidationIssue(
            command,
            "--save-key",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(PotentialAnalysisCommandTest, TrainingReportDirSetterAcceptsEmptyAndRelativePath)
{
    rg::PotentialAnalysisCommand command{ command_test::BuildDataIoServices() };
    command.SetTrainingReportDir(std::filesystem::path{});
    EXPECT_EQ(
        command_test::FindValidationIssue(
            command,
            "--training-report-dir",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);

    command.SetTrainingReportDir(std::filesystem::path{ "reports" });
    EXPECT_EQ(
        command_test::FindValidationIssue(
            command,
            "--training-report-dir",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(PotentialAnalysisCommandTest, TrainingReportEmissionIsSkippedWhenDirectoryIsUnset)
{
    const Eigen::MatrixXd bias_matrix{ Eigen::MatrixXd::Constant(3, 2, 1.0) };
    const std::vector<double> alpha_list{ 0.1, 0.2 };
    rg::PotentialAnalysisCommandOptions options;
    options.training_report_dir.clear();

    const bool emitted{
        rg::detail::EmitTrainingReportIfRequested(
            bias_matrix,
            alpha_list,
            "#alpha_{r}",
            "Deviation with OLS",
            options,
            "alpha_r_bias.pdf")
    };
    EXPECT_FALSE(emitted);
}

TEST(PotentialAnalysisCommandTest, TrainingReportEmissionAndAlphaRGridBehaviorAreDeterministic)
{
    const auto alpha_r_list{ rg::detail::BuildOrderedAlphaRTrainingList() };
    ASSERT_EQ(alpha_r_list.size(), 11u);
    EXPECT_DOUBLE_EQ(alpha_r_list.front(), 0.0);
    EXPECT_DOUBLE_EQ(alpha_r_list.back(), 1.0);

    command_test::ScopedTempDir temp_dir{ "potential_analysis_training_report" };
    const auto report_dir{ temp_dir.path() / "reports" };
    const Eigen::MatrixXd bias_matrix{ Eigen::MatrixXd::Constant(3, 2, 1.0) };
    const std::vector<double> alpha_list{ 0.1, 0.2 };
    rg::PotentialAnalysisCommandOptions options;
    options.training_report_dir = report_dir;

    const bool emitted{
        rg::detail::EmitTrainingReportIfRequested(
            bias_matrix,
            alpha_list,
            "#alpha_{g}",
            "Deviation with Mean",
            options,
            "alpha_g_bias.pdf")
    };

#ifdef HAVE_ROOT
    EXPECT_TRUE(emitted);
    EXPECT_TRUE(std::filesystem::exists(report_dir / "alpha_g_bias.pdf"));
#else
    EXPECT_FALSE(emitted);
    EXPECT_FALSE(std::filesystem::exists(report_dir / "alpha_g_bias.pdf"));
#endif
}
