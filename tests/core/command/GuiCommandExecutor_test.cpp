#include <gtest/gtest.h>

#include <fstream>
#include <type_traits>

#include "CommandTestHelpers.hpp"
#include "CommandValidationAssertions.hpp"
#include <rhbm_gem/gui/GuiCommandExecutor.hpp>

namespace rg = rhbm_gem;

TEST(GuiCommandExecutorTest, MapSimulationRequiresModelFilePath)
{
    const auto result{
        rg::gui::GuiCommandExecutor::ExecuteMapSimulation(rg::gui::MapSimulationRequest{})
    };

    EXPECT_FALSE(result.prepared);
    EXPECT_FALSE(result.executed);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.validation_issues,
            "--model",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}

TEST(GuiCommandExecutorTest, MapSimulationRejectsInvalidBlurringWidthList)
{
    rg::gui::MapSimulationRequest request;
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.blurring_width_list = "-1.0,0.0";

    const auto result{ rg::gui::GuiCommandExecutor::ExecuteMapSimulation(request) };
    EXPECT_FALSE(result.prepared);
    EXPECT_FALSE(result.executed);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.validation_issues,
            "--blurring-width",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

TEST(GuiCommandExecutorTest, PotentialAnalysisValidatesSimulationResolution)
{
    command_test::ScopedTempDir temp_dir{ "gui_executor_potential_analysis" };
    const auto dummy_map_path{ temp_dir.path() / "dummy.map" };
    std::ofstream(dummy_map_path) << "dummy";

    rg::gui::PotentialAnalysisRequest request;
    request.common.database_path = temp_dir.path() / "demo.sqlite";
    request.model_file_path = command_test::TestDataPath("test_model.cif");
    request.map_file_path = dummy_map_path;
    request.simulation_flag = true;
    request.simulated_map_resolution = 0.0;

    const auto result{ rg::gui::GuiCommandExecutor::ExecutePotentialAnalysis(request) };
    EXPECT_FALSE(result.prepared);
    EXPECT_FALSE(result.executed);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.validation_issues,
            "--sim-resolution",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

TEST(GuiCommandExecutorTest, ResultDumpMapPathIsOnlyRequiredForMapPrinter)
{
    command_test::ScopedTempDir temp_dir{ "gui_executor_result_dump" };

    rg::gui::ResultDumpRequest map_request;
    map_request.common.database_path = temp_dir.path() / "demo.sqlite";
    map_request.model_key_tag_list = "demo_model";
    map_request.printer_choice = rg::PrinterType::MAP_VALUE;

    const auto map_result{ rg::gui::GuiCommandExecutor::ExecuteResultDump(map_request) };
    EXPECT_FALSE(map_result.prepared);
    EXPECT_FALSE(map_result.executed);
    EXPECT_NE(
        command_test::FindValidationIssue(
            map_result.validation_issues,
            "--map",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);

    rg::gui::ResultDumpRequest gaus_request;
    gaus_request.common.database_path = temp_dir.path() / "demo.sqlite";
    gaus_request.model_key_tag_list = "demo_model";
    gaus_request.printer_choice = rg::PrinterType::GAUS_ESTIMATES;

    const auto gaus_result{ rg::gui::GuiCommandExecutor::ExecuteResultDump(gaus_request) };
    EXPECT_TRUE(gaus_result.prepared);
    EXPECT_FALSE(gaus_result.executed);
    EXPECT_EQ(
        command_test::FindValidationIssue(
            gaus_result.validation_issues,
            "--map",
            rg::ValidationPhase::Prepare,
            LogLevel::Error),
        nullptr);
}

namespace {

template <typename Type, typename = void>
struct HasRunMapSimulation : std::false_type
{
};

template <typename Type>
struct HasRunMapSimulation<Type, std::void_t<decltype(&Type::RunMapSimulation)>>
    : std::true_type
{
};

template <typename Type, typename = void>
struct HasRunPotentialAnalysis : std::false_type
{
};

template <typename Type>
struct HasRunPotentialAnalysis<Type, std::void_t<decltype(&Type::RunPotentialAnalysis)>>
    : std::true_type
{
};

template <typename Type, typename = void>
struct HasRunResultDump : std::false_type
{
};

template <typename Type>
struct HasRunResultDump<Type, std::void_t<decltype(&Type::RunResultDump)>>
    : std::true_type
{
};

} // namespace

TEST(GuiCommandExecutorTest, ExposesExecuteOnlySurfaceForCommandEntryPoints)
{
    EXPECT_FALSE(HasRunMapSimulation<rg::gui::GuiCommandExecutor>::value);
    EXPECT_FALSE(HasRunPotentialAnalysis<rg::gui::GuiCommandExecutor>::value);
    EXPECT_FALSE(HasRunResultDump<rg::gui::GuiCommandExecutor>::value);

    const auto result{
        rg::gui::GuiCommandExecutor::ExecuteMapSimulation(rg::gui::MapSimulationRequest{})
    };
    EXPECT_FALSE(result.prepared);
    EXPECT_FALSE(result.executed);
    EXPECT_NE(
        command_test::FindValidationIssue(
            result.validation_issues,
            "--model",
            rg::ValidationPhase::Parse,
            LogLevel::Error),
        nullptr);
}
