#include <rhbm_gem/core/command/HRLModelTestCommand.hpp>
#include <rhbm_gem/core/command/CommandOptionBinding.hpp>
#include "workflow/HRLModelTestWorkflowInternal.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <string_view>

namespace rhbm_gem {

namespace {
constexpr std::string_view kTesterFlags{ "-t,--tester" };
constexpr std::string_view kTesterOption{ "--tester" };
constexpr std::string_view kFitMinOption{ "--fit-min" };
constexpr std::string_view kFitMaxOption{ "--fit-max" };
constexpr std::string_view kFitRangeIssue{ "--fit-range" };
constexpr std::string_view kAlphaROption{ "--alpha-r" };
constexpr std::string_view kAlphaGOption{ "--alpha-g" };
} // namespace

HRLModelTestCommand::HRLModelTestCommand(const DataIoServices & data_io_services) :
    CommandWithProfileOptions<HRLModelTestCommandOptions, CommandId::ModelTest>{ data_io_services }
{
}

void HRLModelTestCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    command_cli::AddEnumOption<TesterType>(
        cmd, kTesterFlags,
        [&](TesterType value) { SetTesterChoice(value); },
        "Tester option",
        TesterType::DATA_OUTLIER);
    command_cli::AddScalarOption<double>(
        cmd, kFitMinOption,
        [&](double value) { SetFitRangeMinimum(value); },
        "Minimum fitting range",
        m_options.fit_range_min);
    command_cli::AddScalarOption<double>(
        cmd, kFitMaxOption,
        [&](double value) { SetFitRangeMaximum(value); },
        "Maximum fitting range",
        m_options.fit_range_max);
    command_cli::AddScalarOption<double>(
        cmd, kAlphaROption,
        [&](double value) { SetAlphaR(value); },
        "Alpha value for R",
        m_options.alpha_r);
    command_cli::AddScalarOption<double>(
        cmd, kAlphaGOption,
        [&](double value) { SetAlphaG(value); },
        "Alpha value for G",
        m_options.alpha_g);
}

bool HRLModelTestCommand::ExecuteImpl()
{
    return detail::RunHRLModelTestWorkflow(detail::HRLModelTestWorkflowContext{
        m_options,
        [this](std::string_view stem)
        {
            return BuildOutputPath(stem, "");
        },
    });
}

void HRLModelTestCommand::ValidateOptions()
{
    ResetPrepareIssues(kFitRangeIssue);
    if (m_options.fit_range_min > m_options.fit_range_max)
    {
        AddValidationError(
            kFitRangeIssue,
            "Expected --fit-min <= --fit-max.");
    }
}

void HRLModelTestCommand::SetTesterChoice(TesterType value)
{
    SetValidatedEnumOption(
        m_options.tester_choice,
        value,
        kTesterOption,
        TesterType::BENCHMARK,
        "Tester choice");
}

void HRLModelTestCommand::SetFitRangeMinimum(double value)
{
    SetFiniteNonNegativeScalarOption(
        m_options.fit_range_min,
        value,
        kFitMinOption,
        0.0,
        "Minimum fitting range must be a finite non-negative value.");
}

void HRLModelTestCommand::SetFitRangeMaximum(double value)
{
    SetFiniteNonNegativeScalarOption(
        m_options.fit_range_max,
        value,
        kFitMaxOption,
        1.0,
        "Maximum fitting range must be a finite non-negative value.");
}

void HRLModelTestCommand::SetAlphaR(double value)
{
    SetFinitePositiveScalarOption(
        m_options.alpha_r,
        value,
        kAlphaROption,
        0.1,
        "Alpha-R must be a finite positive value.");
}

void HRLModelTestCommand::SetAlphaG(double value)
{
    SetFinitePositiveScalarOption(
        m_options.alpha_g,
        value,
        kAlphaGOption,
        0.2,
        "Alpha-G must be a finite positive value.");
}

} // namespace rhbm_gem
