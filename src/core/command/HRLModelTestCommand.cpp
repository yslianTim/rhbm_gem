#include "HRLModelTestCommand.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <string_view>

namespace rhbm_gem {

namespace {
constexpr std::string_view kTesterOption{ "--tester" };
constexpr std::string_view kFitMinOption{ "--fit-min" };
constexpr std::string_view kFitMaxOption{ "--fit-max" };
constexpr std::string_view kFitRangeIssue{ "--fit-range" };
constexpr std::string_view kAlphaROption{ "--alpha-r" };
constexpr std::string_view kAlphaGOption{ "--alpha-g" };
} // namespace

HRLModelTestCommand::HRLModelTestCommand() :
    CommandWithProfileOptions<HRLModelTestCommandOptions, CommandId::ModelTest>{}
{
}

void HRLModelTestCommand::ApplyRequest(const HRLModelTestRequest & request)
{
    ApplyCommonRequest(request.common);
    SetTesterChoice(request.tester_choice);
    SetFitRangeMinimum(request.fit_range_min);
    SetFitRangeMaximum(request.fit_range_max);
    SetAlphaR(request.alpha_r);
    SetAlphaG(request.alpha_g);
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
