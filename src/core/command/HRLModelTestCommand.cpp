#include <rhbm_gem/core/HRLModelTestCommand.hpp>
#include <rhbm_gem/core/CommandOptionBinding.hpp>
#include <rhbm_gem/utils/Logger.hpp>

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

HRLModelTestCommand::HRLModelTestCommand() = default;

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
    switch (m_options.tester_choice)
    {
        case TesterType::BENCHMARK:
            RunSimulationTestOnBenchMark();
            break;
        case TesterType::DATA_OUTLIER:
            RunSimulationTestOnDataOutlier();
            break;
        case TesterType::MEMBER_OUTLIER:
            RunSimulationTestOnMemberOutlier();
            break;
        case TesterType::MODEL_ALPHA_DATA:
            RunSimulationTestOnModelAlphaData();
            break;
        case TesterType::MODEL_ALPHA_MEMBER:
            RunSimulationTestOnModelAlphaMember();
            break;
        default:
            Logger::Log(LogLevel::Error,
                        "Invalid tester choice reached execution path: ["
                        + std::to_string(static_cast<int>(m_options.tester_choice)) + "]");
            return false;
    }
    return true;
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
