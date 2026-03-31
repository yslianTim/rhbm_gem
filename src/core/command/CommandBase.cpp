#include "internal/command/CommandBase.hpp"
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>

#include <algorithm>
#include <array>
#include <string>
#include <system_error>

namespace rhbm_gem {

namespace {

constexpr std::array<std::string_view, 2> kValidationPhaseLabels{
    "parse",
    "prepare"
};
constexpr std::string_view kJobsOption{ "--jobs" };
constexpr std::string_view kVerboseOption{ "--verbose" };
constexpr std::string_view kFolderOption{ "--folder" };

std::string BuildIssuePrefix(const ValidationIssue & issue)
{
    const auto phase_index{ static_cast<std::size_t>(issue.phase) };
    const auto phase_label{
        phase_index < kValidationPhaseLabels.size()
            ? kValidationPhaseLabels[phase_index]
            : std::string_view{"unknown"}
    };
    std::string prefix{ "[" + std::string(phase_label) };
    if (issue.auto_corrected)
    {
        prefix += "; auto-corrected";
    }
    prefix += "]";
    return prefix;
}

} // namespace

bool CommandBase::Run()
{
    BeginPreparationPass();
    if (!RunValidationPass())
    {
        return false;
    }

    if (!RunFilesystemPreflight())
    {
        return false;
    }

    m_was_prepared = true;
    const bool executed{ ExecuteImpl() };
    if (!executed)
    {
        Logger::Log(LogLevel::Error,
            "Command execution failed. Aborting command execution.");
    }
    return executed;
}

void CommandBase::InvalidatePreparedState()
{
    m_was_prepared = false;
    ClearValidationIssues(ValidationPhase::Prepare);
}

void CommandBase::CoerceBaseRequest(CommandRequestBase & request)
{
    InvalidatePreparedState();
    const auto raw_verbose_level{ request.verbosity };
    CoerceScalar(
        request.job_count,
        kJobsOption,
        [](int candidate) { return candidate >= 1; },
        1,
        LogLevel::Warning,
        "Thread size must be positive. Using 1 instead.");
    CoerceScalar(
        request.verbosity,
        kVerboseOption,
        [](int candidate)
        {
            return candidate >= static_cast<int>(LogLevel::Error)
                && candidate <= static_cast<int>(LogLevel::Debug);
        },
        static_cast<int>(LogLevel::Info),
        LogLevel::Warning,
        "Invalid verbose level: " + std::to_string(raw_verbose_level)
            + ", using default level 3 [Info]");
    request.output_dir =
        std::filesystem::path(FilePathHelper::EnsureTrailingSlash(request.output_dir));
}

void CommandBase::ReportValidationIssues() const
{
    for (const auto & issue : m_validation_issues)
    {
        Logger::Log(
            issue.level,
            BuildIssuePrefix(issue) + " Option " + issue.option_name + ": " + issue.message);
    }
}

bool CommandBase::HasValidationErrors() const
{
    return std::any_of(
        m_validation_issues.begin(),
        m_validation_issues.end(),
        [](const ValidationIssue & issue)
        {
            return issue.level == LogLevel::Error;
        });
}

void CommandBase::AddValidationError(
    std::string_view option_name,
    const std::string & message,
    ValidationPhase phase)
{
    AddValidationIssue(option_name, phase, LogLevel::Error, message, false);
}

void CommandBase::AddNormalizationWarning(
    std::string_view option_name,
    const std::string & message)
{
    AddValidationIssue(option_name, ValidationPhase::Parse, LogLevel::Warning, message, true);
}

void CommandBase::ClearParseIssues(std::string_view option_name)
{
    ClearValidationIssues(option_name, ValidationPhase::Parse);
}

void CommandBase::ClearPrepareIssues(std::string_view option_name)
{
    ClearValidationIssues(option_name, ValidationPhase::Prepare);
}

void CommandBase::BeginValidationMutation(ValidationPhase phase)
{
    if (phase == ValidationPhase::Parse)
    {
        InvalidatePreparedState();
        return;
    }

    m_was_prepared = false;
}

void CommandBase::RequireNonEmptyText(
    std::string_view field,
    std::string_view option_name,
    std::string_view label,
    ValidationPhase phase)
{
    BeginValidationMutation(phase);
    ClearValidationIssues(option_name, phase);
    if (!field.empty())
    {
        return;
    }

    AddValidationError(
        option_name,
        std::string(label) + " cannot be empty.",
        phase);
}

void CommandBase::RequireCondition(
    bool condition,
    std::string_view option_name,
    const std::string & message,
    ValidationPhase phase)
{
    BeginValidationMutation(phase);
    ClearValidationIssues(option_name, phase);
    if (condition)
    {
        return;
    }

    AddValidationError(option_name, message, phase);
}

void CommandBase::ClearValidationIssues(
    std::string_view option_name,
    std::optional<ValidationPhase> phase)
{
    m_validation_issues.erase(
        std::remove_if(
            m_validation_issues.begin(),
            m_validation_issues.end(),
            [option_name, phase](const ValidationIssue & issue)
            {
                return issue.option_name == option_name
                    && (!phase.has_value() || issue.phase == phase.value());
            }),
        m_validation_issues.end());
}

void CommandBase::ClearValidationIssues(std::optional<ValidationPhase> phase)
{
    m_validation_issues.erase(
        std::remove_if(
            m_validation_issues.begin(),
            m_validation_issues.end(),
            [phase](const ValidationIssue & issue)
            {
                return !phase.has_value() || issue.phase == phase.value();
            }),
        m_validation_issues.end());
}

void CommandBase::ValidateRequiredExistingPath(
    const std::filesystem::path & path,
    std::string_view option_name,
    std::string_view label)
{
    ClearValidationIssues(option_name, ValidationPhase::Parse);
    if (path.empty())
    {
        AddValidationError(
            std::string(option_name),
            std::string(label) + " path is required.",
            ValidationPhase::Parse);
        return;
    }

    std::error_code error_code;
    if (!std::filesystem::exists(path, error_code))
    {
        AddValidationError(
            std::string(option_name),
            std::string(label) + " does not exist: " + path.string(),
            ValidationPhase::Parse);
    }
}

void CommandBase::ValidateOptionalExistingPath(
    const std::filesystem::path & path,
    std::string_view option_name,
    std::string_view label)
{
    ClearValidationIssues(option_name, ValidationPhase::Parse);
    if (path.empty()) return;

    std::error_code error_code;
    if (!std::filesystem::exists(path, error_code))
    {
        AddValidationError(
            std::string(option_name),
            std::string(label) + " does not exist: " + path.string(),
            ValidationPhase::Parse);
    }
}

void CommandBase::ValidateRequiredPath(
    std::filesystem::path & field,
    std::string_view option_name,
    std::string_view label)
{
    InvalidatePreparedState();
    ValidateRequiredExistingPath(field, option_name, label);
}

void CommandBase::ValidateOptionalPath(
    std::filesystem::path & field,
    std::string_view option_name,
    std::string_view label)
{
    InvalidatePreparedState();
    ValidateOptionalExistingPath(field, option_name, label);
}

void CommandBase::AddValidationIssue(
    std::string_view option_name,
    ValidationPhase phase,
    LogLevel level,
    const std::string & message,
    bool auto_corrected)
{
    m_validation_issues.push_back(ValidationIssue{
        std::string(option_name),
        phase,
        level,
        message,
        auto_corrected
    });
}

std::filesystem::path CommandBase::BuildOutputPath(
    std::string_view stem,
    std::string_view extension) const
{
    std::string normalized_extension{ extension };
    if (!normalized_extension.empty() && normalized_extension.front() != '.')
    {
        normalized_extension.insert(normalized_extension.begin(), '.');
    }
    return OutputFolder() / (std::string(stem) + normalized_extension);
}

void CommandBase::BeginPreparationPass()
{
    Logger::SetLogLevel(VerboseLevel());
    ResetRuntimeState();
    m_data_manager.ClearDataObjects();
    InvalidatePreparedState();
}

bool CommandBase::RunValidationPass()
{
    ValidateOptions();
    if (HasValidationErrors())
    {
        ReportValidationIssues();
        return false;
    }

    return true;
}

bool CommandBase::RunFilesystemPreflight()
{
    ClearValidationIssues(kFolderOption, ValidationPhase::Prepare);
    const auto & folder_path{ OutputFolder() };
    if (!folder_path.empty() && !std::filesystem::exists(folder_path))
    {
        std::error_code error_code;
        std::filesystem::create_directories(folder_path, error_code);
        if (error_code)
        {
            AddValidationError(
                kFolderOption,
                "Failed to create output directory '" + folder_path.string()
                    + "': " + error_code.message());
        }
    }

    ReportValidationIssues();
    return !HasValidationErrors();
}

} // namespace rhbm_gem
