#include "CommandBase.hpp"
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

std::string BuildIssuePrefix(const ValidationIssueRecord & issue)
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
    if (!RunValidationPass()) return false;
    if (!RunFilesystemPreflight()) return false;

    m_was_prepared = true;
    const bool executed{ ExecuteImpl() };
    if (!executed)
    {
        Logger::Log(LogLevel::Error,"Command execution failed. Aborting command execution.");
    }
    return executed;
}

void CommandBase::InvalidatePreparedState()
{
    m_was_prepared = false;
    ClearValidationIssues(ValidationPhase::Prepare);
}

void CommandBase::BeginRequestApply(CommandRequestBase & request)
{
    m_was_prepared = false;
    ClearValidationIssues();
    CoerceBaseRequest(request);
}

void CommandBase::CoerceBaseRequest(CommandRequestBase & request)
{
    const auto raw_verbose_level{ request.verbosity };
    CoercePositiveScalar(
        request.job_count,
        "--jobs",
        1,
        LogLevel::Warning,
        "Thread size");
    CoerceScalar(
        request.verbosity,
        "--verbose",
        [](int candidate)
        {
            return candidate >= static_cast<int>(LogLevel::Error)
                && candidate <= static_cast<int>(LogLevel::Debug);
        },
        static_cast<int>(LogLevel::Info),
        LogLevel::Warning,
        "Invalid verbose level: " + std::to_string(raw_verbose_level)
            + ", using default level 3 [Info]");
    request.output_dir = std::filesystem::path(path_helper::EnsureTrailingSlash(request.output_dir));
}

void CommandBase::ReportValidationIssues() const
{
    for (const auto & issue : m_validation_issues)
    {
        Logger::Log(issue.level,
            BuildIssuePrefix(issue) + " Option " + issue.option_name + ": " + issue.message);
    }
}

bool CommandBase::HasValidationErrors() const
{
    return std::any_of(
        m_validation_issues.begin(),
        m_validation_issues.end(),
        [](const ValidationIssueRecord & issue)
        {
            return issue.level == LogLevel::Error;
        });
}

void CommandBase::AddParseError(
    std::string_view option_name,
    const std::string & message)
{
    AddValidationIssue(option_name, ValidationPhase::Parse, LogLevel::Error, message, false);
}

void CommandBase::AddParseNormalizationWarning(
    std::string_view option_name,
    const std::string & message)
{
    AddValidationIssue(option_name, ValidationPhase::Parse, LogLevel::Warning, message, true);
}

void CommandBase::RequirePrepareCondition(
    bool condition,
    std::string_view option_name,
    const std::string & message)
{
    m_was_prepared = false;
    if (condition) return;
    AddValidationIssue(option_name, ValidationPhase::Prepare, LogLevel::Error, message, false);
}

void CommandBase::ClearValidationIssues(std::optional<ValidationPhase> phase)
{
    m_validation_issues.erase(
        std::remove_if(
            m_validation_issues.begin(),
            m_validation_issues.end(),
            [phase](const ValidationIssueRecord & issue)
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
    if (path.empty())
    {
        AddValidationIssue(
            std::string(option_name),
            ValidationPhase::Parse,
            LogLevel::Error,
            std::string(label) + " path is required.",
            false);
        return;
    }

    std::error_code error_code;
    if (!std::filesystem::exists(path, error_code))
    {
        AddValidationIssue(
            std::string(option_name),
            ValidationPhase::Parse,
            LogLevel::Error,
            std::string(label) + " does not exist: " + path.string(),
            false);
    }
}

void CommandBase::ValidateOptionalExistingPath(
    const std::filesystem::path & path,
    std::string_view option_name,
    std::string_view label)
{
    if (path.empty()) return;

    std::error_code error_code;
    if (!std::filesystem::exists(path, error_code))
    {
        AddValidationIssue(
            std::string(option_name),
            ValidationPhase::Parse,
            LogLevel::Error,
            std::string(label) + " does not exist: " + path.string(),
            false);
    }
}

void CommandBase::ValidateRequiredPath(
    const std::filesystem::path & path,
    std::string_view option_name,
    std::string_view label)
{
    InvalidatePreparedState();
    ValidateRequiredExistingPath(path, option_name, label);
}

void CommandBase::ValidateOptionalPath(
    const std::filesystem::path & path,
    std::string_view option_name,
    std::string_view label)
{
    InvalidatePreparedState();
    ValidateOptionalExistingPath(path, option_name, label);
}

void CommandBase::AddValidationIssue(
    std::string_view option_name,
    ValidationPhase phase,
    LogLevel level,
    const std::string & message,
    bool auto_corrected)
{
    m_validation_issues.push_back(ValidationIssueRecord{
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
    InvalidatePreparedState();
}

bool CommandBase::RunValidationPass()
{
    ValidatePreparedRequest();
    if (HasValidationErrors())
    {
        ReportValidationIssues();
        return false;
    }
    return true;
}

bool CommandBase::RunFilesystemPreflight()
{
    const auto & folder_path{ OutputFolder() };
    if (!folder_path.empty() && !std::filesystem::exists(folder_path))
    {
        std::error_code error_code;
        std::filesystem::create_directories(folder_path, error_code);
        if (error_code)
        {
            AddValidationIssue(
                "--folder",
                ValidationPhase::Prepare,
                LogLevel::Error,
                "Failed to create output directory '" + folder_path.string()
                    + "': " + error_code.message(),
                false);
        }
    }

    ReportValidationIssues();
    return !HasValidationErrors();
}

} // namespace rhbm_gem
