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

CommandBase::CommandBase(CommonOptionMask common_options) :
    m_data_manager{},
    m_common_options{ common_options }
{
}

CommandBase::CommandBase(CommonOptionProfile profile) :
    CommandBase{ CommonOptionMaskForProfile(profile) }
{
}

bool CommandBase::Execute()
{
    const bool was_prepared{ m_is_prepared_for_execution };
    if (!was_prepared && !PrepareForExecution())
    {
        return false;
    }

    if (was_prepared)
    {
        Logger::SetLogLevel(VerboseLevel());
    }

    const bool executed{ ExecuteImpl() };
    m_is_prepared_for_execution = false;
    if (!executed)
    {
        Logger::Log(LogLevel::Error,
            "Command execution failed. Aborting command execution.");
    }
    return executed;
}

bool CommandBase::PrepareForExecution()
{
    BeginPreparationPass();
    if (!RunValidationPass())
    {
        return false;
    }

    return RunFilesystemPreflight();
}

void CommandBase::InvalidatePreparedState()
{
    m_is_prepared_for_execution = false;
    ClearValidationIssues(ValidationPhase::Prepare);
}

void CommandBase::SetThreadSize(int value)
{
    SetNormalizedScalarOption(
        m_common_request.thread_size,
        value,
        kJobsOption,
        [](int candidate) { return candidate >= 1; },
        1,
        "Thread size must be positive. Using 1 instead.");
}

void CommandBase::SetVerboseLevel(int value)
{
    SetNormalizedScalarOption(
        m_common_request.verbose_level,
        value,
        kVerboseOption,
        [](int candidate)
        {
            return candidate >= static_cast<int>(LogLevel::Error)
                && candidate <= static_cast<int>(LogLevel::Debug);
        },
        static_cast<int>(LogLevel::Info),
        "Invalid verbose level: " + std::to_string(value) + ", using default level 3 [Info]");
}

void CommandBase::SetDatabasePath(const std::filesystem::path & path)
{
    MutateOptions([&]()
    {
        m_common_request.database_path = path.empty()
            ? std::filesystem::path{}
            : path.lexically_normal();
    });
}

void CommandBase::SetFolderPath(const std::filesystem::path & path)
{
    MutateOptions([&]()
    {
        m_common_request.folder_path =
            std::filesystem::path(FilePathHelper::EnsureTrailingSlash(path));
    });
}

void CommandBase::ApplyCommonRequest(const CommonCommandRequest & request)
{
    SetThreadSize(request.thread_size);
    SetVerboseLevel(request.verbose_level);

    const auto common_options{ GetCommonOptionsMask() };
    if (HasCommonOption(common_options, CommonOption::Database))
    {
        SetDatabasePath(request.database_path);
    }
    if (HasCommonOption(common_options, CommonOption::OutputFolder))
    {
        SetFolderPath(request.folder_path);
    }
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

void CommandBase::ResetParseIssues(std::string_view option_name)
{
    ClearValidationIssues(option_name, ValidationPhase::Parse);
}

void CommandBase::ResetPrepareIssues(std::string_view option_name)
{
    ClearValidationIssues(option_name, ValidationPhase::Prepare);
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

void CommandBase::SetRequiredExistingPathOption(
    std::filesystem::path & field,
    const std::filesystem::path & value,
    std::string_view option_name,
    std::string_view label)
{
    MutateOptions([&]()
    {
        field = value;
        ValidateRequiredExistingPath(field, option_name, label);
    });
}

void CommandBase::SetOptionalExistingPathOption(
    std::filesystem::path & field,
    const std::filesystem::path & value,
    std::string_view option_name,
    std::string_view label)
{
    MutateOptions([&]()
    {
        field = value;
        ValidateOptionalExistingPath(field, option_name, label);
    });
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
        m_is_prepared_for_execution = false;
        return false;
    }

    return true;
}

bool CommandBase::RunFilesystemPreflight()
{
    const auto common_options{ GetCommonOptionsMask() };

    if (HasCommonOption(common_options, CommonOption::OutputFolder))
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
    }

    ReportValidationIssues();
    m_is_prepared_for_execution = !HasValidationErrors();
    return m_is_prepared_for_execution;
}

} // namespace rhbm_gem
