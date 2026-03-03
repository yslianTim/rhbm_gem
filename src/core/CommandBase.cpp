#include "CommandBase.hpp"
#include "CommandOptionBinding.hpp"
#include "FilePathHelper.hpp"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <array>
#include <string>
#include <system_error>

namespace rhbm_gem {

namespace {

constexpr std::array<std::string_view, 3> kValidationPhaseLabels{
    "parse",
    "prepare",
    "runtime"
};

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

bool CommandBase::Execute()
{
    Logger::SetLogLevel(GetOptions().verbose_level);
    if (!m_is_prepared_for_execution && !PrepareForExecution())
    {
        return false;
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

const CommandDescriptor & CommandBase::GetCommandDescriptor() const
{
    return FindCommandDescriptor(GetCommandId());
}

void CommandBase::RegisterCLIOptions(CLI::App * command)
{
    RegisterCLIOptionsBasic(command);
    RegisterDeprecatedCommonAliases(command);
    RegisterCLIOptionsExtend(command);
}

void CommandBase::RegisterCLIOptionsBasic(CLI::App * command)
{
    auto & options{ GetOptions() };
    const auto common_options{ GetCommonOptionMask() };

    if (HasCommonOption(common_options, CommonOption::Threading))
    {
        command_cli::AddScalarOption<int>(
            command,
            "-j,--jobs",
            [&](int value) { SetThreadSize(value); },
            "Number of threads",
            options.thread_size);
    }
    if (HasCommonOption(common_options, CommonOption::Verbose))
    {
        command_cli::AddScalarOption<int>(
            command,
            "-v,--verbose",
            [&](int value) { SetVerboseLevel(value); },
            "Verbose level",
            options.verbose_level);
    }
    if (HasCommonOption(common_options, CommonOption::Database))
    {
        command_cli::AddPathOption(
            command,
            "-d,--database",
            [&](const std::filesystem::path & value) { SetDatabasePath(value); },
            "Database file path",
            options.database_path);
    }
    if (HasCommonOption(common_options, CommonOption::OutputFolder))
    {
        command_cli::AddPathOption(
            command,
            "-o,--folder",
            [&](const std::filesystem::path & value) { SetFolderPath(value); },
            "folder path for output files",
            options.folder_path);
    }
}

void CommandBase::RegisterDeprecatedCommonAliases(CLI::App * command)
{
    const auto deprecated_options{ GetCommandDescriptor().surface.deprecated_hidden_options };
    if (!HasCommonOption(deprecated_options, CommonOption::Database))
    {
        return;
    }

    auto * hidden_group{ command->add_option_group("") };
    command_cli::AddPathOption(
        hidden_group,
        "--database",
        [&](const std::filesystem::path & value)
        {
            m_deprecated_database_option_used = true;
            SetDatabasePath(value);
        },
        "Deprecated hidden database compatibility alias");
}

bool CommandBase::PrepareForExecution()
{
    Logger::SetLogLevel(GetOptions().verbose_level);
    ResetRuntimeState();
    m_data_manager.ClearDataObjects();
    InvalidatePreparedState();

    ClearValidationIssues("--database (deprecated)", ValidationPhase::Prepare);
    if (m_deprecated_database_option_used &&
        !HasCommonOption(GetCommonOptionMask(), CommonOption::Database))
    {
        AddValidationWarning(
            "--database (deprecated)",
            "This command ignores the deprecated hidden '--database' compatibility alias.");
    }

    ValidateOptions();
    if (HasValidationErrors())
    {
        ReportValidationIssues();
        m_is_prepared_for_execution = false;
        return false;
    }

    PrepareFilesystemTargets();
    ReportValidationIssues();
    m_is_prepared_for_execution = !HasValidationErrors();
    return m_is_prepared_for_execution;
}

void CommandBase::InvalidatePreparedState()
{
    m_is_prepared_for_execution = false;
    ClearValidationIssues(ValidationPhase::Prepare);
    ClearValidationIssues(ValidationPhase::Runtime);
}

void CommandBase::SetThreadSize(int value)
{
    auto & options{ GetOptions() };
    SetNormalizedScalarOption(
        options.thread_size,
        value,
        "--jobs",
        [](int candidate) { return candidate >= 1; },
        1,
        "Thread size must be positive. Using 1 instead.");
}

void CommandBase::SetVerboseLevel(int value)
{
    auto & options{ GetOptions() };
    SetNormalizedScalarOption(
        options.verbose_level,
        value,
        "--verbose",
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
    auto & options{ GetOptions() };
    MutateOptions([&]()
    {
        options.database_path = path.empty() ? std::filesystem::path{} : path.lexically_normal();
    });
}

void CommandBase::SetFolderPath(const std::filesystem::path & path)
{
    auto & options{ GetOptions() };
    MutateOptions([&]()
    {
        options.folder_path = std::filesystem::path(FilePathHelper::EnsureTrailingSlash(path));
    });
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

bool CommandBase::HasValidationErrors(std::optional<ValidationPhase> phase) const
{
    return std::any_of(
        m_validation_issues.begin(),
        m_validation_issues.end(),
        [phase](const ValidationIssue & issue)
        {
            return issue.level == LogLevel::Error
                && (!phase.has_value() || issue.phase == phase.value());
        });
}

void CommandBase::AddValidationError(
    const std::string & option_name,
    const std::string & message,
    ValidationPhase phase)
{
    AddValidationIssue(option_name, phase, LogLevel::Error, message, false);
}

void CommandBase::AddValidationWarning(
    const std::string & option_name,
    const std::string & message,
    ValidationPhase phase)
{
    AddValidationIssue(option_name, phase, LogLevel::Warning, message, false);
}

void CommandBase::AddNormalizationWarning(
    const std::string & option_name,
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
    const std::string & option_name,
    ValidationPhase phase,
    LogLevel level,
    const std::string & message,
    bool auto_corrected)
{
    m_validation_issues.push_back(ValidationIssue{
        option_name,
        phase,
        level,
        message,
        auto_corrected
    });
}

void CommandBase::RequireDatabaseManager()
{
    m_data_manager.SetDatabaseManager(GetOptions().database_path);
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
    return GetOptions().folder_path / (std::string(stem) + normalized_extension);
}

void CommandBase::PrepareFilesystemTargets()
{
    const auto common_options{ GetCommonOptionMask() };

    if (HasCommonOption(common_options, CommonOption::Database))
    {
        ClearValidationIssues("--database", ValidationPhase::Prepare);
        const auto parent_path{ GetOptions().database_path.parent_path() };
        if (!parent_path.empty() && !std::filesystem::exists(parent_path))
        {
            std::error_code error_code;
            std::filesystem::create_directories(parent_path, error_code);
            if (error_code)
            {
                AddValidationError(
                    "--database",
                    "Failed to create parent directory '" + parent_path.string()
                        + "': " + error_code.message());
            }
        }
    }

    if (HasCommonOption(common_options, CommonOption::OutputFolder))
    {
        ClearValidationIssues("--folder", ValidationPhase::Prepare);
        const auto & folder_path{ GetOptions().folder_path };
        if (!folder_path.empty() && !std::filesystem::exists(folder_path))
        {
            std::error_code error_code;
            std::filesystem::create_directories(folder_path, error_code);
            if (error_code)
            {
                AddValidationError(
                    "--folder",
                    "Failed to create output directory '" + folder_path.string()
                        + "': " + error_code.message());
            }
        }
    }
}

} // namespace rhbm_gem
