#include "CommandBase.hpp"
#include "FilePathHelper.hpp"

#include <CLI/CLI.hpp>
#include <algorithm>
#include <system_error>

namespace rhbm_gem {

void CommandBase::RegisterCLIOptions(CLI::App * command)
{
    RegisterCLIOptionsBasic(command);
    RegisterCLIOptionsExtend(command);
}

void CommandBase::RegisterCLIOptionsBasic(CLI::App * command)
{
    auto & options{ GetOptions() };
    const auto common_options{ GetCommonOptionMask() };

    if (HasCommonOption(common_options, CommonOption::Threading))
    {
        command->add_option_function<int>("-j,--jobs",
            [&](int value) { SetThreadSize(value); },
            "Number of threads")->default_val(options.thread_size);
    }
    if (HasCommonOption(common_options, CommonOption::Verbose))
    {
        command->add_option_function<int>("-v,--verbose",
            [&](int value) { SetVerboseLevel(value); },
            "Verbose level")->default_val(options.verbose_level);
    }
    if (HasCommonOption(common_options, CommonOption::Database))
    {
        command->add_option_function<std::string>("-d,--database",
            [&](const std::string & value) { SetDatabasePath(value); },
            "Database file path")->default_val(options.database_path.string());
    }
    if (HasCommonOption(common_options, CommonOption::OutputFolder))
    {
        command->add_option_function<std::string>("-o,--folder",
            [&](const std::string & value) { SetFolderPath(value); },
            "folder path for output files")->default_val(options.folder_path.string());
    }
}

void CommandBase::RegisterDeprecatedDatabasePathAlias(CLI::App * command)
{
    auto * hidden_group{ command->add_option_group("") };
    hidden_group->add_option_function<std::string>("--database",
        [&](const std::string & value)
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

    ClearValidationIssuesForOption("--database (deprecated)");
    if (m_deprecated_database_option_used &&
        !HasCommonOption(GetCommonOptionMask(), CommonOption::Database))
    {
        AddValidationWarning("--database (deprecated)",
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

void CommandBase::SetThreadSize(int value)
{
    auto & options{ GetOptions() };
    if (value < 1)
    {
        Logger::Log(LogLevel::Warning,
            "Thread size must be positive. Using 1 instead.");
        options.thread_size = 1;
    }
    else
    {
        options.thread_size = value;
    }
}

void CommandBase::SetVerboseLevel(int value)
{
    auto & options{ GetOptions() };
    if (value < static_cast<int>(LogLevel::Error) || value > static_cast<int>(LogLevel::Debug))
    {
        Logger::Log(LogLevel::Warning,
            "Invalid verbose level: " + std::to_string(value) +
            ", using default level 3 [Info]");
        options.verbose_level = static_cast<int>(LogLevel::Info);
    }
    else
    {
        options.verbose_level = value;
    }
}

void CommandBase::SetDatabasePath(const std::filesystem::path & path)
{
    auto & options{ GetOptions() };
    options.database_path = path.empty() ? std::filesystem::path{} : path.lexically_normal();
}

void CommandBase::SetFolderPath(const std::filesystem::path & path)
{
    auto & options{ GetOptions() };
    options.folder_path = std::filesystem::path(FilePathHelper::EnsureTrailingSlash(path));
}

bool CommandBase::HasValidationErrors() const
{
    return std::any_of(
        m_validation_issues.begin(),
        m_validation_issues.end(),
        [](const ValidationIssue & issue) { return issue.level == LogLevel::Error; });
}

void CommandBase::ReportValidationIssues() const
{
    for (const auto & issue : m_validation_issues)
    {
        Logger::Log(issue.level, "Option " + issue.option_name + ": " + issue.message);
    }
}

void CommandBase::AddValidationError(const std::string & option_name, const std::string & message)
{
    AddValidationIssue(option_name, LogLevel::Error, message);
}

void CommandBase::AddValidationWarning(const std::string & option_name, const std::string & message)
{
    AddValidationIssue(option_name, LogLevel::Warning, message);
}

void CommandBase::ClearValidationIssuesForOption(std::string_view option_name)
{
    m_validation_issues.erase(
        std::remove_if(
            m_validation_issues.begin(),
            m_validation_issues.end(),
            [option_name](const ValidationIssue & issue)
            {
                return issue.option_name == option_name;
            }),
        m_validation_issues.end());
}

void CommandBase::ValidateRequiredExistingPath(
    const std::filesystem::path & path,
    std::string_view option_name,
    std::string_view label)
{
    ClearValidationIssuesForOption(option_name);
    if (path.empty())
    {
        AddValidationError(std::string(option_name), std::string(label) + " path is required.");
        return;
    }

    std::error_code error_code;
    if (!std::filesystem::exists(path, error_code))
    {
        AddValidationError(
            std::string(option_name),
            std::string(label) + " does not exist: " + path.string());
    }
}

void CommandBase::ValidateOptionalExistingPath(
    const std::filesystem::path & path,
    std::string_view option_name,
    std::string_view label)
{
    ClearValidationIssuesForOption(option_name);
    if (path.empty()) return;

    std::error_code error_code;
    if (!std::filesystem::exists(path, error_code))
    {
        AddValidationError(
            std::string(option_name),
            std::string(label) + " does not exist: " + path.string());
    }
}

void CommandBase::AddValidationIssue(
    const std::string & option_name,
    LogLevel level,
    const std::string & message)
{
    ClearValidationIssuesForOption(option_name);
    m_validation_issues.push_back(ValidationIssue{ option_name, level, message });
}

bool CommandBase::EnsurePreparedForExecution()
{
    if (!m_is_prepared_for_execution && !PrepareForExecution())
    {
        return false;
    }

    m_is_prepared_for_execution = false;
    return true;
}

void CommandBase::PrepareFilesystemTargets()
{
    const auto common_options{ GetCommonOptionMask() };

    if (HasCommonOption(common_options, CommonOption::Database))
    {
        ClearValidationIssuesForOption("--database");
        const auto parent_path{ GetOptions().database_path.parent_path() };
        if (!parent_path.empty() && !std::filesystem::exists(parent_path))
        {
            std::error_code error_code;
            std::filesystem::create_directories(parent_path, error_code);
            if (error_code)
            {
                AddValidationError("--database",
                    "Failed to create parent directory '" + parent_path.string()
                    + "': " + error_code.message());
            }
        }
    }

    if (HasCommonOption(common_options, CommonOption::OutputFolder))
    {
        ClearValidationIssuesForOption("--folder");
        const auto & folder_path{ GetOptions().folder_path };
        if (!folder_path.empty() && !std::filesystem::exists(folder_path))
        {
            std::error_code error_code;
            std::filesystem::create_directories(folder_path, error_code);
            if (error_code)
            {
                AddValidationError("--folder",
                    "Failed to create output directory '" + folder_path.string()
                    + "': " + error_code.message());
            }
        }
    }
}

} // namespace rhbm_gem
