#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "DataObjectManager.hpp"
#include "Logger.hpp"

namespace CLI
{
    class App;
}

namespace rhbm_gem {

struct ValidationIssue
{
    std::string option_name;
    LogLevel level;
    std::string message;
};

enum class CommonOption : std::uint8_t
{
    Threading = 0,
    Verbose = 1,
    Database = 2,
    OutputFolder = 3
};

using CommonOptionMask = std::uint32_t;

constexpr CommonOptionMask ToMask(CommonOption option)
{
    return static_cast<CommonOptionMask>(1u << static_cast<std::uint8_t>(option));
}

constexpr CommonOptionMask operator|(CommonOption lhs, CommonOption rhs)
{
    return ToMask(lhs) | ToMask(rhs);
}

constexpr CommonOptionMask operator|(CommonOptionMask lhs, CommonOption rhs)
{
    return lhs | ToMask(rhs);
}

constexpr bool HasCommonOption(CommonOptionMask mask, CommonOption option)
{
    return (mask & ToMask(option)) != 0u;
}

struct CommandOptions
{
    int thread_size{ 1 };
    int verbose_level{ 3 };
    std::filesystem::path database_path{"database.sqlite"};
    std::filesystem::path folder_path{""};
};

class CommandBase
{
public:
    virtual ~CommandBase() = default;
    virtual bool Execute() = 0;
    virtual void RegisterCLIOptionsExtend(::CLI::App * command) = 0;
    virtual const CommandOptions & GetOptions() const = 0;
    virtual CommandOptions & GetOptions() = 0;
    virtual void ValidateOptions() {}
    virtual void ResetRuntimeState() {}
    virtual CommonOptionMask GetCommonOptionMask() const
    {
        return CommonOption::Threading
             | CommonOption::Verbose
             | CommonOption::Database
             | CommonOption::OutputFolder;
    }

    void RegisterCLIOptions(::CLI::App * command);
    bool PrepareForExecution();
    void SetThreadSize(int value);
    void SetVerboseLevel(int value);
    void SetDatabasePath(const std::filesystem::path & path);
    void SetFolderPath(const std::filesystem::path & path);
    bool IsValidateOptions() const { return !HasValidationErrors(); }
    bool HasValidationErrors() const;
    void ReportValidationIssues() const;
    const std::vector<ValidationIssue> & GetValidationIssues() const { return m_validation_issues; }
    DataObjectManager * GetDataManagerPtr() { return &m_data_manager; }
    const DataObjectManager * GetDataManagerPtr() const { return &m_data_manager; }

protected:
    DataObjectManager m_data_manager;
    std::vector<ValidationIssue> m_validation_issues;

    CommandBase() = default;
    void RegisterCLIOptionsBasic(::CLI::App * command);
    void RegisterDeprecatedDatabasePathAlias(::CLI::App * command);
    void AddValidationError(const std::string & option_name, const std::string & message);
    void AddValidationWarning(const std::string & option_name, const std::string & message);
    void ClearValidationIssuesForOption(std::string_view option_name);
    void ValidateRequiredExistingPath(
        const std::filesystem::path & path,
        std::string_view option_name,
        std::string_view label);
    void ValidateOptionalExistingPath(
        const std::filesystem::path & path,
        std::string_view option_name,
        std::string_view label);
    bool EnsurePreparedForExecution();

private:
    bool m_deprecated_database_option_used{ false };
    bool m_is_prepared_for_execution{ false };

    void AddValidationIssue(
        const std::string & option_name,
        LogLevel level,
        const std::string & message);
    void PrepareFilesystemTargets();
};

} // namespace rhbm_gem
