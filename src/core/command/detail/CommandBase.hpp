#pragma once

#include <algorithm>
#include <filesystem>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <vector>

#include <rhbm_gem/core/command/CommandTypes.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include "CommandEnumCatalog.hpp"

namespace rhbm_gem {

template <typename Request>
class CommandBase
{
    static_assert(std::is_base_of_v<CommandRequestBase, Request>,
        "CommandBase<Request> requires Request to derive from CommandRequestBase.");

    struct PendingIssue
    {
        CommandDiagnostic diagnostic;
        LogLevel level;
        bool auto_corrected{ false };
    };

    std::vector<PendingIssue> m_issues;

public:
    virtual ~CommandBase() = default;

    CommandResult ExecuteRequest(const Request & request)
    {
        Request command_request{ request };
        m_issues.clear();
        NormalizeCommonRequestOptions(command_request);
        NormalizeAndValidateRequest(command_request);
        CommandResult result;

        Logger::SetLogLevel(command_request.verbosity);
        ValidatePreparedRequest(command_request);
        if (HasValidationErrors())
        {
            ReportPendingIssues();
            result.succeeded = false;
            result.issues = BuildDiagnostics();
            return result;
        }

        result.succeeded = RunFilesystemPreflight(command_request);
        if (result.succeeded)
        {
            result.succeeded = ExecuteImpl(command_request);
            if (!result.succeeded)
            {
                Logger::Log(LogLevel::Error, "Command execution failed. Aborting command execution.");
            }
        }
        ReportPendingIssues();
        result.issues = BuildDiagnostics();
        return result;
    }

protected:
    CommandBase() = default;
    virtual void NormalizeAndValidateRequest(Request &) {}
    virtual void ValidatePreparedRequest(const Request &) {}
    virtual bool ExecuteImpl(const Request &) = 0;
    void AddParseError(std::string_view option_name, const std::string & message)
    {
        AddPendingIssue(option_name, LogLevel::Error, message, false);
    }
    void AddParseNormalizationWarning(std::string_view option_name, const std::string & message)
    {
        AddPendingIssue(option_name, LogLevel::Warning, message, true);
    }
    void ValidateRequiredPath(
        const std::filesystem::path & path,
        std::string_view option_name,
        std::string_view label)
    {
        if (path.empty())
        {
            AddPendingIssue(option_name, LogLevel::Error,
                std::string(label) + " path is required.", false);
            return;
        }

        std::error_code error_code;
        if (!std::filesystem::exists(path, error_code))
        {
            AddPendingIssue(option_name, LogLevel::Error,
                std::string(label) + " does not exist: " + path.string(), false);
        }
    }
    void ValidateOptionalPath(
        const std::filesystem::path & path,
        std::string_view option_name,
        std::string_view label)
    {
        if (path.empty()) return;

        std::error_code error_code;
        if (!std::filesystem::exists(path, error_code))
        {
            AddPendingIssue(option_name, LogLevel::Error,
                std::string(label) + " does not exist: " + path.string(), false);
        }
    }
    template <typename Container>
    void RequireNonEmptyList(
        const Container & field,
        std::string_view option_name,
        std::string_view label)
    {
        if (!field.empty()) return;
        AddPendingIssue(option_name, LogLevel::Error,
            std::string(label) + " cannot be empty.", false);
    }
    void RequirePrepareCondition(
        bool condition,
        std::string_view option_name,
        const std::string & message)
    {
        if (condition) return;
        AddPendingIssue(option_name, LogLevel::Error, message, false);
    }
    template <typename FieldType>
    void ValidatePositiveScalar(
        FieldType & field,
        std::string_view option_name,
        FieldType fallback_value,
        LogLevel issue_level,
        std::string_view label)
    {
        std::string message{
            std::string(label) + " must be positive. Using "
            + string_helper::ToStringWithPrecision(fallback_value) + " instead."
        };
        ValidateScalar(field, option_name,
            [](const auto candidate) { return numeric_validation::IsPositive(candidate); },
            fallback_value, issue_level, message);
    }
    template <typename FieldType>
    void ValidateFinitePositiveScalar(
        FieldType & field,
        std::string_view option_name,
        FieldType fallback_value,
        LogLevel issue_level,
        std::string_view label)
    {
        std::string message{
            std::string(label) + " must be a finite positive value. Using "
            + string_helper::ToStringWithPrecision(fallback_value) + " instead."
        };
        ValidateScalar(field, option_name,
            [](const auto candidate) { return numeric_validation::IsFinitePositive(candidate); },
            fallback_value, issue_level, message);
    }
    template <typename FieldType>
    void ValidateFiniteNonNegativeScalar(
        FieldType & field,
        std::string_view option_name,
        FieldType fallback_value,
        LogLevel issue_level,
        std::string_view label)
    {
        std::string message{
            std::string(label) + " must be a finite non-negative value. Using "
            + string_helper::ToStringWithPrecision(fallback_value) + " instead."
        };
        ValidateScalar(field, option_name,
            [](const auto candidate) { return numeric_validation::IsFiniteNonNegative(candidate); },
            fallback_value, issue_level, message);
    }
    template <typename FieldType, typename LowerType, typename UpperType>
    void ValidateFiniteExclusiveInclusiveRangeScalar(
        FieldType & field,
        std::string_view option_name,
        LowerType lower,
        UpperType upper,
        FieldType fallback_value,
        LogLevel issue_level,
        std::string_view label)
    {
        std::string message{
            std::string(label) + " must be a finite value within ("
                + string_helper::ToStringWithPrecision(lower) + ", "
                + string_helper::ToStringWithPrecision(upper) + "]. Using "
                + string_helper::ToStringWithPrecision(fallback_value) + " instead."
        };
        ValidateScalar(field, option_name,
            [lower, upper](const auto candidate)
            {
                return numeric_validation::IsFiniteExclusiveInclusiveRange(candidate, lower, upper);
            },
            fallback_value, issue_level, message);
    }
    template <typename FieldType>
    void ValidateEnum(
        FieldType & field,
        std::string_view option_name,
        FieldType fallback_value,
        std::string_view label)
    {
        using UnderlyingType = std::underlying_type_t<FieldType>;
        const auto raw_numeric{ static_cast<UnderlyingType>(field) };
        for (const auto & option : command_internal::CommandEnumTraits<FieldType>::kOptions)
        {
            if (static_cast<UnderlyingType>(option.value) == raw_numeric)
            {
                return;
            }
        }
        field = fallback_value;
        AddPendingIssue(option_name, LogLevel::Error,
            std::string(label) + " must be one of the supported values. Received: "
                + std::to_string(static_cast<long long>(raw_numeric)), false);
    }

private:
    std::vector<CommandDiagnostic> BuildDiagnostics() const
    {
        std::vector<CommandDiagnostic> diagnostics;
        diagnostics.reserve(m_issues.size());
        for (const auto & issue : m_issues)
        {
            diagnostics.push_back(issue.diagnostic);
        }
        return diagnostics;
    }
    bool HasValidationErrors() const
    {
        return std::any_of(
            m_issues.begin(),
            m_issues.end(),
            [](const PendingIssue & issue)
            {
                return issue.level == LogLevel::Error;
            });
    }
    void NormalizeCommonRequestOptions(Request & request)
    {
        ValidatePositiveScalar(request.job_count, "--jobs", 1, LogLevel::Warning, "Thread size");
        ValidateScalar(request.verbosity, "--verbose",
            [](int candidate)
            {
                return candidate >= static_cast<int>(LogLevel::Error)
                    && candidate <= static_cast<int>(LogLevel::Debug);
            },
            static_cast<int>(LogLevel::Info), LogLevel::Warning,
            "Invalid verbose level: " + std::to_string(request.verbosity)
                + ", using default level 3 [Info]");
        request.output_dir =
            std::filesystem::path(path_helper::EnsureTrailingSlash(request.output_dir));
    }
    template <typename FieldType, typename Predicate>
    void ValidateScalar(
        FieldType & field,
        std::string_view option_name,
        Predicate is_valid,
        FieldType fallback_value,
        LogLevel issue_level,
        const std::string & message)
    {
        if (is_valid(field)) return;
        field = fallback_value;
        if (issue_level == LogLevel::Warning)
        {
            AddPendingIssue(option_name, LogLevel::Warning, message, true);
            return;
        }
        AddPendingIssue(option_name, issue_level, message, false);
    }
    bool RunFilesystemPreflight(const Request & request)
    {
        const auto & folder_path{ request.output_dir };
        if (!folder_path.empty() && !std::filesystem::exists(folder_path))
        {
            std::error_code error_code;
            std::filesystem::create_directories(folder_path, error_code);
            if (error_code)
            {
                AddPendingIssue("--folder", LogLevel::Error,
                    "Failed to create output directory '" + folder_path.string()
                    + "': " + error_code.message(), false);
            }
        }
        return !HasValidationErrors();
    }
    void ReportPendingIssues() const
    {
        for (const auto & issue : m_issues)
        {
            std::string prefix{ "[validation" };
            if (issue.auto_corrected) prefix += "; auto-corrected";
            prefix += "]";
            Logger::Log(issue.level, prefix + " Option " + issue.diagnostic.option_name + ": "
                + issue.diagnostic.message);
        }
    }
    void AddPendingIssue(
        std::string_view option_name,
        LogLevel level,
        const std::string & message,
        bool auto_corrected)
    {
        m_issues.push_back(PendingIssue{ CommandDiagnostic{ std::string(option_name), message },
            level, auto_corrected });
    }

};

} // namespace rhbm_gem
