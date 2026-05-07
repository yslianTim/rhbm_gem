#pragma once

#include <algorithm>
#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <vector>

#include <rhbm_gem/core/command/CommandSystem.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include "CommandEnumCatalog.hpp"

namespace rhbm_gem {

enum class ValidationPhase : std::uint8_t
{
    Parse = 0,
    Prepare = 1
};

struct ValidationIssueRecord
{
    std::string option_name;
    ValidationPhase phase;
    LogLevel level;
    std::string message;
    bool auto_corrected{ false };
};

template <typename Request>
class CommandBase
{
    static_assert(std::is_base_of_v<CommandRequestBase, Request>,
        "CommandBase<Request> requires Request to derive from CommandRequestBase.");

    Request m_request{};
    std::vector<ValidationIssueRecord> m_validation_issues;

public:
    using RequestType = Request;

    virtual ~CommandBase() = default;
    void ApplyRequest(const Request & request)
    {
        m_request = request;
        BeginRequestApply();
        NormalizeAndValidateRequest();
    }
    bool Run()
    {
        Logger::SetLogLevel(m_request.verbosity);
        ClearPrepareValidationIssues();
        ValidatePreparedRequest();
        if (HasValidationErrors())
        {
            ReportValidationIssues();
            return false;
        }
        if (!RunFilesystemPreflight()) return false;

        const bool executed{ ExecuteImpl() };
        if (!executed)
        {
            Logger::Log(LogLevel::Error, "Command execution failed. Aborting command execution.");
        }
        return executed;
    }
    bool HasValidationErrors() const
    {
        return std::any_of(
            m_validation_issues.begin(),
            m_validation_issues.end(),
            [](const ValidationIssueRecord & issue)
            {
                return issue.level == LogLevel::Error;
            });
    }
    const std::vector<ValidationIssueRecord> & GetValidationIssues() const { return m_validation_issues; }

protected:
    CommandBase() = default;
    const Request & RequestOptions() const { return m_request; }
    Request & MutableRequest() { return m_request; }
    virtual void NormalizeAndValidateRequest() {}
    virtual void ValidatePreparedRequest() {}
    virtual bool ExecuteImpl() = 0;
    int ThreadSize() const { return m_request.job_count; }
    const std::filesystem::path & OutputFolder() const { return m_request.output_dir; }
    void AddParseError(std::string_view option_name, const std::string & message)
    {
        AddValidationIssue(option_name, ValidationPhase::Parse, LogLevel::Error, message, false);
    }
    void AddParseNormalizationWarning(std::string_view option_name, const std::string & message)
    {
        AddValidationIssue(option_name, ValidationPhase::Parse, LogLevel::Warning, message, true);
    }
    void ValidateRequiredPath(
        const std::filesystem::path & path,
        std::string_view option_name,
        std::string_view label)
    {
        ClearPrepareValidationIssues();
        if (path.empty())
        {
            AddValidationIssue(option_name, ValidationPhase::Parse, LogLevel::Error,
                std::string(label) + " path is required.", false);
            return;
        }

        std::error_code error_code;
        if (!std::filesystem::exists(path, error_code))
        {
            AddValidationIssue(option_name, ValidationPhase::Parse, LogLevel::Error,
                std::string(label) + " does not exist: " + path.string(), false);
        }
    }
    void ValidateOptionalPath(
        const std::filesystem::path & path,
        std::string_view option_name,
        std::string_view label)
    {
        ClearPrepareValidationIssues();
        if (path.empty()) return;

        std::error_code error_code;
        if (!std::filesystem::exists(path, error_code))
        {
            AddValidationIssue(option_name, ValidationPhase::Parse, LogLevel::Error,
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
        AddValidationIssue(option_name, ValidationPhase::Parse, LogLevel::Error,
            std::string(label) + " cannot be empty.", false);
    }
    void RequirePrepareCondition(
        bool condition,
        std::string_view option_name,
        const std::string & message)
    {
        if (condition) return;
        AddValidationIssue(option_name, ValidationPhase::Prepare, LogLevel::Error, message, false);
    }
    template <typename FieldType>
    void CoercePositiveScalar(
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
        CoerceScalar(field, option_name,
            [](const auto candidate) { return numeric_validation::IsPositive(candidate); },
            fallback_value, issue_level, message);
    }
    template <typename FieldType>
    void CoerceFinitePositiveScalar(
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
        CoerceScalar(field, option_name,
            [](const auto candidate) { return numeric_validation::IsFinitePositive(candidate); },
            fallback_value, issue_level, message);
    }
    template <typename FieldType>
    void CoerceFiniteNonNegativeScalar(
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
        CoerceScalar(field, option_name,
            [](const auto candidate) { return numeric_validation::IsFiniteNonNegative(candidate); },
            fallback_value, issue_level, message);
    }
    template <typename FieldType, typename LowerType, typename UpperType>
    void CoerceFiniteExclusiveInclusiveRangeScalar(
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
        CoerceScalar(field, option_name,
            [lower, upper](const auto candidate)
            {
                return numeric_validation::IsFiniteExclusiveInclusiveRange(candidate, lower, upper);
            },
            fallback_value, issue_level, message);
    }
    template <typename FieldType>
    void CoerceEnum(
        FieldType & field,
        std::string_view option_name,
        FieldType fallback_value,
        std::string_view label)
    {
        ClearPrepareValidationIssues();
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
        AddValidationIssue(option_name, ValidationPhase::Parse, LogLevel::Error,
            std::string(label) + " must be one of the supported values. Received: "
                + std::to_string(static_cast<long long>(raw_numeric)), false);
    }

private:
    void BeginRequestApply()
    {
        ClearValidationIssues();
        CoercePositiveScalar(m_request.job_count, "--jobs", 1, LogLevel::Warning, "Thread size");
        CoerceScalar(m_request.verbosity, "--verbose",
            [](int candidate)
            {
                return candidate >= static_cast<int>(LogLevel::Error)
                    && candidate <= static_cast<int>(LogLevel::Debug);
            },
            static_cast<int>(LogLevel::Info), LogLevel::Warning,
            "Invalid verbose level: " + std::to_string(m_request.verbosity)
                + ", using default level 3 [Info]");
        m_request.output_dir =
            std::filesystem::path(path_helper::EnsureTrailingSlash(m_request.output_dir));
    }
    void ClearPrepareValidationIssues()
    {
        ClearValidationIssues(ValidationPhase::Prepare);
    }
    template <typename FieldType, typename Predicate>
    void CoerceScalar(
        FieldType & field,
        std::string_view option_name,
        Predicate is_valid,
        FieldType fallback_value,
        LogLevel issue_level,
        const std::string & message)
    {
        ClearPrepareValidationIssues();
        if (is_valid(field)) return;
        field = fallback_value;
        if (issue_level == LogLevel::Warning)
        {
            AddValidationIssue(option_name, ValidationPhase::Parse, LogLevel::Warning, message, true);
            return;
        }
        AddValidationIssue(option_name, ValidationPhase::Parse, issue_level, message, false);
    }
    bool RunFilesystemPreflight()
    {
        const auto & folder_path{ OutputFolder() };
        if (!folder_path.empty() && !std::filesystem::exists(folder_path))
        {
            std::error_code error_code;
            std::filesystem::create_directories(folder_path, error_code);
            if (error_code)
            {
                AddValidationIssue("--folder", ValidationPhase::Prepare, LogLevel::Error,
                    "Failed to create output directory '" + folder_path.string()
                    + "': " + error_code.message(), false);
            }
        }
        ReportValidationIssues();
        return !HasValidationErrors();
    }
    void ReportValidationIssues() const
    {
        static constexpr std::array<std::string_view, 2> phase_labels{ "parse", "prepare" };
        for (const auto & issue : m_validation_issues)
        {
            const auto phase_index{ static_cast<std::size_t>(issue.phase) };
            const auto phase_label{
                phase_index < phase_labels.size()
                    ? phase_labels[phase_index]
                    : std::string_view{ "unknown" }
            };
            std::string prefix{ "[" + std::string(phase_label) };
            if (issue.auto_corrected) prefix += "; auto-corrected";
            prefix += "]";
            Logger::Log(issue.level, prefix + " Option " + issue.option_name + ": " + issue.message);
        }
    }
    void ClearValidationIssues(std::optional<ValidationPhase> phase = std::nullopt)
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
    void AddValidationIssue(
        std::string_view option_name,
        ValidationPhase phase,
        LogLevel level,
        const std::string & message,
        bool auto_corrected)
    {
        m_validation_issues.push_back(
            ValidationIssueRecord{ std::string(option_name), phase, level, message, auto_corrected }
        );
    }
    
};

} // namespace rhbm_gem
