#pragma once

#include <algorithm>
#include <filesystem>
#include <string>
#include <string_view>
#include <system_error>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <rhbm_gem/core/command/CommandTypes.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include "CommandCatalog.hpp"
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

    struct FieldMetadata
    {
        std::string option_name;
        std::string field_name;
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
        if (HasValidationErrors())
        {
            ReportPendingIssues();
            result.succeeded = false;
            result.issues = BuildDiagnostics();
            return result;
        }
        ValidatePreparedRequest(command_request);
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
    template <typename Owner, typename FieldType>
    void AddParseError(
        const Request & request,
        FieldType Owner::* member,
        const std::string & message)
    {
        static_cast<void>(request);
        const auto metadata{ ResolveFieldMetadata(member) };
        AddParseError(metadata.option_name, message);
    }
    void AddParseNormalizationWarning(std::string_view option_name, const std::string & message)
    {
        AddPendingIssue(option_name, LogLevel::Warning, message, true);
    }
    template <typename Owner, typename FieldType>
    void AddParseNormalizationWarning(
        const Request & request,
        FieldType Owner::* member,
        const std::string & message)
    {
        static_cast<void>(request);
        const auto metadata{ ResolveFieldMetadata(member) };
        AddParseNormalizationWarning(metadata.option_name, message);
    }
    template <typename Owner>
    void RequireExistingPath(
        const Request & request,
        std::filesystem::path Owner::* member)
    {
        const auto metadata{ ResolveFieldMetadata(member) };
        const auto & path{ request.*member };
        if (path.empty())
        {
            AddPendingIssue(metadata.option_name, LogLevel::Error,
                metadata.field_name + " path is required.", false);
            return;
        }

        std::error_code error_code;
        if (!std::filesystem::exists(path, error_code))
        {
            AddPendingIssue(metadata.option_name, LogLevel::Error,
                metadata.field_name + " does not exist: " + path.string(), false);
        }
    }
    template <typename Owner>
    void RequireOptionalExistingPath(
        const Request & request,
        std::filesystem::path Owner::* member)
    {
        const auto metadata{ ResolveFieldMetadata(member) };
        const auto & path{ request.*member };
        if (path.empty()) return;

        std::error_code error_code;
        if (!std::filesystem::exists(path, error_code))
        {
            AddPendingIssue(metadata.option_name, LogLevel::Error,
                metadata.field_name + " does not exist: " + path.string(), false);
        }
    }
    template <typename Owner, typename Container>
    void RequireNonEmptyList(
        const Request & request,
        Container Owner::* member)
    {
        const auto metadata{ ResolveFieldMetadata(member) };
        const auto & field{ request.*member };
        if (!field.empty()) return;
        AddPendingIssue(metadata.option_name, LogLevel::Error,
            metadata.field_name + " cannot be empty.", false);
    }
    void RequirePrepareCondition(
        bool condition,
        const std::string & message)
    {
        if (condition) return;
        AddPendingIssue("request", LogLevel::Error, message, false);
    }
    template <typename Owner, typename FieldType>
    void NormalizePositiveScalar(
        Request & request,
        FieldType Owner::* member,
        FieldType fallback_value)
    {
        const auto metadata{ ResolveFieldMetadata(member) };
        std::string message{
            metadata.field_name + " must be positive. Using "
            + string_helper::ToStringWithPrecision(fallback_value) + " instead."
        };
        NormalizeScalar(request.*member, metadata,
            [](const auto candidate) { return numeric_validation::IsPositive(candidate); },
            fallback_value, message);
    }
    template <typename Owner, typename FieldType>
    void NormalizeFinitePositiveScalar(
        Request & request,
        FieldType Owner::* member,
        FieldType fallback_value)
    {
        const auto metadata{ ResolveFieldMetadata(member) };
        std::string message{
            metadata.field_name + " must be a finite positive value. Using "
            + string_helper::ToStringWithPrecision(fallback_value) + " instead."
        };
        NormalizeScalar(request.*member, metadata,
            [](const auto candidate) { return numeric_validation::IsFinitePositive(candidate); },
            fallback_value, message);
    }
    template <typename Owner, typename FieldType>
    void NormalizeFiniteNonNegativeScalar(
        Request & request,
        FieldType Owner::* member,
        FieldType fallback_value)
    {
        const auto metadata{ ResolveFieldMetadata(member) };
        std::string message{
            metadata.field_name + " must be a finite non-negative value. Using "
            + string_helper::ToStringWithPrecision(fallback_value) + " instead."
        };
        NormalizeScalar(request.*member, metadata,
            [](const auto candidate) { return numeric_validation::IsFiniteNonNegative(candidate); },
            fallback_value, message);
    }
    template <typename Owner, typename FieldType, typename LowerType, typename UpperType>
    void NormalizeFiniteExclusiveInclusiveRangeScalar(
        Request & request,
        FieldType Owner::* member,
        LowerType lower,
        UpperType upper,
        FieldType fallback_value)
    {
        const auto metadata{ ResolveFieldMetadata(member) };
        std::string message{
            metadata.field_name + " must be a finite value within ("
                + string_helper::ToStringWithPrecision(lower) + ", "
                + string_helper::ToStringWithPrecision(upper) + "]. Using "
                + string_helper::ToStringWithPrecision(fallback_value) + " instead."
        };
        NormalizeScalar(request.*member, metadata,
            [lower, upper](const auto candidate)
            {
                return numeric_validation::IsFiniteExclusiveInclusiveRange(candidate, lower, upper);
            },
            fallback_value, message);
    }
    template <typename Owner, typename FieldType>
    void RequirePositiveScalar(const Request & request, FieldType Owner::* member)
    {
        const auto metadata{ ResolveFieldMetadata(member) };
        RequireScalar(request.*member, metadata,
            [](const auto candidate) { return numeric_validation::IsPositive(candidate); },
            metadata.field_name + " must be positive.");
    }
    template <typename Owner, typename FieldType>
    void RequireFinitePositiveScalar(const Request & request, FieldType Owner::* member)
    {
        const auto metadata{ ResolveFieldMetadata(member) };
        RequireScalar(request.*member, metadata,
            [](const auto candidate) { return numeric_validation::IsFinitePositive(candidate); },
            metadata.field_name + " must be a finite positive value.");
    }
    template <typename Owner, typename FieldType>
    void RequireFiniteNonNegativeScalar(
        const Request & request,
        FieldType Owner::* member)
    {
        const auto metadata{ ResolveFieldMetadata(member) };
        RequireScalar(request.*member, metadata,
            [](const auto candidate) { return numeric_validation::IsFiniteNonNegative(candidate); },
            metadata.field_name + " must be a finite non-negative value.");
    }
    template <typename Owner, typename FieldType, typename LowerType, typename UpperType>
    void RequireFiniteExclusiveInclusiveRangeScalar(
        const Request & request,
        FieldType Owner::* member,
        LowerType lower,
        UpperType upper)
    {
        const auto metadata{ ResolveFieldMetadata(member) };
        std::string message{
            metadata.field_name + " must be a finite value within ("
                + string_helper::ToStringWithPrecision(lower) + ", "
                + string_helper::ToStringWithPrecision(upper) + "]."
        };
        RequireScalar(request.*member, metadata,
            [lower, upper](const auto candidate)
            {
                return numeric_validation::IsFiniteExclusiveInclusiveRange(candidate, lower, upper);
            },
            message);
    }
    template <typename Owner, typename FieldType>
    void RequireEnum(
        const Request & request,
        FieldType Owner::* member)
    {
        const auto metadata{ ResolveFieldMetadata(member) };
        const auto & field{ request.*member };
        using UnderlyingType = std::underlying_type_t<FieldType>;
        const auto raw_numeric{ static_cast<UnderlyingType>(field) };
        for (const auto & option : command_internal::CommandEnumTraits<FieldType>::kOptions)
        {
            if (static_cast<UnderlyingType>(option.value) == raw_numeric)
            {
                return;
            }
        }
        AddPendingIssue(metadata.option_name, LogLevel::Error,
            metadata.field_name + " must be one of the supported values. Received: "
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
        NormalizePositiveScalar(request, &CommandRequestBase::job_count, 1);
        const auto verbosity_metadata{ ResolveFieldMetadata(&CommandRequestBase::verbosity) };
        NormalizeScalar(request.verbosity, verbosity_metadata,
            [](int candidate)
            {
                return candidate >= static_cast<int>(LogLevel::Error)
                    && candidate <= static_cast<int>(LogLevel::Debug);
            },
            static_cast<int>(LogLevel::Info),
            "Invalid verbose level: " + std::to_string(request.verbosity)
                + ", using default level 3 [Info]");
        request.output_dir =
            std::filesystem::path(path_helper::EnsureTrailingSlash(request.output_dir));
    }
    template <typename FieldType, typename Predicate>
    void NormalizeScalar(
        FieldType & field,
        const FieldMetadata & metadata,
        Predicate is_valid,
        FieldType fallback_value,
        const std::string & message)
    {
        if (is_valid(field)) return;
        field = fallback_value;
        AddPendingIssue(metadata.option_name, LogLevel::Warning, message, true);
    }
    template <typename FieldType, typename Predicate>
    void RequireScalar(
        const FieldType & field,
        const FieldMetadata & metadata,
        Predicate is_valid,
        const std::string & message)
    {
        if (is_valid(field)) return;
        AddPendingIssue(metadata.option_name, LogLevel::Error, message, false);
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
    template <typename Owner, typename FieldType>
    FieldMetadata ResolveFieldMetadata(FieldType Owner::* member) const
    {
        static_assert(std::is_base_of_v<Owner, Request>,
            "Validated field owner must be the command request type or a base request type.");
        FieldMetadata metadata;
        bool found{ false };
        command_internal::RequestFieldCatalog<Owner>::Visit([&](const auto & field)
        {
            if (found) return;
            if constexpr (std::is_same_v<decltype(field.member), FieldType Owner::*>)
            {
                if (field.member != member) return;
                metadata.option_name = field.cli_flags;
                metadata.field_name = field.field_name;
                found = true;
            }
        });
        if (!found)
        {
            throw std::logic_error(
                "No RequestField metadata registered for validated command field.");
        }
        return metadata;
    }

};

} // namespace rhbm_gem
