#pragma once

#include <cmath>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/core/command/CommandEnumClass.hpp>
#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

namespace rhbm_gem {

class CommandBase
{
public:
    virtual ~CommandBase() = default;
    bool Execute();
    bool PrepareForExecution();
    bool HasValidationErrors() const;
    const std::vector<ValidationIssue> & GetValidationIssues() const { return m_validation_issues; }

protected:
    DataObjectManager m_data_manager;
    std::vector<ValidationIssue> m_validation_issues;

    CommandBase() = default;

    virtual void ValidateOptions() {}
    virtual void ResetRuntimeState() {}
    virtual bool ExecuteImpl() = 0;
    int ThreadSize() const { return m_common_request.thread_size; }
    int VerboseLevel() const { return m_common_request.verbose_level; }
    const std::filesystem::path & OutputFolder() const { return m_common_request.folder_path; }
    const CommonCommandRequest & CommonRequest() const { return m_common_request; }
    void SetThreadSize(int value);
    void SetVerboseLevel(int value);
    void SetFolderPath(const std::filesystem::path & path);
    void ApplyCommonRequest(const CommonCommandRequest & request);

    // Core command extension API used by concrete command setters and validators.
    template <typename Mutator>
    void MutateOptions(Mutator && mutator)
    {
        InvalidatePreparedState();
        std::forward<Mutator>(mutator)();
    }
    template <typename FieldType, typename ValueType>
    void AssignOption(FieldType & field, ValueType && value)
    {
        MutateOptions([&]()
        {
            field = std::forward<ValueType>(value);
        });
    }
    void AddValidationError(
        std::string_view option_name,
        const std::string & message,
        ValidationPhase phase = ValidationPhase::Prepare);
    void AddNormalizationWarning(std::string_view option_name, const std::string & message);
    void ResetParseIssues(std::string_view option_name);
    void ResetPrepareIssues(std::string_view option_name);

    // Base convenience helpers built on top of the core extension API.
    void SetRequiredExistingPathOption(
        std::filesystem::path & field,
        const std::filesystem::path & value,
        std::string_view option_name,
        std::string_view label);
    void SetOptionalExistingPathOption(
        std::filesystem::path & field,
        const std::filesystem::path & value,
        std::string_view option_name,
        std::string_view label);
    // Invalid input is normalized to a fallback value and recorded as a parse warning.
    template <typename FieldType, typename RawType, typename Predicate>
    void SetNormalizedScalarOption(
        FieldType & field,
        RawType raw_value,
        std::string_view option_name,
        Predicate is_valid,
        FieldType fallback_value,
        const std::string & warning_message)
    {
        MutateOptions([&]()
        {
            ResetParseIssues(option_name);
            if (is_valid(raw_value))
            {
                field = static_cast<FieldType>(raw_value);
                return;
            }

            field = fallback_value;
            AddNormalizationWarning(std::string(option_name), warning_message);
        });
    }
    template <typename FieldType, typename RawType>
    void SetFinitePositiveScalarOption(
        FieldType & field,
        RawType raw_value,
        std::string_view option_name,
        FieldType fallback_value,
        const std::string & error_message)
    {
        SetValidatedScalarOption(
            field,
            raw_value,
            option_name,
            [](RawType candidate)
            {
                return std::isfinite(static_cast<long double>(candidate)) && candidate > 0;
            },
            fallback_value,
            error_message);
    }
    template <typename FieldType, typename RawType>
    void SetFiniteNonNegativeScalarOption(
        FieldType & field,
        RawType raw_value,
        std::string_view option_name,
        FieldType fallback_value,
        const std::string & error_message)
    {
        SetValidatedScalarOption(
            field,
            raw_value,
            option_name,
            [](RawType candidate)
            {
                return std::isfinite(static_cast<long double>(candidate)) && candidate >= 0;
            },
            fallback_value,
            error_message);
    }
    template <typename FieldType, typename RawType>
    void SetPositiveScalarOption(
        FieldType & field,
        RawType raw_value,
        std::string_view option_name,
        FieldType fallback_value,
        const std::string & error_message)
    {
        SetValidatedScalarOption(
            field,
            raw_value,
            option_name,
            [](RawType candidate) { return candidate > 0; },
            fallback_value,
            error_message);
    }
    // Invalid enum input falls back to a safe value and is recorded as a parse error.
    template <typename EnumType>
    void SetValidatedEnumOption(
        EnumType & field,
        EnumType raw_value,
        std::string_view option_name,
        EnumType fallback_value,
        std::string_view label)
    {
        MutateOptions([&]()
        {
            ResetParseIssues(option_name);
            using UnderlyingType = std::underlying_type_t<EnumType>;
            const auto raw_numeric{ static_cast<UnderlyingType>(raw_value) };
            if (IsSupportedCommandEnumValue(raw_value))
            {
                field = raw_value;
                return;
            }

            field = fallback_value;
            AddValidationError(
                std::string(option_name),
                std::string(label) + " must be one of the supported values. Received: "
                    + std::to_string(static_cast<long long>(raw_numeric)),
                ValidationPhase::Parse);
        });
    }
    std::filesystem::path BuildOutputPath(
        std::string_view stem,
        std::string_view extension) const;

private:
    CommonCommandRequest m_common_request{};
    bool m_is_prepared_for_execution{ false };
    void InvalidatePreparedState();
    void BeginPreparationPass();
    bool RunValidationPass();
    bool RunFilesystemPreflight();
    void ReportValidationIssues() const;
    void ClearValidationIssues(std::string_view option_name, std::optional<ValidationPhase> phase);
    void ClearValidationIssues(std::optional<ValidationPhase> phase = std::nullopt);
    // Invalid input falls back to a safe value and is recorded as a parse error.
    template <typename FieldType, typename RawType, typename Predicate>
    void SetValidatedScalarOption(
        FieldType & field,
        RawType raw_value,
        std::string_view option_name,
        Predicate is_valid,
        FieldType fallback_value,
        const std::string & error_message)
    {
        MutateOptions([&]()
        {
            ResetParseIssues(option_name);
            if (is_valid(raw_value))
            {
                field = static_cast<FieldType>(raw_value);
                return;
            }

            field = fallback_value;
            AddValidationError(
                std::string(option_name),
                error_message,
                ValidationPhase::Parse);
        });
    }
    void ValidateRequiredExistingPath(
        const std::filesystem::path & path,
        std::string_view option_name,
        std::string_view label);
    void ValidateOptionalExistingPath(
        const std::filesystem::path & path,
        std::string_view option_name,
        std::string_view label);

    void AddValidationIssue(
        std::string_view option_name,
        ValidationPhase phase,
        LogLevel level,
        const std::string & message,
        bool auto_corrected);
};

template <typename Request>
class CommandWithRequest : public CommandBase
{
public:
    CommandWithRequest() :
        CommandBase{}
    {
    }

    void ApplyRequest(const Request & request)
    {
        AssignOption(m_request, request);
        ApplyCommonRequest(m_request.common);
        m_request.common = CommonRequest();
        NormalizeRequest();
    }

protected:
    const Request & RequestOptions() const { return m_request; }
    Request & MutableRequest() { return m_request; }
    virtual void NormalizeRequest() {}

private:
    Request m_request{};
};

} // namespace rhbm_gem
