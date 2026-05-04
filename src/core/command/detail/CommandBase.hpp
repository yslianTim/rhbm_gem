#pragma once

#include <cmath>
#include <filesystem>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem {

template <typename Request>
class CommandWithRequest;

namespace command_validation_detail {

template <typename Type>
std::string ToString(const Type & value)
{
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

template <typename Type>
std::string BuildConstraintMessage(
    std::string_view label,
    std::string_view constraint,
    const Type & fallback_value,
    LogLevel issue_level)
{
    std::string message{ std::string(label) + " must be " + std::string(constraint) };
    if (issue_level == LogLevel::Warning)
    {
        message += ". Using " + ToString(fallback_value) + " instead.";
    }
    return message;
}

template <typename FieldType, typename LowerType, typename UpperType>
std::string BuildExclusiveInclusiveRangeMessage(
    std::string_view label,
    LowerType lower,
    UpperType upper,
    const FieldType & fallback_value,
    LogLevel issue_level)
{
    std::string message{
        std::string(label) + " must be a finite value within ("
            + ToString(lower) + ", " + ToString(upper) + "]"
    };
    if (issue_level == LogLevel::Warning)
    {
        message += ". Using " + ToString(fallback_value) + " instead.";
    }
    return message;
}

} // namespace command_validation_detail

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

class CommandBase
{
    std::vector<ValidationIssueRecord> m_validation_issues;
    bool m_was_prepared{ false };

    template <typename Request>
    friend class CommandWithRequest;

public:
    virtual ~CommandBase() = default;
    bool Run();
    bool WasPrepared() const { return m_was_prepared; }
    bool HasValidationErrors() const;
    const std::vector<ValidationIssueRecord> & GetValidationIssues() const { return m_validation_issues; }

protected:
    CommandBase() = default;
    virtual void ValidatePreparedRequest() {}
    virtual void ResetRuntimeState() {}
    virtual bool ExecuteImpl() = 0;
    int ThreadSize() const { return BaseRequest().job_count; }
    const std::filesystem::path & OutputFolder() const { return BaseRequest().output_dir; }
    void AddParseError(std::string_view option_name, const std::string & message);
    void AddParseNormalizationWarning(std::string_view option_name, const std::string & message);
    void ValidateRequiredPath(
        const std::filesystem::path & path,
        std::string_view option_name,
        std::string_view label);
    void ValidateOptionalPath(
        const std::filesystem::path & path,
        std::string_view option_name,
        std::string_view label);
    template <typename Container>
    void RequireNonEmptyList(
        const Container & field,
        std::string_view option_name,
        std::string_view label)
    {
        if (!field.empty()) return;
        AddValidationIssue(
            option_name,
            ValidationPhase::Parse,
            LogLevel::Error,
            std::string(label) + " cannot be empty.",
            false);
    }
    void RequirePrepareCondition(
        bool condition,
        std::string_view option_name,
        const std::string & message);
    template <typename FieldType>
    void CoercePositiveScalar(
        FieldType & field,
        std::string_view option_name,
        FieldType fallback_value,
        LogLevel issue_level,
        std::string_view label)
    {
        CoerceScalar(
            field,
            option_name,
            [](const auto candidate) { return numeric_validation::IsPositive(candidate); },
            fallback_value,
            issue_level,
            command_validation_detail::BuildConstraintMessage(label, "positive", fallback_value, issue_level));
    }
    template <typename FieldType>
    void CoerceFinitePositiveScalar(
        FieldType & field,
        std::string_view option_name,
        FieldType fallback_value,
        LogLevel issue_level,
        std::string_view label)
    {
        CoerceScalar(
            field,
            option_name,
            [](const auto candidate) { return numeric_validation::IsFinitePositive(candidate); },
            fallback_value,
            issue_level,
            command_validation_detail::BuildConstraintMessage(
                label, "a finite positive value", fallback_value, issue_level));
    }
    template <typename FieldType>
    void CoerceFiniteNonNegativeScalar(
        FieldType & field,
        std::string_view option_name,
        FieldType fallback_value,
        LogLevel issue_level,
        std::string_view label)
    {
        CoerceScalar(
            field,
            option_name,
            [](const auto candidate) { return numeric_validation::IsFiniteNonNegative(candidate); },
            fallback_value,
            issue_level,
            command_validation_detail::BuildConstraintMessage(
                label, "a finite non-negative value", fallback_value, issue_level));
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
        CoerceScalar(
            field,
            option_name,
            [lower, upper](const auto candidate)
            {
                return numeric_validation::IsFiniteExclusiveInclusiveRange(candidate, lower, upper);
            },
            fallback_value,
            issue_level,
            command_validation_detail::BuildExclusiveInclusiveRangeMessage(
                label, lower, upper, fallback_value, issue_level));
    }
    template <typename FieldType>
    void CoerceEnum(
        FieldType & field,
        std::string_view option_name,
        FieldType fallback_value,
        std::string_view label)
    {
        InvalidatePreparedState();
        using UnderlyingType = std::underlying_type_t<FieldType>;
        const auto raw_numeric{ static_cast<UnderlyingType>(field) };
        for (const auto & option : internal::CommandEnumTraits<FieldType>::kOptions)
        {
            if (static_cast<UnderlyingType>(option.value) == raw_numeric)
            {
                return;
            }
        }
        field = fallback_value;
        AddValidationIssue(
            option_name,
            ValidationPhase::Parse,
            LogLevel::Error,
            std::string(label) + " must be one of the supported values. Received: "
                + std::to_string(static_cast<long long>(raw_numeric)), false);
    }
    std::filesystem::path BuildOutputPath(std::string_view stem, std::string_view extension) const;

private:
    virtual CommandRequestBase & BaseRequest() = 0;
    virtual const CommandRequestBase & BaseRequest() const = 0;
    int VerboseLevel() const { return BaseRequest().verbosity; }
    void BeginRequestApply(CommandRequestBase & request);
    void InvalidatePreparedState();
    void CoerceBaseRequest(CommandRequestBase & request);
    template <typename FieldType, typename Predicate>
    void CoerceScalar(
        FieldType & field,
        std::string_view option_name,
        Predicate is_valid,
        FieldType fallback_value,
        LogLevel issue_level,
        const std::string & message)
    {
        InvalidatePreparedState();
        if (is_valid(field)) return;
        field = fallback_value;
        if (issue_level == LogLevel::Warning)
        {
            AddValidationIssue(
                option_name,
                ValidationPhase::Parse,
                LogLevel::Warning,
                message,
                true);
            return;
        }
        AddValidationIssue(option_name, ValidationPhase::Parse, issue_level, message, false);
    }
    void BeginPreparationPass();
    bool RunValidationPass();
    bool RunFilesystemPreflight();
    void ReportValidationIssues() const;
    void ClearValidationIssues(std::optional<ValidationPhase> phase = std::nullopt);
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
    Request m_request{};

public:
    CommandWithRequest() : CommandBase{} {}
    void ApplyRequest(const Request & request)
    {
        m_request = request;
        BeginRequestApply(m_request);
        NormalizeAndValidateRequest();
    }

protected:
    const Request & RequestOptions() const { return m_request; }
    Request & MutableRequest() { return m_request; }
    virtual void NormalizeAndValidateRequest() {}

private:
    CommandRequestBase & BaseRequest() override { return m_request; }
    const CommandRequestBase & BaseRequest() const override { return m_request; }
    
};

} // namespace rhbm_gem
