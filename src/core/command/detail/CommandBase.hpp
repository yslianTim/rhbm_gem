#pragma once

#include <cmath>
#include <filesystem>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include "CommandEnumMetadata.hpp"

namespace rhbm_gem {

namespace command_validation_detail {

template <typename Type>
std::string ToString(const Type & value)
{
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

inline std::string BuildConstraintMessage(
    std::string_view label,
    std::string_view constraint,
    LogLevel issue_level)
{
    std::string message{ std::string(label) + " must be " + std::string(constraint) };
    if (issue_level == LogLevel::Warning)
    {
        message += '.';
    }
    return message;
}

template <typename Type>
std::string BuildConstraintMessage(
    std::string_view label,
    std::string_view constraint,
    const Type & fallback_value,
    LogLevel issue_level)
{
    std::string message{ BuildConstraintMessage(label, constraint, issue_level) };
    if (issue_level == LogLevel::Warning)
    {
        message += " Using " + ToString(fallback_value) + " instead.";
    }
    return message;
}

template <typename LowerType, typename UpperType>
std::string BuildExclusiveInclusiveRangeMessage(
    std::string_view label,
    LowerType lower,
    UpperType upper,
    LogLevel issue_level)
{
    std::string message{
        std::string(label) + " must be a finite value within ("
            + ToString(lower) + ", " + ToString(upper) + "]"
    };
    if (issue_level == LogLevel::Warning)
    {
        message += '.';
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
    std::string message{ BuildExclusiveInclusiveRangeMessage(label, lower, upper, issue_level) };
    if (issue_level == LogLevel::Warning)
    {
        message += " Using " + ToString(fallback_value) + " instead.";
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
    CommandRequestBase * m_base_request{ nullptr };
    std::unique_ptr<DataRepository> m_data_repository;
    std::vector<ValidationIssueRecord> m_validation_issues;
    bool m_was_prepared{ false };

public:
    virtual ~CommandBase() = default;
    bool Run();
    bool WasPrepared() const { return m_was_prepared; }
    bool HasValidationErrors() const;
    const std::vector<ValidationIssueRecord> & GetValidationIssues() const { return m_validation_issues; }

protected:
    CommandBase() = default;
    virtual void ValidateOptions() {}
    virtual void ResetRuntimeState() {}
    virtual bool ExecuteImpl() = 0;
    int ThreadSize() const { return BaseRequest().job_count; }
    int VerboseLevel() const { return BaseRequest().verbosity; }
    const std::filesystem::path & OutputFolder() const { return BaseRequest().output_dir; }
    void BindBaseRequest(CommandRequestBase & request) { m_base_request = &request; }
    void InvalidatePreparedState();
    void CoerceBaseRequest(CommandRequestBase & request);
    void AddValidationError(
        std::string_view option_name,
        const std::string & message,
        ValidationPhase phase = ValidationPhase::Prepare);
    void AddNormalizationWarning(std::string_view option_name, const std::string & message);
    void ClearParseIssues(std::string_view option_name);
    void ClearPrepareIssues(std::string_view option_name);

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
        std::string_view label,
        ValidationPhase phase = ValidationPhase::Parse)
    {
        BeginValidationMutation(phase);
        ClearValidationIssues(option_name, phase);
        if (!field.empty()) return;
        AddValidationError(option_name, std::string(label) + " cannot be empty.", phase);
    }
    void RequireCondition(
        bool condition,
        std::string_view option_name,
        const std::string & message,
        ValidationPhase phase = ValidationPhase::Prepare);
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
        ClearParseIssues(option_name);
        if (is_valid(field)) return;
        field = fallback_value;
        if (issue_level == LogLevel::Warning)
        {
            AddNormalizationWarning(option_name, message);
            return;
        }
        AddValidationIssue(option_name, ValidationPhase::Parse, issue_level, message, false);
    }
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
            command_validation_detail::BuildConstraintMessage(
                label, "positive", fallback_value, issue_level));
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
        ClearParseIssues(option_name);
        using UnderlyingType = std::underlying_type_t<FieldType>;
        const auto raw_numeric{ static_cast<UnderlyingType>(field) };
        if (internal::IsSupportedCommandEnumValue(field)) return;
        field = fallback_value;
        AddValidationError(
            option_name,
            std::string(label) + " must be one of the supported values. Received: "
                + std::to_string(static_cast<long long>(raw_numeric)),
            ValidationPhase::Parse);
    }
    std::filesystem::path BuildOutputPath(std::string_view stem, std::string_view extension) const;
    void OpenDataRepository(const std::filesystem::path & database_path);
    std::shared_ptr<ModelObject> LoadModelFile(
        const std::filesystem::path & filename,
        const std::string & key_tag);
    std::shared_ptr<MapObject> LoadMapFile(
        const std::filesystem::path & filename,
        const std::string & key_tag);
    std::shared_ptr<ModelObject> LoadModelFromRepository(const std::string & key_tag);
    std::shared_ptr<MapObject> LoadMapFromRepository(const std::string & key_tag);
    void SaveModelToRepository(const ModelObject & model_object, const std::string & key_tag);
    void SaveMapToRepository(const MapObject & map_object, const std::string & key_tag);

private:
    CommandRequestBase & BaseRequest() { return *m_base_request; }
    const CommandRequestBase & BaseRequest() const { return *m_base_request; }
    DataRepository & RequireDataRepository();
    void BeginValidationMutation(ValidationPhase phase);
    void BeginPreparationPass();
    bool RunValidationPass();
    bool RunFilesystemPreflight();
    void ReportValidationIssues() const;
    void ClearValidationIssues(std::string_view option_name, std::optional<ValidationPhase> phase);
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
    CommandWithRequest() : CommandBase{}
    {
        BindBaseRequest(m_request);
    }

    void ApplyRequest(const Request & request)
    {
        m_request = request;
        InvalidatePreparedState();
        CoerceBaseRequest(m_request);
        NormalizeRequest();
    }

protected:
    const Request & RequestOptions() const { return m_request; }
    Request & MutableRequest() { return m_request; }
    virtual void NormalizeRequest() {}
    
};

} // namespace rhbm_gem
