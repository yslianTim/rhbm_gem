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

    virtual CommonCommandRequest & MutableCommonRequest() = 0;
    virtual const CommonCommandRequest & CommonRequest() const = 0;
    virtual void ValidateOptions() {}
    virtual void ResetRuntimeState() {}
    virtual bool ExecuteImpl() = 0;
    int ThreadSize() const { return CommonRequest().thread_size; }
    int VerboseLevel() const { return CommonRequest().verbose_level; }
    const std::filesystem::path & OutputFolder() const { return CommonRequest().folder_path; }
    void InvalidatePreparedState();
    void NormalizeCommonRequest(CommonCommandRequest & request);
    void AddValidationError(
        std::string_view option_name,
        const std::string & message,
        ValidationPhase phase = ValidationPhase::Prepare);
    void AddNormalizationWarning(std::string_view option_name, const std::string & message);
    void ClearParseIssues(std::string_view option_name);
    void ClearPrepareIssues(std::string_view option_name);

    void NormalizeRequiredPath(
        std::filesystem::path & field,
        std::string_view option_name,
        std::string_view label);
    void NormalizeOptionalPath(
        std::filesystem::path & field,
        std::string_view option_name,
        std::string_view label);
    template <typename FieldType, typename Predicate>
    void NormalizeScalar(
        FieldType & field,
        std::string_view option_name,
        Predicate is_valid,
        FieldType fallback_value,
        LogLevel issue_level,
        const std::string & message)
    {
        InvalidatePreparedState();
        ClearParseIssues(option_name);
        if (is_valid(field))
        {
            return;
        }

        field = fallback_value;
        if (issue_level == LogLevel::Warning)
        {
            AddNormalizationWarning(option_name, message);
            return;
        }

        AddValidationIssue(option_name, ValidationPhase::Parse, issue_level, message, false);
    }
    template <typename FieldType>
    void NormalizeEnum(
        FieldType & field,
        std::string_view option_name,
        FieldType fallback_value,
        std::string_view label)
    {
        InvalidatePreparedState();
        ClearParseIssues(option_name);
        using UnderlyingType = std::underlying_type_t<FieldType>;
        const auto raw_numeric{ static_cast<UnderlyingType>(field) };
        if (IsSupportedCommandEnumValue(field))
        {
            return;
        }

        field = fallback_value;
        AddValidationError(
            option_name,
            std::string(label) + " must be one of the supported values. Received: "
                + std::to_string(static_cast<long long>(raw_numeric)),
            ValidationPhase::Parse);
    }
    std::filesystem::path BuildOutputPath(
        std::string_view stem,
        std::string_view extension) const;

private:
    bool m_is_prepared_for_execution{ false };
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
public:
    CommandWithRequest() :
        CommandBase{}
    {
    }

    void ApplyRequest(const Request & request)
    {
        m_request = request;
        InvalidatePreparedState();
        NormalizeCommonRequest(m_request.common);
        NormalizeRequest();
    }

protected:
    CommonCommandRequest & MutableCommonRequest() override { return m_request.common; }
    const CommonCommandRequest & CommonRequest() const override { return m_request.common; }
    const Request & RequestOptions() const { return m_request; }
    Request & MutableRequest() { return m_request; }
    virtual void NormalizeRequest() {}

private:
    Request m_request{};
};

} // namespace rhbm_gem
