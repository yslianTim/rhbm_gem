#pragma once

#include <cmath>
#include <filesystem>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include "CommandObjectCache.hpp"
#include "CommandEnumMetadata.hpp"

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

class CommandBase
{
public:
    virtual ~CommandBase() = default;
    bool Run();
    bool WasPrepared() const { return m_was_prepared; }
    bool HasValidationErrors() const;
    const std::vector<ValidationIssueRecord> & GetValidationIssues() const { return m_validation_issues; }

protected:
    command_detail::CommandObjectCache m_runtime_object_cache;
    std::unique_ptr<DataRepository> m_data_repository;
    std::vector<ValidationIssueRecord> m_validation_issues;

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
        std::filesystem::path & field,
        std::string_view option_name,
        std::string_view label);
    void ValidateOptionalPath(
        std::filesystem::path & field,
        std::string_view option_name,
        std::string_view label);
    void RequireNonEmptyText(
        std::string_view field,
        std::string_view option_name,
        std::string_view label,
        ValidationPhase phase = ValidationPhase::Parse);
    template <typename Container>
    void RequireNonEmptyList(
        const Container & field,
        std::string_view option_name,
        std::string_view label,
        ValidationPhase phase = ValidationPhase::Parse)
    {
        BeginValidationMutation(phase);
        ClearValidationIssues(option_name, phase);
        if (!field.empty())
        {
            return;
        }

        AddValidationError(
            option_name,
            std::string(label) + " cannot be empty.",
            phase);
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
        if (internal::IsSupportedCommandEnumValue(field))
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
    void AttachDataRepository(const std::filesystem::path & database_path)
    {
        m_data_repository = std::make_unique<DataRepository>(database_path);
    }
    template <typename TypedDataObject>
    std::shared_ptr<TypedDataObject> LoadInputFile(
        const std::filesystem::path & filename,
        const std::string & key_tag)
    {
        std::unique_ptr<TypedDataObject> data_object;
        if constexpr (std::is_same_v<TypedDataObject, ModelObject>)
        {
            data_object = ReadModel(filename);
        }
        else if constexpr (std::is_same_v<TypedDataObject, MapObject>)
        {
            data_object = ReadMap(filename);
        }
        else
        {
            static_assert(
                std::is_same_v<TypedDataObject, ModelObject>
                    || std::is_same_v<TypedDataObject, MapObject>,
                "LoadInputFile only supports ModelObject or MapObject.");
        }
        data_object->SetKeyTag(key_tag);
        std::shared_ptr<TypedDataObject> shared_object{ std::move(data_object) };
        if constexpr (std::is_same_v<TypedDataObject, ModelObject>)
        {
            m_runtime_object_cache.PutModel(key_tag, shared_object);
        }
        else
        {
            m_runtime_object_cache.PutMap(key_tag, shared_object);
        }
        return shared_object;
    }
    template <typename TypedDataObject>
    std::shared_ptr<TypedDataObject> LoadPersistedObject(const std::string & key_tag)
    {
        if (m_data_repository == nullptr)
        {
            throw std::runtime_error("Database repository is not initialized.");
        }
        std::unique_ptr<TypedDataObject> data_object;
        if constexpr (std::is_same_v<TypedDataObject, ModelObject>)
        {
            data_object = m_data_repository->LoadModel(key_tag);
        }
        else if constexpr (std::is_same_v<TypedDataObject, MapObject>)
        {
            data_object = m_data_repository->LoadMap(key_tag);
        }
        else
        {
            static_assert(
                std::is_same_v<TypedDataObject, ModelObject>
                    || std::is_same_v<TypedDataObject, MapObject>,
                "LoadPersistedObject only supports ModelObject or MapObject.");
        }
        std::shared_ptr<TypedDataObject> shared_object{ std::move(data_object) };
        if constexpr (std::is_same_v<TypedDataObject, ModelObject>)
        {
            m_runtime_object_cache.PutModel(key_tag, shared_object);
        }
        else
        {
            m_runtime_object_cache.PutMap(key_tag, shared_object);
        }
        return shared_object;
    }
    void SaveStoredObject(
        const std::string & key_tag,
        const std::string & persisted_key = "")
    {
        if (m_data_repository == nullptr)
        {
            throw std::runtime_error("Database repository is not initialized.");
        }
        const auto saved_key_tag{ persisted_key.empty() ? key_tag : persisted_key };
        switch (m_runtime_object_cache.GetKind(key_tag))
        {
        case command_detail::CommandObjectCache::ObjectKind::Model:
            m_data_repository->SaveModel(*m_runtime_object_cache.GetModel(key_tag), saved_key_tag);
            return;
        case command_detail::CommandObjectCache::ObjectKind::Map:
            m_data_repository->SaveMap(*m_runtime_object_cache.GetMap(key_tag), saved_key_tag);
            return;
        }
    }

private:
    CommandRequestBase * m_base_request{ nullptr };
    bool m_was_prepared{ false };
    CommandRequestBase & BaseRequest() { return *m_base_request; }
    const CommandRequestBase & BaseRequest() const { return *m_base_request; }
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
public:
    CommandWithRequest() :
        CommandBase{}
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

private:
    Request m_request{};
};

} // namespace rhbm_gem
