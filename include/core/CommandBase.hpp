#pragma once

#include <cstdint>
#include <cmath>
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

#include "BuiltInCommandCatalog.hpp"
#include "CommandMetadata.hpp"
#include "DataObjectManager.hpp"
#include "Logger.hpp"
#include "OptionEnumTraits.hpp"

namespace CLI
{
    class App;
}

namespace rhbm_gem {

enum class ValidationPhase : std::uint8_t
{
    Parse = 0,
    Prepare = 1,
    Runtime = 2
};

struct ValidationIssue
{
    std::string option_name;
    ValidationPhase phase;
    LogLevel level;
    std::string message;
    bool auto_corrected{ false };
};

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
    bool Execute();
    virtual void RegisterCLIOptionsExtend(::CLI::App * command) = 0;
    virtual const CommandOptions & GetOptions() const = 0;
    virtual CommandOptions & GetOptions() = 0;
    virtual CommandId GetCommandId() const = 0;
    virtual void ValidateOptions() {}
    virtual void ResetRuntimeState() {}
    const CommandDescriptor & GetCommandDescriptor() const;
    CommonOptionMask GetCommonOptionMask() const
    {
        return GetCommandDescriptor().surface.common_options;
    }
    DatabaseUsage GetDatabaseUsage() const
    {
        return GetCommandDescriptor().database_usage;
    }
    BindingExposure GetBindingExposure() const
    {
        return GetCommandDescriptor().binding_exposure;
    }

    void RegisterCLIOptions(::CLI::App * command);
    bool PrepareForExecution();
    void SetThreadSize(int value);
    void SetVerboseLevel(int value);
    void SetDatabasePath(const std::filesystem::path & path);
    void SetFolderPath(const std::filesystem::path & path);
    bool IsValidateOptions() const { return !HasValidationErrors(); }
    bool HasValidationErrors(std::optional<ValidationPhase> phase = std::nullopt) const;
    void ReportValidationIssues() const;
    const std::vector<ValidationIssue> & GetValidationIssues() const { return m_validation_issues; }
    DataObjectManager * GetDataManagerPtr() { return &m_data_manager; }
    const DataObjectManager * GetDataManagerPtr() const { return &m_data_manager; }

protected:
    DataObjectManager m_data_manager;
    std::vector<ValidationIssue> m_validation_issues;

    CommandBase() = default;
    void InvalidatePreparedState();
    bool IsPreparedForExecution() const { return m_is_prepared_for_execution; }
    template <typename Mutator>
    void MutateOptions(Mutator && mutator)
    {
        InvalidatePreparedState();
        std::forward<Mutator>(mutator)();
    }
    void RegisterCLIOptionsBasic(::CLI::App * command);
    void RegisterDeprecatedCommonAliases(::CLI::App * command);
    void AddValidationError(
        const std::string & option_name,
        const std::string & message,
        ValidationPhase phase = ValidationPhase::Prepare);
    void AddValidationWarning(
        const std::string & option_name,
        const std::string & message,
        ValidationPhase phase = ValidationPhase::Prepare);
    void AddNormalizationWarning(const std::string & option_name, const std::string & message);
    void ResetParseIssues(std::string_view option_name);
    void ResetPrepareIssues(std::string_view option_name);
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
            for (const auto & entry : EnumOptionTraits<EnumType>::kBindingEntries)
            {
                if (static_cast<UnderlyingType>(entry.value) == raw_numeric)
                {
                    field = raw_value;
                    return;
                }
            }

            field = fallback_value;
            AddValidationError(
                std::string(option_name),
                std::string(label) + " must be one of the supported values. Received: "
                    + std::to_string(static_cast<long long>(raw_numeric)),
                ValidationPhase::Parse);
        });
    }
    void RequireDatabaseManager();
    std::filesystem::path BuildOutputPath(
        std::string_view stem,
        std::string_view extension) const;
    virtual bool ExecuteImpl() = 0;

    template <typename TypedDataObject>
    std::shared_ptr<TypedDataObject> ProcessTypedFile(
        const std::filesystem::path & path,
        std::string_view key_tag,
        std::string_view label)
    {
        auto * data_manager{ GetDataManagerPtr() };
        const auto key{ std::string(key_tag) };
        try
        {
            data_manager->ProcessFile(path, key);
            return data_manager->GetTypedDataObject<TypedDataObject>(key);
        }
        catch (const std::exception & ex)
        {
            throw std::runtime_error(
                "Failed to process " + std::string(label) + " from '" + path.string()
                + "' as " + typeid(TypedDataObject).name() + ": " + ex.what());
        }
    }

    template <typename TypedDataObject>
    std::shared_ptr<TypedDataObject> OptionalProcessTypedFile(
        const std::filesystem::path & path,
        std::string_view key_tag,
        std::string_view label)
    {
        if (path.empty())
        {
            return nullptr;
        }
        return ProcessTypedFile<TypedDataObject>(path, key_tag, label);
    }

    template <typename TypedDataObject>
    std::shared_ptr<TypedDataObject> LoadTypedObject(
        std::string_view key_tag,
        std::string_view label)
    {
        auto * data_manager{ GetDataManagerPtr() };
        const auto key{ std::string(key_tag) };
        try
        {
            data_manager->LoadDataObject(key);
            return data_manager->GetTypedDataObject<TypedDataObject>(key);
        }
        catch (const std::exception & ex)
        {
            throw std::runtime_error(
                "Failed to load " + std::string(label) + " with key tag '"
                + key + "': " + ex.what());
        }
    }

private:
    bool m_deprecated_database_option_used{ false };
    bool m_is_prepared_for_execution{ false };

    void AddValidationIssue(
        const std::string & option_name,
        ValidationPhase phase,
        LogLevel level,
        const std::string & message,
        bool auto_corrected);
    void PrepareFilesystemTargets();
};

} // namespace rhbm_gem
