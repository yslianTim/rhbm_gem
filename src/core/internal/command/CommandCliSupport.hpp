#pragma once

#include <CLI/CLI.hpp>
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/core/command/OptionEnumClass.hpp>

#include <filesystem>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace rhbm_gem::command_cli {

template <typename ValueType, typename Setter>
inline CLI::Option * AddScalarOption(
    CLI::App * command,
    std::string_view flags,
    Setter && setter,
    std::string_view help,
    std::optional<ValueType> default_value = std::nullopt,
    bool required = false)
{
    auto * option{
        command->add_option_function<ValueType>(
            std::string(flags),
            [callback = std::forward<Setter>(setter)](const ValueType & value) mutable
            {
                std::invoke(callback, value);
            },
            std::string(help))
    };
    if (required)
    {
        option->required();
    }
    if (default_value.has_value())
    {
        option->default_val(*default_value);
    }
    return option;
}

template <typename Setter>
inline CLI::Option * AddStringOption(
    CLI::App * command,
    std::string_view flags,
    Setter && setter,
    std::string_view help,
    std::optional<std::string> default_value = std::nullopt,
    bool required = false)
{
    return AddScalarOption<std::string>(
        command,
        flags,
        std::forward<Setter>(setter),
        help,
        std::move(default_value),
        required);
}

template <typename Setter>
inline CLI::Option * AddPathOption(
    CLI::App * command,
    std::string_view flags,
    Setter && setter,
    std::string_view help,
    std::optional<std::filesystem::path> default_value = std::nullopt,
    bool required = false)
{
    auto * option{
        command->add_option_function<std::string>(
            std::string(flags),
            [callback = std::forward<Setter>(setter)](const std::string & value) mutable
            {
                std::invoke(callback, std::filesystem::path{ value });
            },
            std::string(help))
    };
    if (required)
    {
        option->required();
    }
    if (default_value.has_value())
    {
        option->default_val(default_value->string());
    }
    return option;
}

template <typename EnumType, typename Setter>
inline CLI::Option * AddEnumOption(
    CLI::App * command,
    std::string_view flags,
    Setter && setter,
    std::string_view help,
    std::optional<EnumType> default_value = std::nullopt,
    bool required = false)
{
    auto * option{
        AddScalarOption<EnumType>(
            command,
            flags,
            std::forward<Setter>(setter),
            help,
            std::move(default_value),
            required)
    };
    option->transform(CLI::CheckedTransformer(
        rhbm_gem::BuildEnumCLIMap<EnumType>(), CLI::ignore_case));
    return option;
}

inline void BindCommonOptions(
    CLI::App * command,
    CommonCommandRequest & common,
    CommonOptionProfile profile)
{
    const auto common_options{ CommonOptionMaskForProfile(profile) };
    if (HasCommonOption(common_options, CommonOption::Threading))
    {
        AddScalarOption<int>(
            command,
            "-j,--jobs",
            [&](int value) { common.thread_size = value; },
            "Number of threads",
            common.thread_size);
    }
    if (HasCommonOption(common_options, CommonOption::Verbose))
    {
        AddScalarOption<int>(
            command,
            "-v,--verbose",
            [&](int value) { common.verbose_level = value; },
            "Verbose level",
            common.verbose_level);
    }
    if (HasCommonOption(common_options, CommonOption::Database))
    {
        AddPathOption(
            command,
            "-d,--database",
            [&](const std::filesystem::path & value) { common.database_path = value; },
            "Database file path",
            common.database_path);
    }
    if (HasCommonOption(common_options, CommonOption::OutputFolder))
    {
        AddPathOption(
            command,
            "-o,--folder",
            [&](const std::filesystem::path & value) { common.folder_path = value; },
            "folder path for output files",
            common.folder_path);
    }
}

} // namespace rhbm_gem::command_cli
