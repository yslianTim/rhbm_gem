#pragma once

#include <filesystem>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

#include <CLI/CLI.hpp>

#include "OptionEnumTraits.hpp"

namespace rhbm_gem::command_cli {

template <typename ValueType, typename Setter>
CLI::Option * AddScalarOption(
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
CLI::Option * AddStringOption(
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
CLI::Option * AddPathOption(
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
CLI::Option * AddEnumOption(
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
        BuildEnumCLIMap<EnumType>(), CLI::ignore_case));
    return option;
}

} // namespace rhbm_gem::command_cli
