#pragma once

#include <array>
#include <map>
#include <string>
#include <string_view>
#include <type_traits>

#include <rhbm_gem/core/command/CommandTypes.hpp>

namespace rhbm_gem::internal {

template <typename EnumType>
struct CommandEnumEntry
{
    std::string_view token;
    EnumType value;
};

template <typename EnumType, std::size_t AliasCount>
struct CommandEnumDefinition
{
    EnumType value;
    std::string_view binding_token;
    std::array<std::string_view, AliasCount> cli_aliases;
};

template <typename EnumType>
struct CommandEnumTraits;

template <typename EnumType, std::size_t AliasCount, std::size_t OptionCount>
constexpr std::array<CommandEnumEntry<EnumType>, OptionCount> MakeCommandEnumBindingEntries(
    const std::array<CommandEnumDefinition<EnumType, AliasCount>, OptionCount> & options)
{
    std::array<CommandEnumEntry<EnumType>, OptionCount> binding_entries{};
    for (std::size_t index = 0; index < OptionCount; ++index)
    {
        binding_entries[index] = { options[index].binding_token, options[index].value };
    }
    return binding_entries;
}

template <typename EnumType>
constexpr auto GetCommandEnumBindingEntries()
{
    return MakeCommandEnumBindingEntries(CommandEnumTraits<EnumType>::kOptions);
}

template <typename EnumType>
std::map<std::string, EnumType> BuildCommandEnumCliMap()
{
    std::map<std::string, EnumType> option_map;
    for (const auto & option : CommandEnumTraits<EnumType>::kOptions)
    {
        for (const auto alias : option.cli_aliases)
        {
            option_map.emplace(std::string(alias), option.value);
        }
    }
    return option_map;
}

template <typename EnumType>
constexpr bool IsSupportedCommandEnumValue(EnumType value)
{
    using UnderlyingType = std::underlying_type_t<EnumType>;
    const auto raw_numeric{ static_cast<UnderlyingType>(value) };
    for (const auto & option : CommandEnumTraits<EnumType>::kOptions)
    {
        if (static_cast<UnderlyingType>(option.value) == raw_numeric)
        {
            return true;
        }
    }
    return false;
}

template <>
struct CommandEnumTraits<PainterType>
{
    inline static constexpr std::array<CommandEnumDefinition<PainterType, 2>, 5> kOptions{{
        { PainterType::GAUS, "GAUS", { "0", "gaus" } },
        { PainterType::ATOM, "ATOM", { "4", "atom" } },
        { PainterType::MODEL, "MODEL", { "1", "model" } },
        { PainterType::COMPARISON, "COMPARISON", { "2", "comparison" } },
        { PainterType::DEMO, "DEMO", { "3", "demo" } }
    }};
};

template <>
struct CommandEnumTraits<PrinterType>
{
    inline static constexpr std::array<CommandEnumDefinition<PrinterType, 2>, 4> kOptions{{
        { PrinterType::ATOM_POSITION, "ATOM_POSITION", { "0", "atom_pos" } },
        { PrinterType::MAP_VALUE, "MAP_VALUE", { "1", "map" } },
        { PrinterType::GAUS_ESTIMATES, "GAUS_ESTIMATES", { "2", "gaus" } },
        { PrinterType::ATOM_OUTLIER, "ATOM_OUTLIER", { "3", "atom_out" } }
    }};
};

template <>
struct CommandEnumTraits<PotentialModel>
{
    inline static constexpr std::array<CommandEnumDefinition<PotentialModel, 2>, 3> kOptions{{
        { PotentialModel::SINGLE_GAUS, "SINGLE_GAUS", { "0", "single" } },
        { PotentialModel::FIVE_GAUS_CHARGE, "FIVE_GAUS_CHARGE", { "1", "five" } },
        { PotentialModel::SINGLE_GAUS_USER, "SINGLE_GAUS_USER", { "2", "user" } }
    }};
};

template <>
struct CommandEnumTraits<PartialCharge>
{
    inline static constexpr std::array<CommandEnumDefinition<PartialCharge, 2>, 3> kOptions{{
        { PartialCharge::NEUTRAL, "NEUTRAL", { "0", "neutral" } },
        { PartialCharge::PARTIAL, "PARTIAL", { "1", "partial" } },
        { PartialCharge::AMBER, "AMBER", { "2", "amber" } }
    }};
};

template <>
struct CommandEnumTraits<TesterType>
{
    inline static constexpr std::array<CommandEnumDefinition<TesterType, 2>, 6> kOptions{{
        { TesterType::BENCHMARK, "BENCHMARK", { "0", "benchmark" } },
        { TesterType::DATA_OUTLIER, "DATA_OUTLIER", { "1", "data_outlier" } },
        { TesterType::MEMBER_OUTLIER, "MEMBER_OUTLIER", { "2", "member_outlier" } },
        { TesterType::MODEL_ALPHA_DATA, "MODEL_ALPHA_DATA", { "3", "alpha_data" } },
        { TesterType::MODEL_ALPHA_MEMBER, "MODEL_ALPHA_MEMBER", { "4", "alpha_member" } },
        { TesterType::NEIGHBOR_DISTANCE, "NEIGHBOR_DISTANCE", { "5", "neighbor_distance" } }
    }};
};

} // namespace rhbm_gem::internal
