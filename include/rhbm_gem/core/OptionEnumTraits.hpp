#pragma once

#include <array>
#include <map>
#include <string>
#include <string_view>

#include <rhbm_gem/core/OptionEnumClass.hpp>

namespace rhbm_gem {

template <typename EnumType>
struct EnumOptionEntry
{
    std::string_view token;
    EnumType value;
};

template <typename EnumType>
struct EnumOptionTraits;

template <typename EnumType>
std::map<std::string, EnumType> BuildEnumCLIMap()
{
    std::map<std::string, EnumType> option_map;
    for (const auto & entry : EnumOptionTraits<EnumType>::kCliEntries)
    {
        option_map.emplace(std::string(entry.token), entry.value);
    }
    return option_map;
}

template <>
struct EnumOptionTraits<PainterType>
{
    inline static constexpr std::array<EnumOptionEntry<PainterType>, 10> kCliEntries{{
        { "0", PainterType::GAUS },
        { "gaus", PainterType::GAUS },
        { "1", PainterType::MODEL },
        { "model", PainterType::MODEL },
        { "2", PainterType::COMPARISON },
        { "comparison", PainterType::COMPARISON },
        { "3", PainterType::DEMO },
        { "demo", PainterType::DEMO },
        { "4", PainterType::ATOM },
        { "atom", PainterType::ATOM }
    }};
    inline static constexpr std::array<EnumOptionEntry<PainterType>, 5> kBindingEntries{{
        { "GAUS", PainterType::GAUS },
        { "ATOM", PainterType::ATOM },
        { "MODEL", PainterType::MODEL },
        { "COMPARISON", PainterType::COMPARISON },
        { "DEMO", PainterType::DEMO }
    }};
};

template <>
struct EnumOptionTraits<PrinterType>
{
    inline static constexpr std::array<EnumOptionEntry<PrinterType>, 8> kCliEntries{{
        { "0", PrinterType::ATOM_POSITION },
        { "atom_pos", PrinterType::ATOM_POSITION },
        { "1", PrinterType::MAP_VALUE },
        { "map", PrinterType::MAP_VALUE },
        { "2", PrinterType::GAUS_ESTIMATES },
        { "gaus", PrinterType::GAUS_ESTIMATES },
        { "3", PrinterType::ATOM_OUTLIER },
        { "atom_out", PrinterType::ATOM_OUTLIER }
    }};
    inline static constexpr std::array<EnumOptionEntry<PrinterType>, 4> kBindingEntries{{
        { "ATOM_POSITION", PrinterType::ATOM_POSITION },
        { "MAP_VALUE", PrinterType::MAP_VALUE },
        { "GAUS_ESTIMATES", PrinterType::GAUS_ESTIMATES },
        { "ATOM_OUTLIER", PrinterType::ATOM_OUTLIER }
    }};
};

template <>
struct EnumOptionTraits<PotentialModel>
{
    inline static constexpr std::array<EnumOptionEntry<PotentialModel>, 6> kCliEntries{{
        { "0", PotentialModel::SINGLE_GAUS },
        { "single", PotentialModel::SINGLE_GAUS },
        { "1", PotentialModel::FIVE_GAUS_CHARGE },
        { "five", PotentialModel::FIVE_GAUS_CHARGE },
        { "2", PotentialModel::SINGLE_GAUS_USER },
        { "user", PotentialModel::SINGLE_GAUS_USER }
    }};
    inline static constexpr std::array<EnumOptionEntry<PotentialModel>, 3> kBindingEntries{{
        { "SINGLE_GAUS", PotentialModel::SINGLE_GAUS },
        { "FIVE_GAUS_CHARGE", PotentialModel::FIVE_GAUS_CHARGE },
        { "SINGLE_GAUS_USER", PotentialModel::SINGLE_GAUS_USER }
    }};
};

template <>
struct EnumOptionTraits<PartialCharge>
{
    inline static constexpr std::array<EnumOptionEntry<PartialCharge>, 6> kCliEntries{{
        { "0", PartialCharge::NEUTRAL },
        { "neutral", PartialCharge::NEUTRAL },
        { "1", PartialCharge::PARTIAL },
        { "partial", PartialCharge::PARTIAL },
        { "2", PartialCharge::AMBER },
        { "amber", PartialCharge::AMBER }
    }};
    inline static constexpr std::array<EnumOptionEntry<PartialCharge>, 3> kBindingEntries{{
        { "NEUTRAL", PartialCharge::NEUTRAL },
        { "PARTIAL", PartialCharge::PARTIAL },
        { "AMBER", PartialCharge::AMBER }
    }};
};

template <>
struct EnumOptionTraits<TesterType>
{
    inline static constexpr std::array<EnumOptionEntry<TesterType>, 10> kCliEntries{{
        { "0", TesterType::BENCHMARK },
        { "benchmark", TesterType::BENCHMARK },
        { "1", TesterType::DATA_OUTLIER },
        { "data_outlier", TesterType::DATA_OUTLIER },
        { "2", TesterType::MEMBER_OUTLIER },
        { "member_outlier", TesterType::MEMBER_OUTLIER },
        { "3", TesterType::MODEL_ALPHA_DATA },
        { "alpha_data", TesterType::MODEL_ALPHA_DATA },
        { "4", TesterType::MODEL_ALPHA_MEMBER },
        { "alpha_member", TesterType::MODEL_ALPHA_MEMBER }
    }};
    inline static constexpr std::array<EnumOptionEntry<TesterType>, 5> kBindingEntries{{
        { "BENCHMARK", TesterType::BENCHMARK },
        { "DATA_OUTLIER", TesterType::DATA_OUTLIER },
        { "MEMBER_OUTLIER", TesterType::MEMBER_OUTLIER },
        { "MODEL_ALPHA_DATA", TesterType::MODEL_ALPHA_DATA },
        { "MODEL_ALPHA_MEMBER", TesterType::MODEL_ALPHA_MEMBER }
    }};
};

} // namespace rhbm_gem
