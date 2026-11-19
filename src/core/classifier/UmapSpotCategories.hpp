#pragma once

#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>

namespace rhbm_gem::core::classifier {

struct ConfiguredUmapSpot
{
    std::string_view label;
    std::string_view spot;
    std::optional<std::string_view> residue;
};

inline constexpr std::array<ConfiguredUmapSpot, 5> kConfiguredUmapSpots{
    ConfiguredUmapSpot{ "C", "C", std::nullopt },
    ConfiguredUmapSpot{ "CA", "CA", std::nullopt },
    ConfiguredUmapSpot{ "N", "N", std::nullopt },
    ConfiguredUmapSpot{ "O", "O", std::nullopt },
    ConfiguredUmapSpot{ "O (Water)", "O", "HOH" },
};

inline std::optional<std::size_t> FindConfiguredUmapSpotIndex(
    std::string_view residue,
    std::string_view spot)
{
    // Match residue-specific categories before the generic spot categories.
    for (std::size_t index = 0; index < kConfiguredUmapSpots.size(); ++index)
    {
        const auto & configured_spot{ kConfiguredUmapSpots[index] };
        if (configured_spot.residue
            && residue == *configured_spot.residue
            && spot == configured_spot.spot)
        {
            return index;
        }
    }
    for (std::size_t index = 0; index < kConfiguredUmapSpots.size(); ++index)
    {
        const auto & configured_spot{ kConfiguredUmapSpots[index] };
        if (!configured_spot.residue && spot == configured_spot.spot) return index;
    }
    return std::nullopt;
}

inline std::string GetOtherUmapSpotLabel(Element element)
{
    return element == Element::UNK
        ? "Other Unknown"
        : "Side-Chain " + std::string(ChemicalDataHelper::GetLabel(element));
}

} // namespace rhbm_gem::core::classifier
