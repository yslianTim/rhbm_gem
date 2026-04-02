#pragma once

#include <cstddef>

#include "core/painter/StyleSpec.hpp"
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

namespace rhbm_gem {

class AtomStyleCatalog
{
public:
    static StyleSpec GetMainChainElementStyle(size_t member_id);
    static StyleSpec GetMainChainSpotStyle(Spot spot);
    static StyleSpec GetMainChainStructureStyle(Structure structure);

    static short GetMainChainElementColor(size_t member_id)
    {
        return GetMainChainElementStyle(member_id).color;
    }
    static short GetMainChainSpotColor(Spot spot)
    {
        return GetMainChainSpotStyle(spot).color;
    }
    static short GetMainChainStructureColor(Structure structure)
    {
        return GetMainChainStructureStyle(structure).color;
    }
    static short GetMainChainElementSolidMarker(size_t member_id)
    {
        return GetMainChainElementStyle(member_id).solid_marker;
    }
    static short GetMainChainSpotSolidMarker(Spot spot)
    {
        return GetMainChainSpotStyle(spot).solid_marker;
    }
    static short GetMainChainStructureMarker(Structure structure)
    {
        return GetMainChainStructureStyle(structure).solid_marker;
    }
    static short GetMainChainElementOpenMarker(size_t member_id)
    {
        return GetMainChainElementStyle(member_id).open_marker;
    }
    static short GetMainChainSpotOpenMarker(Spot spot)
    {
        return GetMainChainSpotStyle(spot).open_marker;
    }
    static short GetMainChainSpotLineStyle(Spot spot)
    {
        return GetMainChainSpotStyle(spot).line_style;
    }
    static short GetMainChainStructureLineStyle(Structure structure)
    {
        return GetMainChainStructureStyle(structure).line_style;
    }
    static std::string GetMainChainElementLabel(size_t member_id)
    {
        return GetMainChainElementStyle(member_id).label;
    }
    static std::string GetMainChainSpotLabel(Spot spot)
    {
        return GetMainChainSpotStyle(spot).label;
    }
    static std::string GetMainChainStructureLabel(Structure structure)
    {
        return GetMainChainStructureStyle(structure).label;
    }
};

} // namespace rhbm_gem
