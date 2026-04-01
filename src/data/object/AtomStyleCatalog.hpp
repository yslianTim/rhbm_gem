#pragma once

#include <cstddef>
#include <string>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

namespace rhbm_gem {

class AtomStyleCatalog
{
public:
    static short GetMainChainElementColor(size_t member_id);
    static short GetMainChainSpotColor(Spot spot);
    static short GetMainChainStructureColor(Structure structure);
    static short GetNucleotideMainChainElementColor(size_t member_id);
    static short GetMainChainElementSolidMarker(size_t member_id);
    static short GetMainChainSpotSolidMarker(Spot spot);
    static short GetMainChainStructureMarker(Structure structure);
    static short GetMainChainElementOpenMarker(size_t member_id);
    static short GetMainChainSpotOpenMarker(Spot spot);
    static short GetMainChainSpotLineStyle(Spot spot);
    static short GetMainChainStructureLineStyle(Structure structure);
    static const std::string & GetMainChainElementLabel(size_t member_id);
    static const std::string & GetMainChainSpotLabel(Spot spot);
    static const std::string & GetMainChainStructureLabel(Structure structure);
    static const std::string & GetMainChainElementTitle(size_t member_id);
};

} // namespace rhbm_gem
