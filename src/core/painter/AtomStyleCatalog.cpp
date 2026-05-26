#include "core/painter/AtomStyleCatalog.hpp"

#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <string>
#include <vector>

namespace rhbm_gem {

namespace {

struct StyleSpec
{
    short color{ 1 };
    short solid_marker{ 1 };
    short open_marker{ 1 };
    short line_style{ 1 };
    std::string label{ "UNK" };
};

const std::vector<StyleSpec> kMainChainElementStyleList
{
    { 633, 21, 25, 1, "C_{#alpha}" },
    { 881, 20, 24, 1, "C" },
    { 418, 22, 26, 1, "N" },
    { 862, 23, 32, 1, "O" },
    { 111, 33, 27, 1, "P" },
    { 862, 23, 32, 1, "O1" },
    { 862, 23, 32, 1, "O2" },
    { 862, 23, 32, 1, "O5'" },
    { 862, 23, 32, 1, "O4'" },
    { 862, 23, 32, 1, "O3'" },
    { 862, 23, 32, 1, "O2'" },
    { 881, 20, 24, 1, "C5'" },
    { 881, 20, 24, 1, "C4'" },
    { 881, 20, 24, 1, "C3'" },
    { 881, 20, 24, 1, "C2'" },
    { 881, 20, 24, 1, "C1'" }
};

const StyleSpec kUnknownStyle{};

bool IsValidMainChainMemberID(size_t member_id)
{
    if (member_id >= AtomClassifier::GetMainChainMemberCount())
    {
        Logger::Log(
            LogLevel::Error,
            "Invalid main chain member ID: " + std::to_string(member_id));
        return false;
    }
    return true;
}

StyleSpec GetMainChainElementStyle(size_t member_id)
{
    if (IsValidMainChainMemberID(member_id) == false) return kUnknownStyle;
    return kMainChainElementStyleList.at(member_id);
}

StyleSpec GetMainChainSpotStyle(Spot spot)
{
    switch (spot)
    {
    case Spot::CA:
        return StyleSpec{ 633, 21, 25, 1, "C_{#alpha}" };
    case Spot::C:
        return StyleSpec{ 881, 20, 24, 2, "C" };
    case Spot::N:
        return StyleSpec{ 418, 22, 26, 3, "N" };
    case Spot::O:
        return StyleSpec{ 862, 23, 32, 4, "O" };
    default:
        return kUnknownStyle;
    }
}

StyleSpec GetMainChainStructureStyle(Structure structure)
{
    switch (structure)
    {
    case Structure::FREE:
        return StyleSpec{ 633, 25, 25, 1, "N" };
    case Structure::HELX_P:
        return StyleSpec{ 418, 26, 26, 2, "#alpha" };
    case Structure::SHEET:
        return StyleSpec{ 862, 32, 32, 3, "#beta" };
    default:
        return kUnknownStyle;
    }
}

} // namespace

short AtomStyleCatalog::GetMainChainElementColor(size_t member_id)
{
    return GetMainChainElementStyle(member_id).color;
}

short AtomStyleCatalog::GetMainChainSpotColor(Spot spot)
{
    return GetMainChainSpotStyle(spot).color;
}

short AtomStyleCatalog::GetMainChainStructureColor(Structure structure)
{
    return GetMainChainStructureStyle(structure).color;
}

short AtomStyleCatalog::GetMainChainElementSolidMarker(size_t member_id)
{
    return GetMainChainElementStyle(member_id).solid_marker;
}

short AtomStyleCatalog::GetMainChainSpotSolidMarker(Spot spot)
{
    return GetMainChainSpotStyle(spot).solid_marker;
}

short AtomStyleCatalog::GetMainChainStructureMarker(Structure structure)
{
    return GetMainChainStructureStyle(structure).solid_marker;
}

short AtomStyleCatalog::GetMainChainElementOpenMarker(size_t member_id)
{
    return GetMainChainElementStyle(member_id).open_marker;
}

short AtomStyleCatalog::GetMainChainSpotOpenMarker(Spot spot)
{
    return GetMainChainSpotStyle(spot).open_marker;
}

short AtomStyleCatalog::GetMainChainSpotLineStyle(Spot spot)
{
    return GetMainChainSpotStyle(spot).line_style;
}

short AtomStyleCatalog::GetMainChainStructureLineStyle(Structure structure)
{
    return GetMainChainStructureStyle(structure).line_style;
}

std::string AtomStyleCatalog::GetMainChainElementLabel(size_t member_id)
{
    return GetMainChainElementStyle(member_id).label;
}

std::string AtomStyleCatalog::GetMainChainSpotLabel(Spot spot)
{
    return GetMainChainSpotStyle(spot).label;
}

std::string AtomStyleCatalog::GetMainChainStructureLabel(Structure structure)
{
    return GetMainChainStructureStyle(structure).label;
}

} // namespace rhbm_gem
