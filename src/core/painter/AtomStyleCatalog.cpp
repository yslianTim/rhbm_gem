#include "core/painter/AtomStyleCatalog.hpp"

#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <vector>

namespace rhbm_gem {

namespace {

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

} // namespace

StyleSpec AtomStyleCatalog::GetMainChainElementStyle(size_t member_id)
{
    if (IsValidMainChainMemberID(member_id) == false) return kUnknownStyle;
    return kMainChainElementStyleList.at(member_id);
}

StyleSpec AtomStyleCatalog::GetMainChainSpotStyle(Spot spot)
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

StyleSpec AtomStyleCatalog::GetMainChainStructureStyle(Structure structure)
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

} // namespace rhbm_gem
