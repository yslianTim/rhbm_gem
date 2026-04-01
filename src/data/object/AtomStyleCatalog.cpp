#include "data/object/AtomStyleCatalog.hpp"

#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <unordered_map>
#include <vector>

namespace rhbm_gem {

namespace {

const std::unordered_map<Spot, short> kMainChainColorMap
{
    { Spot::CA, 633 },
    { Spot::C, 881 },
    { Spot::N, 418 },
    { Spot::O, 862 }
};

const std::unordered_map<Spot, short> kMainChainSolidMarkerMap
{
    { Spot::CA, 21 },
    { Spot::C, 20 },
    { Spot::N, 22 },
    { Spot::O, 23 }
};

const std::unordered_map<Spot, short> kMainChainOpenMarkerMap
{
    { Spot::CA, 25 },
    { Spot::C, 24 },
    { Spot::N, 26 },
    { Spot::O, 32 }
};

const std::unordered_map<Spot, short> kMainChainLineStyleMap
{
    { Spot::CA, 1 },
    { Spot::C, 2 },
    { Spot::N, 3 },
    { Spot::O, 4 }
};

const std::unordered_map<Spot, std::string> kMainChainSpotLabelMap
{
    { Spot::CA, "C_{#alpha}" },
    { Spot::C, "C" },
    { Spot::N, "N" },
    { Spot::O, "O" }
};

const std::unordered_map<Structure, short> kMainChainStructureColorMap
{
    { Structure::FREE, 633 },
    { Structure::HELX_P, 418 },
    { Structure::SHEET, 862 }
};

const std::unordered_map<Structure, short> kMainChainStructureMarkerMap
{
    { Structure::FREE, 25 },
    { Structure::HELX_P, 26 },
    { Structure::SHEET, 32 }
};

const std::unordered_map<Structure, short> kMainChainStructureLineStyleMap
{
    { Structure::FREE, 1 },
    { Structure::HELX_P, 2 },
    { Structure::SHEET, 3 }
};

const std::unordered_map<Structure, std::string> kMainChainStructureLabelMap
{
    { Structure::FREE, "N" },
    { Structure::HELX_P, "#alpha" },
    { Structure::SHEET, "#beta" }
};

const std::vector<short> kMainChainMemberColorList
{
    633, 881, 418, 862,
    111, 862, 862,
    862, 862, 862, 862,
    881, 881, 881,
    881, 881
};

const std::vector<short> kMainChainMemberSolidMarkerList
{
    21, 20, 22, 23,
    33, 23, 23,
    23, 23, 23, 23,
    20, 20, 20,
    20, 20
};

const std::vector<short> kMainChainMemberOpenMarkerList
{
    25, 24, 26, 32,
    27, 32, 32,
    32, 32, 32, 32,
    24, 24, 24,
    24, 24
};

const std::vector<std::string> kMainChainMemberTitleList
{
    "Alpha Carbon", "Carbonyl Carbon", "Peptide Nitrogen", "Carbonyl Oxygen",
    "Phosphorus", "Phosphate Oxygen 1", "Phosphate Oxygen 2",
    "Sugar Oxygen 5'", "Sugar Oxygen 4'", "Sugar Oxygen 3'", "Sugar Oxygen 2'",
    "Sugar Carbon 5'", "Sugar Carbon 4'", "Sugar Carbon 3'",
    "Sugar Carbon 2'", "Sugar Carbon 1'"
};

const std::vector<std::string> kMainChainMemberLabelList
{
    "C_{#alpha}", "C", "N", "O",
    "P", "O1", "O2",
    "O5'", "O4'", "O3'", "O2'",
    "C5'", "C4'", "C3'",
    "C2'", "C1'"
};

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

template <typename T>
const T & RequireMappedValue(const std::unordered_map<Spot, T> & lookup, Spot spot, const T & fallback)
{
    const auto iter{ lookup.find(spot) };
    if (iter == lookup.end())
    {
        return fallback;
    }
    return iter->second;
}

template <typename T>
const T & RequireMappedValue(
    const std::unordered_map<Structure, T> & lookup,
    Structure structure,
    const T & fallback)
{
    const auto iter{ lookup.find(structure) };
    if (iter == lookup.end())
    {
        return fallback;
    }
    return iter->second;
}

} // namespace

short AtomStyleCatalog::GetMainChainElementColor(size_t member_id)
{
    if (IsValidMainChainMemberID(member_id) == false) return 1;
    return kMainChainMemberColorList.at(member_id);
}

short AtomStyleCatalog::GetMainChainSpotColor(Spot spot)
{
    return RequireMappedValue(kMainChainColorMap, spot, short{ 1 });
}

short AtomStyleCatalog::GetMainChainStructureColor(Structure structure)
{
    return RequireMappedValue(kMainChainStructureColorMap, structure, short{ 1 });
}

short AtomStyleCatalog::GetNucleotideMainChainElementColor(size_t member_id)
{
    member_id += 4;
    if (IsValidMainChainMemberID(member_id) == false) return 1;
    return kMainChainMemberColorList.at(member_id);
}

short AtomStyleCatalog::GetMainChainElementSolidMarker(size_t member_id)
{
    if (IsValidMainChainMemberID(member_id) == false) return 1;
    return kMainChainMemberSolidMarkerList.at(member_id);
}

short AtomStyleCatalog::GetMainChainSpotSolidMarker(Spot spot)
{
    return RequireMappedValue(kMainChainSolidMarkerMap, spot, short{ 1 });
}

short AtomStyleCatalog::GetMainChainStructureMarker(Structure structure)
{
    return RequireMappedValue(kMainChainStructureMarkerMap, structure, short{ 1 });
}

short AtomStyleCatalog::GetMainChainElementOpenMarker(size_t member_id)
{
    if (IsValidMainChainMemberID(member_id) == false) return 1;
    return kMainChainMemberOpenMarkerList.at(member_id);
}

short AtomStyleCatalog::GetMainChainSpotOpenMarker(Spot spot)
{
    return RequireMappedValue(kMainChainOpenMarkerMap, spot, short{ 1 });
}

short AtomStyleCatalog::GetMainChainSpotLineStyle(Spot spot)
{
    return RequireMappedValue(kMainChainLineStyleMap, spot, short{ 1 });
}

short AtomStyleCatalog::GetMainChainStructureLineStyle(Structure structure)
{
    return RequireMappedValue(kMainChainStructureLineStyleMap, structure, short{ 1 });
}

const std::string & AtomStyleCatalog::GetMainChainElementLabel(size_t member_id)
{
    if (IsValidMainChainMemberID(member_id) == false) return kMainChainMemberLabelList.at(0);
    return kMainChainMemberLabelList.at(member_id);
}

const std::string & AtomStyleCatalog::GetMainChainSpotLabel(Spot spot)
{
    static const std::string kFallback{ "UNK" };
    return RequireMappedValue(kMainChainSpotLabelMap, spot, kFallback);
}

const std::string & AtomStyleCatalog::GetMainChainStructureLabel(Structure structure)
{
    static const std::string kFallback{ "UNK" };
    return RequireMappedValue(kMainChainStructureLabelMap, structure, kFallback);
}

const std::string & AtomStyleCatalog::GetMainChainElementTitle(size_t member_id)
{
    if (IsValidMainChainMemberID(member_id) == false) return kMainChainMemberTitleList.at(0);
    return kMainChainMemberTitleList.at(member_id);
}

} // namespace rhbm_gem
