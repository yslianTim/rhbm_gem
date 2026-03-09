#include <rhbm_gem/data/AtomClassifier.hpp>
#include <rhbm_gem/utils/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/KeyPacker.hpp>
#include <rhbm_gem/data/AtomObject.hpp>
#include <rhbm_gem/utils/Logger.hpp>

#include <stdexcept>

namespace rhbm_gem {

const std::unordered_map<Spot, short> AtomClassifier::m_main_chain_color_map
{
    { Spot::CA, 633 }, // kRed+1
    { Spot::C, 881 },  // kViolet+1
    { Spot::N, 418 },  // kGreen+2
    { Spot::O, 862 }   // kAzure+2
};

const std::unordered_map<Spot, short> AtomClassifier::m_main_chain_solid_marker_map
{
    { Spot::CA, 21 },
    { Spot::C, 20 },
    { Spot::N, 22 },
    { Spot::O, 23 }
};

const std::unordered_map<Spot, short> AtomClassifier::m_main_chain_open_marker_map
{
    { Spot::CA, 25 },
    { Spot::C, 24 },
    { Spot::N, 26 },
    { Spot::O, 32 }
};

const std::unordered_map<Spot, short> AtomClassifier::m_main_chain_line_style_map
{
    { Spot::CA, 1 },
    { Spot::C, 2 },
    { Spot::N, 3 },
    { Spot::O, 4 }
};

const std::unordered_map<Spot, std::string> AtomClassifier::m_main_chain_spot_label_map
{
    { Spot::CA, "C_{#alpha}" },
    { Spot::C,  "C" },
    { Spot::N,  "N" },
    { Spot::O,  "O" }
};

const std::unordered_map<Structure, short> AtomClassifier::m_main_chain_structure_color_map
{
    { Structure::FREE,   633 },      // kRed+1
    { Structure::HELX_P, 418 },      // kGreen+2
    { Structure::SHEET,  862 }       // kAzure+2
};

const std::unordered_map<Structure, short> AtomClassifier::m_main_chain_structure_marker_map
{
    { Structure::FREE,   25 },
    { Structure::HELX_P, 26 },
    { Structure::SHEET,  32 }
};

const std::unordered_map<Structure, short> AtomClassifier::m_main_chain_structure_line_style_map
{
    { Structure::FREE,   1 },
    { Structure::HELX_P, 2 },
    { Structure::SHEET,  3 }
};

const std::unordered_map<Structure, std::string> AtomClassifier::m_main_chain_structure_label_map
{
    { Structure::FREE,   "N" },
    { Structure::HELX_P, "#alpha" },
    { Structure::SHEET,  "#beta" }
};

const std::vector<short> AtomClassifier::m_main_chain_member_color_list
{   // [Color defined in ROOT style]
    // kRed+1, kViolet+1, kGreen+2, kAzure+2 
    633, 881, 418, 862,
    111, 862, 862,
    862, 862, 862, 862,
    881, 881, 881,
    881, 881
};

const std::vector<short> AtomClassifier::m_main_chain_member_solid_marker_list
{   // [Marker defined in ROOT style]
    21, 20, 22, 23,
    33, 23, 23,
    23, 23, 23, 23,
    20, 20, 20,
    20, 20
};

const std::vector<short> AtomClassifier::m_main_chain_member_open_marker_list
{   // [Marker defined in ROOT style]
    25, 24, 26, 32,
    27, 32, 32,
    32, 32, 32, 32,
    24, 24, 24,
    24, 24
};

const std::vector<Element> AtomClassifier::m_main_chain_member_element_list
{
    Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
    Element::PHOSPHORUS, Element::OXYGEN, Element::OXYGEN,
    Element::OXYGEN, Element::OXYGEN, Element::OXYGEN, Element::OXYGEN,
    Element::CARBON, Element::CARBON, Element::CARBON,
    Element::CARBON, Element::CARBON
};

const std::vector<Spot> AtomClassifier::m_main_chain_member_spot_list
{
    Spot::CA, Spot::C, Spot::N, Spot::O,
    Spot::P, Spot::OP1, Spot::OP2,
    Spot::O5p, Spot::O4p, Spot::O3p, Spot::O2p,
    Spot::C5p, Spot::C4p, Spot::C3p,
    Spot::C2p, Spot::C1p
};

const std::vector<std::string> AtomClassifier::m_main_chain_member_title_list
{
    "Alpha Carbon", "Carbonyl Carbon", "Peptide Nitrogen", "Carbonyl Oxygen",
    "Phosphorus", "Phosphate Oxygen 1", "Phosphate Oxygen 2",
    "Sugar Oxygen 5'", "Sugar Oxygen 4'", "Sugar Oxygen 3'", "Sugar Oxygen 2'",
    "Sugar Carbon 5'", "Sugar Carbon 4'", "Sugar Carbon 3'",
    "Sugar Carbon 2'", "Sugar Carbon 1'"
};

const std::vector<std::string> AtomClassifier::m_main_chain_member_label_list
{
    "C_{#alpha}", "C", "N", "O",
    "P", "O1", "O2",
    "O5'", "O4'", "O3'", "O2'",
    "C5'", "C4'", "C3'",
    "C2'", "C1'"
};

AtomClassifier::AtomClassifier()
{

}

AtomClassifier::~AtomClassifier()
{

}

bool AtomClassifier::IsMainChainMember(Spot spot, size_t & main_chain_member_id)
{
    for (size_t i = 0; i < m_main_chain_member_count; i++)
    {
        if (m_main_chain_member_spot_list.at(i) == spot)
        {
            main_chain_member_id = i;
            return true;
        }
    }
    return false;
}

size_t AtomClassifier::GetMainChainMemberCount()
{
    return m_main_chain_member_count;
}

bool AtomClassifier::IsValidMainChainMemberID(size_t member_id)
{
    if (member_id >= m_main_chain_member_count)
    {
        Logger::Log(LogLevel::Error,
            "Invalid main chain member ID: " + std::to_string(member_id));
        return false;
    }
    return true;
}

short AtomClassifier::GetMainChainElementColor(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return m_main_chain_member_color_list.at(id);
}

short AtomClassifier::GetMainChainSpotColor(Spot spot)
{
    if (m_main_chain_color_map.find(spot) == m_main_chain_color_map.end())
    {
        return 1;
    }
    return m_main_chain_color_map.at(spot);
}

short AtomClassifier::GetMainChainStructureColor(Structure structure)
{
    if (m_main_chain_structure_color_map.find(structure) ==
        m_main_chain_structure_color_map.end())
    {
        return 1;
    }
    return m_main_chain_structure_color_map.at(structure);
}

short AtomClassifier::GetNucleotideMainChainElementColor(size_t id)
{
    id += m_nucleotide_main_chain_member_start_id;
    if (IsValidMainChainMemberID(id) == false) return 1;
    return m_main_chain_member_color_list.at(id);
}

short AtomClassifier::GetMainChainElementSolidMarker(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return m_main_chain_member_solid_marker_list.at(id);
}

short AtomClassifier::GetMainChainSpotSolidMarker(Spot spot)
{
    if (m_main_chain_solid_marker_map.find(spot) ==
        m_main_chain_solid_marker_map.end())
    {
        return 1;
    }
    return m_main_chain_solid_marker_map.at(spot);
}

short AtomClassifier::GetMainChainStructureMarker(Structure structure)
{
    if (m_main_chain_structure_marker_map.find(structure) ==
        m_main_chain_structure_marker_map.end())
    {
        return 1;
    }
    return m_main_chain_structure_marker_map.at(structure);
}

short AtomClassifier::GetMainChainElementOpenMarker(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return m_main_chain_member_open_marker_list.at(id);
}

short AtomClassifier::GetMainChainSpotOpenMarker(Spot spot)
{
    if (m_main_chain_open_marker_map.find(spot) ==
        m_main_chain_open_marker_map.end())
    {
        return 1;
    }
    return m_main_chain_open_marker_map.at(spot);
}

short AtomClassifier::GetMainChainSpotLineStyle(Spot spot)
{
    if (m_main_chain_line_style_map.find(spot) ==
        m_main_chain_line_style_map.end())
    {
        return 1;
    }
    return m_main_chain_line_style_map.at(spot);
}

short AtomClassifier::GetMainChainStructureLineStyle(Structure structure)
{
    if (m_main_chain_structure_line_style_map.find(structure) ==
        m_main_chain_structure_line_style_map.end())
    {
        return 1;
    }
    return m_main_chain_structure_line_style_map.at(structure);
}

Element AtomClassifier::GetMainChainElement(size_t id)
{
    if (IsValidMainChainMemberID(id) == false)  return Element::UNK;
    return m_main_chain_member_element_list.at(id);
}

Spot AtomClassifier::GetMainChainSpot(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return Spot::UNK;
    return m_main_chain_member_spot_list.at(id);
}

Spot AtomClassifier::GetNucleotideMainChainSpot(size_t id)
{
    id += m_nucleotide_main_chain_member_start_id;
    if (IsValidMainChainMemberID(id) == false) return Spot::UNK;
    return m_main_chain_member_spot_list.at(id);
}

const std::string & AtomClassifier::GetMainChainElementLabel(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return m_main_chain_member_label_list.at(0);
    return m_main_chain_member_label_list.at(id);
}

const std::string & AtomClassifier::GetMainChainSpotLabel(Spot spot)
{
    if (m_main_chain_spot_label_map.find(spot) ==
        m_main_chain_spot_label_map.end())
    {
        static const std::string default_label{ "UNK" };
        return default_label;
    }
    return m_main_chain_spot_label_map.at(spot);
}

const std::string & AtomClassifier::GetMainChainStructureLabel(Structure structure)
{
    if (m_main_chain_structure_label_map.find(structure) ==
        m_main_chain_structure_label_map.end())
    {
        static const std::string default_label{ "UNK" };
        return default_label;
    }
    return m_main_chain_structure_label_map.at(structure);
}

const std::string & AtomClassifier::GetMainChainElementTitle(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return m_main_chain_member_title_list.at(0);
    return m_main_chain_member_title_list.at(id);
}

GroupKey AtomClassifier::GetGroupKeyInClass(
    const AtomObject * atom_object, const std::string & class_key)
{
    if (class_key == ChemicalDataHelper::GetSimpleAtomClassKey())
    {
        return KeyPackerSimpleAtomClass::Pack(
            atom_object->GetAtomKey(),
            atom_object->GetSpecialAtomFlag());
    }
    else if (class_key == ChemicalDataHelper::GetComponentAtomClassKey())
    {
        return KeyPackerComponentAtomClass::Pack(
            atom_object->GetComponentKey(),
            atom_object->GetAtomKey());
    }
    else if (class_key == ChemicalDataHelper::GetStructureAtomClassKey())
    {
        return KeyPackerStructureAtomClass::Pack(
            atom_object->GetStructure(),
            atom_object->GetComponentKey(),
            atom_object->GetAtomKey());
    }
    else
    {
        throw std::runtime_error("Unsupported class key."+ class_key);
    }
}

GroupKey AtomClassifier::GetGroupKeyInClass(ComponentKey component_key, AtomKey atom_key)
{
    return KeyPackerComponentAtomClass::Pack(component_key, atom_key);
}

GroupKey AtomClassifier::GetGroupKeyInClass(
    Structure structure, ComponentKey component_key, AtomKey atom_key)
{
    return KeyPackerStructureAtomClass::Pack(structure, component_key, atom_key);
}

GroupKey AtomClassifier::GetMainChainSimpleAtomClassGroupKey(size_t id) const
{
    if (IsValidMainChainMemberID(id) == false)
    {
        Logger::Log(LogLevel::Warning, "Invalid main chain member ID: " + std::to_string(id));
        return 0;
    }
    auto atom_key{ static_cast<AtomKey>(m_main_chain_member_spot_list.at(id)) };
    return KeyPackerSimpleAtomClass::Pack(atom_key, false);
}

GroupKey AtomClassifier::GetMainChainComponentAtomClassGroupKey(size_t id, Residue residue) const
{
    if (IsValidMainChainMemberID(id) == false)
    {
        Logger::Log(LogLevel::Warning, "Invalid main chain member ID: " + std::to_string(id));
        return 0;
    }
    auto atom_key{ static_cast<AtomKey>(m_main_chain_member_spot_list.at(id)) };
    auto component_key{ static_cast<ComponentKey>(residue) };
    return KeyPackerComponentAtomClass::Pack(component_key, atom_key);
}

GroupKey AtomClassifier::GetMainChainStructureAtomClassGroupKey(
    size_t id, Structure structure, Residue residue) const
{
    if (IsValidMainChainMemberID(id) == false)
    {
        Logger::Log(LogLevel::Warning, "Invalid main chain member ID: " + std::to_string(id));
        return 0;
    }
    auto atom_key{ static_cast<AtomKey>(m_main_chain_member_spot_list.at(id)) };
    auto component_key{ static_cast<ComponentKey>(residue) };
    return KeyPackerStructureAtomClass::Pack(structure, component_key, atom_key);
}

std::vector<GroupKey> AtomClassifier::GetMainChainComponentAtomClassGroupKeyList(
    size_t id) const
{
    if (IsValidMainChainMemberID(id) == false) return {};
    std::vector<GroupKey> group_key_list;
    group_key_list.reserve(ChemicalDataHelper::GetStandardResidueCount());
    for (auto residue_id : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        group_key_list.emplace_back(GetMainChainComponentAtomClassGroupKey(id, residue_id));
    }
    for (auto residue_id : ChemicalDataHelper::GetStandardNucleotideList())
    {
        group_key_list.emplace_back(GetMainChainComponentAtomClassGroupKey(id, residue_id));
    }
    return group_key_list;
}

std::vector<GroupKey> AtomClassifier::GetMainChainStructureAtomClassGroupKeyList(
    size_t id, Structure structure) const
{
    if (IsValidMainChainMemberID(id) == false) return {};
    std::vector<GroupKey> group_key_list;
    group_key_list.reserve(ChemicalDataHelper::GetStandardResidueCount());
    for (auto residue_id : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        group_key_list.emplace_back(GetMainChainStructureAtomClassGroupKey(id, structure, residue_id));
    }
    for (auto residue_id : ChemicalDataHelper::GetStandardNucleotideList())
    {
        group_key_list.emplace_back(GetMainChainStructureAtomClassGroupKey(id, structure, residue_id));
    }
    return group_key_list;
}

GroupKey AtomClassifier::GetNucleotideMainChainComponentAtomClassGroupKey(
    size_t member_id, Residue residue) const
{
    member_id += m_nucleotide_main_chain_member_start_id;
    if (IsValidMainChainMemberID(member_id) == false) return 0;
    auto atom_key{ static_cast<AtomKey>(m_main_chain_member_spot_list.at(member_id)) };
    auto component_key{ static_cast<ComponentKey>(residue) };
    return KeyPackerComponentAtomClass::Pack(component_key, atom_key);
}

std::vector<GroupKey> AtomClassifier::GetNucleotideMainChainComponentAtomClassGroupKeyList(
    Residue component_id) const
{
    std::vector<GroupKey> group_key_list;
    group_key_list.reserve(m_nucleotide_main_chain_member_count);
    for (size_t i = 0; i < m_nucleotide_main_chain_member_count; i++)
    {
        group_key_list.emplace_back(
            GetNucleotideMainChainComponentAtomClassGroupKey(i, component_id));
    }
    return group_key_list;
}

} // namespace rhbm_gem
