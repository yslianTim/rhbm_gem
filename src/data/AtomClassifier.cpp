#include "AtomClassifier.hpp"
#include "ChemicalDataHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "KeyPacker.hpp"
#include "AtomObject.hpp"
#include "Logger.hpp"

#include <stdexcept>

const std::vector<short> AtomClassifier::m_main_chain_member_color_list
{   // [Color defined in ROOT style]
    // kRed+1, kViolet+1, kGreen+2, kAzure+2 
    633, 881, 418, 862,
    111, 862, 862,
    862, 633, 881,
    862, 881, 862,
    881, 862, 881
};

const std::vector<short> AtomClassifier::m_main_chain_member_solid_marker_list
{   // [Marker defined in ROOT style]
    21, 20, 22, 23,
    33, 23, 23,
    23, 21, 20,
    23, 20, 23,
    20, 23, 20
};

const std::vector<short> AtomClassifier::m_main_chain_member_open_marker_list
{   // [Marker defined in ROOT style]
    25, 24, 26, 32,
    27, 32, 32,
    32, 25, 24,
    32, 24, 32,
    24, 32, 24
};

const std::vector<Element> AtomClassifier::m_main_chain_member_element_list
{
    Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN,
    Element::PHOSPHORUS, Element::OXYGEN, Element::OXYGEN,
    Element::OXYGEN, Element::CARBON, Element::CARBON,
    Element::OXYGEN, Element::CARBON, Element::OXYGEN,
    Element::CARBON, Element::OXYGEN, Element::CARBON
};

const std::vector<Spot> AtomClassifier::m_main_chain_member_spot_list
{
    Spot::CA, Spot::C, Spot::N, Spot::O,
    Spot::P, Spot::OP1, Spot::OP2,
    Spot::O5p, Spot::C5p, Spot::C4p,
    Spot::O4p, Spot::C3p, Spot::O3p,
    Spot::C2p, Spot::O2p, Spot::C1p
};

const std::vector<std::string> AtomClassifier::m_main_chain_member_title_list
{
    "Alpha Carbon", "Carbonyl Carbon", "Peptide Nitrogen", "Carbonyl Oxygen",
    "Phosphorus", "Phosphate Oxygen 1", "Phosphate Oxygen 2",
    "Sugar Oxygen 5'", "Sugar Carbon 5'", "Sugar Carbon 4'",
    "Sugar Oxygen 4'", "Sugar Carbon 3'", "Sugar Oxygen 3'",
    "Sugar Carbon 2'", "Sugar Oxygen 2'", "Sugar Carbon 1'"
};

const std::vector<std::string> AtomClassifier::m_main_chain_member_label_list
{
    "C_{#alpha}", "C", "N", "O",
    "P", "O1", "O2",
    "O5'", "C5'", "C4'",
    "O4'", "C3'", "O3'",
    "C2'", "O2'", "C1'"
};

AtomClassifier::AtomClassifier(void)
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

size_t AtomClassifier::GetMainChainMemberCount(void)
{
    return m_main_chain_member_count;
}

bool AtomClassifier::IsValidMainChainMemberID(size_t id)
{
    if (id >= m_main_chain_member_count)
    {
        Logger::Log(LogLevel::Error, "Invalid main chain member ID: " + std::to_string(id));
        return false;
    }
    return true;
}

short AtomClassifier::GetMainChainElementColor(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return m_main_chain_member_color_list.at(id);
}

short AtomClassifier::GetMainChainElementSolidMarker(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return m_main_chain_member_solid_marker_list.at(id);
}

short AtomClassifier::GetMainChainElementOpenMarker(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return m_main_chain_member_open_marker_list.at(id);
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

const std::string & AtomClassifier::GetMainChainElementLabel(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return m_main_chain_member_label_list.at(0);
    return m_main_chain_member_label_list.at(id);
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
    for (auto residue_id : ChemicalDataHelper::GetStandardResidueList())
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
    for (auto residue_id : ChemicalDataHelper::GetStandardResidueList())
    {
        group_key_list.emplace_back(GetMainChainStructureAtomClassGroupKey(id, structure, residue_id));
    }
    return group_key_list;
}
