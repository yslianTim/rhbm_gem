#include "AtomClassifier.hpp"
#include "AtomicInfoHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "KeyPacker.hpp"
#include "AtomObject.hpp"
#include "Logger.hpp"

#include <stdexcept>

const std::vector<short> AtomClassifier::m_main_chain_member_color_list
{   // [Color defined in ROOT style]
    // kRed+1, kViolet+1, kGreen+2, kAzure+2 
    633,       881,      418,      862
};

const std::vector<short> AtomClassifier::m_main_chain_member_solid_marker_list
{   // [Marker defined in ROOT style]
    21, 20, 22, 23
};

const std::vector<short> AtomClassifier::m_main_chain_member_open_marker_list
{   // [Marker defined in ROOT style]
    25, 24, 26, 32
};

const std::vector<Element> AtomClassifier::m_main_chain_member_element_list
{
    Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN
};

const std::vector<Spot> AtomClassifier::m_main_chain_member_spot_list
{
    Spot::CA, Spot::C, Spot::N, Spot::O
};

const std::vector<std::string> AtomClassifier::m_main_chain_member_title_list
{
    "Alpha Carbon",
    "Carbonyl Carbon",
    "Peptide Nitrogen",
    "Carbonyl Oxygen"
};

const std::vector<std::string> AtomClassifier::m_main_chain_member_label_list
{
    "C_{#alpha}",
    "C",
    "N",
    "O"
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
    if (class_key == AtomicInfoHelper::GetSimpleAtomClassKey())
    {
        return KeyPackerSimpleAtomClass::Pack(
            atom_object->GetAtomKey(),
            atom_object->GetSpecialAtomFlag());
    }
    else if (class_key == AtomicInfoHelper::GetComponentAtomClassKey())
    {
        return KeyPackerComponentAtomClass::Pack(
            atom_object->GetComponentKey(),
            atom_object->GetAtomKey());
    }
    else if (class_key == AtomicInfoHelper::GetStructureAtomClassKey())
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
    group_key_list.reserve(AtomicInfoHelper::GetStandardResidueCount());
    for (auto residue_id : AtomicInfoHelper::GetStandardResidueList())
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
    group_key_list.reserve(AtomicInfoHelper::GetStandardResidueCount());
    for (auto residue_id : AtomicInfoHelper::GetStandardResidueList())
    {
        group_key_list.emplace_back(GetMainChainStructureAtomClassGroupKey(id, structure, residue_id));
    }
    return group_key_list;
}
