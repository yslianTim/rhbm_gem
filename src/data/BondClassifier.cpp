#include "BondClassifier.hpp"
#include "AtomicInfoHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "KeyPacker.hpp"
#include "AtomObject.hpp"
#include "Logger.hpp"

#include <stdexcept>

const std::vector<short> BondClassifier::m_main_chain_member_color_list
{   // [Color defined in ROOT style]
    // kRed+1, kViolet+1, kGreen+2, kAzure+2 
    633,       881,      418,      862
};

const std::vector<short> BondClassifier::m_main_chain_member_solid_marker_list
{   // [Marker defined in ROOT style]
    21, 20, 22, 23
};

const std::vector<short> BondClassifier::m_main_chain_member_open_marker_list
{   // [Marker defined in ROOT style]
    25, 24, 26, 32
};

const std::vector<Element> BondClassifier::m_main_chain_member_element_list
{
    Element::CARBON, Element::CARBON, Element::NITROGEN, Element::OXYGEN
};

const std::vector<Spot> BondClassifier::m_main_chain_member_spot_list
{
    Spot::CA, Spot::C, Spot::N, Spot::O
};

const std::vector<std::string> BondClassifier::m_main_chain_member_title_list
{
    "Alpha Carbon",
    "Carbonyl Carbon",
    "Peptide Nitrogen",
    "Carbonyl Oxygen"
};

const std::vector<std::string> BondClassifier::m_main_chain_member_label_list
{
    "C_{#alpha}",
    "C",
    "N",
    "O"
};

BondClassifier::BondClassifier(void)
{

}

BondClassifier::~BondClassifier()
{

}

bool BondClassifier::IsMainChainMember(Spot spot, size_t & main_chain_member_id)
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

size_t BondClassifier::GetMainChainMemberCount(void)
{
    return m_main_chain_member_count;
}

bool BondClassifier::IsValidMainChainMemberID(size_t id)
{
    if (id >= m_main_chain_member_count)
    {
        Logger::Log(LogLevel::Error, "Invalid main chain member ID: " + std::to_string(id));
        return false;
    }
    return true;
}

short BondClassifier::GetMainChainElementColor(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return m_main_chain_member_color_list.at(id);
}

short BondClassifier::GetMainChainElementSolidMarker(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return m_main_chain_member_solid_marker_list.at(id);
}

short BondClassifier::GetMainChainElementOpenMarker(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return m_main_chain_member_open_marker_list.at(id);
}

Element BondClassifier::GetMainChainElement(size_t id)
{
    if (IsValidMainChainMemberID(id) == false)  return Element::UNK;
    return m_main_chain_member_element_list.at(id);
}

Spot BondClassifier::GetMainChainSpot(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return Spot::UNK;
    return m_main_chain_member_spot_list.at(id);
}

const std::string & BondClassifier::GetMainChainElementLabel(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return m_main_chain_member_label_list.at(0);
    return m_main_chain_member_label_list.at(id);
}

const std::string & BondClassifier::GetMainChainElementTitle(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return m_main_chain_member_title_list.at(0);
    return m_main_chain_member_title_list.at(id);
}

uint64_t BondClassifier::GetGroupKeyInClass(
    const BondObject * bond_object, const std::string & class_key)
{
    if (class_key == AtomicInfoHelper::GetElementClassKey())
    {
        //return KeyPackerElementClass::Pack(
        //    bond_object->GetAtomKey(),
        //    bond_object->GetSpecialAtomFlag());
    }
    else if (class_key == AtomicInfoHelper::GetResidueClassKey())
    {
        //return KeyPackerResidueClass::Pack(
        //    bond_object->GetComponentKey(),
        //    bond_object->GetAtomKey());
    }
    else if (class_key == AtomicInfoHelper::GetStructureClassKey())
    {
        //return KeyPackerStructureClass::Pack(
        //    bond_object->GetStructure(),
        //    bond_object->GetComponentKey(),
        //    bond_object->GetAtomKey());
    }
    else
    {
        throw std::runtime_error("Unsupported class key."+ class_key);
    }
}

uint64_t BondClassifier::GetMainChainElementClassGroupKey(size_t id) const
{
    if (IsValidMainChainMemberID(id) == false)
    {
        Logger::Log(LogLevel::Warning, "Invalid main chain member ID: " + std::to_string(id));
        return 0;
    }
    auto atom_key{ static_cast<AtomKey>(m_main_chain_member_spot_list.at(id)) };
    return KeyPackerElementClass::Pack(atom_key, false);
}

uint64_t BondClassifier::GetMainChainResidueClassGroupKey(size_t id, Residue residue) const
{
    if (IsValidMainChainMemberID(id) == false)
    {
        Logger::Log(LogLevel::Warning, "Invalid main chain member ID: " + std::to_string(id));
        return 0;
    }
    auto atom_key{ static_cast<AtomKey>(m_main_chain_member_spot_list.at(id)) };
    auto component_key{ static_cast<ComponentKey>(residue) };
    return KeyPackerResidueClass::Pack(component_key, atom_key);
}

uint64_t BondClassifier::GetMainChainStructureClassGroupKey(
    size_t id, Structure structure, Residue residue) const
{
    if (IsValidMainChainMemberID(id) == false)
    {
        Logger::Log(LogLevel::Warning, "Invalid main chain member ID: " + std::to_string(id));
        return 0;
    }
    auto atom_key{ static_cast<AtomKey>(m_main_chain_member_spot_list.at(id)) };
    auto component_key{ static_cast<ComponentKey>(residue) };
    return KeyPackerStructureClass::Pack(structure, component_key, atom_key);
}

std::vector<uint64_t> BondClassifier::GetMainChainResidueClassGroupKeyList(
    size_t id) const
{
    if (IsValidMainChainMemberID(id) == false) return {};
    std::vector<uint64_t> group_key_list;
    group_key_list.reserve(AtomicInfoHelper::GetStandardResidueCount());
    for (auto residue_id : AtomicInfoHelper::GetStandardResidueList())
    {
        group_key_list.emplace_back(GetMainChainResidueClassGroupKey(id, residue_id));
    }
    return group_key_list;
}

std::vector<uint64_t> BondClassifier::GetMainChainStructureClassGroupKeyList(
    size_t id, Structure structure) const
{
    if (IsValidMainChainMemberID(id) == false) return {};
    std::vector<uint64_t> group_key_list;
    group_key_list.reserve(AtomicInfoHelper::GetStandardResidueCount());
    for (auto residue_id : AtomicInfoHelper::GetStandardResidueList())
    {
        group_key_list.emplace_back(GetMainChainStructureClassGroupKey(id, structure, residue_id));
    }
    return group_key_list;
}
