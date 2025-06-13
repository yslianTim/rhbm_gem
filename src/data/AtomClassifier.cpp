#include "AtomClassifier.hpp"
#include "AtomicInfoHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "KeyPacker.hpp"
#include "AtomObject.hpp"

#include <iostream>

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

const std::vector<Remoteness> AtomClassifier::m_main_chain_member_remoteness_list
{
    Remoteness::ALPHA, Remoteness::NONE, Remoteness::NONE, Remoteness::NONE
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

bool AtomClassifier::IsMainChainMember(
    Element element, Remoteness remoteness, size_t & main_chain_member_id)
{
    for (size_t i = 0; i < m_main_chain_member_count; i++)
    {
        if (m_main_chain_member_element_list.at(i) == element &&
            m_main_chain_member_remoteness_list.at(i) == remoteness)
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
        std::cerr << "Invalid id: " << id << std::endl;
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

Remoteness AtomClassifier::GetMainChainRemoteness(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return Remoteness::UNK;
    return m_main_chain_member_remoteness_list.at(id);
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

uint64_t AtomClassifier::GetGroupKeyInClass(
    const AtomObject * atom_object, const std::string & class_key)
{
    if (class_key == AtomicInfoHelper::GetElementClassKey())
    {
        return KeyPackerElementClass::Pack(
            atom_object->GetElement(),
            atom_object->GetRemoteness(),
            atom_object->GetSpecialAtomFlag());
    }
    else if (class_key == AtomicInfoHelper::GetResidueClassKey())
    {
        return KeyPackerResidueClass::Pack(
            atom_object->GetResidue(),
            atom_object->GetElement(),
            atom_object->GetRemoteness(),
            atom_object->GetBranch(),
            atom_object->GetSpecialAtomFlag());
    }
    else if (class_key == AtomicInfoHelper::GetStructureClassKey())
    {
        return KeyPackerStructureClass::Pack(
            atom_object->GetStructure(),
            atom_object->GetResidue(),
            atom_object->GetElement(),
            atom_object->GetRemoteness(),
            atom_object->GetBranch(),
            atom_object->GetSpecialAtomFlag());
    }
    else
    {
        throw std::runtime_error("Unsupported class key."+ class_key);
    }
}

uint64_t AtomClassifier::GetMainChainElementClassGroupKey(size_t id) const
{
    if (IsValidMainChainMemberID(id) == false) return 0;
    auto element_id{ m_main_chain_member_element_list.at(id) };
    auto remoteness_id{ m_main_chain_member_remoteness_list.at(id) };
    return KeyPackerElementClass::Pack(element_id, remoteness_id, false);
}

uint64_t AtomClassifier::GetMainChainResidueClassGroupKey(size_t id, Residue residue) const
{
    if (IsValidMainChainMemberID(id) == false) return 0;
    auto element_id{ m_main_chain_member_element_list.at(id) };
    auto remoteness_id{ m_main_chain_member_remoteness_list.at(id) };
    return KeyPackerResidueClass::Pack(residue, element_id, remoteness_id, Branch::NONE, false);
}

uint64_t AtomClassifier::GetMainChainStructureClassGroupKey(
    size_t id, Structure structure, Residue residue) const
{
    if (IsValidMainChainMemberID(id) == false) return 0;
    auto element_id{ m_main_chain_member_element_list.at(id) };
    auto remoteness_id{ m_main_chain_member_remoteness_list.at(id) };
    return KeyPackerStructureClass::Pack(structure, residue, element_id, remoteness_id, Branch::NONE, false);
}

std::vector<uint64_t> AtomClassifier::GetMainChainResidueClassGroupKeyList(
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

std::vector<uint64_t> AtomClassifier::GetMainChainStructureClassGroupKeyList(
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