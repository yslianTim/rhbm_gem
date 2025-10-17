#include "BondClassifier.hpp"
#include "ChemicalDataHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "KeyPacker.hpp"
#include "AtomObject.hpp"
#include "BondObject.hpp"
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

const std::vector<Link> BondClassifier::m_main_chain_member_link_list
{
    Link::N_CA, Link::CA_C, Link::C_O, Link::C_N
};

const std::vector<std::string> BondClassifier::m_main_chain_member_title_list
{
    "N-CA Bond",
    "CA-C Bond",
    "C=O Bond",
    "C-N Bond"
};

const std::vector<std::string> BondClassifier::m_main_chain_member_label_list
{
    "N-C_{#alpha}",
    "C_{#alpha}-C",
    "C=O",
    "C-N"
};

BondClassifier::BondClassifier(void)
{

}

BondClassifier::~BondClassifier()
{

}

bool BondClassifier::IsMainChainMember(Link link, size_t & main_chain_member_id)
{
    for (size_t i = 0; i < m_main_chain_member_count; i++)
    {
        if (m_main_chain_member_link_list.at(i) == link)
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

short BondClassifier::GetMainChainMemberColor(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return m_main_chain_member_color_list.at(id);
}

short BondClassifier::GetMainChainMemberSolidMarker(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return m_main_chain_member_solid_marker_list.at(id);
}

short BondClassifier::GetMainChainMemberOpenMarker(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return 1;
    return m_main_chain_member_open_marker_list.at(id);
}

Link BondClassifier::GetMainChainLink(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return Link::UNK;
    return m_main_chain_member_link_list.at(id);
}

const std::string & BondClassifier::GetMainChainMemberLabel(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return m_main_chain_member_label_list.at(0);
    return m_main_chain_member_label_list.at(id);
}

const std::string & BondClassifier::GetMainChainMemberTitle(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return m_main_chain_member_title_list.at(0);
    return m_main_chain_member_title_list.at(id);
}

GroupKey BondClassifier::GetGroupKeyInClass(
    const BondObject * bond_object, const std::string & class_key)
{
    if (class_key == ChemicalDataHelper::GetSimpleBondClassKey())
    {
        return KeyPackerSimpleBondClass::Pack(
            bond_object->GetBondKey(),
            bond_object->GetSpecialBondFlag());
    }
    else if (class_key == ChemicalDataHelper::GetComponentBondClassKey())
    {
        return KeyPackerComponentBondClass::Pack(
            bond_object->GetComponentKey(),
            bond_object->GetBondKey());
    }
    else
    {
        throw std::runtime_error("Unsupported class key."+ class_key);
    }
}

GroupKey BondClassifier::GetGroupKeyInClass(ComponentKey component_key, BondKey bond_key)
{
    return KeyPackerComponentBondClass::Pack(component_key, bond_key);
}

GroupKey BondClassifier::GetMainChainSimpleBondClassGroupKey(size_t id)
{
    if (IsValidMainChainMemberID(id) == false)
    {
        Logger::Log(LogLevel::Warning, "Invalid main chain member ID: " + std::to_string(id));
        return 0;
    }
    auto bond_key{ static_cast<BondKey>(m_main_chain_member_link_list.at(id)) };
    return KeyPackerSimpleBondClass::Pack(bond_key, false);
}

GroupKey BondClassifier::GetMainChainComponentBondClassGroupKey(size_t id, Residue residue) const
{
    if (IsValidMainChainMemberID(id) == false)
    {
        Logger::Log(LogLevel::Warning, "Invalid main chain member ID: " + std::to_string(id));
        return 0;
    }
    auto bond_key{ static_cast<BondKey>(m_main_chain_member_link_list.at(id)) };
    auto component_key{ static_cast<ComponentKey>(residue) };
    return KeyPackerComponentBondClass::Pack(component_key, bond_key);
}

std::vector<GroupKey> BondClassifier::GetMainChainComponentBondClassGroupKeyList(
    size_t id) const
{
    if (IsValidMainChainMemberID(id) == false) return {};
    std::vector<GroupKey> group_key_list;
    group_key_list.reserve(ChemicalDataHelper::GetStandardResidueCount());
    for (auto residue_id : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        group_key_list.emplace_back(GetMainChainComponentBondClassGroupKey(id, residue_id));
    }
    for (auto residue_id : ChemicalDataHelper::GetStandardNucleotideList())
    {
        group_key_list.emplace_back(GetMainChainComponentBondClassGroupKey(id, residue_id));
    }
    return group_key_list;
}
