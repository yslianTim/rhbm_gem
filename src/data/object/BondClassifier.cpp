#include "data/detail/BondClassifier.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/KeyPacker.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <stdexcept>

namespace rhbm_gem {

const std::vector<Link> BondClassifier::m_main_chain_member_link_list
{
    Link::N_CA, Link::CA_C, Link::C_O, Link::C_N
};

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

size_t BondClassifier::GetMainChainMemberCount()
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

Link BondClassifier::GetMainChainLink(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return Link::UNK;
    return m_main_chain_member_link_list.at(id);
}

GroupKey BondClassifier::GetGroupKeyInClass(
    const BondObject * bond_object,
    const std::string & class_key)
{
    if (class_key == ChemicalDataHelper::GetSimpleBondClassKey())
    {
        return KeyPackerSimpleBondClass::Pack(
            bond_object->GetBondKey(),
            bond_object->GetSpecialBondFlag());
    }
    if (class_key == ChemicalDataHelper::GetComponentBondClassKey())
    {
        return KeyPackerComponentBondClass::Pack(
            bond_object->GetAtomObject1()->GetComponentKey(),
            bond_object->GetBondKey());
    }
    throw std::runtime_error("Unsupported class key." + class_key);
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
    return KeyPackerSimpleBondClass::Pack(
        static_cast<BondKey>(m_main_chain_member_link_list.at(id)),
        false);
}

GroupKey BondClassifier::GetMainChainComponentBondClassGroupKey(size_t id, Residue residue)
{
    if (IsValidMainChainMemberID(id) == false)
    {
        Logger::Log(LogLevel::Warning, "Invalid main chain member ID: " + std::to_string(id));
        return 0;
    }
    const auto component_key{ static_cast<ComponentKey>(residue) };
    return KeyPackerComponentBondClass::Pack(
        component_key,
        static_cast<BondKey>(m_main_chain_member_link_list.at(id)));
}

std::vector<GroupKey> BondClassifier::GetMainChainComponentBondClassGroupKeyList(size_t id)
{
    if (IsValidMainChainMemberID(id) == false) return {};
    std::vector<GroupKey> group_key_list;
    group_key_list.reserve(ChemicalDataHelper::GetStandardResidueCount());
    for (const auto residue_id : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        group_key_list.emplace_back(GetMainChainComponentBondClassGroupKey(id, residue_id));
    }
    for (const auto residue_id : ChemicalDataHelper::GetStandardNucleotideList())
    {
        group_key_list.emplace_back(GetMainChainComponentBondClassGroupKey(id, residue_id));
    }
    return group_key_list;
}

} // namespace rhbm_gem
