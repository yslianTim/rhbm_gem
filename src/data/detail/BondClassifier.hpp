#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <rhbm_gem/utils/domain/BondKeySystem.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

namespace rhbm_gem {

class BondObject;
class BondClassifier
{
    static const size_t m_main_chain_member_count{ 4 };
    static const std::vector<Link> m_main_chain_member_link_list;

public:
    static bool IsMainChainMember(Link link, size_t & main_chain_member_id);
    static size_t GetMainChainMemberCount();
    static Link GetMainChainLink(size_t id);

    static GroupKey GetGroupKeyInClass(const BondObject * bond_object, const std::string & class_key);
    static GroupKey GetGroupKeyInClass(ComponentKey component_key, BondKey bond_key);

    static GroupKey GetMainChainSimpleBondClassGroupKey(size_t id);
    static GroupKey GetMainChainComponentBondClassGroupKey(size_t id, Residue residue);
    static std::vector<GroupKey> GetMainChainComponentBondClassGroupKeyList(size_t id);

private:
    static bool IsValidMainChainMemberID(size_t id);
};

} // namespace rhbm_gem
