#pragma once

#include <cstddef>
#include <vector>

#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

namespace rhbm_gem {

class AtomObject;

namespace data_internal {

bool IsMainChainMember(Spot spot, size_t & main_chain_member_id);
size_t GetMainChainMemberCount();
Element GetMainChainElement(size_t member_id);
Spot GetMainChainSpot(size_t member_id);

GroupKey GetGroupKey(const AtomObject * atom_object);
Residue GetResidueFromGroupKey(GroupKey group_key);

GroupKey GetMainChainGroupKey(size_t member_id, Residue residue);
std::vector<GroupKey> GetMainChainGroupKeyList(size_t member_id);

} // namespace data_internal

} // namespace rhbm_gem
