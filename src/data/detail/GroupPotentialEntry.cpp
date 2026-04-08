#include "data/detail/GroupPotentialEntry.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>

namespace rhbm_gem {

template class GroupPotentialEntry<AtomObject>;
template class GroupPotentialEntry<BondObject>;

} // namespace rhbm_gem
