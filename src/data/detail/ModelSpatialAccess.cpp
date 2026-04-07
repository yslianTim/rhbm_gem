#include "data/detail/ModelSpatialAccess.hpp"

#include "data/detail/ModelSpatialCache.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>

namespace rhbm_gem {

std::vector<AtomObject *> ModelSpatialAccess::FindAtomsInRange(
    ModelObject & model_object,
    const AtomObject & center_atom,
    double range)
{
    model_object.EnsureKDTreeRoot();
    if (model_object.m_spatial_cache == nullptr ||
        model_object.m_spatial_cache->kd_tree_root == nullptr)
    {
        return {};
    }

    return KDTreeAlgorithm<AtomObject>::RangeSearch(
        model_object.m_spatial_cache->kd_tree_root.get(),
        const_cast<AtomObject *>(&center_atom),
        range);
}

} // namespace rhbm_gem
