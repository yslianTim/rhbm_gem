#pragma once

#include <memory>

#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>

namespace rhbm_gem {

class AtomObject;

struct ModelSpatialCache
{
    std::unique_ptr<::KDNode<AtomObject>> kd_tree_root;
};

} // namespace rhbm_gem
