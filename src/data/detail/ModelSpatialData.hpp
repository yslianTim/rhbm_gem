#pragma once

#include <memory>
#include <vector>

#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>

namespace rhbm_gem {

class AtomObject;
class ModelObject;

class ModelSpatialData
{
    std::unique_ptr<::KDNode<AtomObject>> m_kd_tree_root;

public:
    ModelSpatialData();
    ~ModelSpatialData();

    static ModelSpatialData & Of(ModelObject & model_object);
    static const ModelSpatialData & Of(const ModelObject & model_object);
    static std::vector<AtomObject *> FindAtomsInRange(
        ModelObject & model_object,
        const AtomObject & center_atom,
        double range);

    void Clear();

private:
    void EnsureKDTreeRoot(ModelObject & model_object);
};

} // namespace rhbm_gem
