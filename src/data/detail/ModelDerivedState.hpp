#pragma once

#include <array>
#include <memory>
#include <optional>
#include <tuple>
#include <vector>

#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>

namespace rhbm_gem {

class AtomObject;
class ModelObject;

class ModelDerivedState
{
    std::unique_ptr<::KDNode<AtomObject>> m_kd_tree_root;
    std::optional<std::array<float, 3>> m_center_of_mass_position;
    std::array<std::optional<std::tuple<double, double>>, 3> m_model_position_range;

public:
    ModelDerivedState();
    ~ModelDerivedState();

    static ModelDerivedState & Of(ModelObject & model_object);
    static const ModelDerivedState & Of(const ModelObject & model_object);

    std::array<float, 3> GetCenterOfMassPosition(ModelObject & model_object);
    std::tuple<double, double> GetModelPositionRange(ModelObject & model_object, int axis);
    std::vector<AtomObject *> FindAtomsInRange(
        ModelObject & model_object,
        const AtomObject & center_atom,
        double range);
    void Clear();

private:
    void EnsureKDTreeRoot(ModelObject & model_object);
};

} // namespace rhbm_gem
