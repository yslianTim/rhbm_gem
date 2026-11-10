#pragma once

#include <vector>

namespace rhbm_gem {
class AtomObject;
class ModelObject;
}

namespace rhbm_gem::core {
class JointProblem;

namespace detail {

struct FittingWorkset
{
    std::vector<AtomObject *> contributors;
    std::vector<bool> target_mask;
};

FittingWorkset MakeJointFittingWorkset(ModelObject &, const JointProblem &);

} // namespace detail
} // namespace rhbm_gem::core
