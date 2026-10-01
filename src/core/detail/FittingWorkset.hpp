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
    std::vector<bool> full_parameter_mask;
};

FittingWorkset MakeTwoStageFittingWorkset(ModelObject &);
FittingWorkset MakeJointFittingWorkset(ModelObject &, const JointProblem &);
std::vector<AtomObject *> CollectFirstStageAtoms(const FittingWorkset &);

} // namespace detail
} // namespace rhbm_gem::core
