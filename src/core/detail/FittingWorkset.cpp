#include "FittingWorkset.hpp"

#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>

namespace rhbm_gem::core::detail {

FittingWorkset MakeTwoStageFittingWorkset(ModelObject & model)
{
    auto contributors = model.GetSelectedAtoms();
    return {
        contributors,
        std::vector<bool>(contributors.size(), true),
        std::vector<bool>(contributors.size(), true)
    };
}

FittingWorkset MakeJointFittingWorkset(ModelObject & model, const JointProblem & problem)
{
    const auto & input = problem.Input();
    if (!input.selection_domain) throw std::invalid_argument("Joint workset requires recorded selection.");
    FittingWorkset workset;
    for (std::size_t i = 0; i < input.atom_ids.size(); ++i)
    {
        auto * atom = model.FindAtomPtr(std::stoi(input.atom_ids[i]));
        if (!atom || std::to_string(atom->GetSerialID()) != input.atom_ids[i])
            throw std::invalid_argument("Unknown joint contributor identity.");
        workset.contributors.push_back(atom);
        const auto & targets = input.selection_domain->target_indices;
        workset.target_mask.push_back(std::binary_search(targets.begin(), targets.end(), i));
        const auto & full = problem.ParameterLayout().full_atoms;
        workset.full_parameter_mask.push_back(std::binary_search(full.begin(), full.end(), i));
    }
    return workset;
}

std::vector<AtomObject *> CollectFirstStageAtoms(const FittingWorkset & workset)
{
    if (workset.full_parameter_mask.size() != workset.contributors.size())
        throw std::invalid_argument("First-stage workset parameter count mismatch.");
    std::vector<AtomObject *> atoms;
    for (std::size_t i = 0; i < workset.contributors.size(); ++i)
        if (workset.full_parameter_mask[i]) atoms.push_back(workset.contributors[i]);
    return atoms;
}

} // namespace rhbm_gem::core::detail
