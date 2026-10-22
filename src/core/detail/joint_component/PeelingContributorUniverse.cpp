#include "PostFitPeeling.hpp"

#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "Problem.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

namespace rhbm_gem::core::detail {

JointPeelingContributorUniverse BuildJointPeelingContributorUniverse(
    const ModelObject & model,const JointProblem & problem,const JointAnalysisResult * result)
{
    JointPeelingContributorUniverse universe;
    if (result && result->halo_median_reconstruction)
        universe.joint_run_id=result->halo_median_reconstruction->joint_run_id;
    const auto & input=problem.OriginalInput();
    std::unordered_map<int,std::size_t> contributor_indices;
    contributor_indices.reserve(input.atom_ids.size());
    for (std::size_t index=0;index<input.atom_ids.size();++index)
        contributor_indices.emplace(std::stoi(input.atom_ids[index]),index);
    const auto * domain=input.selection_domain ? &*input.selection_domain : nullptr;
    universe.contributors.reserve(model.GetNumberOfAtom());
    for (const auto & object:model.GetAtomList())
    {
        JointPeelingContributor contributor;
        contributor.atom_id=object->GetSerialID();
        contributor.position=object->GetPosition();
        contributor.reconstruction_eligible=object->GetElement()!=Element::HYDROGEN;
        if (contributor.reconstruction_eligible && !std::all_of(
                contributor.position.begin(),contributor.position.end(),
                [](double coordinate){return std::isfinite(coordinate);}))
            throw std::invalid_argument("Nonfinite Joint Peeling contributor position.");
        if (const auto found=contributor_indices.find(contributor.atom_id);found!=contributor_indices.end())
        {
            contributor.contributor_index=found->second;
            contributor.target=domain && std::binary_search(domain->target_indices.begin(),
                domain->target_indices.end(),found->second);
        }
        if (universe.joint_run_id.empty() && contributor.target)
        {
            const auto local=AtomLocalPotentialView::For(*object);
            if (local.IsAvailable())
            {
                const auto stage=local.GetStageEstimate(FittingStage::Second);
                if (stage.source.method==EstimateMethod::JointComponents &&
                    stage.source.role==FittingRole::Target)
                    universe.joint_run_id=stage.source.run_id;
            }
        }
        universe.contributors.push_back(std::move(contributor));
    }
    std::sort(universe.contributors.begin(),universe.contributors.end(),
        [](const auto & left,const auto & right){return left.atom_id<right.atom_id;});
    return universe;
}

} // namespace rhbm_gem::core::detail
