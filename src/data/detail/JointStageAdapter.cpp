#include "JointStageAdapter.hpp"
#include "ModelAnalysisData.hpp"
#include "LocalPotentialEntry.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace rhbm_gem::data_internal {
std::map<int, LocalStageEstimate> BuildJointStageEstimates(
    const JointAnalysisResult & result, const std::string & run_id)
{
    std::map<int, LocalStageEstimate> estimates;
    const bool fixed_target_only=result.parameterization_contract==JointTargetOnlyFixedHaloParameterizationContract;
    const bool pure_exclusion=result.parameterization_contract==JointPureHaloSupportExclusionParameterizationContract;
    if(fixed_target_only && (!result.selection_domain ||
        result.selection_domain->contributor_policy!=FixedHaloContributorPolicy || result.layout))
        throw std::invalid_argument("Invalid target-only fixed-halo result contract.");
    if(pure_exclusion && (!result.selection_domain ||
        result.selection_domain->contributor_policy!=PureHaloSupportExclusionContributorPolicy ||
        result.layout || !result.fixed_halos.empty()))
        throw std::invalid_argument("Invalid target-only Pure Exclusion result contract.");
    std::vector<int> identities;
    for (std::size_t i = 0; i < result.atom_ids.size(); ++i)
    {
        const auto & id = result.atom_ids[i];
        const int serial = std::stoi(id);
        if (std::to_string(serial) != id || estimates.contains(serial))
            throw std::invalid_argument("Invalid or duplicate joint atom identity.");
        identities.push_back(serial);
        auto & estimate = estimates[serial];
        estimate.source = {EstimateMethod::JointComponents, id, {}, run_id, FittingRole::NotRecorded};
        if (result.selection_domain)
        {
            const auto & targets = result.selection_domain->target_indices;
            const bool target=std::binary_search(targets.begin(), targets.end(), i);
            estimate.source.role = target ?
                FittingRole::Target : FittingRole::Halo;
        estimate.source.fixed_from_first_stage=fixed_target_only && !target;
        }
        estimate.reason = "missing-component-state";
        estimate.convergence = JointCheckStatus::Unavailable;
    }
    std::set<std::size_t> assigned;
    for (const auto & component : result.components)
    {
        const auto & full=component.layout ? component.layout->full_atoms : component.atoms;
        if (component.state && (component.state->ac.size() != 2 * full.size() ||
            component.state->b.size() != full.size()))
            throw std::invalid_argument("Joint component state/mapping size mismatch.");
        for (std::size_t local = 0; local < component.atoms.size(); ++local)
        {
            const auto global = component.atoms[local];
            if (global >= identities.size() || !assigned.insert(global).second)
                throw std::invalid_argument("Invalid or duplicate joint component mapping.");
            auto & estimate = estimates.at(identities[global]);
            estimate.source.component_id = component.id;
            const auto found=std::find(full.begin(),full.end(),global);
            if(found==full.end())
            {
                estimate.reason="observable-contribution-only";
                estimate.convergence=JointCheckStatus::Unavailable;
                continue;
            }
            if(component.target_evidence)
            {
                if(component.target_evidence->contract==JointTargetPointValidationContractVersion)
                {
                    const auto & points=component.target_evidence->point_atoms;
                    const auto point=std::find_if(points.begin(),points.end(),[&](const auto & value) {
                        return value.atom==global;
                    });
                    if(point==points.end() || point->parameter_validity!=JointCheckStatus::Passed ||
                        point->component_convergence!=JointCheckStatus::Passed)
                    {
                        estimate.reason=point==points.end() ? "missing-point-validation" : point->reason;
                        estimate.convergence=JointCheckStatus::Unavailable; continue;
                    }
                }
                else if(estimate.source.role==FittingRole::Halo)
                {
                    const auto & atoms=component.target_evidence->atoms;
                    const auto identified=std::find_if(atoms.begin(),atoms.end(),[&](const auto & a){return a.atom==global;});
                    if(identified==atoms.end() || identified->status!=JointCheckStatus::Passed)
                    {
                        estimate.reason=identified==atoms.end() ? "missing-identifiability" : identified->reason;
                        estimate.convergence=JointCheckStatus::Unavailable; continue;
                    }
                }
            }
            estimate.convergence = component.state ?
                (component.target_evidence && estimate.source.role==FittingRole::Target ? component.target_runtime_convergence : component.runtime_convergence) :
                JointCheckStatus::Unavailable;
            if (!component.state) { estimate.reason = component.stop_reason; continue; }
            const auto & state = *component.state;
            const auto index=static_cast<std::size_t>(found-full.begin());
            const double a = state.ac[2 * index], c = state.ac[2 * index + 1], b = state.b[index];
            if (!std::isfinite(a) || a < 0 || !std::isfinite(c) || !std::isfinite(b) || b <= 0)
                throw std::invalid_argument("Invalid joint point estimate.");
            estimate.point = GaussianModel3D{a, b, c};
            estimate.reason.clear();
        }
    }
    if(fixed_target_only)
    {
        std::size_t previous{};
        for(std::size_t k=0;k<result.fixed_halos.size();++k)
        {
            const auto & fixed=result.fixed_halos[k];
            if(fixed.contributor_index>=identities.size() ||
                (k && fixed.contributor_index<=previous) ||
                result.atom_ids[fixed.contributor_index]!=fixed.atom_id ||
                std::binary_search(result.selection_domain->target_indices.begin(),
                    result.selection_domain->target_indices.end(),fixed.contributor_index) ||
                !std::isfinite(fixed.amplitude) || fixed.amplitude<0 ||
                !std::isfinite(fixed.width) || fixed.width<=0 || fixed.charge!=0.0 ||
                !IsQualifiedFirstStageFixedHalo(fixed) ||
                fixed.source_contract!=FixedHaloParameterSourceContractVersion)
                throw std::invalid_argument("Invalid fixed FirstStage halo source.");
            previous=fixed.contributor_index;
            auto & estimate=estimates.at(identities[fixed.contributor_index]);
            estimate.point=GaussianModel3D{fixed.amplitude,fixed.width,0.0};
            estimate.reason="fixed-from-first-stage";
            estimate.convergence=JointCheckStatus::NotRun;
        }
        for(std::size_t atom=0;atom<identities.size();++atom)
        {
            if(std::binary_search(result.selection_domain->target_indices.begin(),
                result.selection_domain->target_indices.end(),atom)) continue;
            auto & estimate=estimates.at(identities[atom]);
            if(estimate.point) continue;
            const auto initialized=std::find_if(result.initialization.atoms.begin(),result.initialization.atoms.end(),
                [&](const auto & value){return value.id==result.atom_ids[atom];});
            estimate.reason=initialized!=result.initialization.atoms.end() && !initialized->reason.empty() ?
                "fixed-halo-initialization-unavailable:"+initialized->reason :
                "fixed-halo-initialization-unavailable";
            estimate.convergence=JointCheckStatus::Unavailable;
        }
    }
    return estimates;
}

void ApplyJointStageEstimates(ModelObject & model, const JointAnalysisResult & result, const std::string & run_id)
{
    const auto estimates = BuildJointStageEstimates(result, run_id);
    model.EditAnalysis().ApplySecondStageEstimates(estimates);
}
}
