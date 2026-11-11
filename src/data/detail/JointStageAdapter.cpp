#include "JointStageAdapter.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <vector>

namespace rhbm_gem::data_internal {
std::map<int, LocalStageEstimate> BuildJointStageEstimates(
    const JointAnalysisResult & result, const std::string & run_id)
{
    if(result.parameterization_contract!=JointPureHaloSupportExclusionParameterizationContract ||
        !result.selection_domain ||
        result.selection_domain->contributor_policy!=PureHaloSupportExclusionContributorPolicy)
        throw std::invalid_argument("Invalid target-only Pure Exclusion result contract.");
    std::map<int, LocalStageEstimate> estimates;
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
        const auto & targets = result.selection_domain->target_indices;
        estimate.source.role = std::binary_search(targets.begin(), targets.end(), i) ?
            FittingRole::Target : FittingRole::Halo;
        estimate.reason = "missing-component-state";
        estimate.convergence = JointCheckStatus::Unavailable;
    }
    std::set<std::size_t> assigned;
    for (const auto & component : result.components)
    {
        const auto & full=component.layout ? component.layout->full_atoms : component.atoms;
        if (component.state && (component.state->ac.size() != 2 * full.size() ||
            component.state->tau.size() != full.size()))
            throw std::invalid_argument("Joint component state/mapping size mismatch.");
        std::map<std::size_t,std::size_t> state_indices;
        for (std::size_t k=0;k<full.size();++k) state_indices.emplace(full[k],k);
        std::map<std::size_t,const JointAtomPointValidation *> point_evidence;
        if (component.target_evidence)
            for (const auto & point:component.target_evidence->point_atoms)
                point_evidence.emplace(point.atom,&point);
        for (std::size_t local = 0; local < component.atoms.size(); ++local)
        {
            const auto global = component.atoms[local];
            if (global >= identities.size() || !assigned.insert(global).second)
                throw std::invalid_argument("Invalid or duplicate joint component mapping.");
            auto & estimate = estimates.at(identities[global]);
            estimate.source.component_id = component.id;
            const auto found=state_indices.find(global);
            if(found==state_indices.end())
                throw std::invalid_argument("Joint target component is missing its state mapping.");
            if(component.target_evidence)
            {
                const auto point=point_evidence.find(global);
                if(point==point_evidence.end() || point->second->parameter_validity!=JointCheckStatus::Passed ||
                    point->second->component_convergence!=JointCheckStatus::Passed)
                {
                    estimate.reason=point==point_evidence.end() ? "missing-point-validation" : point->second->reason;
                    estimate.convergence=JointCheckStatus::Unavailable;
                    continue;
                }
            }
            estimate.convergence = component.state ?
                (component.target_evidence ? component.target_runtime_convergence : component.runtime_convergence) :
                JointCheckStatus::Unavailable;
            if (!component.state) { estimate.reason = component.stop_reason; continue; }
            const auto & state = *component.state;
            const auto index=found->second;
            const double a = state.ac[2 * index], c = state.ac[2 * index + 1], tau = state.tau[index];
            if (!std::isfinite(a) || a < 0 || !std::isfinite(c) || !std::isfinite(tau) || tau <= 0)
                throw std::invalid_argument("Invalid joint point estimate.");
            estimate.point = GaussianModel3D{a, tau, c};
            estimate.reason.clear();
        }
    }
    return estimates;
}
}
