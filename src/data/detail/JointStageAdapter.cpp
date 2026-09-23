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
            estimate.source.role = std::binary_search(targets.begin(), targets.end(), i) ?
                FittingRole::Target : FittingRole::Halo;
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
            estimate.convergence = component.state ? component.runtime_convergence : JointCheckStatus::Unavailable;
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
    return estimates;
}

void ApplyJointStageEstimates(ModelObject & model, const JointAnalysisResult & result, const std::string & run_id)
{
    const auto estimates = BuildJointStageEstimates(result, run_id);
    model.EditAnalysis().ApplySecondStageEstimates(estimates);
}
}
