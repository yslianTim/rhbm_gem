#include "JointWidthTrustPolicy.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>

namespace joint_test::width_trust {
namespace {
using Trial=rhbm_gem::core::joint_component::FixedNeighborProfileTrial;
using Atom=Trial::AtomState;
const Atom * Find(const std::unordered_map<std::string,const Atom *> & values,const std::string & id)
{
    const auto found=values.find(id);
    return found==values.end() ? nullptr : found->second;
}

double Minimum(const std::vector<double> & values)
{
    if(values.empty()) return std::numeric_limits<double>::quiet_NaN();
    return *std::min_element(values.begin(),values.end());
}

std::string Evaluate(Policy policy,const Trial & candidate,const Trial & baseline)
{
    if(policy==Policy::Baseline) return {};
    (void)baseline;
    std::unordered_map<std::string,const Atom *> baseline_atoms;
    for(const auto & atom:baseline.tracked_atoms) baseline_atoms.emplace(atom.atom_id,&atom);
    if(candidate.tracked_atoms.empty() || candidate.projected_width_singular_values.empty() ||
        baseline.projected_width_singular_values.empty())
        return "width-trust:evidence-unavailable";
    const double candidate_projected_min=Minimum(candidate.projected_width_singular_values);
    if(!std::isfinite(candidate_projected_min) || !(candidate_projected_min>0))
        return "width-trust:projected-width-singular-value-unavailable";
    double amplification=1.0;
    for(const auto & atom:candidate.tracked_atoms)
    {
        const auto * reference=Find(baseline_atoms,atom.atom_id);
        if(!reference) return "width-trust:evidence-unavailable";
        if(!std::isfinite(atom.local_design_condition) || !(atom.local_design_condition>0) ||
            !std::isfinite(atom.gaussian_column_norm) || !(atom.gaussian_column_norm>0) ||
            !std::isfinite(atom.charge_column_norm) || !(atom.charge_column_norm>0))
            return "width-trust:raw-or-normalized-design-evidence-unavailable:"+atom.atom_id;
        amplification=std::max(amplification,atom.local_design_condition);
        if(policy==Policy::ConditioningAndCancellationAware)
        {
            if(!std::isfinite(atom.ac_cancellation_ratio) || !(atom.ac_cancellation_ratio>0))
                return "width-trust:ac-cancellation-evidence-unavailable:"+atom.atom_id;
            amplification=std::max(amplification,atom.local_design_condition*
                std::max(1.0,atom.ac_cancellation_ratio));
        }
    }
    const double roundoff_budget=std::numeric_limits<double>::epsilon()*amplification*
        std::max(1.0,std::abs(candidate.local_objective_before));
    if(!std::isfinite(candidate.objective_reduction) || candidate.objective_reduction<=roundoff_budget)
        return policy==Policy::ConditioningAndCancellationAware ?
            "width-trust:cancellation-amplified-objective-evidence" :
            "width-trust:conditioning-amplified-objective-evidence";
    return {};
}
}

const char * Name(Policy policy)
{
    switch(policy)
    {
    case Policy::Baseline:return "A-production-baseline";
    case Policy::ConditioningAware:return "B1-conditioning-aware";
    case Policy::ConditioningAndCancellationAware:return "B2-conditioning-and-cancellation-aware";
    }
    return "unknown";
}

std::function<std::string(std::size_t,std::size_t,const Trial &,const Trial &)> Make(Policy policy)
{
    if(policy==Policy::Baseline) return {};
    return [policy](std::size_t,std::size_t,const Trial & candidate,const Trial & baseline) {
        return Evaluate(policy,candidate,baseline);
    };
}

} // namespace joint_test::width_trust
