#include "JointWidthTrustPolicy.hpp"

#include "core/detail/joint_component/StructuralPartition.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace joint_test::width_trust {
namespace {
namespace core=rhbm_gem::core::joint_component;
using Layout=rhbm_gem::JointParameterLayout;
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

Layout ComponentLayout(const Layout & layout,
    const core::ComponentView & view)
{
    Layout out;
    for(const auto atom:layout.full_atoms)
        if(view.LocalAtom(static_cast<Eigen::Index>(atom))>=0) out.full_atoms.push_back(atom);
    for(const auto row:layout.informative_rows)
        if(view.LocalRow(static_cast<Eigen::Index>(row))>=0) out.informative_rows.push_back(row);
    return out;
}

std::string BlockKey(const BlockScope & block)
{ return block.component_id+":"+std::to_string(block.block_index+1); }

void CountRejection(PolicyStatistics & statistics,const BlockScope & block,std::size_t component_index,
    std::size_t block_index,std::size_t sweep_index,const std::string & reason)
{
    ++statistics.in_scope_rejected_trials;
    ++statistics.rejections_by_reason[reason];
    ++statistics.rejections_by_block[BlockKey(block)];
    if(!statistics.first_rejection_sweep)
    {
        statistics.first_rejection_sweep=sweep_index+1;
        statistics.first_rejection_block=block_index+1;
        statistics.first_rejection_component=component_index;
    }
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

std::vector<BlockScope> LocateMutableCoreScope(const rhbm_gem::core::JointProblem & problem,
    const std::vector<std::string> & target_atom_ids,std::size_t core_atoms)
{
    if(target_atom_ids.empty() || core_atoms==0) throw std::invalid_argument("Invalid width-trust target scope.");
    const auto & data=rhbm_gem::core::JointProblemAccess::Get(problem);
    if(!data.input->selection_domain) throw std::invalid_argument("Width-trust targets are not selected targets.");

    std::set<std::string> unresolved(target_atom_ids.begin(),target_atom_ids.end());
    if(unresolved.size()!=target_atom_ids.size()) throw std::invalid_argument("Duplicate width-trust target atom id.");
    std::vector<std::size_t> parent_targets;
    parent_targets.reserve(target_atom_ids.size());
    for(const auto & id:target_atom_ids)
    {
        const auto found=std::find(data.input->atom_ids.begin(),data.input->atom_ids.end(),id);
        if(found==data.input->atom_ids.end()) throw std::invalid_argument("Width-trust target is absent from the problem: "+id);
        const auto index=static_cast<std::size_t>(found-data.input->atom_ids.begin());
        if(std::find(data.input->selection_domain->target_indices.begin(),
            data.input->selection_domain->target_indices.end(),index)==data.input->selection_domain->target_indices.end())
            throw std::invalid_argument("Width-trust target is not a selected target: "+id);
        parent_targets.push_back(index);
    }

    std::vector<BlockScope> scopes;
    for(const auto & view:data.partition.components)
    {
        bool contains_target=false;
        for(const auto target:parent_targets)
            contains_target|=view.LocalAtom(static_cast<Eigen::Index>(target))>=0;
        if(!contains_target) continue;

        const auto layout=ComponentLayout(data.layout,view);
        const auto prepared=core::MakeComponentProblemView(layout,view);
        core::StructuralTopology topology(prepared.domain,prepared.atom_ids,prepared.local_layout);
        const auto partition=topology.BuildCores(core_atoms,false);
        for(const auto target:parent_targets)
        {
            const auto & target_id=data.input->atom_ids[target];
            Eigen::Index local_target=-1;
            for(std::size_t atom=0;atom<prepared.atom_ids.size();++atom)
                if(prepared.atom_ids[atom]==target_id)
                {local_target=static_cast<Eigen::Index>(atom); break;}
            if(local_target<0) continue;
            for(std::size_t block_index=0;block_index<partition.cores.size();++block_index)
            {
                const auto & core_atoms_in_block=partition.cores[block_index].atoms;
                if(std::find(core_atoms_in_block.begin(),core_atoms_in_block.end(),local_target)==core_atoms_in_block.end())
                    continue;
                const auto already=std::find_if(scopes.begin(),scopes.end(),[&](const auto & block) {
                    return block.component_index==static_cast<std::size_t>(view.component_index) &&
                        block.block_index==block_index;
                });
                if(already!=scopes.end())
                {unresolved.erase(target_id); break;}
                BlockScope block;
                block.component_index=static_cast<std::size_t>(view.component_index);
                block.block_index=block_index;
                block.component_id=view.id;
                for(const auto atom:core_atoms_in_block)
                    block.core_atom_ids.push_back(prepared.atom_ids[static_cast<std::size_t>(atom)]);
                scopes.push_back(std::move(block));
                unresolved.erase(target_id);
                break;
            }
        }
    }
    if(!unresolved.empty()) throw std::invalid_argument("Width-trust targets have no mutable core: "+*unresolved.begin());
    std::sort(scopes.begin(),scopes.end(),[](const auto & left,const auto & right) {
        return std::tie(left.component_index,left.block_index)<std::tie(right.component_index,right.block_index);
    });
    return scopes;
}

Evaluation Evaluate(Policy policy,const Trial & candidate,const Trial & baseline)
{
    if(policy==Policy::Baseline) return {Decision::NoIntervention,{}};
    const auto unavailable=[](const std::string & reason) {return Evaluation{Decision::EvidenceUnavailable,reason};};
    if(!candidate.width_trust_evidence_unavailable_reason.empty())
        return unavailable(candidate.width_trust_evidence_unavailable_reason);
    if(!baseline.width_trust_evidence_unavailable_reason.empty())
        return unavailable(baseline.width_trust_evidence_unavailable_reason);
    if(candidate.tracked_atoms.empty() || baseline.tracked_atoms.empty())
        return unavailable("data-fields-incomplete");
    if(candidate.projected_width_singular_values.empty() || baseline.projected_width_singular_values.empty())
        return unavailable("spectrum-unavailable");
    const double candidate_projected_min=Minimum(candidate.projected_width_singular_values);
    const double baseline_projected_min=Minimum(baseline.projected_width_singular_values);
    if(!std::isfinite(candidate_projected_min) || !(candidate_projected_min>0) ||
        !std::isfinite(baseline_projected_min) || !(baseline_projected_min>0))
        return unavailable("spectrum-unavailable");

    std::unordered_map<std::string,const Atom *> baseline_atoms;
    for(const auto & atom:baseline.tracked_atoms) baseline_atoms.emplace(atom.atom_id,&atom);
    double amplification=1.0;
    for(const auto & atom:candidate.tracked_atoms)
    {
        if(!Find(baseline_atoms,atom.atom_id)) return unavailable("data-fields-incomplete");
        if(!std::isfinite(atom.local_design_condition) || !(atom.local_design_condition>0) ||
            !std::isfinite(atom.gaussian_column_norm) || !(atom.gaussian_column_norm>0) ||
            !std::isfinite(atom.charge_column_norm) || !(atom.charge_column_norm>0))
            return unavailable("data-fields-incomplete");
        amplification=std::max(amplification,atom.local_design_condition);
        if(policy==Policy::ConditioningAndCancellationAware)
        {
            if(!std::isfinite(atom.ac_cancellation_ratio) || !(atom.ac_cancellation_ratio>0))
                return unavailable("cancellation-evidence-unavailable");
            amplification=std::max(amplification,atom.local_design_condition*
                std::max(1.0,atom.ac_cancellation_ratio));
        }
    }
    if(!std::isfinite(amplification) || !std::isfinite(candidate.local_objective_before) ||
        !std::isfinite(candidate.objective_reduction))
        return unavailable("data-fields-incomplete");
    const double roundoff_budget=std::numeric_limits<double>::epsilon()*amplification*
        std::max(1.0,std::abs(candidate.local_objective_before));
    if(candidate.objective_reduction<=roundoff_budget)
        return {Decision::Reject,policy==Policy::ConditioningAndCancellationAware ?
            "width-trust:cancellation-amplified-objective-evidence" :
            "width-trust:conditioning-amplified-objective-evidence"};
    return {Decision::Accept,{}};
}

TrustCallback Make(Policy policy,const std::vector<BlockScope> & scope,std::size_t component_index,
    PolicyStatistics * statistics)
{
    if(policy==Policy::Baseline) return {};
    return [policy,scope,component_index,statistics](std::size_t block_index,std::size_t sweep_index,
        const Trial & candidate,const Trial & baseline) {
        const auto scoped=std::find_if(scope.begin(),scope.end(),[&](const auto & block) {
            return block.component_index==component_index && block.block_index==block_index;
        });
        if(scoped==scope.end())
        {
            if(statistics) ++statistics->out_of_scope_evaluated_trials;
            return std::string{};
        }
        if(statistics) ++statistics->in_scope_evaluated_trials;
        const auto evaluation=Evaluate(policy,candidate,baseline);
        if(evaluation.decision==Decision::EvidenceUnavailable)
        {
            if(statistics)
            {
                ++statistics->evidence_unavailable_trials;
                ++statistics->evidence_unavailable_by_reason[evaluation.reason];
            }
            return std::string{};
        }
        if(evaluation.decision!=Decision::Reject) return std::string{};
        if(statistics) CountRejection(*statistics,*scoped,component_index,block_index,sweep_index,evaluation.reason);
        return evaluation.reason;
    };
}

} // namespace joint_test::width_trust
