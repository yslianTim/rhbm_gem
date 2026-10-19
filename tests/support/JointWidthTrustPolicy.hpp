#pragma once

#include "core/detail/joint_component/Problem.hpp"

#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace joint_test::width_trust {

enum class Policy {Baseline,ConditioningAware,ConditioningAndCancellationAware};
enum class Decision {NoIntervention,Accept,Reject,EvidenceUnavailable};

using Trial=rhbm_gem::core::joint_component::FixedNeighborProfileTrial;
using TrustCallback=std::function<std::string(std::size_t,std::size_t,const Trial &,const Trial &)>;

struct BlockScope
{
    std::size_t component_index{},block_index{};
    std::string component_id;
    std::vector<std::string> core_atom_ids;
};

struct Evaluation
{
    Decision decision{Decision::EvidenceUnavailable};
    std::string reason;
};

struct PolicyStatistics
{
    std::size_t in_scope_evaluated_trials{},in_scope_rejected_trials{},out_of_scope_evaluated_trials{},
        out_of_scope_interventions{},evidence_unavailable_trials{};
    std::map<std::string,std::size_t> evidence_unavailable_by_reason,rejections_by_reason,rejections_by_block;
    std::optional<std::size_t> first_rejection_sweep,first_rejection_block,first_rejection_component;
};

const char * Name(Policy policy);
std::vector<BlockScope> LocateMutableCoreScope(const rhbm_gem::core::JointProblem & problem,
    const std::vector<std::string> & target_atom_ids,std::size_t core_atoms);
Evaluation Evaluate(Policy policy,const Trial & candidate,const Trial & baseline);
TrustCallback Make(Policy policy,const std::vector<BlockScope> & scope,std::size_t component_index,
    PolicyStatistics * statistics);

} // namespace joint_test::width_trust
