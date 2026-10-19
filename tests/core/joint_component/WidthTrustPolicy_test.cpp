#include <gtest/gtest.h>

#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/StructuralPartition.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include "support/JointWidthTrustPolicy.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace {
namespace core=rhbm_gem::core;
namespace joint=core::joint_component;
namespace trust=joint_test::width_trust;

core::JointProblem SelectedChain(std::size_t atoms)
{
    auto input=second_stage_test::SyntheticJointWorkload("chain",static_cast<int>(atoms));
    input.selection_domain.emplace();
    input.selection_domain->contributor_policy=std::string(rhbm_gem::PureHaloSupportExclusionContributorPolicy);
    input.selection_domain->target_indices.resize(atoms);
    std::iota(input.selection_domain->target_indices.begin(),input.selection_domain->target_indices.end(),0);
    input.selection_domain->pure_halo_support_exclusion=rhbm_gem::JointPureHaloSupportExclusionProvenance{
        input.observations.size(),input.observations.size(),0,0,0,0};
    return core::JointProblem(std::move(input));
}

trust::Trial Trial(double objective_reduction)
{
    trust::Trial trial;
    trial.local_objective_before=1.0;
    trial.objective_reduction=objective_reduction;
    trial.projected_width_singular_values={1.0};
    trust::Trial::AtomState atom;
    atom.atom_id="target";
    atom.local_design_condition=1.0e6;
    atom.gaussian_column_norm=1.0;
    atom.charge_column_norm=1.0;
    atom.ac_cancellation_ratio=1.0;
    trial.tracked_atoms.push_back(atom);
    return trial;
}

std::string TargetId(const core::JointProblem & problem)
{
    const auto & data=core::JointProblemAccess::Get(problem);
    return data.input->atom_ids.at(data.input->selection_domain->target_indices.at(1));
}

TEST(JointWidthTrustPolicyTest, MutableCoreScopeIsDeterministicAndExcludesNeighborOnlyBlocks)
{
    auto problem=SelectedChain(8);
    const auto target_id=TargetId(problem);
    const auto first=trust::LocateMutableCoreScope(problem,{target_id},2);
    const auto second=trust::LocateMutableCoreScope(problem,{target_id},2);
    ASSERT_EQ(first.size(),1u);
    ASSERT_EQ(second.size(),1u);
    EXPECT_EQ(first.front().component_index,second.front().component_index);
    EXPECT_EQ(first.front().component_id,second.front().component_id);
    EXPECT_EQ(first.front().block_index,second.front().block_index);
    EXPECT_EQ(first.front().core_atom_ids,second.front().core_atom_ids);
    EXPECT_NE(std::find(first.front().core_atom_ids.begin(),first.front().core_atom_ids.end(),target_id),
        first.front().core_atom_ids.end());

    const auto & data=core::JointProblemAccess::Get(problem);
    ASSERT_EQ(data.partition.components.size(),1u);
    const auto & view=data.partition.components.front();
    rhbm_gem::JointParameterLayout layout;
    for(const auto atom:data.layout.full_atoms)
        if(view.LocalAtom(static_cast<Eigen::Index>(atom))>=0) layout.full_atoms.push_back(atom);
    for(const auto row:data.layout.informative_rows)
        if(view.LocalRow(static_cast<Eigen::Index>(row))>=0) layout.informative_rows.push_back(row);
    const auto prepared=joint::MakeComponentProblemView(layout,view);
    joint::StructuralTopology topology(prepared.domain,prepared.atom_ids,prepared.local_layout);
    const auto partition=topology.BuildCores(2,true);
    std::size_t target_local=prepared.atom_ids.size();
    for(std::size_t atom=0;atom<prepared.atom_ids.size();++atom)
        if(prepared.atom_ids[atom]==target_id) target_local=atom;
    ASSERT_LT(target_local,prepared.atom_ids.size());
    const auto neighbor_only=std::find_if(partition.cores.begin(),partition.cores.end(),[&](const auto & block) {
        return std::find(block.atoms.begin(),block.atoms.end(),static_cast<Eigen::Index>(target_local))==block.atoms.end() &&
            std::find(block.neighbor_atoms.begin(),block.neighbor_atoms.end(),static_cast<Eigen::Index>(target_local))!=
                block.neighbor_atoms.end();
    });
    ASSERT_NE(neighbor_only,partition.cores.end());
    const auto neighbor_index=static_cast<std::size_t>(neighbor_only-partition.cores.begin());
    trust::PolicyStatistics statistics;
    const auto callback=trust::Make(trust::Policy::ConditioningAware,first,first.front().component_index,&statistics);
    auto rejected_candidate=Trial(0.0);
    EXPECT_TRUE(callback(neighbor_index,0,rejected_candidate,Trial(1.0)).empty());
    EXPECT_EQ(statistics.out_of_scope_evaluated_trials,1u);
    EXPECT_EQ(statistics.out_of_scope_interventions,0u);
    EXPECT_EQ(statistics.in_scope_rejected_trials,0u);
}

TEST(JointWidthTrustPolicyTest, InScopeBlocksEvaluateAndRecordConditioningRejections)
{
    const auto scope=std::vector<trust::BlockScope>{{.component_index=3,.block_index=4,.component_id="A"}};
    for(const auto [policy,expected_reason]:{
        std::pair{trust::Policy::ConditioningAware,"width-trust:conditioning-amplified-objective-evidence"},
        std::pair{trust::Policy::ConditioningAndCancellationAware,"width-trust:cancellation-amplified-objective-evidence"}})
    {
        trust::PolicyStatistics statistics;
        const auto callback=trust::Make(policy,scope,3,&statistics);
        const auto reason=callback(4,6,Trial(0.0),Trial(1.0));
        EXPECT_EQ(reason,expected_reason);
        EXPECT_EQ(statistics.in_scope_evaluated_trials,1u);
        EXPECT_EQ(statistics.in_scope_rejected_trials,1u);
        EXPECT_EQ(statistics.rejections_by_reason[reason],1u);
        EXPECT_EQ(statistics.rejections_by_block["A:5"],1u);
        ASSERT_TRUE(statistics.first_rejection_sweep);
        ASSERT_TRUE(statistics.first_rejection_block);
        EXPECT_EQ(*statistics.first_rejection_sweep,7u);
        EXPECT_EQ(*statistics.first_rejection_block,5u);
    }

    trust::PolicyStatistics accepted_statistics;
    const auto accepting=trust::Make(trust::Policy::ConditioningAware,scope,3,&accepted_statistics);
    EXPECT_TRUE(accepting(4,0,Trial(1.0),Trial(1.0)).empty());
    EXPECT_EQ(accepted_statistics.in_scope_evaluated_trials,1u);
    EXPECT_EQ(accepted_statistics.in_scope_rejected_trials,0u);
}

TEST(JointWidthTrustPolicyTest, MissingDiagnosticEvidenceIsRecordedWithoutVeto)
{
    const auto scope=std::vector<trust::BlockScope>{{.component_index=0,.block_index=0,.component_id="target"}};
    for(const auto * unavailable : {"capture-event-limit-exhausted","capture-atom-limit-exhausted",
        "row-limit","spectrum-unavailable","candidate-invalid","data-fields-incomplete"})
    {
        trust::PolicyStatistics statistics;
        const auto callback=trust::Make(trust::Policy::ConditioningAndCancellationAware,scope,0,&statistics);
        auto candidate=Trial(0.0);
        candidate.width_trust_evidence_unavailable_reason=unavailable;
        EXPECT_TRUE(callback(0,0,candidate,Trial(1.0)).empty())<<unavailable;
        EXPECT_EQ(statistics.evidence_unavailable_trials,1u)<<unavailable;
        EXPECT_EQ(statistics.evidence_unavailable_by_reason[unavailable],1u)<<unavailable;
        EXPECT_EQ(statistics.in_scope_rejected_trials,0u)<<unavailable;
    }
}

TEST(JointWidthTrustPolicyTest, ExhaustedCaptureBudgetPreservesProductionSearch)
{
    auto problem=SelectedChain(8);
    const auto target_id=TargetId(problem);
    const auto scope=trust::LocateMutableCoreScope(problem,{target_id},2);
    joint::FixedNeighborSearchPolicy baseline_policy;
    baseline_policy.core_atoms=2;
    baseline_policy.maximum_sweeps=1;
    const auto baseline=joint::FitFixedNeighborComponentsImpl(problem,std::vector<double>(8,.55),baseline_policy);

    joint::FixedNeighborWidthTrajectoryCapture capture;
    capture.atom_initial_widths.emplace_back(target_id,.55);
    capture.maximum_events=0;
    capture.maximum_atoms=32;
    capture.maximum_rows=4096;
    capture.policy_evidence_only=true;
    trust::PolicyStatistics statistics;
    joint::FixedNeighborSearchPolicy scoped_policy=baseline_policy;
    scoped_policy.configure_search_for_testing=[&](joint::FixedNeighborPolicy & local) {
        local.capture_local_trajectory=true;
        local.width_trajectory=&capture;
        local.width_trial_trust_for_testing=trust::Make(trust::Policy::ConditioningAware,scope,
            local.component_index,&statistics);
    };
    const auto scoped=joint::FitFixedNeighborComponentsImpl(problem,std::vector<double>(8,.55),scoped_policy);
    ASSERT_EQ(baseline.components.size(),1u);
    ASSERT_EQ(scoped.components.size(),1u);
    ASSERT_TRUE(baseline.components.front().candidate_diagnostic);
    ASSERT_TRUE(scoped.components.front().candidate_diagnostic);
    EXPECT_EQ(scoped.components.front().search_sweeps,baseline.components.front().search_sweeps);
    EXPECT_EQ(scoped.components.front().stop_reason,baseline.components.front().stop_reason);
    EXPECT_EQ(scoped.components.front().candidate_diagnostic->candidate_eta,
        baseline.components.front().candidate_diagnostic->candidate_eta);
    EXPECT_EQ(scoped.components.front().candidate_diagnostic->candidate_ac,
        baseline.components.front().candidate_diagnostic->candidate_ac);
    EXPECT_GT(capture.events_skipped_due_to_event_limit,0u);
    EXPECT_GT(statistics.evidence_unavailable_by_reason["capture-event-limit-exhausted"],0u);
    EXPECT_EQ(statistics.in_scope_rejected_trials,0u);
    EXPECT_EQ(statistics.out_of_scope_interventions,0u);
}

TEST(JointWidthTrustPolicyTest, BoundedAssemblyUsesProjectedWidthStationarity)
{
    auto problem=SelectedChain(3);
    const double eta_max=std::log(.4);
    joint::FixedNeighborSearchPolicy policy;
    policy.maximum_sweeps=200;
    policy.eta_upper_bound=eta_max;
    policy.configure_search_for_testing=[&](joint::FixedNeighborPolicy & local) {
        local.eta_upper_bound=eta_max;
    };
    const auto fit=joint::FitFixedNeighborComponentsImpl(problem,std::vector<double>(3,.4),policy);
    const auto stationarity=std::find_if(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="width-stationarity" &&
            check.scope==rhbm_gem::JointEvidenceScope::AssembledGlobal;
    });
    ASSERT_NE(stationarity,fit.evidence.end());
    EXPECT_EQ(stationarity->reason,"independent-bound-aware-projected-width-kkt");
    ASSERT_TRUE(stationarity->value.has_value());
    const auto bound_feasibility=std::find_if(fit.evidence.begin(),fit.evidence.end(),[](const auto & check) {
        return check.name=="width-bound-feasibility" &&
            check.scope==rhbm_gem::JointEvidenceScope::AssembledGlobal;
    });
    ASSERT_NE(bound_feasibility,fit.evidence.end());
    EXPECT_EQ(bound_feasibility->status,rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(fit.assembled_state.has_value());
    double projected_inf{},raw_inf{};
    bool active_bound{};
    for(std::size_t atom=0;atom<fit.assembled_state->log_b.size();++atom)
    {
        const double eta=fit.assembled_state->log_b[atom];
        const double gradient=fit.assembled_state->width_gradient[atom];
        projected_inf=std::max(projected_inf,std::abs(eta-std::min(eta-gradient,eta_max)));
        raw_inf=std::max(raw_inf,std::abs(gradient));
        active_bound|=eta>=eta_max-1e-12;
    }
    EXPECT_TRUE(active_bound);
    EXPECT_GT(raw_inf,policy.convergence_criteria.global_width_gradient);
    EXPECT_NEAR(*stationarity->value,projected_inf,1e-12);
    EXPECT_LE(*stationarity->value,policy.convergence_criteria.global_width_gradient);
    EXPECT_EQ(stationarity->status,rhbm_gem::JointCheckStatus::Passed);
}

TEST(JointWidthTrustPolicyTest, BaselinePolicyCreatesNoAcceptanceCallback)
{
    trust::PolicyStatistics statistics;
    const auto callback=trust::Make(trust::Policy::Baseline,{},0,&statistics);
    EXPECT_FALSE(callback);
}

} // namespace
