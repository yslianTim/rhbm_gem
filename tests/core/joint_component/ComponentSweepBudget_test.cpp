#include <algorithm>
#include <gtest/gtest.h>
#include "core/detail/joint_component/ComponentSweepBudget.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointBcdDiagnostic.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include <limits>

namespace {
namespace n=rhbm_gem::core::joint_component;
void Observe(n::ComponentSweepBudget & budget,double ac,double width,bool recovery=false)
{budget.Observe(ac*1e-7,width*1e-9,recovery);}
void ImprovingWindow(n::ComponentSweepBudget & budget)
{
    for(int i=0;i<20;++i) Observe(budget,i<10 ? 6.0 : 3.0,2.0);
}
rhbm_gem::core::JointProblem CoupledPair()
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids={"1","2"}; input.support.resize(2);
    for(int z=-5;z<=5;++z) for(int y=-5;y<=5;++y) for(int x=-5;x<=7;++x)
    {
        const double first=.25*(x*x+y*y+z*z),second=std::pow(.5*x-1.0,2)+.25*(y*y+z*z);
        if(first>6.25 && second>6.25) continue;
        const auto row=input.observations.size();
        double value=1e-5*std::sin(.19*x+.17*y+.13*z);
        for(std::size_t atom=0;atom<2;++atom)
        {
            const double square=atom ? second : first;
            if(square>6.25) continue;
            input.support[atom].push_back({row,square});
            const auto basis=n::EvaluateKernel(square,.5,2.5);
            value+=2*basis.gaussian+.2*basis.charge;
        }
        input.row_ids.push_back(std::to_string(row)); input.observations.push_back(value);
    }
    input.selection_domain.emplace(); input.selection_domain->target_indices={0,1};
    input.selection_domain->contributor_policy=std::string(rhbm_gem::PureHaloSupportExclusionContributorPolicy);
    input.selection_domain->pure_halo_support_exclusion=rhbm_gem::JointPureHaloSupportExclusionProvenance{
        input.observations.size(),input.observations.size(),0,0,0,0};
    return rhbm_gem::core::JointProblem(std::move(input));
}
}

TEST(JointComponentSweepBudgetTest, NominalSweepsAreAvailableWithoutSpeculativeExtension)
{
    n::ComponentSweepBudget budget(200,300,n::PracticalJointConvergenceCriteria);
    for(int i=0;i<200;++i)
    {ASSERT_TRUE(budget.AllowsNextSweep()); Observe(budget,5.0,2.0);}
    EXPECT_FALSE(budget.AllowsNextSweep());
    EXPECT_FALSE(budget.HardLimitReached());
}

TEST(JointComponentSweepBudgetTest, NearbyImprovementExtendsOnlyToTheHardLimit)
{
    n::ComponentSweepBudget budget(20,30,n::PracticalJointConvergenceCriteria);
    ImprovingWindow(budget);
    for(int i=20;i<30;++i)
    {ASSERT_TRUE(budget.AllowsNextSweep()); Observe(budget,3.0,2.0);}
    EXPECT_FALSE(budget.AllowsNextSweep());
    EXPECT_TRUE(budget.HardLimitReached());
}

TEST(JointComponentSweepBudgetTest, ReboundsAndOddEvenOscillationDoNotEraseObservedProgress)
{
    for(bool oscillation:{false,true})
    {
        n::ComponentSweepBudget budget(20,30,n::PracticalJointConvergenceCriteria);
        for(int i=0;i<20;++i)
        {
            double ac=i<10 ? 6.0 : 3.0;
            if((oscillation && i%2) || i==19) ac=20.0;
            Observe(budget,ac,2.0);
        }
        EXPECT_TRUE(budget.AllowsNextSweep());
    }
}

TEST(JointComponentSweepBudgetTest, WidthProgressUsesItsOwnNormalizedThreshold)
{
    n::ComponentSweepBudget budget(20,30,n::PracticalJointConvergenceCriteria);
    for(int i=0;i<20;++i) Observe(budget,.5,i<10 ? 6.0 : 3.0);
    EXPECT_TRUE(budget.AllowsNextSweep());
    n::ComponentSweepBudget stalled_width(20,30,n::PracticalJointConvergenceCriteria);
    for(int i=0;i<20;++i) Observe(stalled_width,i<10 ? 6.0 : 3.0,8.0);
    EXPECT_FALSE(stalled_width.AllowsNextSweep());
}

TEST(JointComponentSweepBudgetTest, DistantEndpointsDoNotConsumeAnExtension)
{
    n::ComponentSweepBudget budget(20,30,n::PracticalJointConvergenceCriteria);
    for(int i=0;i<20;++i) Observe(budget,i<10 ? 60.0 : 30.0,2.0);
    EXPECT_FALSE(budget.AllowsNextSweep());
}

TEST(JointComponentSweepBudgetTest, NonfiniteAndNegativeMetricsCannotAuthorizeExtension)
{
    for(double bad:{std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity(),-1.0})
        for(bool width:{false,true})
        {
            n::ComponentSweepBudget budget(20,30,n::PracticalJointConvergenceCriteria);
            for(int i=0;i<20;++i)
                Observe(budget,i==15 && !width ? bad : i<10 ? 6.0 : 3.0,
                    i==15 && width ? bad : 2.0);
            EXPECT_FALSE(budget.AllowsNextSweep());
        }
}

TEST(JointComponentSweepBudgetTest, UnresolvedRecoveryBlocksExtensionButRecoveredHistoryDoesNot)
{
    for(bool pending:{false,true})
    {
        n::ComponentSweepBudget budget(20,30,n::PracticalJointConvergenceCriteria);
        for(int i=0;i<20;++i) Observe(budget,i<10 ? 6.0 : 3.0,2.0,i==0 || (pending && i==19));
        EXPECT_EQ(budget.AllowsNextSweep(),!pending);
    }
}

TEST(JointComponentSweepBudgetTest, IndependentComponentsCannotShareBudgetHistory)
{
    n::ComponentSweepBudget improving(20,30,n::PracticalJointConvergenceCriteria);
    n::ComponentSweepBudget stalled(20,30,n::PracticalJointConvergenceCriteria);
    ImprovingWindow(improving);
    for(int i=0;i<20;++i) Observe(stalled,3.0,2.0);
    EXPECT_TRUE(improving.AllowsNextSweep());
    EXPECT_FALSE(stalled.AllowsNextSweep());
}

TEST(JointComponentSweepBudgetTest, FixedAndZeroBudgetsDoNotExtend)
{
    n::ComponentSweepBudget fixed(0,20,n::PracticalJointConvergenceCriteria);
    ImprovingWindow(fixed);
    EXPECT_FALSE(fixed.AllowsNextSweep());
    EXPECT_FALSE(fixed.HardLimitReached());
    n::ComponentSweepBudget empty(0,0,n::PracticalJointConvergenceCriteria);
    EXPECT_FALSE(empty.AllowsNextSweep());
}

TEST(JointComponentSweepBudgetTest, EasyQualifiedSearchStopsWithoutExtraWorkOrParameterChanges)
{
    const rhbm_gem::core::JointProblem problem(second_stage_test::SyntheticJointWorkload("chain",2));
    const n::Vector eta=n::Vector::Constant(2,std::log(.55));
    n::BcdSearchPolicy policy;
    policy.maximum_sweeps=30; policy.convergence_criteria=n::PracticalJointConvergenceCriteria;
    policy.eta_upper_bound=std::log(2.2);
    const auto fixed=second_stage_test::SearchBcdFixture(problem,eta,policy);
    policy.nominal_sweeps=20;
    const auto adaptive=second_stage_test::SearchBcdFixture(problem,eta,policy);
    ASSERT_TRUE(adaptive.search_converged);
    EXPECT_LT(adaptive.sweep_count,20u);
    EXPECT_EQ(adaptive.sweep_count,fixed.sweep_count);
    EXPECT_EQ(adaptive.state.eta,fixed.state.eta);
    EXPECT_EQ(adaptive.state.beta,fixed.state.beta);
    EXPECT_EQ(adaptive.state.prediction,fixed.state.prediction);
    EXPECT_DOUBLE_EQ(adaptive.state.objective,fixed.state.objective);
    ASSERT_TRUE(adaptive.certification)<<adaptive.reason;
    EXPECT_TRUE(adaptive.certification->assessment.independent_replay.passed);
    EXPECT_TRUE(adaptive.certification->qualified());
    EXPECT_TRUE(adaptive.accepted_state_integrity);
}

TEST(JointComponentSweepBudgetTest, ExtensionPreservesTheOriginalPathAndQualifiesProductionAssembly)
{
    const auto problem=CoupledPair();
    n::BcdExecutionPolicy policy; policy.eta_upper_bound=std::log(2.2);
    n::BcdSearchPolicy search;
    search.core_atoms=1; search.maximum_sweeps=40;
    search.eta_upper_bound=std::log(1.0);
    search.component_id="search-only-id";
    search.production_progress=[](const auto &) {ADD_FAILURE()<<"Search callback escaped into component execution";};
    policy.search_policy_for_testing=&search;
    const auto fixed=n::FitBcdComponents(problem,{1.5,1.5},policy);
    search.nominal_sweeps=30;
    const auto adaptive=n::FitBcdComponents(problem,{1.5,1.5},policy);
    ASSERT_EQ(fixed.components.size(),1u);
    ASSERT_EQ(adaptive.components.size(),1u);
    const auto & fixed_component=fixed.components.front();
    const auto & adaptive_component=adaptive.components.front();
    ASSERT_GT(adaptive_component.search_sweeps,30u);
    EXPECT_EQ(adaptive_component.search_sweeps,fixed_component.search_sweeps);
    ASSERT_TRUE(adaptive_component.state && fixed_component.state);
    ASSERT_TRUE(adaptive_component.candidate_diagnostic);
    EXPECT_EQ(adaptive_component.candidate_diagnostic->convergence_criteria_id,"practical-v2");
    ASSERT_FALSE(adaptive_component.local_conditioning.empty());
    EXPECT_EQ(adaptive_component.local_conditioning.front().block_id,adaptive_component.id+"/block-1");
    EXPECT_EQ(adaptive_component.state->log_tau,fixed_component.state->log_tau);
    EXPECT_EQ(adaptive_component.state->ac,fixed_component.state->ac);
    ASSERT_TRUE(adaptive.prediction && fixed.prediction);
    EXPECT_EQ(adaptive.prediction,fixed.prediction);
    for(const auto * name:{"accepted-state-integrity","endpoint-verification","independent-replay","numerical-qualification"})
    {
        const auto check=std::find_if(adaptive_component.evidence.begin(),adaptive_component.evidence.end(),
            [&](const auto & item) {return item.name==name;});
        ASSERT_NE(check,adaptive_component.evidence.end())<<name;
        EXPECT_EQ(check->status,rhbm_gem::JointCheckStatus::Passed)<<name;
    }
    EXPECT_EQ(adaptive.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(adaptive.TargetRuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(adaptive.assembled_state); ASSERT_TRUE(fixed.assembled_state);
    EXPECT_EQ(adaptive.assembled_state->ac,fixed.assembled_state->ac);
    EXPECT_EQ(adaptive.assembled_state->tau,fixed.assembled_state->tau);
}

TEST(JointComponentSweepBudgetTest, UnqualifiedHardLimitAndMemoryRejectionCannotPublishComponentState)
{
    const auto problem=CoupledPair();
    n::BcdExecutionPolicy policy; policy.eta_upper_bound=std::log(2.2);
    n::BcdSearchPolicy search;
    search.core_atoms=1; search.nominal_sweeps=30; search.maximum_sweeps=31;
    policy.search_policy_for_testing=&search;
    const auto limited=n::FitBcdComponents(problem,{1.5,1.5},policy);
    ASSERT_EQ(limited.components.size(),1u);
    EXPECT_EQ(limited.components[0].search_sweeps,31u);
    EXPECT_EQ(limited.components[0].stop_reason,"block-hard-sweep-budget");
    EXPECT_FALSE(limited.components[0].state); EXPECT_FALSE(limited.assembled_state);
    EXPECT_NE(limited.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    search.core_atoms=n::ProductionBcdCoreAtoms; search.maximum_sweeps=40;
    search.global_state_memory_budget_for_testing=1;
    const auto memory=n::FitBcdComponents(problem,{1.5,1.5},policy);
    EXPECT_EQ(memory.components[0].search_sweeps,0u);
    EXPECT_TRUE(memory.components[0].stop_reason.starts_with("resource-limit:"));
    EXPECT_FALSE(memory.components[0].state); EXPECT_FALSE(memory.assembled_state);
}

TEST(JointComponentSweepBudgetTest, MultipleComponentsExtendAndStopIndependently)
{
    auto input=CoupledPair().Input();
    const auto single=second_stage_test::SyntheticJointWorkload("chain",1);
    const auto offset=input.observations.size();
    input.atom_ids.push_back("3"); input.support.emplace_back();
    for(const auto & support:single.support[0])
        input.support.back().push_back({offset+support.row,support.squared_distance});
    for(std::size_t row=0;row<single.observations.size();++row)
    {
        input.row_ids.push_back("separate/"+single.row_ids[row]);
        input.observations.push_back(single.observations[row]);
    }
    input.selection_domain->target_indices={0,1,2};
    input.selection_domain->pure_halo_support_exclusion=rhbm_gem::JointPureHaloSupportExclusionProvenance{
        input.observations.size(),input.observations.size(),0,0,0,0};
    const rhbm_gem::core::JointProblem problem(std::move(input));
    n::BcdExecutionPolicy policy; policy.eta_upper_bound=std::log(2.2);
    n::BcdSearchPolicy search;
    search.core_atoms=1; search.nominal_sweeps=30; search.maximum_sweeps=40;
    policy.search_policy_for_testing=&search;
    const auto fit=n::FitBcdComponents(problem,{1.5,1.5,1.5},policy);
    ASSERT_EQ(fit.components.size(),2u);
    EXPECT_GT(fit.components[0].search_sweeps,30u);
    EXPECT_LT(fit.components[1].search_sweeps,30u);
    EXPECT_EQ(fit.components[0].RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(fit.components[1].RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_EQ(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    EXPECT_TRUE(fit.assembled_state);
}
