#include <gtest/gtest.h>

#include "support/DynamicHaloMean.hpp"
#include "core/detail/joint_component/Numerics.hpp"

#include <algorithm>

using namespace joint_test::dynamic_halo;

namespace {
rhbm_gem::core::JointProblemInput Input()
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids={"z-target","a-target","halo"};
    input.row_ids={"r0","r1","r2","r3"};
    input.observations.assign(4,0.0);
    input.support={{{0,0.1},{1,0.2}},{{1,0.3},{2,0.35}},{{0,0.4},{1,0.5},{2,0.6},{3,0.7}}};
    input.selection_domain.emplace();
    input.selection_domain->target_indices={0,1};
    return input;
}
}

TEST(DynamicHaloMeanTest, AssignsOneOwnerAndExpandsRowsToOwnedHaloSupport)
{
    auto input=Input();
    const auto result=BuildOwnership(input,{{"block-a",0,{1}},{"block-z",0,{0}}});
    ASSERT_EQ(result.halos.size(),1u);
    EXPECT_EQ(result.blocks[result.halos[0].owner_block].id,"block-a");
    EXPECT_EQ(result.halos[0].owner_overlap,2u);
    EXPECT_EQ(result.halos[0].intersecting_blocks.size(),2u);
    EXPECT_EQ(result.halos[0].overlapping_targets.size(),2u);
    EXPECT_EQ(result.block_dynamic_rows[result.halos[0].owner_block],
        (std::vector<std::size_t>{0,1,2,3}));
}

TEST(DynamicHaloMeanTest, TieBreakDoesNotDependOnBlockTraversalOrder)
{
    const auto input=Input();
    const auto forward=BuildOwnership(input,{{"z-block",0,{0}},{"a-block",0,{1}}});
    const auto reverse=BuildOwnership(input,{{"a-block",0,{1}},{"z-block",0,{0}}});
    EXPECT_EQ(forward.blocks[forward.halos[0].owner_block].id,"a-block");
    EXPECT_EQ(reverse.blocks[reverse.halos[0].owner_block].id,"a-block");
}

TEST(DynamicHaloMeanTest, FlagsHaloThatCouplesSeparateComponents)
{
    auto input=Input();
    const auto result=BuildOwnership(input,{{"left",0,{0}},{"right",1,{1}}});
    EXPECT_EQ(result.cross_component_halos,1u);
}

TEST(DynamicHaloMeanTest, ComputesPhysicalMeansAndArithmeticWidthChainRule)
{
    const std::array values{Parameters{2.0,0.4,-0.3},Parameters{4.0,0.8,0.5}};
    const auto ab=Mean(values,false),abc=Mean(values,true);
    EXPECT_DOUBLE_EQ(ab.a,3.0); EXPECT_DOUBLE_EQ(ab.b,0.6); EXPECT_DOUBLE_EQ(ab.c,0.0);
    EXPECT_DOUBLE_EQ(abc.c,0.1);
    EXPECT_DOUBLE_EQ(MeanWidthLogDerivative(0.8,2),0.4);
}

TEST(DynamicHaloMeanTest, AugmentedDesignAndCoupledWidthJacobianMatchReplayFiniteDifferences)
{
    auto input=Input();
    const auto ownership=BuildOwnership(input,{{"whole",0,{0,1}}});
    std::vector<Parameters> targets{{0.0,0.45,-0.3},{2.0,0.75,0.5},{0.0,1.0,0.0}};
    constexpr double epsilon=1e-5;
    for(const auto mode:{Mode::AB,Mode::ABC})
    {
        const auto design=BuildAugmentedDesign(input,ownership,0,targets,mode);
        Eigen::VectorXd beta(4); beta<<targets[0].a,targets[0].c,targets[1].a,targets[1].c;
        const auto jacobian=CoupledWidthJacobian(input,ownership,0,targets,beta,mode);
        ASSERT_EQ(design.values.rows(),static_cast<Eigen::Index>(design.rows.size()));
        ASSERT_EQ(jacobian.cols(),2);

        // The A=0 column remains present at its active nonnegative boundary.
        Eigen::VectorXd positive=beta; positive(0)+=epsilon;
        const auto amplitude_difference=(design.values*positive-design.values*beta)/epsilon;
        EXPECT_TRUE(amplitude_difference.isApprox(design.values.col(0),1e-9));
        Eigen::VectorXd charged=beta; charged(3)+=epsilon;
        const auto charge_difference=(design.values*charged-design.values*beta)/epsilon;
        EXPECT_TRUE(charge_difference.isApprox(design.values.col(3),1e-9));

        for(Eigen::Index target=0;target<2;++target)
        {
            auto plus=targets,minus=targets;
            plus[static_cast<std::size_t>(target)].b*=std::exp(epsilon);
            minus[static_cast<std::size_t>(target)].b*=std::exp(-epsilon);
            const auto plus_design=BuildAugmentedDesign(input,ownership,0,plus,mode);
            const auto minus_design=BuildAugmentedDesign(input,ownership,0,minus,mode);
            const auto finite_difference=(plus_design.values*beta-minus_design.values*beta)/(2*epsilon);
            EXPECT_TRUE(finite_difference.isApprox(jacobian.col(target),2e-6))
                << "mode=" << (mode==Mode::AB?"AB":"ABC") << " target=" << target
                << " finite_difference=" << finite_difference.transpose()
                << " analytic=" << jacobian.col(target).transpose();
        }
    }
}

TEST(DynamicHaloMeanTest, NoHaloDesignMatchesTargetOnlyKernelColumns)
{
    auto input=Input();
    input.support[2].clear();
    const auto ownership=BuildOwnership(input,{{"whole",0,{0,1}}});
    ASSERT_TRUE(ownership.halos.empty());
    const std::vector<Parameters> targets{{1.0,0.45,-0.3},{2.0,0.75,0.5},{0.0,1.0,0.0}};
    const auto design=BuildAugmentedDesign(input,ownership,0,targets,Mode::ABC);
    EXPECT_EQ(design.rows,(std::vector<std::size_t>{0,1,2}));
    for(std::size_t slot=0;slot<2;++slot)
        for(const auto & support:input.support[slot])
        {
            const auto basis=rhbm_gem::core::joint_component::EvaluateKernel(
                support.squared_distance,targets[slot].b,2.5);
            const auto row=std::find(design.rows.begin(),design.rows.end(),support.row)-design.rows.begin();
            EXPECT_DOUBLE_EQ(design.values(static_cast<Eigen::Index>(row),static_cast<Eigen::Index>(2*slot)),basis.gaussian);
            EXPECT_DOUBLE_EQ(design.values(static_cast<Eigen::Index>(row),static_cast<Eigen::Index>(2*slot+1)),basis.charge);
        }
}

TEST(DynamicHaloMeanTest, SharedHaloAddsItsRemoteSupportRowsToOnlyItsOwnersBlock)
{
    const auto input=Input();
    const auto ownership=BuildOwnership(input,{{"left",0,{0}},{"right",0,{1}}});
    ASSERT_EQ(ownership.halos.size(),1u);
    const auto owner=ownership.halos.front().owner_block;
    EXPECT_EQ(ownership.blocks[owner].id,"left");
    const auto & rows=ownership.block_dynamic_rows[owner];
    EXPECT_EQ(rows,(std::vector<std::size_t>{0,1,2,3}));
}

TEST(DynamicHaloMeanTest, OneRowHaloAndMixedSignChargeRemainValidInAbcMode)
{
    auto input=Input();
    input.support[2]={{1,0.4}};
    const auto ownership=BuildOwnership(input,{{"whole",0,{0,1}}});
    std::vector<Parameters> targets{{1.0,0.5,-0.5},{1.0,0.7,0.5},{0.0,1.0,0.0}};
    const auto derived=DeriveHalos(input,ownership,targets,Mode::ABC);
    EXPECT_DOUBLE_EQ(derived[2].c,0.0);
    const auto design=BuildAugmentedDesign(input,ownership,0,targets,Mode::ABC);
    EXPECT_EQ(design.rows.size(),3u);
    Eigen::VectorXd beta(4); beta<<targets[0].a,targets[0].c,targets[1].a,targets[1].c;
    const auto jacobian=CoupledWidthJacobian(input,ownership,0,targets,beta,Mode::ABC);
    EXPECT_TRUE(jacobian.allFinite());
}

TEST(DynamicHaloMeanTest, RejectedCandidateRollsBackTargetAndDerivedHaloTogether)
{
    auto input=Input();
    const auto ownership=BuildOwnership(input,{{"whole",0,{0,1}}});
    std::vector<Parameters> current{{1.0,0.5,0.0},{1.0,0.6,0.0},{0.0,1.0,0.0}};
    input.observations=ReplayPrediction(input,DeriveHalos(input,ownership,current,Mode::ABC));
    auto accepted=BuildState(input,ownership,current,Mode::ABC,1.0);
    auto proposal=current;
    proposal[0].a+=0.2;
    auto candidate=BuildState(input,ownership,proposal,Mode::ABC,1.0);
    ASSERT_GT(candidate.objective,accepted.objective);
    const auto before=accepted;
    EXPECT_FALSE(CommitMonotone(accepted,std::move(candidate)));
    EXPECT_DOUBLE_EQ(accepted.objective,before.objective);
    EXPECT_EQ(accepted.prediction,before.prediction);
    EXPECT_DOUBLE_EQ(accepted.contributors[2].a,before.contributors[2].a);
}

TEST(DynamicHaloMeanTest, FullReplayRebuildsHaloStateFromTargets)
{
    auto input=Input();
    const auto ownership=BuildOwnership(input,{{"whole",0,{0,1}}});
    const std::vector<Parameters> targets{{0.7,0.4,-0.2},{1.2,0.8,0.2},{9.0,1.0,9.0}};
    const auto state=BuildState(input,ownership,targets,Mode::ABC,2.0);
    const auto independently_replayed=ReplayPrediction(input,state.contributors);
    EXPECT_EQ(state.prediction,independently_replayed);
    EXPECT_DOUBLE_EQ(state.contributors[2].a,0.95);
    EXPECT_DOUBLE_EQ(state.contributors[2].b,0.6);
    EXPECT_DOUBLE_EQ(state.contributors[2].c,0.0);
}

TEST(DynamicHaloMeanTest, StaticMeanUsesInitialTargetStateAndDoesNotEnterLocalDesign)
{
    const auto input=Input();
    const auto ownership=BuildOwnership(input,{{"whole",0,{0,1}}});
    std::vector<Parameters> initial{{0.0,0.4,0.0},{0.0,0.8,0.0},{0.0,1.0,0.0}};
    const auto fixed=DeriveHalos(input,ownership,initial,Mode::StaticABC);
    EXPECT_DOUBLE_EQ(fixed[2].b,0.6);
    auto updated=initial;
    updated[0].b=0.5;
    updated[0].a=1.0;
    const auto design=BuildAugmentedDesign(input,ownership,0,updated,Mode::StaticABC);
    const auto dynamic_design=BuildAugmentedDesign(input,ownership,0,updated,Mode::ABC);
    EXPECT_LT(design.values.norm(),dynamic_design.values.norm());
    EXPECT_DOUBLE_EQ(fixed[2].a,0.0);
    EXPECT_DOUBLE_EQ(fixed[2].b,0.6);
    EXPECT_DOUBLE_EQ(fixed[2].c,0.0);
}
