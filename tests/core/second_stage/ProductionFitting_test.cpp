#include <gtest/gtest.h>
#include "support/SecondStageTestSupport.hpp"
#include "core/detail/second_stage/FixedPointRecovery.hpp"
#include "core/detail/second_stage/CandidateEvaluation.hpp"
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <cmath>

namespace {
namespace d = rhbm_gem::core::detail;
namespace rg = rhbm_gem;

second_stage_test::JointPolishFixture Fixture()
{
    auto fixture{ second_stage_test::BuildJointPolishFixture(
        { {6.0, 0.55, 0.12}, {7.0, 0.52, -0.2} },
        { {6.0, 0.5, 0.1}, {7.0, 0.5, -0.1} }) };
    for (auto & atom : fixture.context.atom_list)
    {
        // Nonzero observational noise avoids a zero-variance solver boundary.
        for (std::size_t i = 0; i < atom.raw_sampling_entries.size(); ++i)
            atom.raw_sampling_entries[i].response += 0.001 * std::sin(1.0 + static_cast<double>(i));
        atom.refit_design = d::PreparedLocalGaussianDesign(atom.raw_sampling_entries, 0.0, 0.7);
    }
    return fixture;
}

d::FixedPointOperatorEvidence QualifiedEvidence(const d::FitState & state)
{
    d::FixedPointOperatorEvidence evidence;
    for (const auto & atom : state) evidence.state.push_back(atom.mdpde.GetModel());
    evidence.shape_available_atom_mask.assign(state.size(), 1);
    evidence.offset_available_atom_mask.assign(state.size(), 1);
    evidence.shape_solves.resize(state.size());
    for (std::size_t i = 0; i < state.size(); ++i)
    {
        evidence.shape_solves[i].status = rg::RHBMEstimationStatus::SUCCESS;
        evidence.offset_solves.emplace(d::ClusterKey{i}, d::JointOffsetSolveResult{d::JointOffsetSolveStatus::Converged});
    }
    return evidence;
}
}

TEST(ProductionFittingTest, NominalQualificationBelongsToEveryEndpoint)
{
    auto fixture{ Fixture() };
    auto evidence{ QualifiedEvidence(fixture.state) };
    EXPECT_TRUE(d::AssessNominalOperator(evidence, fixture.state).certificate.StrictOperatorPassed());
    evidence.shape_solves[1].status = rg::RHBMEstimationStatus::MAX_ITERATIONS_REACHED;
    EXPECT_FALSE(d::IsNominalOperatorSolverQualified(evidence));
    EXPECT_FALSE(d::QualifiedNominalResidualMeanSquare(evidence, fixture.state));
    evidence.shape_solves[1].status = rg::RHBMEstimationStatus::SUCCESS;
    evidence.offset_solves.at({1}).status = d::JointOffsetSolveStatus::IrlsMaximumIterationsReached;
    EXPECT_FALSE(d::IsNominalOperatorSolverQualified(evidence));
    evidence.offset_solves.at({1}).status = d::JointOffsetSolveStatus::Converged;
    evidence.shape_available_atom_mask[1] = 0;
    EXPECT_FALSE(d::AssessNominalOperator(evidence, fixture.state).certificate.StrictOperatorPassed());
}


TEST(ProductionFittingTest, NominalEvaluationIsReadOnlyAndUsesUnrestrictedSolveStatus)
{
    auto fixture{ Fixture() };
    const auto before{ fixture.state };
    const auto options{ second_stage_test::MakeSecondStageOptions() };
    const std::vector<d::ClusterKey> keys{{0}, {1}};
    const auto evidence{ d::EvaluateNominalOperator(fixture.context, keys, fixture.state, options, {1.0, 1.0}) };
    d::ClusterSolverWorkspaceMap workspaces;
    for (const auto & key : keys) workspaces.try_emplace(key);
    d::SuspiciousBlockActivity activity{{0,0}, {1,0}, {0,0}};
    const auto proposal{ d::BuildIterationProposal(fixture.context, keys, fixture.state, options,
        {1.0,1.0}, activity, workspaces) };
    for (std::size_t i = 0; i < before.size(); ++i)
    {
        second_stage_test::ExpectGaussianModelsNear(before[i].mdpde.GetModel(), fixture.state[i].mdpde.GetModel(), 0.0);
        second_stage_test::ExpectGaussianModelsNear(evidence.state[i], proposal.fixed_point_operator.state[i], 0.0);
        EXPECT_EQ(evidence.shape_solves[i].status, proposal.fixed_point_operator.shape_solves[i].status);
        EXPECT_EQ(evidence.shape_solves[i].diagnostics.iterations, proposal.fixed_point_operator.shape_solves[i].diagnostics.iterations);
    }
}
