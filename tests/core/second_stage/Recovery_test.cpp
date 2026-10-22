#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "support/SecondStageTestSupport.hpp"
#include "core/detail/second_stage/IterationProcess.hpp"
#include "core/detail/second_stage/Quarantine.hpp"
#include "core/detail/second_stage/SuspiciousUpdate.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

namespace {

namespace detail = rhbm_gem::core::detail;
namespace rg = rhbm_gem;
using rhbm_gem::FittingStage;

using second_stage_test::BuildDefenseModel;
using second_stage_test::BuildNearCollinearDefenseModel;
using second_stage_test::BuildSeparatedRollbackDefenseModel;
using second_stage_test::CalculateSelectedAtomResponseMeanSquaredError;
using second_stage_test::ExpectGaussianModelsNear;
using second_stage_test::ExpectSelectedAtomEstimatesAreFinite;
using second_stage_test::GetEstimateModel;
using second_stage_test::MakeGaussianResult;
using second_stage_test::MakeSecondStageOptions;


std::unique_ptr<rg::ModelObject> BuildSeparatedSystemBuildFailureDefenseModel()
{
    auto model{ BuildSeparatedRollbackDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    auto raw_sampling_entries{
        rg::AtomLocalPotentialView::For(*atom)
            .GetRawSamplingEntries(false)
    };
    raw_sampling_entries.front().response =
        std::numeric_limits<double>::quiet_NaN();
    auto analysis{ model->EditAnalysis() };
    analysis.SetAtomLocalRawSamplingEntries(*atom, std::move(raw_sampling_entries));
    return model;
}


std::unique_ptr<rg::ModelObject> BuildNonFiniteJointOffsetDefenseModel()
{
    auto model{
        BuildDefenseModel(
            { std::array<double, 3>{ 0.0, 0.0, 0.0 } },
            { Spot::O },
            { Element::OXYGEN },
            { rg::GaussianModel3D{ 8.0, 0.5, -0.1 } },
            rg::GaussianModel3D{ 7.0, 0.5, 0.0 })
    };
    auto * atom{ model->GetSelectedAtoms().front() };
    auto raw_sampling_entries{
        rg::AtomLocalPotentialView::For(*atom)
            .GetRawSamplingEntries(false)
    };
    raw_sampling_entries.front().response =
        std::numeric_limits<double>::quiet_NaN();
    auto analysis{ model->EditAnalysis() };
    analysis.SetAtomLocalRawSamplingEntries(*atom, std::move(raw_sampling_entries));
    return model;
}


} // namespace


TEST(EstimatorSecondStageDefenseTest, RunSecondStageIterationsFallsBackWhenJointOffsetSamplesAreNonFinite)
{
    auto model{ BuildNonFiniteJointOffsetDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    const auto previous_model{ GetEstimateModel(*atom) };
    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*model, MakeSecondStageOptions());

    ExpectGaussianModelsNear(GetEstimateModel(*atom), previous_model, 1.0e-12);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, SystemBuildFailureDoesNotBlockRemoteCluster)
{
    auto model{ BuildSeparatedSystemBuildFailureDefenseModel() };
    const auto initial_remote_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4)
    };

    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*model, MakeSecondStageOptions());

    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4),
        initial_remote_error);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}


TEST(EstimatorSecondStageDefenseTest, PersistentQuarantineReasonRequiresStableReasonAndReleasesOnDomainRetry)
{
    const detail::QuarantineTarget target{
        detail::QuarantineTargetKind::ShapeAtom,
        { 0 }
    };
    detail::QuarantineFailureStateMap state_by_target;
    const auto observe = [&](
        detail::SuspiciousGaussianReason reason,
        std::size_t recovery_revision)
    {
        return detail::UpdateQuarantineFailureState(
            {
                {
                    target,
                    detail::StabilizationTerminalFailure{
                        detail::StabilizationTerminalReason::GuardInfeasible,
                        reason
                    }
                }
            },
            {},
            {},
            recovery_revision,
            state_by_target);
    };

    EXPECT_TRUE(observe(detail::SuspiciousGaussianReason::WidthGrowth, 1)
        .entered_target_list.empty());
    EXPECT_NE(
        state_by_target.at(target).lifecycle,
        detail::QuarantineLifecycle::Frozen);
    EXPECT_TRUE(observe(detail::SuspiciousGaussianReason::WidthGrowth, 2)
        .entered_target_list.empty());
    EXPECT_TRUE(observe(
        detail::SuspiciousGaussianReason::AmplitudeOffsetCompensation,
        3).entered_target_list.empty());
    EXPECT_EQ(state_by_target.at(target).stable_iteration_count, 1U);
    for (std::size_t recovery_revision = 4;
        recovery_revision < 7;
        recovery_revision++)
    {
        EXPECT_TRUE(observe(
            detail::SuspiciousGaussianReason::AmplitudeOffsetCompensation,
            recovery_revision).entered_target_list.empty());
    }
    const auto entered{
        observe(
            detail::SuspiciousGaussianReason::AmplitudeOffsetCompensation,
            7)
    };
    ASSERT_EQ(entered.entered_target_list, (std::vector{ target }));
    ASSERT_EQ(
        state_by_target.at(target).lifecycle,
        detail::QuarantineLifecycle::Frozen);
    EXPECT_EQ(
        state_by_target.at(target).last_recovery_revision,
        7U);

    const auto released{
        detail::UpdateQuarantineFailureState(
            {},
            { target },
            { target },
            9,
            state_by_target)
    };
    EXPECT_EQ(released.released_target_list, (std::vector{ target }));
    EXPECT_TRUE(state_by_target.empty());

    state_by_target.emplace(
        target,
        detail::QuarantineFailureState{
            detail::StabilizationTerminalFailure{
                detail::StabilizationTerminalReason::GuardInfeasible,
                detail::SuspiciousGaussianReason::WidthGrowth
            },
            detail::kPersistentQuarantineFailureIterationLimit,
            9,
            detail::QuarantineLifecycle::Frozen
        });
    for (std::size_t retry = 1;
        retry <= 3;
        retry++)
    {
        const auto failed{
            detail::UpdateQuarantineFailureState(
                {
                    {
                        target,
                        detail::StabilizationTerminalFailure{
                            detail::StabilizationTerminalReason::GuardInfeasible,
                            detail::SuspiciousGaussianReason::WidthGrowth
                        }
                    }
                },
                { target },
                {},
                9 + retry,
                state_by_target)
        };
        EXPECT_EQ(failed.failed_retry_target_list, (std::vector{ target }));
        EXPECT_EQ(state_by_target.at(target).last_recovery_revision, 9 + retry);
    }
    EXPECT_EQ(
        state_by_target.at(target).lifecycle,
        detail::QuarantineLifecycle::Frozen);

    const detail::QuarantineTarget offset_target{
        detail::QuarantineTargetKind::OffsetAtom, { 0 }
    };
    const detail::QuarantineTarget other_offset_target{
        detail::QuarantineTargetKind::OffsetAtom, { 1 }
    };
    auto quarantined{ state_by_target.at(target) };
    quarantined.lifecycle = detail::QuarantineLifecycle::Frozen;
    for (const auto & released_target : { target, offset_target })
    {
        state_by_target = {
            { target, quarantined },
            { offset_target, quarantined },
            { other_offset_target, quarantined }
        };
        const auto independent_release{
            detail::UpdateQuarantineFailureState(
                {}, { released_target }, { released_target }, 20, state_by_target)
        };
        EXPECT_EQ(independent_release.released_target_list,
            (std::vector{ released_target }));
        EXPECT_FALSE(state_by_target.contains(released_target));
        EXPECT_EQ(state_by_target.size(), 2U);
        for (const auto & [remaining_target, state] : state_by_target)
        {
            EXPECT_NE(remaining_target, released_target);
            EXPECT_EQ(state.lifecycle, detail::QuarantineLifecycle::Frozen);
        }
    }
}
