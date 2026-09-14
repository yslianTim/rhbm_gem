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

rg::GaussianModel3D MakeGaussianWithCenterSignal(
    double center_signal,
    double width,
    double offset = 0.0)
{
    const auto amplitude{
        center_signal * std::pow(2.0 * std::acos(-1.0) * width * width, 1.5)
    };
    return rg::GaussianModel3D{ amplitude, width, offset };
}

LocalPotentialSampleList BuildSuspiciousGuardSamples(
    const rg::GaussianModel3D & previous_model,
    const std::vector<double> & radius_list,
    const std::vector<std::vector<double>> & zero_offset_response_list_by_radius)
{
    if (radius_list.size() != zero_offset_response_list_by_radius.size())
    {
        throw std::invalid_argument(
            "Suspicious guard sample input sizes are inconsistent.");
    }

    LocalPotentialSampleList sample_list;
    for (std::size_t radius_index = 0;
        radius_index < radius_list.size();
        radius_index++)
    {
        const auto radius{ radius_list.at(radius_index) };
        const auto previous_offset_response{
            previous_model.GetOffset() *
                previous_model.OffsetBasisAtDistance(radius)
        };
        for (const auto zero_offset_response :
            zero_offset_response_list_by_radius.at(radius_index))
        {
            SamplingPoint point;
            point.distance = radius;
            point.position = {
                radius,
                0.0,
                0.0
            };
            sample_list.emplace_back(LocalPotentialSample{
                zero_offset_response + previous_offset_response,
                point
            });
        }
    }
    return sample_list;
}

detail::SuspiciousGaussianReason EvaluateSuspiciousPostRefitUpdateForTest(
    const LocalPotentialSampleList & sample_entries,
    const rg::GaussianModel3D & previous_model,
    const rg::GaussianModel3D & candidate_model)
{
    const auto previous_baseline{
        detail::BuildPreviousSuspiciousProfileBaseline(
            sample_entries, previous_model)
    };
    return detail::AssessSuspiciousGaussianUpdate(
        sample_entries,
        candidate_model,
        previous_baseline,
        detail::SuspiciousUpdateMode::PostRefit).reason;
}

detail::SuspiciousGaussianReason EvaluateSuspiciousOffsetUpdateForTest(
    const LocalPotentialSampleList & sample_entries,
    const rg::GaussianModel3D & previous_model,
    const rg::GaussianModel3D & candidate_model)
{
    return detail::AssessSuspiciousGaussianUpdate(
        sample_entries,
        candidate_model,
        detail::BuildPreviousSuspiciousProfileBaseline(
            sample_entries, previous_model),
        detail::SuspiciousUpdateMode::OffsetOnly).reason;
}

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

std::unique_ptr<rg::ModelObject> BuildSeparatedLocalRefitFallbackDefenseModel()
{
    auto model{
        BuildDefenseModel(
            {
                std::array<double, 3>{ 0.0, 0.0, 0.0 },
                std::array<double, 3>{ 1.0e-4, 0.0, 0.0 },
                std::array<double, 3>{ 10.0, 0.0, 0.0 },
                std::array<double, 3>{ 10.0001, 0.0, 0.0 }
            },
            { Spot::C, Spot::C, Spot::C, Spot::C },
            { Element::CARBON, Element::CARBON, Element::CARBON, Element::CARBON },
            {
                rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
                rg::GaussianModel3D{ 5.5, 0.55, -0.15 },
                rg::GaussianModel3D{ 6.2, 0.55, 0.18 },
                rg::GaussianModel3D{ 5.7, 0.55, -0.12 }
            },
            rg::GaussianModel3D{ 5.8, 0.55, 0.0 })
    };
    auto * atom{ model->GetSelectedAtoms().front() };
    auto raw_sampling_entries{
        rg::AtomLocalPotentialView::For(*atom)
            .GetRawSamplingEntries(false)
    };
    LocalPotentialSampleList fallback_sampling_entries{
        raw_sampling_entries.at(0),
        raw_sampling_entries.at(6)
    };
    auto analysis{ model->EditAnalysis() };
    analysis.SetAtomLocalRawSamplingEntries(
        *atom, std::move(fallback_sampling_entries));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildSeparatedEmptyJointOffsetDefenseModel()
{
    auto model{
        BuildDefenseModel(
            {
                std::array<double, 3>{ 0.0, 0.0, 0.0 },
                std::array<double, 3>{ 10.0, 0.0, 0.0 },
                std::array<double, 3>{ 10.0001, 0.0, 0.0 }
            },
            { Spot::C, Spot::N, Spot::CA },
            { Element::CARBON, Element::NITROGEN, Element::CARBON },
            {
                rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
                rg::GaussianModel3D{ 6.2, 0.55, 0.18 },
                rg::GaussianModel3D{ 5.7, 0.55, -0.12 }
            },
            rg::GaussianModel3D{ 5.8, 0.55, 0.0 })
    };
    auto analysis{ model->EditAnalysis() };
    analysis.SetAtomLocalRawSamplingEntries(
        *model->GetSelectedAtoms().front(), {});
    return model;
}

std::unique_ptr<rg::ModelObject> BuildPostRefitRollbackChainDefenseModel()
{
    auto model{
        BuildDefenseModel(
            {
                std::array<double, 3>{ 0.0, 0.0, 0.0 },
                std::array<double, 3>{ 2.0, 0.0, 0.0 },
                std::array<double, 3>{ 4.0, 0.0, 0.0 },
                std::array<double, 3>{ 6.0, 0.0, 0.0 },
                std::array<double, 3>{ 8.0, 0.0, 0.0 }
            },
            { Spot::C, Spot::O, Spot::N, Spot::CA, Spot::C },
            {
                Element::CARBON,
                Element::OXYGEN,
                Element::NITROGEN,
                Element::CARBON,
                Element::CARBON
            },
            {
                rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
                rg::GaussianModel3D{ 5.5, 0.55, -0.15 },
                rg::GaussianModel3D{ 6.2, 0.55, 0.18 },
                rg::GaussianModel3D{ 5.7, 0.55, -0.12 },
                rg::GaussianModel3D{ 9.0, 0.55, 0.10 }
            },
            rg::GaussianModel3D{ 5.8, 0.55, 0.0 })
    };
    auto analysis{ model->EditAnalysis() };
    analysis.SetAtomLocalGaussianResult(
        FittingStage::Second,
        *model->GetSelectedAtoms().front(),
        MakeGaussianResult(rg::GaussianModel3D{ 1.0e300, 0.55, 0.0 }));
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

std::unique_ptr<rg::ModelObject> BuildFiniteNonphysicalProfileDefenseModel()
{
    const rg::GaussianModel3D initial_model{ 6.0, 0.55, 0.0 };
    auto model{
        BuildDefenseModel(
            { std::array<double, 3>{ 0.0, 0.0, 0.0 } },
            { Spot::O },
            { Element::OXYGEN },
            { initial_model },
            initial_model)
    };
    auto * atom{ model->GetSelectedAtoms().front() };
    auto raw_sampling_entries{
        rg::AtomLocalPotentialView::For(*atom)
            .GetRawSamplingEntries(false)
    };
    for (auto & sample : raw_sampling_entries)
    {
        const auto distance{ sample.point.distance };
        const auto outer_bias{ distance > 0.2 ? 12.0 : 8.0 };
        sample.response = initial_model.SignalAtDistance(distance) + outer_bias;
    }
    auto analysis{ model->EditAnalysis() };
    analysis.SetAtomLocalRawSamplingEntries(*atom, std::move(raw_sampling_entries));
    return model;
}

} // namespace

TEST(EstimatorSecondStageDefenseTest, SuspiciousEvaluatorReportsInvalidAndNonFiniteReasons)
{
    const auto previous_model{ MakeGaussianWithCenterSignal(0.1, 1.0) };
    const auto sample_list{
        BuildSuspiciousGuardSamples(
            previous_model,
            { 0.0 },
            { { 0.1 } })
    };

    EXPECT_EQ(
        EvaluateSuspiciousPostRefitUpdateForTest(
            sample_list,
            previous_model,
            rg::GaussianModel3D{ -1.0, 1.0, 0.0 }),
        detail::SuspiciousGaussianReason::InvalidModel);
    const auto previous_baseline{
        detail::BuildPreviousSuspiciousProfileBaseline(
            sample_list,
            previous_model)
    };
    const auto invalid_assessment{
        detail::AssessSuspiciousGaussianUpdate(
            sample_list,
            rg::GaussianModel3D{ -1.0, 1.0, 0.0 },
            previous_baseline,
            detail::SuspiciousUpdateMode::PostRefit)
    };
    EXPECT_EQ(
        invalid_assessment.reason,
        detail::SuspiciousGaussianReason::InvalidModel);
    EXPECT_TRUE(std::isinf(invalid_assessment.normalized_margin));

    auto non_finite_sample_list{ sample_list };
    non_finite_sample_list.front().response =
        std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(
        EvaluateSuspiciousPostRefitUpdateForTest(
            non_finite_sample_list,
            previous_model,
            previous_model),
        detail::SuspiciousGaussianReason::NonFiniteResponse);
}

TEST(EstimatorSecondStageDefenseTest, OffsetOnlyEvaluatorAppliesMagnitudeButSkipsWidthGuard)
{
    const auto previous_model{ MakeGaussianWithCenterSignal(0.1, 1.0) };
    const auto sample_list{
        BuildSuspiciousGuardSamples(
            previous_model,
            { 0.0 },
            { { 0.1 } })
    };
    const auto large_offset_model{
        previous_model.WithOffset(
            1.0 / previous_model.OffsetBasisAtDistance(0.0))
    };
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            sample_list,
            previous_model,
            large_offset_model),
        detail::SuspiciousGaussianReason::OffsetMagnitude);
    const auto large_offset_assessment{
        detail::AssessSuspiciousGaussianUpdate(
            sample_list,
            large_offset_model,
            detail::BuildPreviousSuspiciousProfileBaseline(
                sample_list,
                previous_model),
            detail::SuspiciousUpdateMode::OffsetOnly)
    };
    EXPECT_GT(large_offset_assessment.normalized_margin, 0.0);

    const auto wide_model{ MakeGaussianWithCenterSignal(0.1, 2.0) };
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            sample_list,
            previous_model,
            wide_model),
        detail::SuspiciousGaussianReason::None);
    const auto safe_offset_assessment{
        detail::AssessSuspiciousGaussianUpdate(
            sample_list,
            previous_model,
            detail::BuildPreviousSuspiciousProfileBaseline(
                sample_list,
                previous_model),
            detail::SuspiciousUpdateMode::OffsetOnly)
    };
    EXPECT_LE(safe_offset_assessment.normalized_margin, 0.0);
    EXPECT_EQ(
        EvaluateSuspiciousPostRefitUpdateForTest(
            sample_list,
            previous_model,
            wide_model),
        detail::SuspiciousGaussianReason::WidthGrowth);
}

TEST(EstimatorSecondStageDefenseTest, OffsetOnlyEvaluatorAcceptsUnchangedShapeOutsideProfileRange)
{
    const auto previous_model{ MakeGaussianWithCenterSignal(0.1, 1.0) };
    const auto sample_list{
        BuildSuspiciousGuardSamples(
            previous_model,
            { 0.0, 0.1 },
            { { 0.1 }, { previous_model.ResponseAtDistance(0.1) } })
    };
    const auto fallback_model{ previous_model.WithOffset(1.0e-4) };

    EXPECT_EQ(
        EvaluateSuspiciousPostRefitUpdateForTest(
            sample_list,
            previous_model,
            fallback_model),
        detail::SuspiciousGaussianReason::WidthGrowth);
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            sample_list,
            previous_model,
            fallback_model),
        detail::SuspiciousGaussianReason::None);
}

TEST(EstimatorSecondStageDefenseTest, CenterSignFlipRequiresPositiveSignalNoiseAndEffectSizeThresholds)
{
    const auto previous_model{ MakeGaussianWithCenterSignal(0.01, 1.0) };
    const auto noisy_sample_list{
        BuildSuspiciousGuardSamples(
            previous_model,
            { 0.0, 0.01, 0.02 },
            {
                { 0.9, 1.0, 1.1 },
                { 0.9, 1.0, 1.1 },
                { 0.9, 1.0, 1.1 }
            })
    };
    const auto candidate_with_center_offset = [&](double response)
    {
        return previous_model.WithOffset(
            response / previous_model.OffsetBasisAtDistance(0.0));
    };

    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            noisy_sample_list,
            previous_model,
            candidate_with_center_offset(1.3)),
        detail::SuspiciousGaussianReason::None);
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            noisy_sample_list,
            previous_model,
            candidate_with_center_offset(1.6)),
        detail::SuspiciousGaussianReason::CenterSignFlip);

    const auto zero_mad_sample_list{
        BuildSuspiciousGuardSamples(
            previous_model,
            { 0.0, 0.01, 0.02 },
            { { 1.0 }, { 1.0 }, { 1.0 } })
    };
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            zero_mad_sample_list,
            previous_model,
            candidate_with_center_offset(1.2)),
        detail::SuspiciousGaussianReason::None);
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            zero_mad_sample_list,
            previous_model,
            candidate_with_center_offset(1.3)),
        detail::SuspiciousGaussianReason::CenterSignFlip);

    const auto low_snr_sample_list{
        BuildSuspiciousGuardSamples(
            previous_model,
            { 0.0, 0.01, 0.02 },
            {
                { 0.0, 0.1, 0.2 },
                { 0.0, 0.1, 0.2 },
                { 0.0, 0.1, 0.2 }
            })
    };
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            low_snr_sample_list,
            previous_model,
            candidate_with_center_offset(0.3)),
        detail::SuspiciousGaussianReason::None);

    const auto negative_profile_samples{
        BuildSuspiciousGuardSamples(
            previous_model,
            { 0.0, 0.01, 0.02 },
            { { -1.0 }, { -1.0 }, { -1.0 } })
    };
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            negative_profile_samples,
            previous_model,
            candidate_with_center_offset(-1.5)),
        detail::SuspiciousGaussianReason::None);
}

TEST(EstimatorSecondStageDefenseTest, RadialReboundUsesResidualNoiseAndExcursionCount)
{
    const auto previous_model{ MakeGaussianWithCenterSignal(1.0e-6, 0.25) };
    const auto candidate_model{
        previous_model.WithOffset(
            0.6 / previous_model.OffsetBasisAtDistance(0.0))
    };
    const auto noisy_sample_list{
        BuildSuspiciousGuardSamples(
            previous_model,
            { 0.0, 0.5, 1.0 },
            {
                { 0.8, 1.0, 1.2 },
                { 0.8, 1.0, 1.2 },
                { 0.8, 1.0, 1.2 }
            })
    };
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            noisy_sample_list,
            previous_model,
            candidate_model),
        detail::SuspiciousGaussianReason::None);

    const auto low_noise_sample_list{
        BuildSuspiciousGuardSamples(
            previous_model,
            { 0.0, 0.5, 1.0 },
            {
                { 0.95, 1.0, 1.05 },
                { 0.95, 1.0, 1.05 },
                { 0.95, 1.0, 1.05 }
            })
    };
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            low_noise_sample_list,
            previous_model,
            candidate_model),
        detail::SuspiciousGaussianReason::RadialRebound);

    const auto excursion_model{ MakeGaussianWithCenterSignal(1.0, 1.0) };
    const std::vector<double> radius_list{ 0.0, 0.1, 0.2, 0.3, 0.4 };
    std::vector<double> innermost_response_list(10, 1.0);
    const auto one_excursion_samples{
        BuildSuspiciousGuardSamples(
            excursion_model,
            radius_list,
            {
                innermost_response_list,
                { 0.6 },
                { 0.9 },
                { 0.5 },
                { 0.6 }
            })
    };
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            one_excursion_samples,
            excursion_model,
            excursion_model),
        detail::SuspiciousGaussianReason::None);

    const auto two_excursion_samples{
        BuildSuspiciousGuardSamples(
            excursion_model,
            radius_list,
            {
                innermost_response_list,
                { 0.6 },
                { 0.9 },
                { 0.5 },
                { 0.8 }
            })
    };
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            two_excursion_samples,
            excursion_model,
            excursion_model),
        detail::SuspiciousGaussianReason::RadialRebound);
}

TEST(EstimatorSecondStageDefenseTest, WidthAndCompensationRemainActiveWithoutTrustedRadialShape)
{
    const auto previous_model{ MakeGaussianWithCenterSignal(0.1, 1.0) };
    const auto short_range_samples{
        BuildSuspiciousGuardSamples(
            previous_model,
            { 0.0, 0.5 },
            { { 0.1 }, { 0.1 } })
    };
    const auto range_wide_model{ MakeGaussianWithCenterSignal(0.1, 1.4) };
    EXPECT_EQ(
        EvaluateSuspiciousPostRefitUpdateForTest(
            short_range_samples,
            previous_model,
            range_wide_model),
        detail::SuspiciousGaussianReason::WidthGrowth);

    const auto previous_compensation_model{
        MakeGaussianWithCenterSignal(0.1, 1.0, 10.0)
    };
    const auto compensation_samples{
        BuildSuspiciousGuardSamples(
            previous_compensation_model,
            { 0.0 },
            { { 0.1 } })
    };
    const auto candidate_compensation_model{
        MakeGaussianWithCenterSignal(0.4, 1.4, 10.0)
    };
    EXPECT_EQ(
        EvaluateSuspiciousPostRefitUpdateForTest(
            compensation_samples,
            previous_compensation_model,
            candidate_compensation_model),
        detail::SuspiciousGaussianReason::AmplitudeOffsetCompensation);

    const auto same_direction_model{
        MakeGaussianWithCenterSignal(0.4, 0.8, 10.0)
    };
    EXPECT_EQ(
        EvaluateSuspiciousPostRefitUpdateForTest(
            compensation_samples,
            previous_compensation_model,
            same_direction_model),
        detail::SuspiciousGaussianReason::None);

    const auto insufficient_signal_model{
        MakeGaussianWithCenterSignal(0.2, 1.4, 10.0)
    };
    EXPECT_EQ(
        EvaluateSuspiciousPostRefitUpdateForTest(
            compensation_samples,
            previous_compensation_model,
            insufficient_signal_model),
        detail::SuspiciousGaussianReason::None);
}

TEST(EstimatorSecondStageDefenseTest, SuspiciousFailureMaskRequiresHardFailureOrSuspiciousFixedEndpoint)
{
    for (const char hard : std::array<char, 2>{ 0, 1 })
    for (const char shape : std::array<char, 2>{ 0, 1 })
    for (const char offset : std::array<char, 2>{ 0, 1 })
    for (const bool suspicious : { false, true })
    {
        const detail::SuspiciousBlockActivity activity{ { shape }, { offset }, { hard } };
        const std::vector<detail::SuspiciousGaussianAssessment> assessments{
            { suspicious ? detail::SuspiciousGaussianReason::WidthGrowth : detail::SuspiciousGaussianReason::None } };
        EXPECT_EQ(detail::BuildSuspiciousFailureAtomMask(activity, assessments),
            (detail::SuspiciousUpdateMask{ static_cast<char>(hard || ((shape || offset) && suspicious)) }));
    }
    EXPECT_TRUE(detail::BuildSuspiciousFailureAtomMask({}, {}).empty());
    for (std::size_t mask = 0; mask < 3; ++mask)
    for (const std::size_t wrong_size : { 0U, 2U })
    {
        detail::SuspiciousBlockActivity activity{ { 0 }, { 0 }, { 0 } };
        std::array<detail::SuspiciousUpdateMask *, 3> masks{
            &activity.shape_fixed_atom_mask, &activity.offset_fixed_atom_mask, &activity.hard_failure_atom_mask };
        masks.at(mask)->resize(wrong_size);
        try
        {
            detail::BuildSuspiciousFailureAtomMask(activity, std::vector<detail::SuspiciousGaussianAssessment>(1));
            FAIL() << "Expected an inconsistent-size exception";
        }
        catch (const std::invalid_argument & error)
        {
            EXPECT_STREQ(error.what(), "Suspicious failure activity and assessment sizes are inconsistent.");
        }
    }
}

TEST(EstimatorSecondStageDefenseTest, PostRefitSuspiciousLongChainKeepsTerminalBlocksFinite)
{
    auto model{ BuildPostRefitRollbackChainDefenseModel() };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    std::vector<rg::GaussianModel3D> previous_model_list;
    previous_model_list.reserve(selected_atoms.size());
    for (const auto * atom : selected_atoms)
    {
        previous_model_list.emplace_back(GetEstimateModel(*atom));
    }

    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*model, MakeSecondStageOptions());

    for (std::size_t i = 0; i < selected_atoms.size(); i++)
    {
        const auto fitted_model{ GetEstimateModel(*selected_atoms.at(i)) };
        ExpectGaussianModelsNear(fitted_model, previous_model_list.at(i), 1.0e-12);
    }
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

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

TEST(EstimatorSecondStageDefenseTest, LocalRefitFallbackDoesNotFreezeSameChemicalKeyAtoms)
{
    auto model{ BuildSeparatedLocalRefitFallbackDefenseModel() };
    const auto previous_fallback_model{
        GetEstimateModel(*model->GetSelectedAtoms().front())
    };
    const auto initial_remote_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4)
    };

    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*model, MakeSecondStageOptions());

    const auto fallback_model{
        GetEstimateModel(*model->GetSelectedAtoms().front())
    };
    EXPECT_DOUBLE_EQ(
        fallback_model.GetAmplitude(),
        previous_fallback_model.GetAmplitude());
    EXPECT_DOUBLE_EQ(
        fallback_model.GetWidth(),
        previous_fallback_model.GetWidth());
    EXPECT_NE(
        fallback_model.GetOffset(),
        previous_fallback_model.GetOffset());
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

TEST(EstimatorSecondStageDefenseTest, PersistentEmptySystemDoesNotBlockRemoteCluster)
{
    auto model{ BuildSeparatedEmptyJointOffsetDefenseModel() };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    const auto previous_empty_model{ GetEstimateModel(*selected_atoms.front()) };
    const auto initial_remote_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 1, 3)
    };

    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*model, MakeSecondStageOptions());

    ExpectGaussianModelsNear(
        GetEstimateModel(*selected_atoms.front()),
        previous_empty_model,
        1.0e-12);
    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 1, 3),
        initial_remote_error);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageIterationsRejectsTerminalFiniteNonphysicalProfile)
{
    auto model{ BuildFiniteNonphysicalProfileDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    const auto previous_model{ GetEstimateModel(*atom) };

    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*model, MakeSecondStageOptions());

    const auto fitted_model{ GetEstimateModel(*atom) };
    EXPECT_GT(fitted_model.GetAmplitude(), 0.0);
    EXPECT_GT(fitted_model.GetWidth(), 0.0);
    ExpectGaussianModelsNear(fitted_model, previous_model, 1.0e-12);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageIterationsAppliesCollinearRidgeGuard)
{
    auto model{ BuildNearCollinearDefenseModel() };
    const auto initial_error{ CalculateSelectedAtomResponseMeanSquaredError(*model) };

    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*model, MakeSecondStageOptions());

    const auto fitted_error{ CalculateSelectedAtomResponseMeanSquaredError(*model) };
    const auto tolerance{ 1.0e-3 * std::max(initial_error, 1.0) };
    EXPECT_LE(fitted_error, initial_error + tolerance);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageIterationsUpdatesHealthyVariablesAcrossClusters)
{
    auto model{ BuildSeparatedRollbackDefenseModel() };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    const std::array<rg::GaussianModel3D, 2> previous_left_model_list{
        GetEstimateModel(*selected_atoms.at(0)),
        GetEstimateModel(*selected_atoms.at(1))
    };
    const auto initial_right_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4)
    };

    model->EditAnalysis().CopyLocalFittingStageResult(FittingStage::Second, FittingStage::First);
    detail::RunSecondStageIterations(*model, MakeSecondStageOptions());

    for (std::size_t i = 0; i < previous_left_model_list.size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*selected_atoms.at(i)),
            previous_left_model_list.at(i),
            1.0e-12);
    }
    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4),
        initial_right_error);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}
