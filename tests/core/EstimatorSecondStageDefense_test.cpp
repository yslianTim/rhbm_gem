#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <sstream>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "core/detail/GaussianModelOperations.hpp"
#include "core/detail/PreparedLocalGaussianFit.hpp"
#include "core/detail/SecondStageFitting.hpp"
#include "core/detail/CouplingGraph.hpp"
#include "core/detail/JointFitting.hpp"
#include "core/detail/CandidateSelection.hpp"
#include "core/detail/Diagnosis.hpp"
#include "core/detail/IterationProcess.hpp"
#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace {
namespace alg = rhbm_gem::algorithm;
namespace audit_detail = rhbm_gem::core::detail;
namespace change_detail = rhbm_gem::core::detail;
namespace conditioning_detail = rhbm_gem::core::detail;
namespace coupling_detail = rhbm_gem::core::detail;
namespace backtracking_detail = rhbm_gem::core::detail;
namespace residual_detail = rhbm_gem::core::detail;
namespace median_detail = rhbm_gem::core::detail;
namespace health_detail = rhbm_gem::core::detail;
namespace offset_detail = rhbm_gem::core::detail;
namespace polish_detail = rhbm_gem::core::detail;
namespace seed_detail = rhbm_gem::core::detail;
using rhbm_gem::FittingStage;
namespace trust_detail = rhbm_gem::core::detail;
namespace rt = rhbm_gem::core;
namespace rg = rhbm_gem;

static_assert(!std::is_default_constructible_v<audit_detail::ClusterHealth>);

rt::FitOptions MakeSecondStageOptions()
{
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;
    return options;
}

double Distance(
    const std::array<double, 3> & lhs,
    const std::array<double, 3> & rhs)
{
    const auto dx{ lhs.at(0) - rhs.at(0) };
    const auto dy{ lhs.at(1) - rhs.at(1) };
    const auto dz{ lhs.at(2) - rhs.at(2) };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void AddCouplingGraphSample(
    coupling_detail::CouplingGraphBuilder & builder,
    coupling_detail::SampleRef sample_id,
    std::vector<coupling_detail::GraphParticipant> participant_list)
{
    builder.AddSample(sample_id, participant_list);
}

bool HasCouplingNeighbor(
    const coupling_detail::GraphTopology & topology,
    std::size_t atom_index,
    std::size_t neighbor_index)
{
    const auto & neighbor_index_list{ topology.adjacency_list.at(atom_index) };
    return std::find(
        neighbor_index_list.begin(),
        neighbor_index_list.end(),
        neighbor_index) != neighbor_index_list.end();
}

std::unique_ptr<rg::AtomObject> MakeAtom(
    int serial_id,
    Spot spot,
    Element element,
    const std::array<double, 3> & position)
{
    auto atom{ std::make_unique<rg::AtomObject>() };
    atom->SetSerialID(serial_id);
    atom->SetChainID("A");
    atom->SetSequenceID(serial_id);
    atom->SetComponentKey(1);
    atom->SetAtomKey(static_cast<AtomKey>(spot));
    atom->SetElement(element);
    atom->SetSpot(spot);
    atom->SetPosition(position);
    return atom;
}

rg::LocalGaussianResult MakeGaussianResult(const rg::GaussianModel3D & model)
{
    rg::LocalGaussianResult result;
    result.ols = rg::GaussianModel3DWithUncertainty{
        model,
        rg::GaussianModel3DUncertainty{}
    };
    result.mdpde = rg::GaussianModel3DWithUncertainty{
        model,
        rg::GaussianModel3DUncertainty{}
    };
    return result;
}

LocalPotentialSampleList BuildSamples(
    const rg::AtomObject & target_atom,
    const std::vector<rg::AtomObject *> & atom_list,
    const std::vector<rg::GaussianModel3D> & truth_model_list)
{
    const std::array<std::array<double, 3>, 6> direction_list{
        std::array<double, 3>{ 1.0, 0.0, 0.0 },
        std::array<double, 3>{ -1.0, 0.0, 0.0 },
        std::array<double, 3>{ 0.0, 1.0, 0.0 },
        std::array<double, 3>{ 0.0, -1.0, 0.0 },
        std::array<double, 3>{ 0.0, 0.0, 1.0 },
        std::array<double, 3>{ 0.0, 0.0, -1.0 }
    };
    const std::array<double, 4> radius_list{ 0.15, 0.35, 0.65, 0.95 };

    LocalPotentialSampleList sample_list;
    const auto target_position{ target_atom.GetPosition() };
    for (const auto radius : radius_list)
    {
        for (const auto & direction : direction_list)
        {
            SamplingPoint point;
            point.distance = radius;
            for (std::size_t axis = 0; axis < point.position.size(); axis++)
            {
                point.position.at(axis) =
                    target_position.at(axis) + radius * direction.at(axis);
            }

            double response{ 0.0 };
            for (std::size_t i = 0; i < atom_list.size(); i++)
            {
                response += truth_model_list.at(i).ResponseAtDistance(
                    Distance(point.position, atom_list.at(i)->GetPosition()));
            }
            sample_list.emplace_back(LocalPotentialSample{ response, point });
        }
    }
    return sample_list;
}

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

std::pair<offset_detail::SecondStageContext,
    offset_detail::SecondStageModelSnapshot>
BuildJointOffsetEstimationFixture(
    const std::vector<rg::GaussianModel3D> & model_list,
    const std::vector<double> & target_offset_list)
{
    if (model_list.size() != target_offset_list.size())
    {
        throw std::invalid_argument(
            "Joint offset fixture input sizes are inconsistent.");
    }

    offset_detail::SecondStageContext context;
    context.atom_list.resize(model_list.size());
    for (std::size_t atom_index = 0;
        atom_index < model_list.size();
        atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
        SamplingPoint point;
        point.distance = 0.35;
        atom_context.raw_sampling_entries.emplace_back(LocalPotentialSample{
            model_list.at(atom_index).SignalAtDistance(point.distance) +
                target_offset_list.at(atom_index) *
                    model_list.at(atom_index).OffsetBasisAtDistance(point.distance),
            point
        });
        atom_context.neighbor_atom_sample_offset_list = { 0, 0 };
    }
    return {
        std::move(context),
        offset_detail::SecondStageModelSnapshot{
            model_list
        }
    };
}

struct JointPolishFixture
{
    polish_detail::SecondStageContext context{};
    polish_detail::FitState state{};
    std::vector<polish_detail::SampleRef>
        sample_ref_list{};
};

JointPolishFixture BuildJointPolishFixture(
    const std::vector<rg::GaussianModel3D> & base_model_list,
    const std::vector<rg::GaussianModel3D> & target_model_list)
{
    if (base_model_list.size() != target_model_list.size() ||
        base_model_list.empty())
    {
        throw std::invalid_argument(
            "Joint polish fixture input sizes are inconsistent.");
    }

    constexpr std::array<double, 5> distance_list{
        0.0,
        0.15,
        0.30,
        0.45,
        0.60
    };
    JointPolishFixture fixture;
    fixture.context.atom_list.resize(base_model_list.size());
    fixture.state.reserve(base_model_list.size());
    for (std::size_t atom_index = 0;
        atom_index < base_model_list.size();
        atom_index++)
    {
        auto & atom_context{ fixture.context.at(atom_index) };
        atom_context.neighbor_atom_sample_offset_list.assign(
            distance_list.size() + 1,
            0);
        for (const auto distance : distance_list)
        {
            atom_context.raw_sampling_entries.emplace_back(
                LocalPotentialSample{
                    target_model_list.at(atom_index).ResponseAtDistance(distance),
                    SamplingPoint{ distance }
                });
            fixture.sample_ref_list.emplace_back(
                polish_detail::SampleRef{
                    atom_index,
                    atom_context.raw_sampling_entries.size() - 1
                });
        }
        fixture.state.emplace_back(MakeGaussianResult(
            base_model_list.at(atom_index)));

    }
    return fixture;
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

audit_detail::SuspiciousGaussianReason EvaluateSuspiciousPostRefitUpdateForTest(
    const LocalPotentialSampleList & sample_entries,
    const rg::GaussianModel3D & previous_model,
    const rg::GaussianModel3D & candidate_model,
    const rt::FitOptions & options)
{
    const auto previous_baseline{
        audit_detail::BuildPreviousSuspiciousProfileBaseline(
            sample_entries, previous_model, options)
    };
    return audit_detail::AssessSuspiciousGaussianUpdate(
        sample_entries,
        candidate_model,
        options,
        previous_baseline,
        audit_detail::SuspiciousUpdateMode::PostRefit).reason;
}

audit_detail::SuspiciousGaussianReason EvaluateSuspiciousOffsetUpdateForTest(
    const LocalPotentialSampleList & sample_entries,
    const rg::GaussianModel3D & previous_model,
    const rg::GaussianModel3D & candidate_model,
    const rt::FitOptions & options)
{
    return audit_detail::AssessSuspiciousGaussianUpdate(
        sample_entries,
        candidate_model,
        options,
        audit_detail::BuildPreviousSuspiciousProfileBaseline(
            sample_entries, previous_model, options),
        audit_detail::SuspiciousUpdateMode::OffsetOnly).reason;
}

TEST(EstimatorSecondStageDefenseTest, SuspiciousEvaluatorReportsInvalidAndNonFiniteReasons)
{
    const auto options{ MakeSecondStageOptions() };
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
            rg::GaussianModel3D{ -1.0, 1.0, 0.0 },
            options),
        audit_detail::SuspiciousGaussianReason::InvalidModel);
    const auto previous_baseline{
        audit_detail::BuildPreviousSuspiciousProfileBaseline(
            sample_list,
            previous_model,
            options)
    };
    const auto invalid_assessment{
        audit_detail::AssessSuspiciousGaussianUpdate(
            sample_list,
            rg::GaussianModel3D{ -1.0, 1.0, 0.0 },
            options,
            previous_baseline,
            audit_detail::SuspiciousUpdateMode::PostRefit)
    };
    EXPECT_EQ(
        invalid_assessment.reason,
        audit_detail::SuspiciousGaussianReason::InvalidModel);
    EXPECT_TRUE(std::isinf(invalid_assessment.normalized_margin));

    auto non_finite_sample_list{ sample_list };
    non_finite_sample_list.front().response =
        std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(
        EvaluateSuspiciousPostRefitUpdateForTest(
            non_finite_sample_list,
            previous_model,
            previous_model,
            options),
        audit_detail::SuspiciousGaussianReason::NonFiniteResponse);
}

TEST(EstimatorSecondStageDefenseTest, OffsetOnlyEvaluatorAppliesMagnitudeButSkipsWidthGuard)
{
    const auto options{ MakeSecondStageOptions() };
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
            large_offset_model,
            options),
        audit_detail::SuspiciousGaussianReason::OffsetMagnitude);
    const auto large_offset_assessment{
        audit_detail::AssessSuspiciousGaussianUpdate(
            sample_list,
            large_offset_model,
            options,
            audit_detail::BuildPreviousSuspiciousProfileBaseline(
                sample_list,
                previous_model,
                options),
            audit_detail::SuspiciousUpdateMode::OffsetOnly)
    };
    EXPECT_GT(large_offset_assessment.normalized_margin, 0.0);

    const auto wide_model{ MakeGaussianWithCenterSignal(0.1, 2.0) };
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            sample_list,
            previous_model,
            wide_model,
            options),
        audit_detail::SuspiciousGaussianReason::None);
    const auto safe_offset_assessment{
        audit_detail::AssessSuspiciousGaussianUpdate(
            sample_list,
            previous_model,
            options,
            audit_detail::BuildPreviousSuspiciousProfileBaseline(
                sample_list,
                previous_model,
                options),
            audit_detail::SuspiciousUpdateMode::OffsetOnly)
    };
    EXPECT_LE(safe_offset_assessment.normalized_margin, 0.0);
    EXPECT_EQ(
        EvaluateSuspiciousPostRefitUpdateForTest(
            sample_list,
            previous_model,
            wide_model,
            options),
        audit_detail::SuspiciousGaussianReason::WidthGrowth);
}

TEST(EstimatorSecondStageDefenseTest, OffsetOnlyEvaluatorAcceptsUnchangedShapeOutsideProfileRange)
{
    const auto options{ MakeSecondStageOptions() };
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
            fallback_model,
            options),
        audit_detail::SuspiciousGaussianReason::WidthGrowth);
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            sample_list,
            previous_model,
            fallback_model,
            options),
        audit_detail::SuspiciousGaussianReason::None);
}

TEST(EstimatorSecondStageDefenseTest, CenterSignFlipRequiresPositiveSignalNoiseAndEffectSizeThresholds)
{
    auto options{ MakeSecondStageOptions() };
    options.distance_max = 0.02;
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
            candidate_with_center_offset(1.3),
            options),
        audit_detail::SuspiciousGaussianReason::None);
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            noisy_sample_list,
            previous_model,
            candidate_with_center_offset(1.6),
            options),
        audit_detail::SuspiciousGaussianReason::CenterSignFlip);

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
            candidate_with_center_offset(1.2),
            options),
        audit_detail::SuspiciousGaussianReason::None);
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            zero_mad_sample_list,
            previous_model,
            candidate_with_center_offset(1.3),
            options),
        audit_detail::SuspiciousGaussianReason::CenterSignFlip);

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
            candidate_with_center_offset(0.3),
            options),
        audit_detail::SuspiciousGaussianReason::None);

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
            candidate_with_center_offset(-1.5),
            options),
        audit_detail::SuspiciousGaussianReason::None);
}

TEST(EstimatorSecondStageDefenseTest, RadialReboundUsesResidualNoiseAndExcursionCount)
{
    auto options{ MakeSecondStageOptions() };
    options.distance_max = 4.0;
    const auto previous_model{ MakeGaussianWithCenterSignal(1.0e-6, 1.0) };
    const auto candidate_model{
        previous_model.WithOffset(
            0.6 / previous_model.OffsetBasisAtDistance(0.0))
    };
    const auto noisy_sample_list{
        BuildSuspiciousGuardSamples(
            previous_model,
            { 0.0, 2.0, 4.0 },
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
            candidate_model,
            options),
        audit_detail::SuspiciousGaussianReason::None);

    const auto low_noise_sample_list{
        BuildSuspiciousGuardSamples(
            previous_model,
            { 0.0, 2.0, 4.0 },
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
            candidate_model,
            options),
        audit_detail::SuspiciousGaussianReason::RadialRebound);

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
    options.distance_max = 0.41;
    EXPECT_EQ(
        EvaluateSuspiciousOffsetUpdateForTest(
            one_excursion_samples,
            excursion_model,
            excursion_model,
            options),
        audit_detail::SuspiciousGaussianReason::None);

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
            excursion_model,
            options),
        audit_detail::SuspiciousGaussianReason::RadialRebound);
}

TEST(EstimatorSecondStageDefenseTest, WidthAndCompensationRemainActiveWithoutTrustedRadialShape)
{
    auto options{ MakeSecondStageOptions() };
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
            range_wide_model,
            options),
        audit_detail::SuspiciousGaussianReason::WidthGrowth);

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
            candidate_compensation_model,
            options),
        audit_detail::SuspiciousGaussianReason::AmplitudeOffsetCompensation);

    const auto same_direction_model{
        MakeGaussianWithCenterSignal(0.4, 0.8, 10.0)
    };
    EXPECT_EQ(
        EvaluateSuspiciousPostRefitUpdateForTest(
            compensation_samples,
            previous_compensation_model,
            same_direction_model,
            options),
        audit_detail::SuspiciousGaussianReason::None);

    const auto insufficient_signal_model{
        MakeGaussianWithCenterSignal(0.2, 1.4, 10.0)
    };
    EXPECT_EQ(
        EvaluateSuspiciousPostRefitUpdateForTest(
            compensation_samples,
            previous_compensation_model,
            insufficient_signal_model,
            options),
        audit_detail::SuspiciousGaussianReason::None);
}

std::unique_ptr<rg::ModelObject> BuildDefenseModel(
    const std::vector<std::array<double, 3>> & position_list,
    const std::vector<Spot> & spot_list,
    const std::vector<Element> & element_list,
    const std::vector<rg::GaussianModel3D> & truth_model_list,
    const rg::GaussianModel3D & initial_model)
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    for (std::size_t i = 0; i < position_list.size(); i++)
    {
        atom_list.emplace_back(MakeAtom(
            static_cast<int>(i + 1),
            spot_list.at(i),
            element_list.at(i),
            position_list.at(i)));
    }

    auto model{ std::make_unique<rg::ModelObject>(std::move(atom_list)) };
    model->SelectAllAtoms();
    auto analysis{ model->EditAnalysis() };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    for (auto * atom : selected_atoms)
    {
        analysis.SetAtomLocalAlphaR(FittingStage::Second, *atom, 0.0);
        analysis.SetAtomLocalGaussianResult(
            FittingStage::Second,
            *atom,
            MakeGaussianResult(initial_model));
        analysis.SetAtomLocalRawSamplingEntries(
            *atom, BuildSamples(*atom, selected_atoms, truth_model_list));
    }
    return model;
}

std::unique_ptr<rg::ModelObject> BuildUnselectedContributorDefenseModel(
    const std::array<rg::GaussianModel3D, 2> & selected_seed_model_list,
    const std::array<rg::GaussianModel3D, 7> & truth_model_list,
    bool use_alternate_group_keys = false,
    bool shared_cluster = false,
    bool shared_contributor = false,
    bool use_alternate_residue_keys = false)
{
    std::vector<std::array<double, 3>> position_list{
        { 0.0, 0.0, 0.0 },
        { 7.0, 0.0, 0.0 },
        { 0.60, 0.0, 0.0 },
        { 7.60, 0.0, 0.0 },
        { 0.90, 0.0, 0.0 },
        { 15.0, 0.0, 0.0 },
        { 22.0, 0.0, 0.0 }
    };
    if (shared_cluster)
    {
        // Selected sample responses overlap, independently of residue identity.
        position_list.at(1) = { 0.8, 0.0, 0.0 };
        position_list.at(3) = { 1.4, 0.0, 0.0 };
    }
    if (shared_contributor)
    {
        // Contributor 3 affects both targets, but selected responses do not overlap.
        position_list.at(1) = { 4.0, 0.0, 0.0 };
        position_list.at(2) = { 2.0, 0.0, 0.0 };
        position_list.at(3) = { 4.6, 0.0, 0.0 };
    }
    auto spot_list = std::vector<Spot>{
        Spot::C,
        Spot::C,
        Spot::C,
        Spot::C,
        Spot::C,
        Spot::O,
        Spot::O
    };
    auto element_list = std::vector<Element>{
        Element::CARBON,
        Element::CARBON,
        Element::CARBON,
        Element::CARBON,
        Element::HYDROGEN,
        Element::OXYGEN,
        Element::OXYGEN
    };
    if (use_alternate_group_keys)
    {
        spot_list.at(0) = Spot::N;
        spot_list.at(1) = Spot::O;
        spot_list.at(2) = Spot::N;
        spot_list.at(3) = Spot::O;
    }
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    for (std::size_t i = 0; i < position_list.size(); i++)
    {
        atom_list.emplace_back(MakeAtom(
            static_cast<int>(i + 1),
            spot_list.at(i),
            element_list.at(i),
            position_list.at(i)));
        if (use_alternate_residue_keys)
        {
            atom_list.back()->SetChainID("relabeled");
            atom_list.back()->SetSequenceID(42);
            atom_list.back()->SetResidue(Residue::GLY);
        }
    }
    auto model{ std::make_unique<rg::ModelObject>(std::move(atom_list)) };
    model->SelectAllAtoms();
    for (int serial_id = 3; serial_id <= 7; serial_id++)
    {
        model->SetAtomSelected(serial_id, false);
    }

    std::vector<rg::AtomObject *> all_atoms;
    std::vector<rg::GaussianModel3D> truth_models;
    for (std::size_t atom_index = 0;
        atom_index < model->GetAtomList().size(); atom_index++)
    {
        all_atoms.emplace_back(model->GetAtomList().at(atom_index).get());
        truth_models.emplace_back(truth_model_list.at(atom_index));
    }

    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    for (std::size_t atom_index = 0;
        atom_index < selected_atoms.size();
        atom_index++)
    {
        auto * atom{ selected_atoms.at(atom_index) };
        analysis.SetAtomLocalAlphaR(FittingStage::Second, *atom, 0.0);
        analysis.SetAtomLocalGaussianResult(
            FittingStage::Second,
            *atom,
            MakeGaussianResult(selected_seed_model_list.at(atom_index)));
        analysis.SetAtomLocalRawSamplingEntries(
            *atom, BuildSamples(*atom, all_atoms, truth_models));
    }
    return model;
}

std::unique_ptr<rg::ModelObject> BuildNearCollinearDefenseModel(
    double intensity_scale = 1.0)
{
    return BuildDefenseModel(
        {
            std::array<double, 3>{ 0.0, 0.0, 0.0 },
            std::array<double, 3>{ 1.0e-4, 0.0, 0.0 }
        },
        { Spot::C, Spot::O },
        { Element::CARBON, Element::OXYGEN },
        {
            rg::GaussianModel3D{ 6.0 * intensity_scale, 0.55, 0.20 * intensity_scale },
            rg::GaussianModel3D{ 5.5 * intensity_scale, 0.55, -0.15 * intensity_scale }
        },
        rg::GaussianModel3D{ 5.75 * intensity_scale, 0.55, 0.0 });
}

std::unique_ptr<rg::ModelObject> BuildJointPolishDefenseModel()
{
    auto model{ BuildDefenseModel(
        {
            std::array<double, 3>{ 0.0, 0.0, 0.0 },
            std::array<double, 3>{ 0.8, 0.0, 0.0 }
        },
        { Spot::C, Spot::O },
        { Element::CARBON, Element::OXYGEN },
        {
            rg::GaussianModel3D{ 6.0, 0.45, 0.20 },
            rg::GaussianModel3D{ 4.5, 0.70, -0.12 }
        },
        rg::GaussianModel3D{ 5.25, 0.60, 0.02 }) };
    auto analysis{ model->EditAnalysis() };
    const auto & atom_list{ model->GetSelectedAtoms() };
    // Keep the raw joint-offset update inside the trust region now that a
    // noise-scale rebound no longer supplies an incidental rollback.
    analysis.SetAtomLocalGaussianResult(
        FittingStage::Second,
        *atom_list.at(0),
        MakeGaussianResult(rg::GaussianModel3D{ 5.8, 0.47, 0.15 }));
    analysis.SetAtomLocalGaussianResult(
        FittingStage::Second,
        *atom_list.at(1),
        MakeGaussianResult(rg::GaussianModel3D{ 4.7, 0.68, -0.08 }));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildIndependentOffsetDefenseModel(
    double intensity_scale = 1.0,
    bool alternate_keys = false)
{
    auto model{ BuildDefenseModel(
        {
            std::array<double, 3>{ 0.0, 0.0, 0.0 },
            std::array<double, 3>{ 0.8, 0.0, 0.0 }
        },
        alternate_keys ? std::vector<Spot>{ Spot::N, Spot::O } :
            std::vector<Spot>{ Spot::C, Spot::C },
        { Element::CARBON, Element::CARBON },
        {
            rg::GaussianModel3D{
                6.0 * intensity_scale,
                0.48,
                0.05 * intensity_scale
            },
            rg::GaussianModel3D{
                4.8 * intensity_scale,
                0.62,
                -0.08 * intensity_scale
            }
        },
        rg::GaussianModel3D{
            5.3 * intensity_scale,
            0.57,
            0.02 * intensity_scale
        }) };

    const std::array<double, 2> initial_offset_list{
        0.03 * intensity_scale,
        0.07 * intensity_scale
    };
    auto analysis{ model->EditAnalysis() };
    const auto & atom_list{ model->GetSelectedAtoms() };
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        const auto current{
            rg::AtomLocalPotentialView::For(*atom_list.at(i))
                .GetGaussianResult(FittingStage::Second).mdpde.GetModel()
        };
        analysis.SetAtomLocalGaussianResult(
            FittingStage::Second,
            *atom_list.at(i),
            MakeGaussianResult(current.WithOffset(initial_offset_list.at(i))));
    }
    return model;
}

void MakeAtomSamplesSuspicious(rg::ModelObject & model, std::size_t atom_index)
{
    const auto & selected_atoms{ model.GetSelectedAtoms() };
    auto * target_atom{ selected_atoms.at(atom_index) };
    const auto target_position{ target_atom->GetPosition() };

    auto raw_sampling_entries{
        rg::AtomLocalPotentialView::For(*target_atom)
            .GetRawSamplingEntries(false)
    };
    raw_sampling_entries.resize(256);
    raw_sampling_entries.front().response = 0.0;
    raw_sampling_entries.front().point.distance = 0.0;
    raw_sampling_entries.front().point.position = target_position;
    for (std::size_t i = 1; i < raw_sampling_entries.size(); i++)
    {
        auto & sample{ raw_sampling_entries.at(i) };
        const auto response_scale{
            0.5 + 0.5 * static_cast<double>(i) /
                static_cast<double>(raw_sampling_entries.size())
        };
        sample.response = std::numeric_limits<double>::max() * response_scale;
        sample.point.position = target_position;
        sample.point.position.at(0) += 100.0;
        sample.point.distance = 100.0;
    }

    auto analysis{ model.EditAnalysis() };
    analysis.SetAtomLocalRawSamplingEntries(
        *target_atom, std::move(raw_sampling_entries));
}

std::unique_ptr<rg::ModelObject> BuildSeparatedRollbackDefenseModel()
{
    auto model{
        BuildDefenseModel(
            {
                std::array<double, 3>{ 0.0, 0.0, 0.0 },
                std::array<double, 3>{ 1.0e-4, 0.0, 0.0 },
                std::array<double, 3>{ 10.0, 0.0, 0.0 },
                std::array<double, 3>{ 10.0001, 0.0, 0.0 }
            },
            { Spot::C, Spot::O, Spot::N, Spot::CA },
            { Element::CARBON, Element::OXYGEN, Element::NITROGEN, Element::CARBON },
            {
                rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
                rg::GaussianModel3D{ 5.5, 0.55, -0.15 },
                rg::GaussianModel3D{ 6.2, 0.55, 0.18 },
                rg::GaussianModel3D{ 5.7, 0.55, -0.12 }
            },
            rg::GaussianModel3D{ 5.8, 0.55, 0.0 })
    };
    MakeAtomSamplesSuspicious(*model, 0);
    return model;
}

std::unique_ptr<rg::ModelObject> BuildBoundaryComponentConflictDefenseModel(
    double intensity_scale = 1.0)
{
    std::vector<std::array<double, 3>> position_list;
    std::vector<Spot> spot_list;
    std::vector<Element> element_list;
    std::vector<rg::GaussianModel3D> truth_model_list;
    for (std::size_t i = 0; i < 103; i++)
    {
        // Break repeated edge-weight ties while retaining a connected 101-atom chain.
        const auto atom_position{ static_cast<double>(i) };
        const auto x_position{
            i < 101 ? 0.60 * atom_position + 0.002 * atom_position * atom_position :
                110.0 + 0.45 * static_cast<double>(i - 101)
        };
        position_list.push_back({ x_position, 0.0, 0.0 });
        spot_list.push_back(i % 2 == 0 ? Spot::C : Spot::O);
        element_list.push_back(
            i % 2 == 0 ? Element::CARBON : Element::OXYGEN);
        const auto truth_model{
            i >= 101 ? rg::GaussianModel3D{ 7.5, 0.70, 0.10 } :
            i % 2 == 0 ?
                rg::GaussianModel3D{ 10.0, 0.85, 0.25 } :
                rg::GaussianModel3D{ 2.0, 0.40, -0.15 }
        };
        truth_model_list.emplace_back(rg::GaussianModel3D{
            truth_model.GetAmplitude() * intensity_scale,
            truth_model.GetWidth(),
            truth_model.GetOffset() * intensity_scale
        });
    }
    return BuildDefenseModel(
        position_list,
        spot_list,
        element_list,
        truth_model_list,
        rg::GaussianModel3D{ 5.5 * intensity_scale, 0.55, 0.0 });
}

std::unique_ptr<rg::ModelObject> BuildBoundaryJointCorrectionDefenseModel(
    double intensity_scale = 1.0)
{
    return BuildDefenseModel(
        {
            std::array<double, 3>{ 0.0, 0.0, 0.0 },
            std::array<double, 3>{ 3.4, 0.0, 0.0 }
        },
        { Spot::C, Spot::O },
        { Element::CARBON, Element::OXYGEN },
        {
            rg::GaussianModel3D{
                10.0 * intensity_scale,
                0.65,
                0.2 * intensity_scale
            },
            rg::GaussianModel3D{
                2.0 * intensity_scale,
                0.65,
                -0.1 * intensity_scale
            }
        },
        rg::GaussianModel3D{
            5.5 * intensity_scale,
            0.6,
            0.0
        });
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

double CalculateSelectedAtomResponseMeanSquaredError(
    const rg::ModelObject & model,
    std::size_t target_begin,
    std::size_t target_end)
{
    double squared_error_sum{ 0.0 };
    std::size_t sample_count{ 0 };
    const auto & selected_atoms{ model.GetSelectedAtoms() };
    for (std::size_t target_index = target_begin; target_index < target_end; target_index++)
    {
        const auto * atom{ selected_atoms.at(target_index) };
        const auto local_view{ rg::AtomLocalPotentialView::For(*atom) };
        for (const auto & sample : local_view.GetRawSamplingEntries(false))
        {
            double fitted_response{ 0.0 };
            for (const auto * fitted_atom : selected_atoms)
            {
                const auto fitted_view{
                    rg::AtomLocalPotentialView::For(*fitted_atom)
                };
                fitted_response += fitted_view.GetEstimateMDPDE(
                    FittingStage::Second).ResponseAtDistance(
                    Distance(sample.point.position, fitted_atom->GetPosition()));
            }
            const auto residual{ sample.response - fitted_response };
            squared_error_sum += residual * residual;
            sample_count++;
        }
    }
    return squared_error_sum / static_cast<double>(sample_count);
}

double CalculateSelectedAtomResponseMeanSquaredError(const rg::ModelObject & model)
{
    return CalculateSelectedAtomResponseMeanSquaredError(
        model,
        0,
        model.GetSelectedAtomCount());
}

rg::GaussianModel3D GetEstimateModel(const rg::AtomObject & atom)
{
    return rg::AtomLocalPotentialView::For(atom).GetEstimateMDPDE(
        FittingStage::Second);
}

void ExpectGaussianModelsNear(
    const rg::GaussianModel3D & actual,
    const rg::GaussianModel3D & expected,
    double tolerance)
{
    EXPECT_NEAR(actual.GetAmplitude(), expected.GetAmplitude(), tolerance);
    EXPECT_NEAR(actual.GetWidth(), expected.GetWidth(), tolerance);
    EXPECT_NEAR(actual.GetOffset(), expected.GetOffset(), tolerance);
}

void ExpectSelectedAtomEstimatesAreFinite(const rg::ModelObject & model)
{
    for (const auto * atom : model.GetSelectedAtoms())
    {
        const auto estimate{
            rg::AtomLocalPotentialView::For(*atom).GetEstimateMDPDE(
                FittingStage::Second)
        };
        EXPECT_TRUE(std::isfinite(estimate.GetAmplitude()));
        EXPECT_TRUE(std::isfinite(estimate.GetWidth()));
        EXPECT_TRUE(std::isfinite(estimate.GetOffset()));
    }
}

void ExpectPeelingSamplingEntriesMatchFinalModels(
    const rg::ModelObject & model)
{
    constexpr double neighbor_atom_search_range{ 5.0 };
    constexpr double neighbor_contribution_distance_max{ 2.5 };
    const auto & selected_atoms{ model.GetSelectedAtoms() };
    for (const auto * target_atom : selected_atoms)
    {
        const auto target_view{
            rg::AtomLocalPotentialView::For(*target_atom)
        };
        const auto raw_sampling_entries{
            target_view.GetRawSamplingEntries(false)
        };
        const auto peeling_sampling_entries{
            target_view.GetPeelingSamplingEntries(false)
        };
        ASSERT_EQ(peeling_sampling_entries.size(), raw_sampling_entries.size());
        for (std::size_t sample_index = 0;
            sample_index < raw_sampling_entries.size();
            sample_index++)
        {
            const auto & raw_sample{ raw_sampling_entries.at(sample_index) };
            const auto & peeling_sample{
                peeling_sampling_entries.at(sample_index)
            };
            auto expected_response{ raw_sample.response };
            for (const auto * neighbor_atom : selected_atoms)
            {
                if (neighbor_atom == target_atom) continue;
                if (Distance(
                        target_atom->GetPosition(),
                        neighbor_atom->GetPosition()) >
                    neighbor_atom_search_range)
                {
                    continue;
                }
                const auto sample_distance{
                    Distance(
                        raw_sample.point.position,
                        neighbor_atom->GetPosition())
                };
                if (sample_distance > neighbor_contribution_distance_max)
                {
                    continue;
                }
                expected_response -= GetEstimateModel(*neighbor_atom)
                    .ResponseAtDistance(sample_distance);
            }

            EXPECT_DOUBLE_EQ(
                peeling_sample.response,
                expected_response);
            EXPECT_DOUBLE_EQ(
                peeling_sample.point.distance,
                raw_sample.point.distance);
            EXPECT_EQ(
                peeling_sample.point.position,
                raw_sample.point.position);
            EXPECT_EQ(
                peeling_sample.point.is_selected,
                raw_sample.point.is_selected);
        }
    }
}

} // namespace

TEST(EstimatorSecondStageDefenseTest, SeedSelectionUsesLocalMdpdeThenGlobalMedian)
{
    const auto make_candidate = [](double amplitude)
    {
        return rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ amplitude, 0.5, -0.2 },
            rg::GaussianModel3DUncertainty{ 0.1, 0.02, 0.03 }
        };
    };
    const auto invalid_candidate{
        rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 0.0, 0.0, 0.0 },
            rg::GaussianModel3DUncertainty{}
        }
    };
    auto local_mdpde{ make_candidate(1.0) };
    std::optional<rg::GaussianModel3D> global_median{
        make_candidate(4.0).GetModel()
    };

    const auto expect_source = [&](seed_detail::SecondStageSeedSource source)
    {
        const auto selection{
            seed_detail::SelectSecondStageSeed(local_mdpde, global_median)
        };
        ASSERT_TRUE(selection.has_value());
        EXPECT_EQ(selection->source, source);
    };
    expect_source(seed_detail::SecondStageSeedSource::LocalMdpde);
    local_mdpde = invalid_candidate;
    expect_source(seed_detail::SecondStageSeedSource::GlobalMedian);
    global_median = invalid_candidate.GetModel();
    EXPECT_FALSE(seed_detail::SelectSecondStageSeed(
        local_mdpde,
        global_median).has_value());
}

TEST(EstimatorSecondStageDefenseTest, SeedSelectionReturnsCompleteSourceModelAndUncertainty)
{
    auto local_mdpde{
        rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 6.0, 0.55, -0.2 },
            rg::GaussianModel3DUncertainty{ 0.1, 0.02, 0.03 }
        }
    };
    const std::optional<rg::GaussianModel3D> global_median{
        rg::GaussianModel3D{ 4.0, 0.65, 0.4 }
    };

    const auto local_selection{
        seed_detail::SelectSecondStageSeed(local_mdpde, global_median)
    };
    ASSERT_TRUE(local_selection.has_value());
    EXPECT_EQ(local_selection->source, seed_detail::SecondStageSeedSource::LocalMdpde);
    EXPECT_DOUBLE_EQ(local_selection->model.GetModel().GetAmplitude(), 6.0);
    EXPECT_DOUBLE_EQ(local_selection->model.GetModel().GetWidth(), 0.55);
    EXPECT_DOUBLE_EQ(local_selection->model.GetModel().GetOffset(), -0.2);
    EXPECT_DOUBLE_EQ(
        local_selection->model.GetStandardDeviationModel().GetAmplitude(),
        0.1);
    EXPECT_DOUBLE_EQ(
        local_selection->model.GetStandardDeviationModel().GetWidth(),
        0.02);
    EXPECT_DOUBLE_EQ(
        local_selection->model.GetStandardDeviationModel().GetOffset(),
        0.03);

    local_mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 0.0, 0.0, 0.0 },
        rg::GaussianModel3DUncertainty{ 0.4, 0.5, 0.6 }
    };
    const auto median_selection{
        seed_detail::SelectSecondStageSeed(local_mdpde, global_median)
    };
    ASSERT_TRUE(median_selection.has_value());
    EXPECT_EQ(median_selection->source, seed_detail::SecondStageSeedSource::GlobalMedian);
    EXPECT_DOUBLE_EQ(median_selection->model.GetModel().GetAmplitude(), 4.0);
    EXPECT_DOUBLE_EQ(median_selection->model.GetModel().GetWidth(), 0.65);
    EXPECT_DOUBLE_EQ(median_selection->model.GetModel().GetOffset(), 0.4);
    EXPECT_DOUBLE_EQ(
        median_selection->model.GetStandardDeviationModel().GetAmplitude(),
        0.0);
    EXPECT_DOUBLE_EQ(
        median_selection->model.GetStandardDeviationModel().GetWidth(),
        0.0);
    EXPECT_DOUBLE_EQ(
        median_selection->model.GetStandardDeviationModel().GetOffset(),
        0.0);
}

TEST(EstimatorSecondStageDefenseTest, AuditObjectiveKeepsEarlierBestOnTie)
{
    const audit_detail::ObjectiveTolerance tolerance{
        1.0e-10,
        1.0e-8
    };
    EXPECT_TRUE(audit_detail::IsBetterAuditObjective(
        0.8, 1.0, tolerance));
    EXPECT_FALSE(audit_detail::IsBetterAuditObjective(
        1.2, 1.0, tolerance));
    EXPECT_FALSE(audit_detail::IsBetterAuditObjective(
        1.0 - 0.5e-8, 1.0, tolerance));
}

TEST(EstimatorSecondStageDefenseTest, BestAuditStateUpdateUsesPrecomputedObjective)
{
    audit_detail::BestAuditState audit_state;
    const audit_detail::FitState initial_state;
    const auto initial_objective{
        audit_detail::BuildObjectiveBreakdown(2.0, 0.0, 1.0)
    };
    const auto tied_objective{
        audit_detail::BuildObjectiveBreakdown(3.0, 0.0, 0.0)
    };
    const auto worse_objective{
        audit_detail::BuildObjectiveBreakdown(3.5, 0.0, 0.0)
    };
    const auto improved_objective{
        audit_detail::BuildObjectiveBreakdown(1.0, 0.0, 0.0)
    };
    ASSERT_TRUE(initial_objective.has_value());
    ASSERT_TRUE(tied_objective.has_value());
    ASSERT_TRUE(worse_objective.has_value());
    ASSERT_TRUE(improved_objective.has_value());

    EXPECT_TRUE(audit_detail::TryUpdateBestAuditState(
        initial_state,
        false,
        0,
        *initial_objective,
        audit_state));
    ASSERT_TRUE(audit_state.has_value());
    EXPECT_DOUBLE_EQ(
        audit_state->objective.GetTotalObjective(),
        initial_objective->GetTotalObjective());
    EXPECT_FALSE(audit_state->uses_polish);
    EXPECT_EQ(audit_state->source_iteration, 0U);

    EXPECT_FALSE(audit_detail::TryUpdateBestAuditState(
        initial_state,
        false,
        0,
        *tied_objective,
        audit_state));
    EXPECT_FALSE(audit_detail::TryUpdateBestAuditState(
        initial_state,
        false,
        0,
        *worse_objective,
        audit_state));

    EXPECT_TRUE(audit_detail::TryUpdateBestAuditState(
        initial_state,
        true,
        7,
        *improved_objective,
        audit_state));
    ASSERT_TRUE(audit_state.has_value());
    EXPECT_DOUBLE_EQ(
        audit_state->objective.GetTotalObjective(),
        improved_objective->GetTotalObjective());
    EXPECT_TRUE(audit_state->uses_polish);
    EXPECT_EQ(audit_state->source_iteration, 7U);
    // A background refresh must rescore both retained best and previous under
    // the same fixed domain; an old zero score must not block a new exact fit.
    audit_detail::SecondStageContext context;
    context.atom_list.resize(1);
    context.at(0).neighbor_atom_sample_offset_list = { 0, 0, 0, 0 };
    const audit_detail::FitState seed{ MakeGaussianResult({ 4.0, 0.5, 0.0 }) };
    const audit_detail::FitState earlier_best{ MakeGaussianResult({ 6.0, 0.5, 0.0 }) };
    const audit_detail::FitState previous{ MakeGaussianResult({ 5.0, 0.5, 0.0 }) };
    for (const double distance : { 0.15, 0.35, 0.60 })
    {
        context.at(0).raw_sampling_entries.emplace_back(LocalPotentialSample{
            2.0 * previous.front().mdpde.GetModel().ResponseAtDistance(distance),
            SamplingPoint{ distance } });
        context.at(0).unselected_distance_list_by_sample.push_back({ distance });
    }
    context.frozen_background = audit_detail::BuildFrozenBackground(context, seed, { { 0 } });
    ASSERT_TRUE(context.frozen_background);
    const auto old_snapshot{ audit_detail::BuildSecondStageModelSnapshot(context, earlier_best) };
    const auto domain{ audit_detail::BuildObjectiveDomain(context, old_snapshot, { { 0 } }, 0.0, 1.0) };
    const auto old_score{ audit_detail::EvaluateAuditObjective(
        domain, audit_detail::SnapshotResidualEvaluator{ context, old_snapshot }) };
    ASSERT_TRUE(old_score.has_value());
    audit_state.reset();
    ASSERT_TRUE(audit_detail::TryUpdateBestAuditState(earlier_best, true, 3, *old_score, audit_state));
    context.frozen_background = audit_detail::BuildFrozenBackground(context, previous, { { 0 } });
    ASSERT_TRUE(context.frozen_background);
    audit_detail::ReevaluateBestAuditState(context, domain, audit_state);
    ASSERT_TRUE(audit_state.has_value());
    EXPECT_GT(audit_state->objective.GetTotalObjective(), old_score->GetTotalObjective());
    EXPECT_EQ(audit_state->source_iteration, 3U);
    EXPECT_TRUE(audit_state->uses_polish);
    ExpectGaussianModelsNear(audit_state->state.front().mdpde.GetModel(), earlier_best.front().mdpde.GetModel(), 0.0);
    const auto refreshed_snapshot{ audit_detail::BuildSecondStageModelSnapshot(context, previous) };
    const auto refreshed_score{ audit_detail::EvaluateAuditObjective(
        domain, audit_detail::SnapshotResidualEvaluator{ context, refreshed_snapshot }) };
    ASSERT_TRUE(refreshed_score.has_value());
    EXPECT_TRUE(audit_detail::TryUpdateBestAuditState(previous, false, 4, *refreshed_score, audit_state));
    EXPECT_EQ(audit_state->source_iteration, 4U);
    EXPECT_NEAR(audit_state->objective.GetTotalObjective(), 0.0, 1.0e-12);

}

TEST(EstimatorSecondStageDefenseTest, AuditObjectiveProgressGuardChecksPreviousAndBest)
{
    const audit_detail::ObjectiveTolerance tolerance{
        1.0e-8,
        1.0e-3
    };
    const audit_detail::ObjectiveBreakdown best_one{ 1.0, 0.0, 0.0 };
    const audit_detail::ObjectiveBreakdown best_below{ 0.99, 0.0, 0.0 };
    EXPECT_TRUE(audit_detail::IsAuditObjectiveAcceptableForProgress(
        1.0005, 1.0, &best_one, tolerance));
    EXPECT_FALSE(audit_detail::IsAuditObjectiveAcceptableForProgress(
        1.002, 1.0, &best_one, tolerance));
    EXPECT_FALSE(audit_detail::IsAuditObjectiveAcceptableForProgress(
        1.0, 1.0, &best_below, tolerance));
}

TEST(EstimatorSecondStageDefenseTest, AuditToleranceUsesAbsolutePlusRelativeReference)
{
    const audit_detail::ObjectiveTolerance tolerance{
        1.0e-8,
        1.0e-3
    };
    EXPECT_TRUE(audit_detail::IsAuditObjectiveAcceptableForProgress(
        1.0e-8,
        0.0,
        nullptr,
        tolerance));
    EXPECT_TRUE(audit_detail::IsAuditObjectiveAcceptableForProgress(
        -2.0 + 1.0e-8 + 2.0e-3,
        -2.0,
        nullptr,
        tolerance));
    const auto boundary{ 2.0 + 1.0e-8 + 2.0e-3 };
    EXPECT_TRUE(audit_detail::IsAuditObjectiveAcceptableForProgress(
        boundary,
        2.0,
        nullptr,
        tolerance));
    EXPECT_FALSE(audit_detail::IsAuditObjectiveAcceptableForProgress(
        boundary + 1.0e-9,
        2.0,
        nullptr,
        tolerance));
}

TEST(EstimatorSecondStageDefenseTest, ScientificObjectiveUsesFitTailAndOffsetOnly)
{
    const auto objective{
        audit_detail::BuildObjectiveBreakdown(
            0.4,
            2.0,
            0.18)
    };
    const auto previous{
        audit_detail::BuildObjectiveBreakdown(
            1.4,
            0.0,
            0.0)
    };

    ASSERT_TRUE(objective.has_value());
    ASSERT_TRUE(previous.has_value());
    const auto empty_tail{
        audit_detail::BuildObjectiveBreakdown(
            0.4,
            0.0,
            0.0)
    };
    ASSERT_TRUE(empty_tail.has_value());
    EXPECT_DOUBLE_EQ(empty_tail->tail_validation_loss, 0.0);
    EXPECT_DOUBLE_EQ(empty_tail->GetTailValidationPenalty(), 0.0);
    EXPECT_DOUBLE_EQ(objective->tail_validation_loss, 2.0);
    EXPECT_DOUBLE_EQ(objective->GetTailValidationPenalty(), 0.5);
    EXPECT_DOUBLE_EQ(objective->offset_plausibility_penalty, 0.18);
    EXPECT_DOUBLE_EQ(objective->GetTotalObjective(), 1.08);
    EXPECT_DOUBLE_EQ(
        objective->GetTotalObjective(),
        objective->fit_range_residual_objective +
            objective->GetTailValidationPenalty() +
            objective->offset_plausibility_penalty);
    EXPECT_TRUE(audit_detail::IsBetterAuditObjective(
        objective->GetTotalObjective(),
        previous->GetTotalObjective(),
        audit_detail::ObjectiveTolerance{ 1.0e-10, 1.0e-8 }));
    EXPECT_FALSE(
        audit_detail::BuildObjectiveBreakdown(
            std::numeric_limits<double>::infinity(),
            2.0,
            0.18).has_value());
    EXPECT_FALSE(
        audit_detail::BuildObjectiveBreakdown(
            std::numeric_limits<double>::max(),
            0.0,
            std::numeric_limits<double>::max()).has_value());
}

TEST(EstimatorSecondStageDefenseTest, GlobalObjectiveWeightsClustersByAtomCount)
{
    EXPECT_DOUBLE_EQ(
        audit_detail::CalculateClusterAtomWeight(1, 4),
        0.25);
    EXPECT_DOUBLE_EQ(
        audit_detail::CalculateClusterAtomWeight(3, 4),
        0.75);
    EXPECT_DOUBLE_EQ(
        0.25 * 2.0 + 0.75 * 6.0,
        5.0);
}

TEST(EstimatorSecondStageDefenseTest, ObjectiveClusterStateLifecycleReconcilesPartition)
{
    const audit_detail::ClusterKey existing_key{ 0 };
    const audit_detail::ClusterKey new_key{ 1 };
    coupling_detail::CouplingGraphPartition partition;
    partition.sample_id_list_by_key.emplace(
        existing_key,
        std::vector<coupling_detail::SampleRef>{ { 0, 0 } });
    partition.sample_id_list_by_key.emplace(
        new_key,
        std::vector<coupling_detail::SampleRef>{ { 1, 0 } });

    const audit_detail::ObjectiveBreakdown previous_breakdown{
        1.0, 2.0, 0.0
    };
    const audit_detail::ObjectiveBreakdown existing_best{
        0.5, 0.5, 0.0
    };
    audit_detail::ObjectiveByKey previous_objective_by_key;
    previous_objective_by_key.emplace(existing_key, previous_breakdown);
    previous_objective_by_key.emplace(new_key, previous_breakdown);

    audit_detail::ClusterObjectiveStateMap state_by_key;
    state_by_key.emplace(
        existing_key,
        audit_detail::ClusterObjectiveState{
            existing_best,
            0.25
        });

    audit_detail::ReconcileClusterObjectiveState(
        previous_objective_by_key,
        state_by_key);

    ASSERT_EQ(state_by_key.size(), 2U);
    ASSERT_TRUE(state_by_key.at(existing_key).best_objective.has_value());
    EXPECT_DOUBLE_EQ(
        state_by_key.at(existing_key).best_objective->GetTotalObjective(),
        existing_best.GetTotalObjective());
    ASSERT_TRUE(state_by_key.at(new_key).best_objective.has_value());
    EXPECT_DOUBLE_EQ(
        state_by_key.at(new_key).best_objective->GetTotalObjective(),
        previous_breakdown.GetTotalObjective());
    EXPECT_DOUBLE_EQ(
        state_by_key.at(new_key).best_maximum_transformed_change,
        0.0);
}

TEST(EstimatorSecondStageDefenseTest, TrustRegionStateReconcilesShrinksGrowsAndSaturates)
{
    using Action = trust_detail::TrustRegionRadiusAction;
    trust_detail::ObjectiveAttemptDiagnostic accepted_diagnostic;
    accepted_diagnostic.accepted_factor = 0.5;
    accepted_diagnostic.trust_region_radius = 1.0;
    accepted_diagnostic.trust_region_step_norm = 1.0;
    accepted_diagnostic.previous_objective =
        trust_detail::BuildObjectiveBreakdown(2.0, 0.0, 0.0);
    accepted_diagnostic.candidate_objective =
        trust_detail::BuildObjectiveBreakdown(1.999, 0.0, 0.0);

    EXPECT_EQ(
        trust_detail::DetermineAcceptedTrustRegionRadiusAction(
            0.5, accepted_diagnostic),
        Action::Keep);

    accepted_diagnostic.candidate_objective =
        trust_detail::BuildObjectiveBreakdown(1.997, 0.0, 0.0);
    EXPECT_EQ(
        trust_detail::DetermineAcceptedTrustRegionRadiusAction(
            0.5, accepted_diagnostic),
        Action::Grow);

    accepted_diagnostic.accepted_factor = 0.25;
    EXPECT_EQ(
        trust_detail::DetermineAcceptedTrustRegionRadiusAction(
            0.5, accepted_diagnostic),
        Action::Shrink);

    accepted_diagnostic.accepted_factor = 0.5;
    EXPECT_EQ(
        trust_detail::DetermineAcceptedTrustRegionRadiusAction(
            0.5, accepted_diagnostic),
        Action::Grow);

#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
    trust_detail::TrustModelShadowDiagnostic shadow{
        .status = trust_detail::TrustModelPredictionStatus::Available,
        .rho = 0.20,
        .boundary_utilization = 1.0,
        .current_action = Action::Grow
    };
    EXPECT_EQ(
        trust_detail::DetermineTrustModelShadowAction(shadow),
        Action::Shrink);
    shadow.rho = 0.50;
    EXPECT_EQ(
        trust_detail::DetermineTrustModelShadowAction(shadow),
        Action::Keep);
    shadow.rho = 0.90;
    shadow.boundary_utilization = 0.79;
    EXPECT_EQ(
        trust_detail::DetermineTrustModelShadowAction(shadow),
        Action::Keep);
    shadow.boundary_utilization = 0.80;
    EXPECT_EQ(
        trust_detail::DetermineTrustModelShadowAction(shadow),
        Action::Grow);
    shadow.objective_backtracked = true;
    EXPECT_EQ(
        trust_detail::DetermineTrustModelShadowAction(shadow),
        Action::Shrink);
    shadow.objective_backtracked = false;
    shadow.status = trust_detail::TrustModelPredictionStatus::NonmaterialPrediction;
    shadow.current_action = Action::Grow;
    EXPECT_EQ(
        trust_detail::DetermineTrustModelShadowAction(shadow),
        Action::Grow);
#endif

    trust_detail::TrustRegionStateSet state;
    const trust_detail::ClusterKey key{ 0 };
    state.Reconcile({ key });
    EXPECT_DOUBLE_EQ(state.GetRadius(key), 1.0);

    for (const auto expected : { 0.5, 0.25, 0.125, 0.0625 })
    {
        const auto update{
            state.ApplyRadiusUpdates({}, { key }, {}, {})
        };
        EXPECT_EQ(
            update.changed_key_list,
            std::vector<trust_detail::ClusterKey>{ key });
        EXPECT_TRUE(update.saturated_key_list.empty());
        EXPECT_DOUBLE_EQ(state.GetRadius(key), expected);
    }
    const auto saturated{
        state.ApplyRadiusUpdates({}, { key }, {}, {})
    };
    EXPECT_TRUE(saturated.changed_key_list.empty());
    EXPECT_EQ(
        saturated.saturated_key_list,
        std::vector<trust_detail::ClusterKey>{ key });

    state.ApplyRadiusUpdates({ key }, {}, {}, {});
    EXPECT_DOUBLE_EQ(state.GetRadius(key), 0.125);

    const trust_detail::ClusterKey replacement_key{ 1 };
    state.Reconcile({ replacement_key });
    EXPECT_DOUBLE_EQ(state.GetRadius(replacement_key), 1.0);
}

#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
TEST(EstimatorSecondStageDefenseTest, TrustModelShadowUsesFrozenIrlsDirectionalPrediction)
{
    using Action = trust_detail::TrustRegionRadiusAction;
    const rg::GaussianModel3D previous_model{ 5.4, 0.52, 0.05 };
    const rg::GaussianModel3D target_model{ 6.0, 0.55, 0.10 };
    auto fixture{
        BuildJointPolishFixture( { previous_model }, { target_model })
    };
    constexpr double unselected_distance{ 0.25 };
    auto & atom_context{ fixture.context.at(0) };
    atom_context.unselected_distance_list_by_sample.assign(atom_context.raw_sampling_entries.size(),
        { unselected_distance });
    for (auto & sample : atom_context.raw_sampling_entries)
        sample.response += previous_model.ResponseAtDistance(unselected_distance);
    const trust_detail::ClusterKey key{ 0 };
    fixture.context.frozen_background = trust_detail::BuildFrozenBackground(fixture.context, fixture.state, { key });
    ASSERT_TRUE(fixture.context.frozen_background);
    const auto previous_snapshot{
        trust_detail::BuildSecondStageModelSnapshot(
            fixture.context,
            fixture.state)
    };
    const auto objective_domain{
        trust_detail::BuildObjectiveDomain(
            fixture.context,
            previous_snapshot,
            { key },
            0.0,
            1.0)
    };
    const auto residual_baseline{
        trust_detail::BuildResidualBaseline(fixture.context, fixture.state)
    };
    const auto previous_objective{
        trust_detail::EvaluateAuditObjective(
            objective_domain,
            residual_baseline)
    };
    ASSERT_TRUE(previous_objective.has_value());

    const auto previous_coordinates{
        previous_model.ToTransformedCoordinates()
    };
    const auto target_coordinates{
        target_model.ToTransformedCoordinates()
    };
    ASSERT_TRUE(previous_coordinates.has_value());
    ASSERT_TRUE(target_coordinates.has_value());
    const auto candidate_model{
        rg::GaussianModel3D::FromTransformedCoordinates(
            *previous_coordinates + 0.05 *
                (*target_coordinates - *previous_coordinates))
    };
    ASSERT_TRUE(candidate_model.has_value());
    auto candidate_state{ fixture.state };
    candidate_state.at(0) = MakeGaussianResult(*candidate_model);
    const auto candidate_patch{
        trust_detail::FitStatePatch::FromState(candidate_state, key)
    };
    const auto candidate_snapshot{
        trust_detail::BuildSecondStageModelSnapshot(
            fixture.context,
            candidate_state)
    };
    const auto candidate_objective{
        trust_detail::EvaluateAuditObjective(
            objective_domain,
            trust_detail::SnapshotResidualEvaluator{
                fixture.context,
                candidate_snapshot
            })
    };
    ASSERT_TRUE(candidate_objective.has_value());

    const auto diagnostic{
        trust_detail::EvaluateTrustModelShadow(
            fixture.context,
            residual_baseline,
            fixture.state,
            candidate_patch,
            key,
            fixture.sample_ref_list,
            objective_domain,
            previous_objective,
            candidate_objective,
            1.0,
            Action::Keep,
            trust_detail::TrustModelCandidateSource::Base,
            false)
    };
    EXPECT_EQ(
        diagnostic.status,
        trust_detail::TrustModelPredictionStatus::Available);
    ASSERT_TRUE(diagnostic.actual_reduction.has_value());
    ASSERT_TRUE(diagnostic.predicted_reduction.has_value());
    ASSERT_TRUE(diagnostic.rho.has_value());
    EXPECT_GT(*diagnostic.actual_reduction, 0.0);
    EXPECT_GT(*diagnostic.predicted_reduction, 0.0);
    EXPECT_NEAR(*diagnostic.rho, 1.0, 0.10);
    EXPECT_EQ(candidate_state.size(), 1U);

    constexpr double intensity_scale{ 1.0e4 };
    auto scaled_fixture{ fixture };
    auto scaled_background{ std::make_shared<trust_detail::FrozenBackground>(*fixture.context.frozen_background) };
    for (auto & responses : scaled_background->response_by_atom)
        for (auto & response : responses) response *= intensity_scale;
    for (auto & model : scaled_background->model_by_atom)
        model = { model.GetAmplitude() * intensity_scale, model.GetWidth(), model.GetOffset() * intensity_scale };
    scaled_fixture.context.frozen_background = std::move(scaled_background);

    for (auto & sample : scaled_fixture.context.at(0).raw_sampling_entries)
    {
        sample.response *= intensity_scale;
    }
    const auto scale_model = [](const rg::GaussianModel3D & model)
    {
        return rg::GaussianModel3D{
            model.GetAmplitude() * intensity_scale,
            model.GetWidth(),
            model.GetOffset() * intensity_scale
        };
    };
    scaled_fixture.state.at(0) = MakeGaussianResult(scale_model(previous_model));
    auto scaled_candidate_state{ scaled_fixture.state };
    scaled_candidate_state.at(0) = MakeGaussianResult(scale_model(*candidate_model));
    const auto scaled_previous_snapshot{
        trust_detail::BuildSecondStageModelSnapshot(
            scaled_fixture.context,
            scaled_fixture.state)
    };
    const auto scaled_domain{
        trust_detail::BuildObjectiveDomain(
            scaled_fixture.context,
            scaled_previous_snapshot,
            { key },
            0.0,
            1.0)
    };
    const auto scaled_baseline{
        trust_detail::BuildResidualBaseline(
            scaled_fixture.context,
            scaled_fixture.state)
    };
    const auto scaled_previous_objective{
        trust_detail::EvaluateAuditObjective(scaled_domain, scaled_baseline)
    };
    const auto scaled_candidate_snapshot{
        trust_detail::BuildSecondStageModelSnapshot(
            scaled_fixture.context,
            scaled_candidate_state)
    };
    const auto scaled_candidate_objective{
        trust_detail::EvaluateAuditObjective(
            scaled_domain,
            trust_detail::SnapshotResidualEvaluator{
                scaled_fixture.context,
                scaled_candidate_snapshot
            })
    };
    ASSERT_TRUE(scaled_previous_objective.has_value());
    ASSERT_TRUE(scaled_candidate_objective.has_value());
    const auto scaled_diagnostic{
        trust_detail::EvaluateTrustModelShadow(
            scaled_fixture.context,
            scaled_baseline,
            scaled_fixture.state,
            trust_detail::FitStatePatch::FromState(
                scaled_candidate_state,
                key),
            key,
            scaled_fixture.sample_ref_list,
            scaled_domain,
            scaled_previous_objective,
            scaled_candidate_objective,
            1.0,
            Action::Keep,
            trust_detail::TrustModelCandidateSource::Base,
            false)
    };
    ASSERT_TRUE(scaled_diagnostic.predicted_reduction.has_value());
    ASSERT_TRUE(scaled_diagnostic.rho.has_value());
    EXPECT_NEAR(
        *scaled_diagnostic.predicted_reduction,
        *diagnostic.predicted_reduction,
        1.0e-6);
    EXPECT_NEAR(*scaled_diagnostic.rho, *diagnostic.rho, 1.0e-6);

    auto previous_with_penalty{ *previous_objective };
    auto candidate_with_penalty{ *candidate_objective };
    previous_with_penalty.offset_plausibility_penalty = 0.20;
    candidate_with_penalty.offset_plausibility_penalty = 0.05;
    const auto penalty_diagnostic{
        trust_detail::EvaluateTrustModelShadow(
            fixture.context,
            residual_baseline,
            fixture.state,
            candidate_patch,
            key,
            fixture.sample_ref_list,
            objective_domain,
            previous_with_penalty,
            candidate_with_penalty,
            1.0,
            Action::Keep,
            trust_detail::TrustModelCandidateSource::Polish,
            false)
    };
    ASSERT_TRUE(penalty_diagnostic.predicted_reduction.has_value());
    EXPECT_NEAR(
        *penalty_diagnostic.predicted_reduction -
            *diagnostic.predicted_reduction,
        0.15,
        1.0e-12);

    const auto tail_domain{
        trust_detail::BuildObjectiveDomain(
            fixture.context,
            previous_snapshot,
            { key },
            0.0,
            0.30)
    };
    const auto tail_previous_objective{
        trust_detail::EvaluateAuditObjective(
            tail_domain,
            residual_baseline)
    };
    const auto tail_candidate_objective{
        trust_detail::EvaluateAuditObjective(
            tail_domain,
            trust_detail::SnapshotResidualEvaluator{
                fixture.context,
                candidate_snapshot
            })
    };
    ASSERT_TRUE(tail_previous_objective.has_value());
    ASSERT_TRUE(tail_candidate_objective.has_value());
    const auto tail_diagnostic{
        trust_detail::EvaluateTrustModelShadow(
            fixture.context,
            residual_baseline,
            fixture.state,
            candidate_patch,
            key,
            fixture.sample_ref_list,
            tail_domain,
            tail_previous_objective,
            tail_candidate_objective,
            1.0,
            Action::Keep,
            trust_detail::TrustModelCandidateSource::Base,
            false)
    };
    EXPECT_EQ(
        tail_diagnostic.status,
        trust_detail::TrustModelPredictionStatus::Available);
    ASSERT_TRUE(tail_diagnostic.rho.has_value());
    EXPECT_NEAR(*tail_diagnostic.rho, 1.0, 0.10);

    auto nonmaterial_prediction_domain{ objective_domain };
    nonmaterial_prediction_domain.cluster_by_key.at(key).scale->fit *= 1.0e6;
    const auto nonmaterial_prediction{
        trust_detail::EvaluateTrustModelShadow(
            fixture.context,
            residual_baseline,
            fixture.state,
            candidate_patch,
            key,
            fixture.sample_ref_list,
            nonmaterial_prediction_domain,
            previous_objective,
            candidate_objective,
            1.0,
            Action::Keep,
            trust_detail::TrustModelCandidateSource::Base,
            false)
    };
    EXPECT_EQ(
        nonmaterial_prediction.status,
        trust_detail::TrustModelPredictionStatus::NonmaterialPrediction);
    EXPECT_FALSE(nonmaterial_prediction.rho.has_value());

    auto nonfinite_baseline{ residual_baseline };
    nonfinite_baseline.sample_list.at(0).at(0)->residual =
        std::numeric_limits<double>::infinity();
    const auto nonfinite_prediction{
        trust_detail::EvaluateTrustModelShadow(
            fixture.context,
            nonfinite_baseline,
            fixture.state,
            candidate_patch,
            key,
            fixture.sample_ref_list,
            objective_domain,
            previous_objective,
            candidate_objective,
            1.0,
            Action::Keep,
            trust_detail::TrustModelCandidateSource::Base,
            false)
    };
    EXPECT_EQ(
        nonfinite_prediction.status,
        trust_detail::TrustModelPredictionStatus::Nonfinite);
    EXPECT_FALSE(nonfinite_prediction.rho.has_value());

    const auto opposite_model{
        rg::GaussianModel3D::FromTransformedCoordinates(
            *previous_coordinates - 0.05 *
                (*target_coordinates - *previous_coordinates))
    };
    ASSERT_TRUE(opposite_model.has_value());
    auto opposite_state{ fixture.state };
    opposite_state.at(0) = MakeGaussianResult(*opposite_model);
    const auto opposite_patch{
        trust_detail::FitStatePatch::FromState(opposite_state, key)
    };
    const auto opposite_snapshot{
        trust_detail::BuildSecondStageModelSnapshot(
            fixture.context,
            opposite_state)
    };
    const auto opposite_objective{
        trust_detail::EvaluateAuditObjective(
            objective_domain,
            trust_detail::SnapshotResidualEvaluator{
                fixture.context,
                opposite_snapshot
            })
    };
    ASSERT_TRUE(opposite_objective.has_value());
    const auto nonpositive_prediction{
        trust_detail::EvaluateTrustModelShadow(
            fixture.context,
            residual_baseline,
            fixture.state,
            opposite_patch,
            key,
            fixture.sample_ref_list,
            objective_domain,
            previous_objective,
            opposite_objective,
            1.0,
            Action::Keep,
            trust_detail::TrustModelCandidateSource::Base,
            false)
    };
    EXPECT_EQ(
        nonpositive_prediction.status,
        trust_detail::TrustModelPredictionStatus::NonpositivePrediction);
    EXPECT_FALSE(nonpositive_prediction.rho.has_value());

    const auto nonmaterial{
        trust_detail::EvaluateTrustModelShadow(
            fixture.context,
            residual_baseline,
            fixture.state,
            trust_detail::FitStatePatch::FromState(fixture.state, key),
            key,
            fixture.sample_ref_list,
            objective_domain,
            previous_objective,
            previous_objective,
            1.0,
            Action::Keep,
            trust_detail::TrustModelCandidateSource::Base,
            false)
    };
    EXPECT_EQ(
        nonmaterial.status,
        trust_detail::TrustModelPredictionStatus::NonmaterialStep);
    EXPECT_FALSE(nonmaterial.rho.has_value());
}

#endif

TEST(EstimatorSecondStageDefenseTest, ExhaustedRejectionsAreExcludedFromRadiusShrink)
{
    using Key = trust_detail::ClusterKey;
    const Key grow_key{ 0 };
    const Key shrink_key{ 1 };

    trust_detail::TrustRegionStateSet state;
    state.Reconcile({ grow_key, shrink_key });
    const auto update{ state.ApplyRadiusUpdates(
        { grow_key },
        {},
        { shrink_key },
        { shrink_key }) };

    EXPECT_DOUBLE_EQ(state.GetRadius(grow_key), 2.0);
    EXPECT_DOUBLE_EQ(state.GetRadius(shrink_key), 1.0);
    EXPECT_TRUE(update.changed_key_list.empty());
    EXPECT_TRUE(update.saturated_key_list.empty());
}

TEST(EstimatorSecondStageDefenseTest, TerminalRejectionShrinksRadiusOncePerControllerUpdate)
{
    using Key = trust_detail::ClusterKey;
    const Key first_key{ 0 };
    const Key second_key{ 1 };
    trust_detail::TrustRegionStateSet state;
    state.Reconcile({ first_key, second_key });

    const auto update{ state.ApplyRadiusUpdates(
        {},
        {},
        { first_key, second_key },
        {}) };

    EXPECT_EQ(
        update.changed_key_list,
        (std::vector<Key>{ first_key, second_key }));
    EXPECT_DOUBLE_EQ(state.GetRadius(first_key), 0.5);
    EXPECT_DOUBLE_EQ(state.GetRadius(second_key), 0.5);
}

TEST(EstimatorSecondStageDefenseTest, JointFittingConditioningDetectsJointDependence)
{
    Eigen::SparseMatrix<double> design_matrix{ 3, 3 };
    const std::vector<Eigen::Triplet<double>> entries{
        { 0, 0, 1.0 }, { 0, 2, 1.0 },
        { 1, 1, 1.0 }, { 1, 2, 1.0 },
        { 2, 0, 1.0 }, { 2, 1, 1.0 }, { 2, 2, 2.0 }
    };
    design_matrix.setFromTriplets(entries.begin(), entries.end());

    const Eigen::MatrixXd dense{ design_matrix };
    for (Eigen::Index left = 0; left < dense.cols(); left++)
    {
        for (Eigen::Index right = left + 1; right < dense.cols(); right++)
        {
            const auto overlap{
                std::abs(dense.col(left).dot(dense.col(right))) /
                (dense.col(left).norm() * dense.col(right).norm())
            };
            EXPECT_LT(overlap, 0.98);
        }
    }

    const auto diagnostics{
        conditioning_detail::EvaluateJointFittingConditioning(
            design_matrix,
            1.0e-8)
    };
    EXPECT_TRUE(diagnostics.guard_required);
    EXPECT_LE(diagnostics.pivot_ratio, 1.0e-8);
}

TEST(EstimatorSecondStageDefenseTest, JointFittingConditioningKeepsIndependentColumns)
{
    Eigen::SparseMatrix<double> design_matrix{ 3, 3 };
    design_matrix.setIdentity();

    const auto diagnostics{
        conditioning_detail::EvaluateJointFittingConditioning(
            design_matrix,
            1.0e-8)
    };
    EXPECT_FALSE(diagnostics.guard_required);
    EXPECT_NEAR(diagnostics.pivot_ratio, 1.0, 1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetEstimatorPreservesIndividualRidgeAnchors)
{
    const std::vector<rg::GaussianModel3D> models{
        { 6.0, 0.55, 1.0 }, { 6.0, 0.55, 6.0 } };
    auto fixture{ BuildJointOffsetEstimationFixture(models, { 1.0, 6.0 }) };
    // Identical columns constrain only the sum; ridge must retain each atom's own seed.
    for (std::size_t atom = 0; atom < 2; atom++)
    {
        fixture.first.at(atom).raw_sampling_entries.front().response +=
            models.at(1 - atom).ResponseAtDistance(0.35);
        fixture.first.at(atom).neighbor_atom_sample_list = { { 1 - atom, 0.35 } };
        fixture.first.at(atom).neighbor_atom_sample_offset_list = { 0, 1 };
    }
    alg::WeightedRidgeSolver solver;
    const auto result{ offset_detail::EstimateJointOffsets(
        fixture.first, { 0, 1 }, fixture.second, { 1.0, 1.0 }, solver, false) };
    ASSERT_EQ(result.status, offset_detail::JointOffsetSolveStatus::Converged);
    ASSERT_EQ(result.offset.size(), 2);
    EXPECT_NEAR(result.offset(0), 1.0, 1.0e-10);
    EXPECT_NEAR(result.offset(1), 6.0, 1.0e-10);
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetEstimatorMapsPermutedAtomColumns)
{
    const std::vector<double> target_offsets{ 1.15, 3.9, 3.2, -0.2 };
    auto fixture{ BuildJointOffsetEstimationFixture(
        { { 6.0, 0.55, 1.0 }, { 7.0, 0.60, 4.0 },
            { 8.0, 0.65, 3.0 }, { 4.0, 0.50, -0.2 } },
        target_offsets) };
    // Each active target sees two active neighbors and one fixed, outside-cluster atom.
    for (std::size_t atom_index = 0; atom_index < 3; atom_index++)
    {
        auto & atom{ fixture.first.at(atom_index) };
        for (std::size_t neighbor_index = 0; neighbor_index < 4; neighbor_index++)
        {
            if (neighbor_index == atom_index) continue;
            const auto distance{ 0.7 + 0.2 * static_cast<double>(neighbor_index) +
                0.1 * static_cast<double>(atom_index) };
            atom.neighbor_atom_sample_list.push_back({ neighbor_index, distance });
            atom.raw_sampling_entries.front().response += fixture.second.node.at(neighbor_index)
                .WithOffset(target_offsets.at(neighbor_index)).ResponseAtDistance(distance);
        }
        atom.neighbor_atom_sample_offset_list = { 0, atom.neighbor_atom_sample_list.size() };
    }
    alg::WeightedRidgeSolver solver;
    const auto original{ offset_detail::EstimateJointOffsets(
        fixture.first, { 0, 1, 2 }, fixture.second, { 1.0, 2.0, 4.0, 9.0 }, solver, false) };
    const auto reordered{ offset_detail::EstimateJointOffsets(
        fixture.first, { 1, 2, 0 }, fixture.second, { 1.0, 2.0, 4.0, 9.0 }, solver, false) };
    ASSERT_EQ(original.status, offset_detail::JointOffsetSolveStatus::Converged);
    ASSERT_EQ(reordered.status, offset_detail::JointOffsetSolveStatus::Converged);
    ASSERT_EQ(original.offset.size(), 3);
    ASSERT_EQ(reordered.offset.size(), 3);
    EXPECT_NEAR(reordered.offset(0), original.offset(1), 1.0e-12);
    EXPECT_NEAR(reordered.offset(1), original.offset(2), 1.0e-12);
    EXPECT_NEAR(reordered.offset(2), original.offset(0), 1.0e-12);
    for (auto & atom : fixture.first.atom_list)
    {
        std::ranges::reverse(atom.neighbor_atom_sample_list);
    }
    const auto reversed_neighbors{ offset_detail::EstimateJointOffsets(
        fixture.first, { 0, 1, 2 }, fixture.second, { 1.0, 2.0, 4.0, 9.0 }, solver, false) };
    ASSERT_EQ(reversed_neighbors.status, original.status);
    ASSERT_EQ(reversed_neighbors.offset.size(), original.offset.size());
    for (Eigen::Index column = 0; column < original.offset.size(); column++)
    {
        EXPECT_NEAR(reversed_neighbors.offset(column), original.offset(column), 1.0e-12);
        const auto atom_index{ static_cast<std::size_t>(column) };
        EXPECT_LT(std::abs(original.offset(column) - target_offsets.at(atom_index)),
            std::abs(fixture.second.node.at(atom_index).GetOffset() - target_offsets.at(atom_index)));
    }
}

TEST(EstimatorSecondStageDefenseTest, LocalRefitHealthTracksSolverQualification)
{
    EXPECT_TRUE(health_detail::IsLocalRefitStatusSolverQualified(
        rg::RHBMEstimationStatus::SUCCESS));
    EXPECT_FALSE(health_detail::IsLocalRefitStatusSolverQualified(
        rg::RHBMEstimationStatus::MAX_ITERATIONS_REACHED));
    for (const auto status : {
        rg::RHBMEstimationStatus::NUMERICAL_FALLBACK,
        rg::RHBMEstimationStatus::INSUFFICIENT_DATA,
        rg::RHBMEstimationStatus::SINGLE_MEMBER })
    {
        EXPECT_FALSE(health_detail::IsLocalRefitStatusSolverQualified(status));
    }
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetHealthSeparatesHardFailureFromSolverQualification)
{
    using Status = offset_detail::JointOffsetSolveStatus;

    EXPECT_TRUE(health_detail::ClusterHealth{ Status::Converged }.IsSolverQualified());
    EXPECT_FALSE(offset_detail::IsJointOffsetSolveHardFailure(Status::Converged));

    for (const auto status : {
        Status::IrlsObjectiveDeteriorated,
        Status::IrlsMaximumIterationsReached })
    {
        EXPECT_FALSE(health_detail::ClusterHealth{ status }.IsSolverQualified());
        EXPECT_FALSE(offset_detail::IsJointOffsetSolveHardFailure(status));
    }

    for (const auto status : {
        Status::SystemBuildFailed,
        Status::EmptySystem,
        Status::InitialSolveFailed,
        Status::IrlsSolveFailed })
    {
        EXPECT_FALSE(health_detail::ClusterHealth{ status }.IsSolverQualified());
        EXPECT_TRUE(offset_detail::IsJointOffsetSolveHardFailure(status));
    }
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetEstimatorFitsIndependentAtomOffsets)
{
    auto fixture{
        BuildJointOffsetEstimationFixture(
            {
                rg::GaussianModel3D{ 6.0, 0.55, 0.0 },
                rg::GaussianModel3D{ 7.0, 0.60, 4.0 }
            },
            { 1.0, 3.0 })
    };
    alg::WeightedRidgeSolver solver;
    const auto result{
        offset_detail::EstimateJointOffsets(
            fixture.first,
            { 0, 1 },
            fixture.second,
            { 1.0, 1.0 },
            solver,
            false)
    };

    EXPECT_EQ(
        result.status,
        offset_detail::JointOffsetSolveStatus::Converged);
    ASSERT_EQ(result.offset.size(), 2);
    EXPECT_NEAR(result.offset(0), 1.0, 0.01);
    EXPECT_NEAR(result.offset(1), 3.0, 0.01);
    EXPECT_GT(result.offset(1) - result.offset(0), 1.9);
    for (const auto multiplier : { 0.0, 0.25, std::numeric_limits<double>::quiet_NaN() })
    {
        const auto clamped{ offset_detail::EstimateJointOffsets(
            fixture.first, { 0, 1 }, fixture.second,
            { multiplier, multiplier }, solver, false) };
        ASSERT_EQ(clamped.status, result.status);
        ASSERT_EQ(clamped.offset.size(), result.offset.size());
        EXPECT_DOUBLE_EQ(clamped.offset(0), result.offset(0));
        EXPECT_DOUBLE_EQ(clamped.offset(1), result.offset(1));
    }
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetEstimatorKeepsFrozenBackgroundInRhs)
{
    auto fixture{
        BuildJointOffsetEstimationFixture(
            {
                rg::GaussianModel3D{ 6.0, 0.55, 1.0 },
                rg::GaussianModel3D{ 7.0, 0.60, 3.0 }
            },
            { 1.0, 3.0 })
    };
    alg::WeightedRidgeSolver solver;
    const auto result{
        offset_detail::EstimateJointOffsets(
            fixture.first,
            { 0, 1 },
            fixture.second,
            { 1.0, 1.0 },
            solver,
            false)
    };

    EXPECT_EQ(
        result.status,
        offset_detail::JointOffsetSolveStatus::Converged);
    ASSERT_EQ(result.offset.size(), 2);
    EXPECT_NEAR(result.offset(0), 1.0, 1.0e-5);
    EXPECT_NEAR(result.offset(1), 3.0, 1.0e-5);
    EXPECT_GT(std::abs(result.offset(0) - result.offset(1)), 1.0);


    offset_detail::FitState state;
    for (const auto & gaussian : fixture.second.node) state.emplace_back(MakeGaussianResult(gaussian));
    for (auto & atom : fixture.first.atom_list)
        atom.unselected_distance_list_by_sample.assign(atom.raw_sampling_entries.size(), { 0.4 });
    fixture.first.frozen_background = offset_detail::BuildFrozenBackground(fixture.first, state, { { 0, 1 } });
    ASSERT_TRUE(fixture.first.frozen_background);
    for (std::size_t node = 0; node < state.size(); node++)
        for (std::size_t row = 0; row < fixture.first.at(node).raw_sampling_entries.size(); row++)
            fixture.first.at(node).raw_sampling_entries.at(row).response +=
                fixture.first.frozen_background->response_by_atom.at(node).at(row);
    fixture.second = offset_detail::BuildSecondStageModelSnapshot(fixture.first, state);
    rg::algorithm::WeightedRidgeSolver background_solver;
    const auto with_background{ offset_detail::EstimateJointOffsets(fixture.first, { 0, 1 },
        fixture.second, { 1.0, 1.0 }, background_solver, false) };
    ASSERT_EQ(with_background.status, offset_detail::JointOffsetSolveStatus::Converged);
    ASSERT_EQ(with_background.offset.size(), 2);
    EXPECT_NEAR(with_background.offset(0), result.offset(0), 1.0e-12);
    EXPECT_NEAR(with_background.offset(1), result.offset(1), 1.0e-12);

}

TEST(EstimatorSecondStageDefenseTest, JointOffsetEstimatorReportsBuildAndEmptyFailures)
{
    auto empty_fixture{
        BuildJointOffsetEstimationFixture(
            { rg::GaussianModel3D{ 6.0, 0.55, 2.0 } },
            { 2.0 })
    };
    empty_fixture.first.at(0).raw_sampling_entries.clear();
    empty_fixture.first.at(0).neighbor_atom_sample_offset_list.clear();
    alg::WeightedRidgeSolver empty_solver;
    const auto empty_result{
        offset_detail::EstimateJointOffsets(
            empty_fixture.first,
            { 0 },
            empty_fixture.second,
            { 1.0 },
            empty_solver,
            false)
    };
    EXPECT_EQ(
        empty_result.status,
        offset_detail::JointOffsetSolveStatus::EmptySystem);
    ASSERT_EQ(empty_result.offset.size(), 1);
    EXPECT_DOUBLE_EQ(empty_result.offset(0), 2.0);

    auto invalid_fixture{
        BuildJointOffsetEstimationFixture(
            { rg::GaussianModel3D{ 6.0, 0.55, 2.0 } },
            { 2.0 })
    };
    invalid_fixture.first.at(0).raw_sampling_entries.at(0).response =
        std::numeric_limits<double>::infinity();
    alg::WeightedRidgeSolver invalid_solver;
    const auto invalid_result{
        offset_detail::EstimateJointOffsets(
            invalid_fixture.first,
            { 0 },
            invalid_fixture.second,
            { 1.0 },
            invalid_solver,
            false)
    };
    EXPECT_EQ(
        invalid_result.status,
        offset_detail::JointOffsetSolveStatus::SystemBuildFailed);
    ASSERT_EQ(invalid_result.offset.size(), 1);
    EXPECT_DOUBLE_EQ(invalid_result.offset(0), 2.0);

    auto negligible_basis_fixture{ BuildJointOffsetEstimationFixture(
        { { 6.0, 0.55, 2.0 } }, { 2.0 }) };
    negligible_basis_fixture.first.at(0).raw_sampling_entries.front().point.distance = 1.0e20;
    const auto negligible_basis_result{ offset_detail::EstimateJointOffsets(
        negligible_basis_fixture.first, { 0 }, negligible_basis_fixture.second,
        { 1.0 }, empty_solver, false) };
    EXPECT_EQ(negligible_basis_result.status, offset_detail::JointOffsetSolveStatus::EmptySystem);
    ASSERT_EQ(negligible_basis_result.offset.size(), 1);
    EXPECT_DOUBLE_EQ(negligible_basis_result.offset(0), 2.0);

    for (const bool active_neighbor : { false, true })
    {
        auto non_finite_neighbor_fixture{ BuildJointOffsetEstimationFixture(
            { { 6.0, 0.55, 2.0 }, { 7.0, 0.60, 3.0 } }, { 2.0, 3.0 }) };
        auto & target{ non_finite_neighbor_fixture.first.at(0) };
        target.neighbor_atom_sample_list = { { 1, std::numeric_limits<double>::quiet_NaN() } };
        target.neighbor_atom_sample_offset_list = { 0, 1 };
        const auto non_finite_neighbor_result{ offset_detail::EstimateJointOffsets(
            non_finite_neighbor_fixture.first,
            active_neighbor ? offset_detail::ClusterKey{ 0, 1 } : offset_detail::ClusterKey{ 0 },
            non_finite_neighbor_fixture.second, { 1.0, 1.0 }, invalid_solver, false) };
        EXPECT_EQ(non_finite_neighbor_result.status, offset_detail::JointOffsetSolveStatus::SystemBuildFailed);
        ASSERT_EQ(non_finite_neighbor_result.offset.size(), active_neighbor ? 2 : 1);
        EXPECT_DOUBLE_EQ(non_finite_neighbor_result.offset(0), 2.0);
    }
}

TEST(EstimatorSecondStageDefenseTest, TransformedChangeIsIntensityScaleInvariant)
{
    const rg::GaussianModel3D previous{ 8.0, 0.50, -0.10 };
    const rg::GaussianModel3D current{ 8.8, 0.55, -0.12 };
    const auto base_change{
        change_detail::CalculateTransformedChange(current, previous)
    };

    for (const auto scale : { 1.0e-2, 1.0e2 })
    {
        const auto scaled_change{
            change_detail::CalculateTransformedChange(
                rg::GaussianModel3D{
                    current.GetAmplitude() * scale,
                    current.GetWidth(),
                    current.GetOffset() * scale
                },
                rg::GaussianModel3D{
                    previous.GetAmplitude() * scale,
                    previous.GetWidth(),
                    previous.GetOffset() * scale
                })
        };
        ASSERT_EQ(base_change.size(), scaled_change.size());
        for (std::size_t i = 0; i < base_change.size(); i++)
        {
            EXPECT_NEAR(
                base_change.at(i),
                scaled_change.at(i),
                1.0e-12);
        }
    }
}

TEST(EstimatorSecondStageDefenseTest, AdaptiveTopologyRebuildUsesDriftAndIntervalTriggers)
{
    const audit_detail::FitState reference_state{
        MakeGaussianResult(rg::GaussianModel3D{ 8.0, 0.50, 0.10 })
    };
    auto small_drift_state{ reference_state };
    auto large_drift_state{ reference_state };
    const auto reference_coordinates{
        reference_state.at(0).mdpde.GetModel().ToTransformedCoordinates()
    };
    ASSERT_TRUE(reference_coordinates.has_value());
    auto small_coordinates{ *reference_coordinates };
    auto large_coordinates{ *reference_coordinates };
    small_coordinates(static_cast<Eigen::Index>(
        rg::GaussianModel3D::LogWidthCoordinateIndex())) += 0.099;
    large_coordinates(static_cast<Eigen::Index>(
        rg::GaussianModel3D::LogWidthCoordinateIndex())) += 0.101;
    const auto small_model{
        rg::GaussianModel3D::FromTransformedCoordinates(small_coordinates)
    };
    const auto large_model{
        rg::GaussianModel3D::FromTransformedCoordinates(large_coordinates)
    };
    ASSERT_TRUE(small_model.has_value());
    ASSERT_TRUE(large_model.has_value());
    small_drift_state.at(0) = MakeGaussianResult(*small_model);
    large_drift_state.at(0) = MakeGaussianResult(*large_model);

    const auto none{
        audit_detail::EvaluateAdaptiveTopologyRebuildTrigger(
            small_drift_state,
            reference_state,
            { 0 },
            2)
    };
    EXPECT_EQ(none.trigger, audit_detail::AdaptiveTopologyRebuildTrigger::None);
    EXPECT_NEAR(none.maximum_transformed_drift, 0.099, 1.0e-12);

    const auto interval{
        audit_detail::EvaluateAdaptiveTopologyRebuildTrigger(
            small_drift_state,
            reference_state,
            { 0 },
            3)
    };
    EXPECT_EQ(
        interval.trigger,
        audit_detail::AdaptiveTopologyRebuildTrigger::Interval);

    const auto drift{
        audit_detail::EvaluateAdaptiveTopologyRebuildTrigger(
            large_drift_state,
            reference_state,
            { 0 },
            1)
    };
    EXPECT_EQ(drift.trigger, audit_detail::AdaptiveTopologyRebuildTrigger::Drift);
    EXPECT_NEAR(drift.maximum_transformed_drift, 0.101, 1.0e-12);


    audit_detail::SecondStageContext context;
    context.atom_list.resize(2);
    audit_detail::FitState partition_state{
        MakeGaussianResult({ 4.0, 0.4, 1.0 }), MakeGaussianResult({ 8.0, 0.8, 3.0 }) };
    for (auto & atom : context.atom_list)
    {
        atom.raw_sampling_entries.resize(1);
        atom.unselected_distance_list_by_sample = { { 0.3 } };
    }
    const auto separate{ audit_detail::BuildFrozenBackground(context, partition_state, { { 0 }, { 1 } }) };
    const auto merged{ audit_detail::BuildFrozenBackground(context, partition_state, { { 0, 1 } }) };
    ASSERT_TRUE(separate);
    ASSERT_TRUE(merged);
    ExpectGaussianModelsNear(separate->model_by_atom.at(0), partition_state.at(0).mdpde.GetModel(), 0.0);
    ExpectGaussianModelsNear(separate->model_by_atom.at(1), partition_state.at(1).mdpde.GetModel(), 0.0);
    for (const auto & model : merged->model_by_atom)
        ExpectGaussianModelsNear(model, { 6.0, 0.6, 2.0 }, 1.0e-12);
    EXPECT_NE(separate->response_by_atom, merged->response_by_atom);
    context.frozen_background = separate;
    context.at(1).unselected_distance_list_by_sample.front().front() = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(audit_detail::BuildFrozenBackground(context, partition_state, { { 0, 1 } }));
    EXPECT_EQ(context.frozen_background, separate);
    context.at(1).unselected_distance_list_by_sample.front().front() = 0.3;
    EXPECT_FALSE(audit_detail::BuildFrozenBackground(context, partition_state, { { 0 } }));
    EXPECT_FALSE(audit_detail::BuildFrozenBackground(context, partition_state, { { 0 }, { 0, 1 } }));
    partition_state.front().mdpde = MakeGaussianResult({ -1.0, 0.4, 1.0 }).mdpde;
    EXPECT_FALSE(audit_detail::BuildFrozenBackground(context, partition_state, { { 0, 1 } }));
    ExpectGaussianModelsNear(separate->model_by_atom.at(0), { 4.0, 0.4, 1.0 }, 0.0);

}

TEST(EstimatorSecondStageDefenseTest, AdaptiveTopologyDriftTriggerIsIntensityScaleInvariant)
{
    constexpr double intensity_scale{ 100.0 };
    const audit_detail::FitState reference_state{
        MakeGaussianResult(rg::GaussianModel3D{ 8.0, 0.50, 0.10 })
    };
    const audit_detail::FitState accepted_state{
        MakeGaussianResult(rg::GaussianModel3D{ 9.0, 0.60, 0.15 })
    };
    const audit_detail::FitState scaled_reference_state{
        MakeGaussianResult(rg::GaussianModel3D{
            8.0 * intensity_scale,
            0.50,
            0.10 * intensity_scale })
    };
    const audit_detail::FitState scaled_accepted_state{
        MakeGaussianResult(rg::GaussianModel3D{
            9.0 * intensity_scale,
            0.60,
            0.15 * intensity_scale })
    };

    const auto base{
        audit_detail::EvaluateAdaptiveTopologyRebuildTrigger(
            accepted_state,
            reference_state,
            { 0 },
            1)
    };
    const auto scaled{
        audit_detail::EvaluateAdaptiveTopologyRebuildTrigger(
            scaled_accepted_state,
            scaled_reference_state,
            { 0 },
            1)
    };
    EXPECT_EQ(base.trigger, scaled.trigger);
    EXPECT_NEAR(
        base.maximum_transformed_drift,
        scaled.maximum_transformed_drift,
        1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, JointPolishJacobianMatchesFiniteDifference)
{
    constexpr double step{ 1.0e-6 };
    const std::array<rg::GaussianModel3D, 2> model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
        rg::GaussianModel3D{ 5.5, 0.70, -0.15 }
    };
    const std::array<double, 4> distance_list{
        0.0,
        5.0e-6,
        1.0e-5,
        0.35
    };
    for (const auto & model : model_list)
    {
        const auto transformed{ model.ToTransformedCoordinates() };
        ASSERT_TRUE(transformed.has_value());
        const auto invariants{
            polish_detail::BuildTransformedModelInvariants(model)
        };
        ASSERT_TRUE(invariants.has_value());
        for (const auto distance : distance_list)
        {
            const auto evaluation{
                polish_detail::EvaluateTransformedJacobian(
                    *invariants,
                    distance)
            };
            ASSERT_TRUE(evaluation.has_value());
            EXPECT_TRUE(evaluation->allFinite());

            for (Eigen::Index parameter_index = 0;
                parameter_index < transformed->size();
                parameter_index++)
            {
                auto lower{ *transformed };
                auto upper{ *transformed };
                lower(parameter_index) -= step;
                upper(parameter_index) += step;
                const auto lower_model{
                    rg::GaussianModel3D::FromTransformedCoordinates(lower)
                };
                const auto upper_model{
                    rg::GaussianModel3D::FromTransformedCoordinates(upper)
                };
                ASSERT_TRUE(lower_model.has_value());
                ASSERT_TRUE(upper_model.has_value());
                const auto finite_difference{
                    (upper_model->ResponseAtDistance(distance) -
                        lower_model->ResponseAtDistance(distance)) /
                    (2.0 * step)
                };
                const auto analytic{
                    (*evaluation)(parameter_index)
                };
                const auto tolerance{
                    1.0e-6 * std::max(std::abs(finite_difference), 1.0)
                };
                EXPECT_NEAR(analytic, finite_difference, tolerance);
            }
        }
    }
}

TEST(EstimatorSecondStageDefenseTest, GaussianParameterMedianUsesValidComponentMedians)
{
    const std::vector<rg::GaussianModel3D> models{
        { 1.0, 0.6, 1.0 }, { 9.0, 0.2, 3.0 }, { 5.0, 0.4, 2.0 }, {} };
    const auto odd{ median_detail::BuildGaussianParameterMedian(models) };
    ASSERT_TRUE(odd.has_value());
    ExpectGaussianModelsNear(*odd, { 5.0, 0.4, 2.0 }, 1.0e-12);
    const auto even{ median_detail::BuildGaussianParameterMedian(
        { { 4.0, 0.5, -2.0 }, { 8.0, 0.9, 2.0 } }) };
    ASSERT_TRUE(even.has_value());
    ExpectGaussianModelsNear(*even, { 6.0, 0.7, 0.0 }, 1.0e-12);
    EXPECT_FALSE(median_detail::BuildGaussianParameterMedian({}).has_value());
    EXPECT_FALSE(median_detail::BuildGaussianParameterMedian({ { } }).has_value());
}

TEST(EstimatorSecondStageDefenseTest,
    DampedModelsInterpolateIndividualPhysicalOffsets)
{
    const std::vector<rg::GaussianModel3D> previous_model_list{
        rg::GaussianModel3D{ 4.0, 0.40, 0.1 },
        rg::GaussianModel3D{ 6.0, 0.60, 0.9 },
        rg::GaussianModel3D{ 8.0, 0.80, -1.0 }
    };
    const std::vector<rg::GaussianModel3D> raw_model_list{
        rg::GaussianModel3D{ 5.0, 0.50, 0.5 },
        rg::GaussianModel3D{ 9.0, 0.90, 0.7 },
        rg::GaussianModel3D{ 7.0, 0.70, 2.0 }
    };

    for (const auto damping : std::array<double, 3>{ 0.0, 0.25, 1.0 })
    {
        const auto candidate_model_list{
            median_detail::BuildDampedModelList(
                previous_model_list,
                raw_model_list,
                damping)
        };
        ASSERT_TRUE(candidate_model_list.has_value());
        for (std::size_t atom_position = 0;
            atom_position < candidate_model_list->size();
            atom_position++)
        {
            EXPECT_NEAR(candidate_model_list->at(atom_position).GetOffset(),
                std::lerp(previous_model_list.at(atom_position).GetOffset(),
                    raw_model_list.at(atom_position).GetOffset(), damping), 1.0e-12);
            const auto previous_coordinates{
                previous_model_list.at(atom_position).ToTransformedCoordinates()
            };
            const auto raw_coordinates{
                raw_model_list.at(atom_position).ToTransformedCoordinates()
            };
            const auto candidate_coordinates{
                candidate_model_list->at(atom_position).ToTransformedCoordinates()
            };
            ASSERT_TRUE(previous_coordinates.has_value());
            ASSERT_TRUE(raw_coordinates.has_value());
            ASSERT_TRUE(candidate_coordinates.has_value());
            for (const auto parameter_index : std::array<int, 2>{
                rg::GaussianModel3D::LogPeakHeightCoordinateIndex(),
                rg::GaussianModel3D::LogWidthCoordinateIndex() })
            {
                const auto eigen_index{
                    static_cast<Eigen::Index>(parameter_index)
                };
                EXPECT_NEAR(
                    (*candidate_coordinates)(eigen_index),
                    (*previous_coordinates)(eigen_index) +
                        damping * ((*raw_coordinates)(eigen_index) -
                            (*previous_coordinates)(eigen_index)),
                    1.0e-12);
            }
        }
    }
}

TEST(EstimatorSecondStageDefenseTest, GaussianParameterMedianIsIntensityScaleInvariant)
{
    constexpr double scale{ 100.0 };
    const std::vector<rg::GaussianModel3D> models{
        { 3.0, 0.4, -0.2 }, { 5.0, 0.6, 0.1 }, { 7.0, 0.8, 0.4 } };
    std::vector<rg::GaussianModel3D> scaled_models;
    for (const auto & model : models)
    {
        scaled_models.emplace_back(scale * model.GetAmplitude(), model.GetWidth(),
            scale * model.GetOffset());
    }
    const auto median{ median_detail::BuildGaussianParameterMedian(models) };
    const auto scaled{ median_detail::BuildGaussianParameterMedian(scaled_models) };
    ASSERT_TRUE(median.has_value());
    ASSERT_TRUE(scaled.has_value());
    EXPECT_DOUBLE_EQ(scale * median->GetAmplitude(), scaled->GetAmplitude());
    EXPECT_DOUBLE_EQ(median->GetWidth(), scaled->GetWidth());
    EXPECT_DOUBLE_EQ(scale * median->GetOffset(), scaled->GetOffset());
}

TEST(EstimatorSecondStageDefenseTest,
    IndividualRefitUsesOwnModelForTargetOffsetResponse)
{
    const rg::GaussianModel3D offset_model{ 3.0, 0.70, 0.2 };
    const rg::GaussianModel3D truth_shape{ 6.0, 0.55, 0.0 };
    LocalPotentialSampleList sample_list;
    for (const auto distance : std::array<double, 5>{
        0.0, 0.15, 0.30, 0.45, 0.60 })
    {
        sample_list.emplace_back(LocalPotentialSample{
            truth_shape.SignalAtDistance(distance) +
                offset_model.GetOffset() *
                    offset_model.OffsetBasisAtDistance(distance),
            SamplingPoint{ distance }
        });
    }

    const auto result{
        rt::EstimateLocalGaussian(
            sample_list,
            0.0,
            MakeSecondStageOptions(),
            offset_model)
    };
    EXPECT_NEAR(
        result.mdpde.GetModel().GetAmplitude(),
        truth_shape.GetAmplitude(),
        1.0e-4);
    EXPECT_NEAR(
        result.mdpde.GetModel().GetWidth(),
        truth_shape.GetWidth(),
        1.0e-6);
    EXPECT_DOUBLE_EQ(
        result.mdpde.GetModel().GetOffset(),
        offset_model.GetOffset());
}

TEST(EstimatorSecondStageDefenseTest,
    JointPolishParameterizationKeepsClusterBackgroundFrozen)
{
    const std::vector<rg::GaussianModel3D> models{
        { 6.0, 0.55, 1.0 }, { 7.0, 0.60, 4.0 }, { 8.0, 0.65, 3.0 } };
    polish_detail::SecondStageContext context;
    context.atom_list.resize(models.size());
    polish_detail::FitState state;
    for (std::size_t i = 0; i < models.size(); i++)
    {
        state.emplace_back(MakeGaussianResult(models.at(i)));
        context.at(i).raw_sampling_entries.resize(1);
        context.at(i).unselected_distance_list_by_sample = { { 0.35 } };
    }
    context.frozen_background = polish_detail::BuildFrozenBackground(context, state, { { 0, 1, 2 } });
    ASSERT_TRUE(context.frozen_background);
    ExpectGaussianModelsNear(context.frozen_background->model_by_atom.front(), { 7.0, 0.60, 3.0 }, 1.0e-12);
    const auto parameterization{ polish_detail::JointPolishParameterization::Build( models) };
    ASSERT_TRUE(parameterization.has_value());
    EXPECT_EQ(parameterization->AtomCount(), 3U);
    EXPECT_EQ(parameterization->ParameterCount(), 9);
    EXPECT_NE(parameterization->OffsetColumn(0), parameterization->OffsetColumn(2));
    const auto frozen{ context.frozen_background };
    auto direction{ Eigen::VectorXd::Zero(9).eval() };
    direction(parameterization->OffsetColumn(0)) = 6.0;
    const auto candidate{ parameterization->DecodeModels(direction, 1.0) };
    ASSERT_TRUE(candidate.has_value());
    EXPECT_EQ(context.frozen_background, frozen);
    ExpectGaussianModelsNear(frozen->model_by_atom.front(), { 7.0, 0.60, 3.0 }, 1.0e-12);

    // A fixed selected member still contributes equally to the even median.
    const auto even{ polish_detail::BuildFrozenBackground(context, state, { { 0, 1 }, { 2 } }) };
    ASSERT_TRUE(even);
    ExpectGaussianModelsNear(even->model_by_atom.front(), { 6.5, 0.575, 2.5 }, 1.0e-12);
    const auto fixed{ polish_detail::JointPolishParameterization::BuildActiveSet( models, { 1, 0, 0 }, { 1, 0, 0 }) };
    ASSERT_TRUE(fixed.has_value());
    EXPECT_EQ(fixed->ParameterCount(), 3);
    constexpr double step{ 1.0e-6 };
    auto perturbation{ Eigen::VectorXd::Zero(3).eval() };
    perturbation(fixed->OffsetColumn(0)) = step;
    const auto upper{ fixed->DecodeModels(perturbation, 1.0) };
    const auto lower{ fixed->DecodeModels(-perturbation, 1.0) };
    ASSERT_TRUE(upper.has_value());
    ASSERT_TRUE(lower.has_value());
    const auto response{ polish_detail::EvaluatePhysicalOffsetResponse(models.front(), 0.35) };
    ASSERT_TRUE(response.has_value());
    const double background{ frozen->response_by_atom.front().front() };
    EXPECT_NEAR(((upper->front().ResponseAtDistance(0.35) + background) -
        (lower->front().ResponseAtDistance(0.35) + background)) / (2.0 * step),
        response->offset_jacobian, 1.0e-8);

}

TEST(EstimatorSecondStageDefenseTest, JointPolishSeedsAndDecodesIndividualOffsets)
{
    const std::vector<rg::GaussianModel3D> models{
        { 6.0, 0.55, 1.0 }, { 7.0, 0.60, 4.0 }, { 8.0, 0.65, 3.0 },
        { 9.0, 0.70, 2.0 }, { 10.0, 0.75, 6.0 } };
    const auto parameterization{ polish_detail::JointPolishParameterization::Build(models) };
    ASSERT_TRUE(parameterization.has_value());
    EXPECT_EQ(parameterization->ParameterCount(), 15);
    const auto seed{ parameterization->DecodeSeedModels() };
    ASSERT_TRUE(seed.has_value());
    auto direction{ Eigen::VectorXd::Zero(parameterization->ParameterCount()).eval() };
    direction(parameterization->OffsetColumn(0)) = -1.0;
    direction(parameterization->OffsetColumn(1)) = 2.0;
    const auto candidate{ parameterization->DecodeModels(direction, 0.5) };
    const auto zero{ parameterization->DecodeModels(direction, 0.0) };
    ASSERT_TRUE(candidate.has_value());
    ASSERT_TRUE(zero.has_value());
    for (std::size_t atom = 0; atom < models.size(); atom++)
    {
        ExpectGaussianModelsNear(seed->at(atom), models.at(atom), 1.0e-12);
        ExpectGaussianModelsNear(zero->at(atom), models.at(atom), 1.0e-12);
        const auto delta{ atom == 0 ? -0.5 : atom == 1 ? 1.0 : 0.0 };
        EXPECT_DOUBLE_EQ(candidate->at(atom).GetOffset(), models.at(atom).GetOffset() + delta);
        for (std::size_t other = atom + 1; other < models.size(); other++)
            EXPECT_NE(parameterization->OffsetColumn(atom), parameterization->OffsetColumn(other));
    }
}

TEST(
    EstimatorSecondStageDefenseTest,
    ActiveSetJointPolishKeepsInactiveCoordinatesAndIndependentOffsets)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.10 },
        rg::GaussianModel3D{ 7.0, 0.60, 0.20 }
    };
    const auto parameterization{
        polish_detail::JointPolishParameterization::BuildActiveSet(
            base_model_list,
            { 1, 0, 1 },
            { 1, 1, 1 })
    };
    ASSERT_TRUE(parameterization.has_value());
    EXPECT_TRUE(parameterization->HasShapeColumn(0));
    EXPECT_FALSE(parameterization->HasShapeColumn(1));
    EXPECT_TRUE(parameterization->HasShapeColumn(2));
    EXPECT_EQ(parameterization->ParameterCount(), 7);

    Eigen::VectorXd direction{
        Eigen::VectorXd::Zero(parameterization->ParameterCount())
    };
    direction(parameterization->ShapeColumn(0, 0)) = 0.2;
    direction(parameterization->ShapeColumn(2, 1)) = -0.1;
    direction(parameterization->OffsetColumn(0)) = 0.4;
    direction(parameterization->OffsetColumn(2)) = -0.2;
    const auto candidate_model_list{
        parameterization->DecodeModels(direction, 1.0)
    };
    ASSERT_TRUE(candidate_model_list.has_value());
    EXPECT_NE(
        candidate_model_list->at(0).GetAmplitude(),
        base_model_list.at(0).GetAmplitude());
    EXPECT_DOUBLE_EQ(
        candidate_model_list->at(1).GetAmplitude(),
        base_model_list.at(1).GetAmplitude());
    EXPECT_DOUBLE_EQ(
        candidate_model_list->at(1).GetWidth(),
        base_model_list.at(1).GetWidth());
    EXPECT_NE(
        candidate_model_list->at(2).GetWidth(),
        base_model_list.at(2).GetWidth());
    EXPECT_NE(
        candidate_model_list->at(0).GetOffset(),
        candidate_model_list->at(1).GetOffset());

    const auto fixed_offset_parameterization{
        polish_detail::JointPolishParameterization::BuildActiveSet(
            base_model_list,
            { 1, 0, 1 },
            { 1, 1, 0 })
    };
    ASSERT_TRUE(fixed_offset_parameterization.has_value());
    EXPECT_TRUE(fixed_offset_parameterization->HasOffsetColumn(0));
    EXPECT_FALSE(fixed_offset_parameterization->HasOffsetColumn(2));
    EXPECT_EQ(fixed_offset_parameterization->ParameterCount(), 6);
    const auto fixed_offset_models{
        fixed_offset_parameterization->DecodeModels(
            Eigen::VectorXd::Zero(fixed_offset_parameterization->ParameterCount()),
            1.0)
    };
    ASSERT_TRUE(fixed_offset_models.has_value());
    EXPECT_DOUBLE_EQ(
        fixed_offset_models->at(2).GetOffset(),
        base_model_list.at(2).GetOffset());

    auto fixture{ BuildJointPolishFixture(
        { { 6.0, 0.5, 0.1 }, { 7.0, 0.6, 0.3 } },
        { { 6.4, 0.55, 0.2 }, { 7.2, 0.62, 0.4 } }) };
    for (auto & atom : fixture.context.atom_list)
        atom.unselected_distance_list_by_sample.assign(atom.raw_sampling_entries.size(), { 0.35 });
    fixture.context.frozen_background = polish_detail::BuildFrozenBackground(fixture.context, fixture.state, { { 0 }, { 1 } });
    ASSERT_TRUE(fixture.context.frozen_background);
    const auto frozen{ fixture.context.frozen_background };
    for (std::size_t node = 0; node < fixture.state.size(); node++)
        for (std::size_t row = 0; row < fixture.context.at(node).raw_sampling_entries.size(); row++)
            fixture.context.at(node).raw_sampling_entries.at(row).response += frozen->response_by_atom.at(node).at(row);
    const polish_detail::FitStatePatch empty_patch;
    const polish_detail::FitStateView endpoint{ fixture.state, empty_patch };
    for (const bool freeze_second_shape : { false, true })
    {
        rg::algorithm::WeightedRidgeSolver solver;
        const auto correction{ polish_detail::BuildBoundaryJointCorrection(fixture.context, endpoint,
            freeze_second_shape ? polish_detail::ClusterKey{ 0 } : polish_detail::ClusterKey{ 0, 1 },
            { 0, 1 }, fixture.sample_ref_list, { 1.0, 1.0 }, { { { 0, 1 }, 100.0 } }, solver) };
        ASSERT_EQ(correction.status, polish_detail::BoundaryJointCorrectionStatus::CandidateReady);
        ASSERT_TRUE(correction.patch.has_value());
        EXPECT_EQ(correction.parameter_count, freeze_second_shape ? 4U : 6U);
        EXPECT_EQ(correction.patch->atom_index_list, (polish_detail::ClusterKey{ 0, 1 }));
        EXPECT_EQ(fixture.context.frozen_background, frozen);
        ExpectGaussianModelsNear(frozen->model_by_atom.at(0), { 6.0, 0.5, 0.1 }, 1.0e-12);
        ExpectGaussianModelsNear(frozen->model_by_atom.at(1), { 7.0, 0.6, 0.3 }, 1.0e-12);
    }

}

TEST(
    EstimatorSecondStageDefenseTest,
    BoundaryJointCorrectionUsesIndependentActiveCoordinatesAndPerClusterTrust)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.10 },
        rg::GaussianModel3D{ 7.0, 0.60, 0.20 }
    };
    const std::vector<rg::GaussianModel3D> target_model_list{
        rg::GaussianModel3D{ 6.8, 0.60, 0.30 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.20 },
        rg::GaussianModel3D{ 6.4, 0.56, -0.15 }
    };
    auto fixture{
        BuildJointPolishFixture(
            base_model_list,
            target_model_list)
    };
    const polish_detail::FitStatePatch endpoint_patch;
    const polish_detail::FitStateView endpoint_state{
        fixture.state,
        endpoint_patch
    };
    alg::WeightedRidgeSolver solver;
    const auto result{
        polish_detail::BuildBoundaryJointCorrection(
            fixture.context,
            endpoint_state,
            { 0, 2 },
            { 0, 1, 2 },
            fixture.sample_ref_list,
            { 1.0, 1.0, 1.0 },
            {
                { { 0, 1 }, 4.0 },
                { { 2 }, 0.5 }
            },
            solver)
    };
    ASSERT_EQ(
        result.status,
        polish_detail::BoundaryJointCorrectionStatus::CandidateReady);
    ASSERT_TRUE(result.patch.has_value());
    EXPECT_EQ(result.parameter_count, 7U);
    EXPECT_LE(result.maximum_normalized_trust_step, 1.0 + 1.0e-12);
    ASSERT_EQ(result.patch->mdpde_list.size(), 3U);
    const auto & offset_only_candidate{
        result.patch->mdpde_list.at(1).GetModel()
    };
    EXPECT_DOUBLE_EQ(
        offset_only_candidate.GetAmplitude(),
        base_model_list.at(1).GetAmplitude());
    EXPECT_DOUBLE_EQ(
        offset_only_candidate.GetWidth(),
        base_model_list.at(1).GetWidth());
    EXPECT_NE(
        result.patch->mdpde_list.at(0).GetModel().GetOffset(),
        offset_only_candidate.GetOffset());
    EXPECT_DOUBLE_EQ(
        result.patch->mdpde_list.at(1)
            .GetStandardDeviationModel().GetAmplitude(),
        fixture.state.at(1).mdpde
            .GetStandardDeviationModel().GetAmplitude());

    // Atom 1 has only an offset column: the ridge floor must not affect shape columns.
    for (const auto multiplier : { 0.25, std::numeric_limits<double>::quiet_NaN() })
    {
        const auto clamped{ polish_detail::BuildBoundaryJointCorrection(
            fixture.context, endpoint_state, { 0, 2 }, { 0, 1, 2 },
            fixture.sample_ref_list, { 1.0, multiplier, 1.0 },
            { { { 0, 1 }, 4.0 }, { { 2 }, 0.5 } }, solver) };
        ASSERT_EQ(clamped.status, result.status);
        ASSERT_TRUE(clamped.patch.has_value());
        ASSERT_EQ(clamped.patch->mdpde_list.size(), result.patch->mdpde_list.size());
        EXPECT_EQ(clamped.parameter_count, result.parameter_count);
        for (std::size_t atom_index = 0; atom_index < result.patch->mdpde_list.size(); atom_index++)
        {
            ExpectGaussianModelsNear(clamped.patch->mdpde_list.at(atom_index).GetModel(),
                result.patch->mdpde_list.at(atom_index).GetModel(), 0.0);
        }
    }

    alg::WeightedRidgeSolver invalid_solver;
    EXPECT_EQ(
        polish_detail::BuildBoundaryJointCorrection(
            fixture.context,
            endpoint_state,
            {},
            {},
            fixture.sample_ref_list,
            { 1.0, 1.0, 1.0 },
            { { { 0, 1, 2 }, 4.0 } },
            invalid_solver).status,
        polish_detail::BoundaryJointCorrectionStatus::InvalidInput);

    alg::WeightedRidgeSolver small_step_solver;
    EXPECT_EQ(
        polish_detail::BuildBoundaryJointCorrection(
            fixture.context,
            endpoint_state,
            { 0, 2 },
            { 0, 1, 2 },
            fixture.sample_ref_list,
            { 1.0, 1.0, 1.0 },
            {
                { { 0, 1 }, 1.0e-8 },
                { { 2 }, 1.0e-8 }
            },
            small_step_solver).status,
        polish_detail::BoundaryJointCorrectionStatus::NoMaterialChange);

    auto outside_patch{ polish_detail::FitStatePatch::FromState(fixture.state, { 0 }) };
    outside_patch.mdpde_list.front() = MakeGaussianResult({ 10.0, 0.55, 0.10 }).mdpde;
    const polish_detail::FitStateView outside_state{ fixture.state, outside_patch };
    alg::WeightedRidgeSolver unavailable_solver;
    EXPECT_EQ(
        polish_detail::BuildBoundaryJointCorrection(
            fixture.context,
            outside_state,
            { 0, 2 },
            { 0, 1, 2 },
            fixture.sample_ref_list,
            { 1.0, 1.0, 1.0 },
            {
                { { 0, 1 }, 1.0e-8 },
                { { 2 }, 1.0e-8 }
            },
            unavailable_solver).status,
        polish_detail::BoundaryJointCorrectionStatus::TrustRegionUnavailable);

    auto non_finite_fixture{ fixture };
    non_finite_fixture.context.at(0).raw_sampling_entries.at(0).response =
        std::numeric_limits<double>::infinity();
    alg::WeightedRidgeSolver non_finite_solver;
    EXPECT_EQ(
        polish_detail::BuildBoundaryJointCorrection(
            non_finite_fixture.context,
            endpoint_state,
            { 0, 2 },
            { 0, 1, 2 },
            non_finite_fixture.sample_ref_list,
            { 1.0, 1.0, 1.0 },
            {
                { { 0, 1 }, 4.0 },
                { { 2 }, 0.5 }
            },
            non_finite_solver).status,
        polish_detail::BoundaryJointCorrectionStatus::SystemBuildFailed);

    const std::vector<rg::GaussianModel3D> stationary_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.0 },
        rg::GaussianModel3D{ 4.5, 0.70, 0.0 },
        rg::GaussianModel3D{ 7.0, 0.60, 0.0 }
    };
    auto stationary_fixture{
        BuildJointPolishFixture(
            stationary_model_list,
            stationary_model_list)
    };
    const polish_detail::FitStateView stationary_endpoint_state{
        stationary_fixture.state,
        endpoint_patch
    };
    alg::WeightedRidgeSolver stationary_solver;
    EXPECT_EQ(
        polish_detail::BuildBoundaryJointCorrection(
            stationary_fixture.context,
            stationary_endpoint_state,
            { 0, 2 },
            { 0, 1, 2 },
            stationary_fixture.sample_ref_list,
            { 1.0, 1.0, 1.0 },
            {
                { { 0, 1 }, 4.0 },
                { { 2 }, 0.5 }
            },
            stationary_solver).status,
        polish_detail::BoundaryJointCorrectionStatus::NoMaterialChange);
}

TEST(
    EstimatorSecondStageDefenseTest,
    JointPolishProposalKeepsIndependentOffsets)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.10 }
    };
    const std::vector<rg::GaussianModel3D> target_model_list{
        rg::GaussianModel3D{ 6.5, 0.60, 0.30 },
        rg::GaussianModel3D{ 4.0, 0.65, -0.20 }
    };
    auto fixture{
        BuildJointPolishFixture(
            base_model_list,
            target_model_list)
    };
    const polish_detail::FitStatePatch base_patch;
    const polish_detail::FitStateView base_state_view{ fixture.state, base_patch };
    alg::WeightedRidgeSolver proposal_solver;
    const auto proposal{
        polish_detail::BuildJointPolishProposal(
            fixture.context,
            base_state_view,
            polish_detail::ClusterKey{ 0, 1 },
            fixture.sample_ref_list,
            { 1.0, 1.0 },
            proposal_solver,
            4.0)
    };
    ASSERT_TRUE(proposal.has_value());
    EXPECT_EQ(proposal->patch.atom_index_list,
        (polish_detail::ClusterKey{ 0, 1 }));
    ASSERT_EQ(proposal->patch.mdpde_list.size(), 2U);
    EXPECT_NE(
        proposal->patch.mdpde_list.at(0).GetModel().GetOffset(),
        proposal->patch.mdpde_list.at(1).GetModel().GetOffset());
    bool shape_changed{ false };
    bool offset_changed{ false };
    for (std::size_t atom_index = 0; atom_index < base_model_list.size(); atom_index++)
    {
        const auto & candidate{
            proposal->patch.mdpde_list.at(atom_index).GetModel()
        };
        const auto & base{ base_model_list.at(atom_index) };
        shape_changed = shape_changed ||
            std::abs(candidate.GetAmplitude() - base.GetAmplitude()) > 1.0e-8 ||
            std::abs(candidate.GetWidth() - base.GetWidth()) > 1.0e-8;
        offset_changed = offset_changed ||
            std::abs(candidate.GetOffset() - base.GetOffset()) > 1.0e-8;
    }
    EXPECT_TRUE(shape_changed);
    EXPECT_TRUE(offset_changed);
    EXPECT_LE(proposal->step_norm, 4.0 + 1.0e-12);
    EXPECT_DOUBLE_EQ(
        proposal->patch.mdpde_list.at(0)
            .GetStandardDeviationModel().GetAmplitude(),
        0.0);
}

TEST(
    EstimatorSecondStageDefenseTest,
    JointPolishProposalRejectsUnchangedAndRespectsSmallTrustRegion)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.10 }
    };
    const std::vector<rg::GaussianModel3D> target_model_list{
        rg::GaussianModel3D{ 6.5, 0.60, 0.30 },
        rg::GaussianModel3D{ 4.0, 0.65, -0.20 }
    };
    auto fixture{
        BuildJointPolishFixture(
            base_model_list,
            target_model_list)
    };
    const polish_detail::FitStatePatch base_patch;
    const polish_detail::FitStateView base_state_view{ fixture.state, base_patch };
    const auto key{ polish_detail::ClusterKey{ 0, 1 } };

    const std::vector<rg::GaussianModel3D> unchanged_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, 0.10 }
    };
    auto unchanged_fixture{
        BuildJointPolishFixture(
            unchanged_model_list,
            unchanged_model_list)
    };
    const polish_detail::FitStatePatch unchanged_patch;
    const polish_detail::FitStateView unchanged_state_view{
        unchanged_fixture.state,
        unchanged_patch
    };
    alg::WeightedRidgeSolver unchanged_solver;
    EXPECT_FALSE(
        polish_detail::BuildJointPolishProposal(
            unchanged_fixture.context,
            unchanged_state_view,
            key,
            unchanged_fixture.sample_ref_list,
            { 1.0, 1.0 },
            unchanged_solver,
            4.0).has_value());

    alg::WeightedRidgeSolver trust_region_solver;
    const auto small_proposal{
        polish_detail::BuildJointPolishProposal(
            fixture.context,
            base_state_view,
            key,
            fixture.sample_ref_list,
            { 1.0, 1.0 },
            trust_region_solver,
            0.01)
    };
    ASSERT_TRUE(small_proposal.has_value());
    EXPECT_LE(small_proposal->step_norm, 0.01 + 1.0e-12);
}

TEST(
    EstimatorSecondStageDefenseTest,
    JointPolishProposalUsesFitStateViewBaseForTrustRegionOrigin)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.10 }
    };
    const std::vector<rg::GaussianModel3D> target_model_list{
        rg::GaussianModel3D{ 6.5, 0.60, 0.30 },
        rg::GaussianModel3D{ 4.0, 0.65, -0.20 }
    };
    auto fixture{
        BuildJointPolishFixture(
            base_model_list,
            target_model_list)
    };
    auto patch{
        polish_detail::FitStatePatch::FromState(
            fixture.state,
            polish_detail::ClusterKey{ 0, 1 })
    };
    patch.mdpde_list.at(0) = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 10.0, 0.55, 0.10 },
        rg::GaussianModel3DUncertainty{}
    };
    const polish_detail::FitStateView base_state_view{ fixture.state, patch };
    EXPECT_DOUBLE_EQ(base_state_view.GetBaseModel(0).GetAmplitude(), 6.0);
    EXPECT_DOUBLE_EQ(base_state_view.GetModel(0).GetAmplitude(), 10.0);
    const auto patched_parameterization{
        polish_detail::JointPolishParameterization::Build(
            std::vector<rg::GaussianModel3D>{
                base_state_view.GetModel(0),
                base_state_view.GetModel(1) })
    };
    ASSERT_TRUE(patched_parameterization.has_value());
    const auto patched_seed{ patched_parameterization->DecodeSeedModels() };
    ASSERT_TRUE(patched_seed.has_value());
    EXPECT_DOUBLE_EQ(patched_seed->at(0).GetAmplitude(), 10.0);

    alg::WeightedRidgeSolver solver;
    EXPECT_FALSE(
        polish_detail::BuildJointPolishProposal(
            fixture.context,
            base_state_view,
            polish_detail::ClusterKey{ 0, 1 },
            fixture.sample_ref_list,
            { 1.0, 1.0 },
            solver,
            0.2).has_value());
}

TEST(EstimatorSecondStageDefenseTest, PhysicalOffsetJacobianMatchesFiniteDifference)
{
    constexpr double step{ 1.0e-6 };
    const std::array<rg::GaussianModel3D, 2> model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
        rg::GaussianModel3D{ 5.5, 0.70, -0.15 }
    };
    const std::array<double, 4> distance_list{
        0.0,
        5.0e-6,
        1.0e-5,
        0.35
    };
    for (const auto & model : model_list)
    {
        const auto transformed{ model.ToTransformedCoordinates() };
        ASSERT_TRUE(transformed.has_value());
        for (const auto distance : distance_list)
        {
            const auto evaluation{
                polish_detail::EvaluatePhysicalOffsetResponse(
                    model,
                    distance)
            };
            ASSERT_TRUE(evaluation.has_value());
            EXPECT_TRUE(std::isfinite(evaluation->response));
            EXPECT_TRUE(evaluation->shape_jacobian.allFinite());
            EXPECT_TRUE(std::isfinite(evaluation->offset_jacobian));

            for (Eigen::Index parameter_index = 0;
                parameter_index < evaluation->shape_jacobian.size();
                parameter_index++)
            {
                auto lower{ *transformed };
                auto upper{ *transformed };
                lower(parameter_index) -= step;
                upper(parameter_index) += step;
                lower(static_cast<Eigen::Index>(
                    rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex())) = 0.0;
                upper(static_cast<Eigen::Index>(
                    rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex())) = 0.0;
                const auto lower_shape{
                    rg::GaussianModel3D::FromTransformedCoordinates(lower)
                };
                const auto upper_shape{
                    rg::GaussianModel3D::FromTransformedCoordinates(upper)
                };
                ASSERT_TRUE(lower_shape.has_value());
                ASSERT_TRUE(upper_shape.has_value());
                const auto finite_difference{
                    (upper_shape->WithOffset(model.GetOffset())
                            .ResponseAtDistance(distance) -
                        lower_shape->WithOffset(model.GetOffset())
                            .ResponseAtDistance(distance)) /
                    (2.0 * step)
                };
                const auto tolerance{
                    1.0e-6 * std::max(std::abs(finite_difference), 1.0)
                };
                EXPECT_NEAR(
                    evaluation->shape_jacobian(parameter_index),
                    finite_difference,
                    tolerance);
            }

            const auto offset_finite_difference{
                (model.WithOffset(model.GetOffset() + step)
                        .ResponseAtDistance(distance) -
                    model.WithOffset(model.GetOffset() - step)
                        .ResponseAtDistance(distance)) /
                (2.0 * step)
            };
            EXPECT_NEAR(
                evaluation->offset_jacobian,
                offset_finite_difference,
                1.0e-8 * std::max(std::abs(offset_finite_difference), 1.0));
        }
    }
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphNormalizesFullJacobianEnergy)
{
    Eigen::Vector3d jacobian;
    jacobian << 1.0, 2.0, 3.0;
    coupling_detail::CouplingGraphBuilder builder{ 2 };
    AddCouplingGraphSample(builder, { 0, 0 }, { { 0, jacobian }, { 1, jacobian } });
    AddCouplingGraphSample(builder, { 0, 1 }, { { 0, 2.0 * jacobian }, { 1, 2.0 * jacobian } });
    const auto topology{ builder.BuildTopology() };
    EXPECT_TRUE(HasCouplingNeighbor(topology, 0, 1));
    EXPECT_NEAR(topology.summary.weight_median, 1.0, 1.0e-12);
    EXPECT_NEAR(topology.summary.weight_percentile_95, 1.0, 1.0e-12);
    EXPECT_NEAR(topology.summary.weight_maximum, 1.0, 1.0e-12);

    coupling_detail::CouplingGraphBuilder scaled_builder{ 2 };
    AddCouplingGraphSample(
        scaled_builder,
        { 0, 0 },
        { { 0, 5.0 * jacobian }, { 1, jacobian } });
    AddCouplingGraphSample(
        scaled_builder,
        { 0, 1 },
        { { 0, 10.0 * jacobian }, { 1, 2.0 * jacobian } });
    const auto scaled_topology{
        scaled_builder.BuildTopology()
    };
    EXPECT_TRUE(HasCouplingNeighbor(scaled_topology, 0, 1));
    EXPECT_NEAR(
        scaled_topology.summary.weight_median,
        topology.summary.weight_median,
        1.0e-12);

    coupling_detail::CouplingGraphBuilder tiny_builder{ 2 };
    AddCouplingGraphSample(
        tiny_builder,
        { 0, 0 },
        { { 0, 1.0e-100 * jacobian }, { 1, 1.0e-100 * jacobian } });
    const auto tiny_topology{
        tiny_builder.BuildTopology()
    };
    EXPECT_TRUE(HasCouplingNeighbor(tiny_topology, 0, 1));
    EXPECT_NEAR(
        tiny_topology.summary.weight_median,
        topology.summary.weight_median,
        1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphNormalizesDuplicateParticipants)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    coupling_detail::CouplingGraphBuilder builder{ 2 };
    AddCouplingGraphSample(
        builder,
        { 0, 0 },
        {
            { 1, unit },
            { 0, unit },
            { 1, unit }
        });
    AddCouplingGraphSample(builder, { 0, 1 }, { { 0, unit } });

    const auto topology{ builder.BuildTopology() };
    ASSERT_EQ(topology.sample_dependency_list.size(), 2U);
    EXPECT_EQ(
        topology.sample_dependency_list.front().contributor_atom_index_list,
        (std::vector<std::size_t>{ 0, 1 }));
    EXPECT_TRUE(HasCouplingNeighbor(topology, 0, 1));
    ASSERT_EQ(topology.retained_edge_list.size(), 1U);
    EXPECT_NEAR(
        topology.retained_edge_list.front().weight,
        1.0 / std::sqrt(2.0),
        1.0e-12);
    EXPECT_EQ(topology.summary.component_count, 1U);
    EXPECT_EQ(topology.summary.maximum_component_size, 2U);
    EXPECT_DOUBLE_EQ(topology.summary.maximum_component_ratio, 1.0);
    EXPECT_DOUBLE_EQ(topology.summary.configured_minimum_weight, 0.05);
    EXPECT_EQ(topology.atom_cutoff_summary.maximum_atom_count_limit, 100U);
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphPropagatesInvalidDuplicateJacobian)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    const auto invalid{
        Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN())
    };
    coupling_detail::CouplingGraphBuilder builder{ 2 };
    AddCouplingGraphSample(
        builder,
        { 0, 0 },
        {
            { 0, unit },
            { 1, unit },
            { 1, invalid }
        });

    const auto topology{ builder.BuildTopology() };
    EXPECT_FALSE(topology.summary.uses_weighted_graph);
    EXPECT_TRUE(topology.summary.threshold_sensitivity_list.empty());
    ASSERT_EQ(topology.sample_dependency_list.size(), 1U);
    EXPECT_EQ(
        topology.sample_dependency_list.front().contributor_atom_index_list,
        (std::vector<std::size_t>{ 0, 1 }));
    EXPECT_TRUE(HasCouplingNeighbor(topology, 0, 1));
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphSummaryUsesOnlySelectedSampleConnectivity)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    coupling_detail::CouplingGraphBuilder builder{ 2 };
    AddCouplingGraphSample(builder, { 0, 0 }, { { 0, unit } });
    AddCouplingGraphSample(builder, { 1, 0 }, { { 1, unit } });

    const auto topology{
        builder.BuildTopology()
    };
    EXPECT_TRUE(topology.adjacency_list.at(0).empty());
    EXPECT_TRUE(topology.adjacency_list.at(1).empty());
    EXPECT_EQ(topology.summary.component_count, 2U);
    EXPECT_EQ(topology.summary.maximum_component_size, 1U);
    EXPECT_DOUBLE_EQ(topology.summary.maximum_component_ratio, 0.5);

    constexpr std::size_t selected_count{ 11 };
    std::vector<std::unique_ptr<rg::AtomObject>> atoms;
    coupling_detail::SecondStageContext context;
    coupling_detail::FitState state;
    context.atom_list.resize(selected_count);
    const rg::GaussianModel3D model{ 6.0, 0.55, 0.10 };
    for (std::size_t i = 0; i < selected_count; i++)
    {
        atoms.emplace_back(MakeAtom(static_cast<int>(i + 1), Spot::C,
            Element::CARBON, { static_cast<double>(i), 0.0, 0.0 }));
        atoms.back()->SetSequenceID(1);
        context.at(i).atom = atoms.back().get();
        state.emplace_back(MakeGaussianResult(model));
        context.at(i).raw_sampling_entries = {
            { model.ResponseAtDistance(0.2) + model.ResponseAtDistance(0.3), SamplingPoint{ 0.2 } } };
        context.at(i).neighbor_atom_sample_offset_list = { 0, 0 };
        context.at(i).unselected_distance_list_by_sample = { { 0.3 } };
    }
    const auto background_topology{ coupling_detail::BuildSecondStageGraphTopology(context, state, true) };
    EXPECT_EQ(background_topology.adjacency_list.size(), selected_count);
    for (const auto & neighbors : background_topology.adjacency_list) EXPECT_TRUE(neighbors.empty());
    EXPECT_EQ(background_topology.summary.component_count, selected_count);
    EXPECT_EQ(background_topology.atom_cutoff_summary.maximum_atom_count, 1U);
    const auto partition{ coupling_detail::BuildGraphPartition(
        background_topology, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }) };
    EXPECT_EQ(partition.sample_id_list_by_key.size(), selected_count);
    for (std::size_t i = 0; i < atoms.size(); i++)
    {
        atoms.at(i)->SetChainID(i % 2 == 0 ? "B" : "C");
        atoms.at(i)->SetSequenceID(static_cast<int>(100 - i));
    }
    const auto relabeled{ coupling_detail::BuildSecondStageGraphTopology(context, state, true) };
    EXPECT_EQ(relabeled.adjacency_list, background_topology.adjacency_list);
    EXPECT_EQ(relabeled.summary.component_count, selected_count);
    EXPECT_EQ(coupling_detail::BuildGraphPartition(
        relabeled, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }).sample_id_list_by_key,
        partition.sample_id_list_by_key);
    EXPECT_THROW(coupling_detail::BuildSecondStageGraphTopology(context, {}, true), std::invalid_argument);

}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphCutsWeakAndCancelledEdges)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    coupling_detail::CouplingGraphBuilder weak_builder{ 2 };
    AddCouplingGraphSample(weak_builder, { 0, 0 }, { { 0, unit }, { 1, unit } });
    AddCouplingGraphSample(weak_builder, { 0, 1 }, { { 0, 10.0 * unit } });
    AddCouplingGraphSample(weak_builder, { 1, 0 }, { { 1, 10.0 * unit } });
    const auto weak_topology{
        weak_builder.BuildTopology()
    };
    EXPECT_FALSE(HasCouplingNeighbor(weak_topology, 0, 1));
    EXPECT_EQ(weak_topology.summary.candidate_edge_count, 1U);
    EXPECT_EQ(weak_topology.summary.cut_edge_count, 1U);

    coupling_detail::CouplingGraphBuilder cancelled_builder{ 2 };
    AddCouplingGraphSample(cancelled_builder, { 0, 0 }, { { 0, unit }, { 1, unit } });
    AddCouplingGraphSample(cancelled_builder, { 0, 1 }, { { 0, unit }, { 1, -unit } });
    const auto cancelled_topology{
        cancelled_builder.BuildTopology()
    };
    EXPECT_FALSE(HasCouplingNeighbor(cancelled_topology, 0, 1));
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphAdaptiveHysteresisAddsAndRemovesEdges)
{
    const auto build_topology = [](
        double edge_weight,
        const coupling_detail::GraphTopology * previous_topology)
    {
        const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
        const auto self_scale{ std::sqrt(1.0 / edge_weight - 1.0) };
        coupling_detail::CouplingGraphBuilder builder{ 2 };
        AddCouplingGraphSample(
            builder,
            { 0, 0 },
            { { 0, unit }, { 1, unit } });
        AddCouplingGraphSample(
            builder,
            { 0, 1 },
            { { 0, self_scale * unit } });
        AddCouplingGraphSample(
            builder,
            { 1, 0 },
            { { 1, self_scale * unit } });
        coupling_detail::CouplingGraphOptions options;
        options.minimum_weight = 0.06;
        options.retained_edge_minimum_weight = 0.04;
        return builder.BuildTopology(
            options,
            previous_topology);
    };

    coupling_detail::GraphTopology absent_previous;
    absent_previous.adjacency_list.resize(2);
    const auto absent_midpoint{ build_topology(0.05, &absent_previous) };
    EXPECT_FALSE(HasCouplingNeighbor(absent_midpoint, 0, 1));
    const auto added{ build_topology(0.061, &absent_previous) };
    EXPECT_TRUE(HasCouplingNeighbor(added, 0, 1));
    const auto cutoff_previous{ coupling_detail::ApplyGraphAtomCutoff(added, 1) };
    EXPECT_FALSE(HasCouplingNeighbor(build_topology(0.05, &cutoff_previous), 0, 1));
    const auto retained_midpoint{ build_topology(0.05, &added) };
    EXPECT_TRUE(HasCouplingNeighbor(retained_midpoint, 0, 1));
    const auto removed{ build_topology(0.039, &retained_midpoint) };
    EXPECT_FALSE(HasCouplingNeighbor(removed, 0, 1));
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphReportsThresholdSensitivity)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    coupling_detail::CouplingGraphBuilder builder{ 7 };
    const std::array<double, 3> edge_weight_list{ 0.06, 0.12, 0.25 };
    for (std::size_t edge_index = 0; edge_index < edge_weight_list.size(); edge_index++)
    {
        const auto left_index{ 2 * edge_index };
        const auto right_index{ left_index + 1 };
        const auto self_scale{
            std::sqrt(1.0 / edge_weight_list.at(edge_index) - 1.0)
        };
        AddCouplingGraphSample(
            builder,
            { edge_index, 0 },
            { { left_index, unit }, { right_index, unit } });
        AddCouplingGraphSample(
            builder,
            { edge_index, 1 },
            { { left_index, self_scale * unit } });
        AddCouplingGraphSample(
            builder,
            { edge_index, 2 },
            { { right_index, self_scale * unit } });
    }

    const std::vector<double> threshold_list{ 0.05, 0.075, 0.10, 0.15, 0.20, 0.30 };
    coupling_detail::CouplingGraphOptions options;
    options.sensitivity_minimum_weight_list = threshold_list;
    options.maximum_atom_count = 1;
    const auto topology{
        builder.BuildTopology(options)
    };
    ASSERT_EQ(topology.retained_edge_list.size(), edge_weight_list.size());
    for (std::size_t edge_index = 0;
        edge_index < topology.retained_edge_list.size();
        edge_index++)
    {
        const auto & edge{ topology.retained_edge_list.at(edge_index) };
        EXPECT_EQ(edge.left_atom_index, 2 * edge_index);
        EXPECT_EQ(edge.right_atom_index, 2 * edge_index + 1);
    }
    ASSERT_EQ(topology.summary.threshold_sensitivity_list.size(), threshold_list.size());

    const std::array<std::size_t, 6> retained_edge_count_list{ 3, 2, 2, 1, 1, 0 };
    for (std::size_t i = 0; i < threshold_list.size(); i++)
    {
        const auto & sensitivity{ topology.summary.threshold_sensitivity_list.at(i) };
        const auto retained_edge_count{ retained_edge_count_list.at(i) };
        EXPECT_DOUBLE_EQ(sensitivity.minimum_weight, threshold_list.at(i));
        EXPECT_EQ(sensitivity.retained_edge_count, retained_edge_count);
        EXPECT_EQ(sensitivity.cut_edge_count, 3U - retained_edge_count);
        EXPECT_EQ(sensitivity.component_count, 7U - retained_edge_count);
        EXPECT_EQ(sensitivity.maximum_component_size, retained_edge_count == 0 ? 1U : 2U);
        EXPECT_NEAR(
            sensitivity.maximum_component_ratio,
            static_cast<double>(sensitivity.maximum_component_size) / 7.0,
            1.0e-12);
    }

    const auto & formal_threshold{ topology.summary.threshold_sensitivity_list.front() };
    EXPECT_EQ(formal_threshold.retained_edge_count, topology.summary.retained_edge_count);
    EXPECT_EQ(
        topology.summary.candidate_edge_count - formal_threshold.retained_edge_count,
        topology.summary.cut_edge_count);
    const auto formal_partition{
        coupling_detail::BuildGraphPartition(
            topology,
            { 0, 1, 2, 3, 4, 5, 6 })
    };
    EXPECT_EQ(formal_threshold.component_count, 4U);
    EXPECT_EQ(formal_partition.sample_id_list_by_key.size(), 7U);
    EXPECT_EQ(formal_threshold.maximum_component_size, 2U);
    EXPECT_EQ(topology.summary.component_count, 7U);
    EXPECT_EQ(topology.summary.maximum_component_size, 1U);
    EXPECT_EQ(topology.atom_cutoff_summary.cut_edge_count, 3U);
    EXPECT_NEAR(topology.summary.maximum_component_ratio, 1.0 / 7.0, 1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, CouplingPartitionCutsWeakBridgeAndDuplicatesBoundarySample)
{
    coupling_detail::GraphTopology topology;
    topology.adjacency_list.resize(3);
    topology.adjacency_list.at(0).push_back(1);
    topology.adjacency_list.at(1).push_back(0);
    topology.sample_dependency_list = {
        { { 0, 0 }, { 0, 1 } },
        { { 1, 0 }, { 1, 2 } },
        { { 0, 1 }, { 1, 2 } }
    };

    const auto partition{
        coupling_detail::BuildGraphPartition(topology, { 0, 1, 2 })
    };
    ASSERT_EQ(partition.sample_id_list_by_key.size(), 2U);
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 0, 1 }), 1U);
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 2 }), 1U);
    ASSERT_EQ(partition.boundary_sample_dependency_list.size(), 2U);
    EXPECT_EQ(
        partition.boundary_sample_dependency_list.front().sample_id,
        (coupling_detail::SampleRef{ 0, 1 }));
    EXPECT_EQ(
        partition.boundary_sample_dependency_list.front().cluster_key_list,
        (std::vector<audit_detail::ClusterKey>{ { 0, 1 }, { 2 } }));
    EXPECT_EQ(
        partition.boundary_sample_dependency_list.front()
            .contributor_atom_index_list,
        (std::vector<std::size_t>{ 1, 2 }));
    EXPECT_EQ(
        partition.boundary_sample_dependency_list.back().sample_id,
        (coupling_detail::SampleRef{ 1, 0 }));
    auto contributor_changed_partition{ partition };
    contributor_changed_partition.boundary_sample_dependency_list.front()
        .contributor_atom_index_list = { 1 };
    EXPECT_NE(
        contributor_changed_partition.boundary_sample_dependency_list,
        partition.boundary_sample_dependency_list);
    EXPECT_EQ(partition.sample_id_list_by_key.at({ 0, 1 }).size(), 3U);
    EXPECT_EQ(partition.sample_id_list_by_key.at({ 2 }).size(), 2U);

    const auto key_list{
        coupling_detail::BuildGraphClusterKeyList(partition)
    };
    EXPECT_EQ(key_list, (std::vector<std::vector<std::size_t>>{ { 0, 1 }, { 2 } }));
    const auto affected_sample_list{
        coupling_detail::BuildGraphAffectedSampleUnion(
            partition,
            key_list)
    };
    EXPECT_EQ(affected_sample_list.size(), 3U);
    EXPECT_EQ(affected_sample_list.at(0).atom_index, 0U);
    EXPECT_EQ(affected_sample_list.at(0).sample_index, 0U);
    EXPECT_EQ(affected_sample_list.at(1).atom_index, 0U);
    EXPECT_EQ(affected_sample_list.at(1).sample_index, 1U);
    EXPECT_EQ(affected_sample_list.at(2).atom_index, 1U);
    EXPECT_EQ(affected_sample_list.at(2).sample_index, 0U);

    const auto inactive_partition{
        coupling_detail::BuildGraphPartition(topology, { 2, 0 })
    };
    EXPECT_EQ(inactive_partition.sample_id_list_by_key.count({ 0 }), 1U);
    EXPECT_EQ(inactive_partition.sample_id_list_by_key.count({ 2 }), 1U);
    EXPECT_TRUE(inactive_partition.boundary_sample_dependency_list.empty());
}

TEST(EstimatorSecondStageDefenseTest, BoundaryReconciliationComponentsUseAcceptedSharedSamples)
{
    const audit_detail::ClusterKey key_a{ 0 };
    const audit_detail::ClusterKey key_b{ 1 };
    const audit_detail::ClusterKey key_c{ 2 };
    const audit_detail::ClusterKey key_d{ 3 };
    const audit_detail::ClusterKey key_e{ 4 };
    const audit_detail::ClusterKey key_f{ 5 };
    const coupling_detail::SampleRef sample_ab{ 0, 0 };
    const coupling_detail::SampleRef sample_bc{ 1, 0 };
    const coupling_detail::SampleRef sample_de{ 3, 0 };
    coupling_detail::CouplingGraphPartition partition;
    partition.sample_id_list_by_key = {
        { key_a, { sample_ab } },
        { key_b, { sample_ab, sample_bc } },
        { key_c, { sample_bc } },
        { key_d, { sample_de } },
        { key_e, { sample_de } },
        { key_f, {} }
    };
    partition.boundary_sample_dependency_list = {
        { sample_ab, { key_a, key_b }, { 0, 1 } },
        { sample_bc, { key_b, key_c }, { 1, 2 } },
        { sample_de, { key_d, key_e }, { 3, 4 } }
    };
    coupling_detail::SecondStageContext context;
    context.atom_list.resize(6);
    const auto component_list{
        coupling_detail::BuildBoundaryReconciliationComponents(
            context,
            partition,
            { key_f, key_e, key_c, key_a, key_d, key_b })
    };
    ASSERT_EQ(component_list.size(), 2U);
    EXPECT_EQ(
        component_list.at(0).key_list,
        (std::vector<audit_detail::ClusterKey>{ key_a, key_b, key_c }));
    EXPECT_EQ(
        component_list.at(0).affected_sample_ref_list,
        (std::vector<coupling_detail::SampleRef>{ sample_ab, sample_bc }));
    EXPECT_EQ(component_list.at(0).boundary_sample_count, 2U);
    EXPECT_EQ(
        component_list.at(0).interface_atom_index_list,
        (std::vector<std::size_t>{ 0, 1, 2 }));
    EXPECT_TRUE(component_list.at(0).halo_atom_index_list.empty());
    EXPECT_EQ(
        component_list.at(1).key_list,
        (std::vector<audit_detail::ClusterKey>{ key_d, key_e }));
    EXPECT_EQ(
        component_list.at(1).affected_sample_ref_list,
        (std::vector<coupling_detail::SampleRef>{ sample_de }));
    EXPECT_EQ(component_list.at(1).boundary_sample_count, 1U);
    EXPECT_EQ(
        component_list.at(1).interface_atom_index_list,
        (std::vector<std::size_t>{ 3, 4 }));
    EXPECT_TRUE(component_list.at(1).halo_atom_index_list.empty());
    for (const auto & component : component_list)
    {
        const auto expanded{ coupling_detail::ExpandBoundaryReconciliationHalo(context, component, 0) };
        EXPECT_EQ(expanded.halo_atom_index_list, component.interface_atom_index_list);
        EXPECT_EQ(expanded.interface_atom_index_list, component.interface_atom_index_list);
        EXPECT_EQ(expanded.key_list, component.key_list);
        EXPECT_EQ(expanded.affected_sample_ref_list, component.affected_sample_ref_list);
    }

    EXPECT_TRUE(
        coupling_detail::BuildBoundaryReconciliationComponents(
            context,
            partition,
            { key_c, key_a }).empty());
    EXPECT_EQ(
        coupling_detail::BuildBoundaryReconciliationComponents(
            context,
            partition,
            { key_c, key_b, key_a }),
        std::vector<coupling_detail::BoundaryReconciliationComponent>{
            component_list.front()
        });
}

TEST(EstimatorSecondStageDefenseTest, BoundaryPhysicalHaloRejectsOutOfRangeAtoms)
{
    const audit_detail::ClusterKey key_a{ 0 };
    const audit_detail::ClusterKey key_b{ 1 };
    coupling_detail::CouplingGraphPartition partition;
    partition.sample_id_list_by_key = {
        { key_a, {} },
        { key_b, {} }
    };
    partition.boundary_sample_dependency_list = {
        { { 0, 0 }, { key_a, key_b }, { 0, 1 } }
    };
    coupling_detail::SecondStageContext context;
    context.atom_list.resize(1);

    EXPECT_THROW(
        coupling_detail::BuildBoundaryReconciliationComponents(
            context,
            partition,
            { key_a, key_b }),
        std::invalid_argument);

    partition.boundary_sample_dependency_list.front()
        .contributor_atom_index_list = { 0 };
    EXPECT_THROW(
        coupling_detail::BuildBoundaryReconciliationComponents(
            context,
            partition,
            { key_a, key_b }),
        std::invalid_argument);

    const coupling_detail::BoundaryReconciliationComponent invalid_interface{
        .key_list = { { 0, 1 } },
        .affected_sample_ref_list = {},
        .interface_atom_index_list = { 1 },
        .halo_atom_index_list = {},
            .boundary_sample_count = 0
    };
    EXPECT_THROW(
        coupling_detail::ExpandBoundaryReconciliationHalo(
            context,
            invalid_interface,
            0),
        std::invalid_argument);

    auto invalid_component{ invalid_interface };
    invalid_component.interface_atom_index_list = { 0 };
    EXPECT_THROW(
        coupling_detail::ExpandBoundaryReconciliationHalo(
            context,
            invalid_component,
            0),
        std::invalid_argument);
}

TEST(EstimatorSecondStageDefenseTest, BoundaryHaloExpandsPhysicalParticipantsByHop)
{
    coupling_detail::SecondStageContext context;
    context.atom_list.resize(5);
    for (auto & atom_context : context.atom_list)
    {
        atom_context.raw_sampling_entries.resize(1);
        atom_context.neighbor_atom_sample_offset_list = { 0, 0 };
    }
    context.at(0).neighbor_atom_sample_list = {
        { 1, 0.5 },
        { 3, 0.5 }
    };
    context.at(0).neighbor_atom_sample_offset_list = { 0, 2 };
    context.at(1).neighbor_atom_sample_list = {
        { 2, 0.5 }
    };
    context.at(1).neighbor_atom_sample_offset_list = { 0, 1 };
    context.at(2).neighbor_atom_sample_list = {
        { 3, 0.5 },
        { 4, 0.5 }
    };
    context.at(2).neighbor_atom_sample_offset_list = { 0, 2 };
    context.at(2).unselected_distance_list_by_sample = { { 0.5 } };

    const coupling_detail::BoundaryReconciliationComponent component{
        .key_list = { { 0, 1 }, { 2, 3 } },
        .affected_sample_ref_list = { { 0, 0 }, { 1, 0 }, { 2, 0 } },
        .interface_atom_index_list = { 0 },
        .halo_atom_index_list = {},
        .boundary_sample_count = 1
    };
    const auto depth_zero{
        coupling_detail::ExpandBoundaryReconciliationHalo(context, component, 0)
    };
    EXPECT_EQ(depth_zero.interface_atom_index_list, (std::vector<std::size_t>{ 0 }));
    EXPECT_EQ(depth_zero.halo_atom_index_list, (std::vector<std::size_t>{ 0 }));

    const auto depth_one{
        coupling_detail::ExpandBoundaryReconciliationHalo(context, component, 1)
    };
    EXPECT_EQ(depth_one.interface_atom_index_list, (std::vector<std::size_t>{ 0 }));
    EXPECT_EQ(
        depth_one.halo_atom_index_list,
        (std::vector<std::size_t>{ 0, 1, 3 }));

    const auto fixed_point{
        coupling_detail::ExpandBoundaryReconciliationHalo(context, component, 10)
    };
    EXPECT_EQ(
        fixed_point.halo_atom_index_list,
        (std::vector<std::size_t>{ 0, 1, 2, 3 }));
}

TEST(EstimatorSecondStageDefenseTest, UncutDependencyPolishMergesWholeActiveClusters)
{
    coupling_detail::GraphTopology topology;
    topology.adjacency_list.resize(6);
    topology.sample_dependency_list = {
        { { 0, 0 }, { 1, 2 } },
        { { 2, 0 }, { 2, 3 } },
        { { 3, 0 }, { 3, 4 } }
    };
    coupling_detail::CouplingGraphPartition partition;
    partition.sample_id_list_by_key = {
        { { 0, 1 }, { { 0, 0 } } },
        { { 2 }, { { 2, 0 } } },
        { { 3 }, { { 3, 0 } } },
        { { 5 }, { { 5, 0 } } }
    };
    const std::vector<audit_detail::ClusterKey> owner_key_by_atom_index{
        { 0, 1 }, { 0, 1 }, { 2 }, { 3 }, {}, { 5 }
    };
    const auto component_list{
        coupling_detail::BuildUncutDependencyPolishComponents(
            topology,
            partition,
            owner_key_by_atom_index)
    };
    ASSERT_EQ(component_list.size(), 1U);
    EXPECT_EQ(
        component_list.front().key_list,
        (std::vector<audit_detail::ClusterKey>{ { 0, 1 }, { 2 }, { 3 } }));
    EXPECT_EQ(
        component_list.front().atom_index_list,
        (std::vector<std::size_t>{ 0, 1, 2, 3 }));
    EXPECT_EQ(
        component_list.front().affected_sample_ref_list,
        (std::vector<coupling_detail::SampleRef>{ { 0, 0 }, { 2, 0 }, { 3, 0 } }));
}

TEST(EstimatorSecondStageDefenseTest, DependencyPolishDefaultsAndIterationValidation)
{
    const rt::FitOptions defaults;
    EXPECT_EQ(defaults.second_stage_boundary_halo_depth, 1U);
    EXPECT_TRUE(defaults.enable_second_stage_dependency_polish);
    EXPECT_EQ(defaults.second_stage_dependency_polish_max_iterations, 10U);

    auto model{ BuildJointPolishDefenseModel() };
    std::vector<rg::GaussianModel3D> original_model_list;
    for (const auto * atom : model->GetSelectedAtoms())
    {
        original_model_list.emplace_back(GetEstimateModel(*atom));
    }
    auto options{ MakeSecondStageOptions() };
    options.second_stage_dependency_polish_max_iterations = 0;
    EXPECT_THROW(
        rt::RunSecondStageLocalFitting(*model, options),
        std::invalid_argument);
    for (std::size_t atom_index = 0;
        atom_index < model->GetSelectedAtoms().size();
        atom_index++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*model->GetSelectedAtoms().at(atom_index)),
            original_model_list.at(atom_index),
            0.0);
    }
}

TEST(EstimatorSecondStageDefenseTest, FinalDependencyPolishImprovesUncutComponent)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        { 6.0, 0.55, 0.0 },
        { 4.5, 0.70, 0.0 }
    };
    const std::vector<rg::GaussianModel3D> target_model_list{
        { 6.5, 0.60, 0.15 },
        { 4.0, 0.65, -0.15 }
    };
    auto fixture{
        BuildJointPolishFixture(
            base_model_list,
            target_model_list)
    };
    coupling_detail::GraphTopology topology;
    topology.adjacency_list.resize(2);
    for (const auto & sample_ref : fixture.sample_ref_list)
    {
        topology.sample_dependency_list.emplace_back(
            coupling_detail::GraphSampleDependency{
                sample_ref,
                sample_ref == fixture.sample_ref_list.front() ?
                    std::vector<std::size_t>{ 0, 1 } :
                    std::vector<std::size_t>{ sample_ref.atom_index }
            });
    }
    const auto partition{
        coupling_detail::BuildGraphPartition(topology, { 0, 1 })
    };
    const auto base_snapshot{
        polish_detail::BuildSecondStageModelSnapshot(
            fixture.context,
            fixture.state)
    };
    const auto objective_domain{
        audit_detail::BuildObjectiveDomain(
            fixture.context,
            base_snapshot,
            coupling_detail::BuildGraphClusterKeyList(partition),
            0.0,
            1.0)
    };
    polish_detail::TrustRegionStateSet trust_region_state;
    trust_region_state.Reconcile(
        coupling_detail::BuildGraphClusterKeyList(partition));
    polish_detail::ClusterSolverWorkspaceMap solver_workspace_by_key;
    polish_detail::BoundaryJointCorrectionWorkspaceMap correction_workspace_by_key;
    polish_detail::PerformanceCounters performance_counters{
        true,
        fixture.context,
        solver_workspace_by_key,
        correction_workspace_by_key
    };
    auto options{ MakeSecondStageOptions() };
    const polish_detail::SuspiciousBlockActivity all_active{
        std::vector<char>(fixture.context.size(), 0),
        std::vector<char>(fixture.context.size(), 0),
        std::vector<char>(fixture.context.size(), 0)
    };
    const auto polish_result{
        polish_detail::RunFinalDependencyPolish(
            fixture.context,
            options,
            topology,
            partition,
            objective_domain,
            all_active,
            trust_region_state,
            fixture.state,
            correction_workspace_by_key,
            performance_counters)
    };
    ASSERT_TRUE(polish_result.accepted);
    ASSERT_TRUE(polish_result.objective.has_value());
    ASSERT_TRUE(polish_result.diagnostic.objective_before.has_value());
    ASSERT_TRUE(polish_result.diagnostic.objective_after.has_value());
    EXPECT_LT(
        *polish_result.diagnostic.objective_after,
        *polish_result.diagnostic.objective_before);
    ASSERT_EQ(polish_result.diagnostic.component_list.size(), 1U);
    EXPECT_GE(polish_result.diagnostic.component_list.front().round_count, 1U);
    EXPECT_LE(
        polish_result.diagnostic.component_list.front().round_count,
        options.second_stage_dependency_polish_max_iterations);
    EXPECT_EQ(polish_result.diagnostic.component_list.front().parameter_count, 6U);
    EXPECT_NE(
        polish_result.state.at(0).mdpde.GetModel().GetOffset(),
        polish_result.state.at(1).mdpde.GetModel().GetOffset());
    EXPECT_NE(
        polish_result.state.at(0).mdpde.GetModel().GetAmplitude(),
        fixture.state.at(0).mdpde.GetModel().GetAmplitude());
}

TEST(EstimatorSecondStageDefenseTest, CouplingAtomCutoffBoundsComponentsAndPreservesDependencies)
{
    for (const std::size_t atom_count : { 100U, 101U, 102U })
    {
        coupling_detail::GraphTopology topology;
        topology.adjacency_list.resize(atom_count);
        // The 102-atom case also has an isolated atom, with no forced closure.
        const std::size_t first_connected_atom{ atom_count == 102 ? 1U : 0U };
        for (std::size_t atom_index = first_connected_atom; atom_index + 1 < atom_count; atom_index++)
        {
            topology.retained_edge_list.push_back({
                atom_index, atom_index + 1, 1.0 - 0.001 * static_cast<double>(atom_index) });
        }
        std::vector<std::size_t> active_index_list(atom_count);
        std::iota(active_index_list.begin(), active_index_list.end(), 0U);
        topology.sample_dependency_list = { { { 0, 0 }, active_index_list } };
        const auto capped_topology{ coupling_detail::ApplyGraphAtomCutoff(topology, 100) };
        EXPECT_EQ(capped_topology.atom_cutoff_summary.atom_count, atom_count);
        EXPECT_EQ(capped_topology.atom_cutoff_summary.maximum_atom_count, 100U);
        EXPECT_EQ(capped_topology.atom_cutoff_summary.cluster_count, atom_count - 99);
        EXPECT_EQ(capped_topology.atom_cutoff_summary.cut_edge_count, atom_count == 100 ? 0U : 1U);
        EXPECT_EQ(capped_topology.retained_edge_list.size(), topology.retained_edge_list.size());
        ASSERT_EQ(capped_topology.sample_dependency_list.size(), 1U);
        EXPECT_EQ(capped_topology.sample_dependency_list.front().contributor_atom_index_list,
            active_index_list);
        const auto partition{ coupling_detail::BuildGraphPartition(capped_topology, active_index_list) };
        EXPECT_EQ(partition.boundary_sample_dependency_list.size(), atom_count == 100 ? 0U : 1U);
        EXPECT_EQ(partition.sample_id_list_by_key.size(), atom_count - 99);
        std::vector<coupling_detail::ClusterKey> owner_key_by_atom_index(atom_count);
        for (const auto & [key, sample_id_list] : partition.sample_id_list_by_key)
        {
            EXPECT_EQ(sample_id_list.size(), 1U);
            EXPECT_LE(key.size(), 100U);
            for (const auto atom_index : key) owner_key_by_atom_index.at(atom_index) = key;
        }
        if (first_connected_atom == 1)
            EXPECT_EQ(partition.sample_id_list_by_key.count({ 0 }), 1U);
        const auto polish_components{ coupling_detail::BuildUncutDependencyPolishComponents(
            capped_topology, partition, owner_key_by_atom_index) };
        ASSERT_EQ(polish_components.size(), 1U);
        EXPECT_EQ(polish_components.front().atom_index_list.size(), atom_count);

        std::reverse(active_index_list.begin(), active_index_list.end());
        const auto reversed_partition{ coupling_detail::BuildGraphPartition(capped_topology, active_index_list) };
        EXPECT_EQ(partition.sample_id_list_by_key, reversed_partition.sample_id_list_by_key);
        ASSERT_EQ(partition.boundary_sample_dependency_list.size(), reversed_partition.boundary_sample_dependency_list.size());
        if (!partition.boundary_sample_dependency_list.empty())
        {
            auto expected{ partition.boundary_sample_dependency_list.front() };
            auto actual{ reversed_partition.boundary_sample_dependency_list.front() };
            std::sort(expected.cluster_key_list.begin(), expected.cluster_key_list.end());
            std::sort(actual.cluster_key_list.begin(), actual.cluster_key_list.end());
            EXPECT_EQ(expected, actual);
        }
    }

    const auto empty{ coupling_detail::CouplingGraphBuilder{ 0 }.BuildTopology() };
    EXPECT_EQ(empty.atom_cutoff_summary.atom_count, 0U);
    EXPECT_EQ(empty.atom_cutoff_summary.cluster_count, 0U);
    EXPECT_EQ(empty.atom_cutoff_summary.maximum_atom_count, 0U);
    EXPECT_TRUE(coupling_detail::BuildGraphPartition(empty, {}).sample_id_list_by_key.empty());
    EXPECT_THROW(coupling_detail::ApplyGraphAtomCutoff({}, 0), std::invalid_argument);
    coupling_detail::CouplingGraphOptions invalid_options;
    invalid_options.maximum_atom_count = 0;
    EXPECT_THROW(coupling_detail::CouplingGraphBuilder{ 0 }.BuildTopology(invalid_options),
        std::invalid_argument);
}

TEST(EstimatorSecondStageDefenseTest, CouplingAtomCutoffPrioritizesStrongEdgesAndStableTies)
{
    coupling_detail::GraphTopology topology;
    topology.adjacency_list.resize(3);
    topology.retained_edge_list = { { 1, 2, 0.80 }, { 0, 1, 0.90 }, { 0, 2, 0.70 } };
    const auto capped_topology{ coupling_detail::ApplyGraphAtomCutoff(topology, 2) };
    const auto partition{ coupling_detail::BuildGraphPartition(capped_topology, { 2, 1, 0 }) };
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 0, 1 }), 1U);
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 2 }), 1U);
    EXPECT_TRUE(HasCouplingNeighbor(capped_topology, 0, 1));
    EXPECT_FALSE(HasCouplingNeighbor(capped_topology, 1, 2));

    const auto singleton_topology{ coupling_detail::ApplyGraphAtomCutoff(topology, 1) };
    EXPECT_EQ(singleton_topology.atom_cutoff_summary.cluster_count, 3U);
    EXPECT_EQ(singleton_topology.atom_cutoff_summary.cut_edge_count, 3U);
    const auto whole_topology{ coupling_detail::ApplyGraphAtomCutoff(topology, 3) };
    EXPECT_EQ(whole_topology.atom_cutoff_summary.cluster_count, 1U);
    // All three internal edges survive, not only the union-find spanning tree.
    for (const auto & neighbors : whole_topology.adjacency_list) EXPECT_EQ(neighbors.size(), 2U);

    for (auto & edge : topology.retained_edge_list) edge.weight = 0.8;
    for (std::size_t order = 0; order < 3; order++)
    {
        std::rotate(topology.retained_edge_list.begin(),
            topology.retained_edge_list.begin() + 1, topology.retained_edge_list.end());
        for (auto & edge : topology.retained_edge_list)
            std::swap(edge.left_atom_index, edge.right_atom_index);
        const auto tied_topology{ coupling_detail::ApplyGraphAtomCutoff(topology, 2) };
        EXPECT_EQ(coupling_detail::BuildGraphPartition(tied_topology, { 0, 1, 2 }).sample_id_list_by_key,
            partition.sample_id_list_by_key);
    }

    topology.retained_edge_list.front().weight = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(coupling_detail::ApplyGraphAtomCutoff(topology, 2), std::invalid_argument);
    topology.retained_edge_list.front().weight = 1.0;
    topology.retained_edge_list.front().left_atom_index = 3;
    EXPECT_THROW(coupling_detail::ApplyGraphAtomCutoff(topology, 2), std::invalid_argument);
}

TEST(EstimatorSecondStageDefenseTest, CouplingPartitionKeepsStrongChainAndBinaryFallback)
{
    coupling_detail::GraphTopology strong_topology;
    strong_topology.adjacency_list = {
        { 1 },
        { 0, 2 },
        { 1 }
    };
    const auto strong_partition{
        coupling_detail::BuildGraphPartition(
            strong_topology,
            { 0, 1, 2 })
    };
    EXPECT_EQ(strong_partition.sample_id_list_by_key.count({ 0, 1, 2 }), 1U);

    coupling_detail::CouplingGraphBuilder builder{ 2 };
    const auto invalid{
        Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN())
    };
    AddCouplingGraphSample(
        builder,
        { 0, 0 },
        { { 0, Eigen::Vector3d::Ones() }, { 1, invalid } });
    coupling_detail::CouplingGraphOptions fallback_options;
    fallback_options.sensitivity_minimum_weight_list = { 0.05, 0.10 };
    const auto binary_topology{
        builder.BuildTopology(fallback_options)
    };
    EXPECT_FALSE(binary_topology.summary.uses_weighted_graph);
    EXPECT_TRUE(binary_topology.summary.threshold_sensitivity_list.empty());
    const auto binary_partition{
        coupling_detail::BuildGraphPartition(
            binary_topology,
            { 0, 1 })
    };
    EXPECT_EQ(binary_partition.sample_id_list_by_key.count({ 0, 1 }), 1U);
    const auto capped_binary_topology{
        coupling_detail::ApplyGraphAtomCutoff(
            binary_topology,
            1)
    };
    const auto capped_binary_partition{
        coupling_detail::BuildGraphPartition(
            capped_binary_topology,
            { 0, 1 })
    };
    EXPECT_EQ(capped_binary_partition.sample_id_list_by_key.size(), 2U);

    coupling_detail::CouplingGraphBuilder overflow_builder{ 2 };
    const auto huge{ Eigen::Vector3d::Constant(1.0e200) };
    AddCouplingGraphSample(
        overflow_builder,
        { 0, 0 },
        { { 0, huge }, { 1, huge } });
    fallback_options.maximum_atom_count = 1;
    const auto overflow_topology{
        overflow_builder.BuildTopology(fallback_options)
    };
    EXPECT_FALSE(overflow_topology.summary.uses_weighted_graph);
    EXPECT_EQ(overflow_topology.atom_cutoff_summary.cluster_count, 2U);
    EXPECT_EQ(overflow_topology.atom_cutoff_summary.maximum_atom_count_limit, 1U);
    EXPECT_FALSE(HasCouplingNeighbor(overflow_topology, 0, 1));
}

TEST(EstimatorSecondStageDefenseTest, TransformedDampingIsIntensityScaleInvariant)
{
    const rg::GaussianModel3D previous{ 8.0, 0.50, -0.10 };
    const rg::GaussianModel3D current{ 9.0, 0.60, -0.15 };
    constexpr double damping{ 0.25 };

    const auto damp = [&](const rg::GaussianModel3D & lhs,
                          const rg::GaussianModel3D & rhs)
    {
        const auto lhs_coordinates{ lhs.ToTransformedCoordinates() };
        const auto rhs_coordinates{ rhs.ToTransformedCoordinates() };
        EXPECT_TRUE(lhs_coordinates.has_value());
        EXPECT_TRUE(rhs_coordinates.has_value());
        return rg::GaussianModel3D::FromTransformedCoordinates(
            *lhs_coordinates + damping * (*rhs_coordinates - *lhs_coordinates));
    };

    const auto base{ damp(previous, current) };
    ASSERT_TRUE(base.has_value());
    for (const auto scale : { 1.0e-2, 1.0e2 })
    {
        const auto scaled{
            damp(
                rg::GaussianModel3D{
                    previous.GetAmplitude() * scale,
                    previous.GetWidth(),
                    previous.GetOffset() * scale
                },
                rg::GaussianModel3D{
                    current.GetAmplitude() * scale,
                    current.GetWidth(),
                    current.GetOffset() * scale
                })
        };
        ASSERT_TRUE(scaled.has_value());
        EXPECT_NEAR(base->GetAmplitude() * scale, scaled->GetAmplitude(), 1.0e-10);
        EXPECT_NEAR(base->GetWidth(), scaled->GetWidth(), 1.0e-12);
        EXPECT_NEAR(base->GetOffset() * scale, scaled->GetOffset(), 1.0e-12);
    }
}

TEST(EstimatorSecondStageDefenseTest,
    BacktrackingWorkspaceGeneratesCandidatesAndMergesProvenance)
{
    const std::vector<rg::GaussianModel3D> previous_model_list{
        rg::GaussianModel3D{ 8.0, 0.50, -0.10 },
        rg::GaussianModel3D{ 10.0, 0.60, 0.20 }
    };
    const std::vector<rg::GaussianModel3D> endpoint_model_list{
        rg::GaussianModel3D{ 12.0, 0.75, 0.40 },
        rg::GaussianModel3D{ 14.0, 0.90, 0.80 }
    };
    const std::vector<rg::GaussianModel3DUncertainty> endpoint_uncertainty_list{
        rg::GaussianModel3DUncertainty{ 0.10, 0.02, 0.03 },
        rg::GaussianModel3DUncertainty{ 0.20, 0.04, 0.05 }
    };

    backtracking_detail::FitState previous_state;
    backtracking_detail::FitState endpoint_state;
    previous_state.resize(previous_model_list.size());
    endpoint_state.resize(endpoint_model_list.size());
    for (std::size_t atom_index = 0;
        atom_index < previous_model_list.size();
        atom_index++)
    {
        previous_state.at(atom_index).mdpde =
            rg::GaussianModel3DWithUncertainty{
                previous_model_list.at(atom_index),
                rg::GaussianModel3DUncertainty{ 0.01, 0.02, 0.03 }
            };
        endpoint_state.at(atom_index).mdpde =
            rg::GaussianModel3DWithUncertainty{
                endpoint_model_list.at(atom_index),
                endpoint_uncertainty_list.at(atom_index)
            };
    }

    const auto endpoint_patch{
        backtracking_detail::FitStatePatch::FromState(
            endpoint_state,
            std::vector<std::size_t>{ 1, 0 })
    };
    backtracking_detail::BacktrackingWorkspace workspace{
        previous_state,
        endpoint_patch,
        1.0e-4
    };
    const auto step{ workspace.BuildNextCandidate() };
    ASSERT_EQ(
        step.status,
        backtracking_detail::BacktrackingStepStatus::CandidateReady);
    EXPECT_DOUBLE_EQ(step.factor, 0.5);
    EXPECT_EQ(step.trial_number, 2U);
    const auto & candidate_patch{ workspace.GetCandidatePatch() };
    EXPECT_EQ(
        candidate_patch.atom_index_list,
        (std::vector<std::size_t>{ 0, 1 }));
    EXPECT_DOUBLE_EQ(
        candidate_patch.mdpde_list.at(0)
            .GetStandardDeviationModel().GetAmplitude(),
        endpoint_uncertainty_list.at(0).GetAmplitude());
    EXPECT_DOUBLE_EQ(
        candidate_patch.mdpde_list.at(1)
            .GetStandardDeviationModel().GetWidth(),
        endpoint_uncertainty_list.at(1).GetWidth());

    auto candidate_state{ previous_state };
    candidate_patch.ApplyTo(candidate_state);
    for (const auto atom_index : candidate_patch.atom_index_list)
    {
        ExpectGaussianModelsNear(
            candidate_state.at(atom_index).mdpde.GetModel(),
            candidate_patch.mdpde_list.at(atom_index).GetModel(),
            1.0e-12);
    }
    EXPECT_DOUBLE_EQ(candidate_state.at(0).mdpde.GetModel().GetOffset(),
        std::lerp(-0.10, 0.40, 0.5));
    EXPECT_DOUBLE_EQ(candidate_state.at(1).mdpde.GetModel().GetOffset(),
        std::lerp(0.20, 0.80, 0.5));

    const auto merged_provenance{
        workspace.BuildCandidatePolishProvenance(
            std::vector<char>{ 0, 1 },
            std::vector<char>{ 1, 0 })
    };
    EXPECT_EQ(merged_provenance, (std::vector<char>{ 1, 0 }));

    backtracking_detail::SecondStageContext median_context;
    median_context.atom_list.resize(3);
    backtracking_detail::FitState median_previous;
    backtracking_detail::FitState median_endpoint;
    const std::array previous_offsets{ 0.0, 10.0, 20.0 };
    const std::array endpoint_offsets{ 30.0, 5.0, 0.0 };
    for (std::size_t node = 0; node < 3; node++)
    {
        median_context.at(node).raw_sampling_entries.resize(1);
        median_context.at(node).unselected_distance_list_by_sample = { { 0.3 } };
        median_previous.emplace_back(MakeGaussianResult({ 8.0 + static_cast<double>(node), 0.5, previous_offsets.at(node) }));
        median_endpoint.emplace_back(MakeGaussianResult({ 17.0 - static_cast<double>(node), 0.7, endpoint_offsets.at(node) }));
    }
    median_context.frozen_background = backtracking_detail::BuildFrozenBackground(median_context, median_previous, { { 0, 1, 2 } });
    ASSERT_TRUE(median_context.frozen_background);
    const auto frozen{ median_context.frozen_background };
    backtracking_detail::BacktrackingWorkspace median_workspace{ median_previous,
        backtracking_detail::FitStatePatch::FromState(median_endpoint, { 0, 1, 2 }), 1.0e-4 };
    ASSERT_EQ(median_workspace.BuildNextCandidate().status,
        backtracking_detail::BacktrackingStepStatus::CandidateReady);
    EXPECT_EQ(median_workspace.GetCandidatePatch().atom_index_list.size(), 3U);
    EXPECT_EQ(median_context.frozen_background, frozen);
    ExpectGaussianModelsNear(frozen->model_by_atom.front(), { 9.0, 0.5, 10.0 }, 1.0e-12);
    const auto next{ backtracking_detail::BuildFrozenBackground(median_context, median_endpoint, { { 0, 1, 2 } }) };
    ASSERT_TRUE(next);
    ExpectGaussianModelsNear(next->model_by_atom.front(), { 16.0, 0.7, 5.0 }, 1.0e-12);
    EXPECT_NE(next->response_by_atom, frozen->response_by_atom);

}

TEST(EstimatorSecondStageDefenseTest,
    BacktrackingWorkspaceStopsWhenChangeBecomesNonmaterial)
{
    backtracking_detail::SecondStageContext context;
    context.atom_list.resize(1);
    backtracking_detail::FitState previous_state(1);
    previous_state.at(0).mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 8.0, 0.50, -0.10 },
        rg::GaussianModel3DUncertainty{ 0.1, 0.02, 0.03 }
    };

    backtracking_detail::FitState endpoint_state(1);
    endpoint_state.at(0).mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 12.0, 0.75, 0.40 },
        rg::GaussianModel3DUncertainty{ 0.2, 0.04, 0.05 }
    };
    const auto endpoint_patch{
        backtracking_detail::FitStatePatch::FromState(
            endpoint_state,
            std::vector<std::size_t>{ 0 })
    };
    backtracking_detail::BacktrackingWorkspace change_exhausted_workspace{
        previous_state,
        endpoint_patch,
        1.0e6
    };
    const auto change_exhausted_step{
        change_exhausted_workspace.BuildNextCandidate()
    };
    EXPECT_EQ(
        change_exhausted_step.status,
        backtracking_detail::BacktrackingStepStatus::Exhausted);
    EXPECT_EQ(change_exhausted_step.trial_number, 1U);
}

TEST(EstimatorSecondStageDefenseTest,
    ResidualBaselineAndOverlayPreserveFrozenBackground)
{
    residual_detail::SecondStageContext context;
    context.atom_list.resize(1);
    context.at(0).neighbor_atom_sample_offset_list = { 0, 0, 0 };
    const rg::GaussianModel3D previous_model{ 8.0, 0.50, -0.10 };
    const rg::GaussianModel3D candidate_model{ 10.0, 0.60, 0.20 };
    constexpr double contributor_distance{ 0.25 };
    for (const auto distance : { 0.15, 0.45 })
    {
        context.at(0).raw_sampling_entries.emplace_back(LocalPotentialSample{
            candidate_model.ResponseAtDistance(distance) + previous_model.ResponseAtDistance(contributor_distance),
            SamplingPoint{ distance } });
        context.at(0).unselected_distance_list_by_sample.push_back({ contributor_distance });
    }
    const residual_detail::FitState previous_state{ MakeGaussianResult(previous_model) };
    context.frozen_background = residual_detail::BuildFrozenBackground(context, previous_state, { { 0 } });
    ASSERT_TRUE(context.frozen_background);
    const auto baseline{ residual_detail::BuildResidualBaseline(context, previous_state) };
    ASSERT_TRUE(baseline.sample_list.at(0).at(1).has_value());
    EXPECT_NEAR(baseline.sample_list.at(0).at(1)->adjusted_response, candidate_model.ResponseAtDistance(0.45), 1.0e-12);
    residual_detail::FitStatePatch patch;
    patch.atom_index_list = { 0 };
    patch.mdpde_list = { MakeGaussianResult(candidate_model).mdpde };
    const residual_detail::FitStateView candidate_view{ previous_state, patch };
    const residual_detail::CandidateEvaluationOverlay overlay{ context, baseline, candidate_view };
    const residual_detail::SampleRef sample_ref{ 0, 1 };
    const auto candidate_snapshot{ residual_detail::BuildSecondStageModelSnapshot(
        context, residual_detail::BuildFittedGaussianSnapshot(candidate_view)) };
    const auto direct{ residual_detail::EvaluateResidualSample(context, sample_ref, candidate_snapshot) };
    const auto overlaid{ overlay(sample_ref) };
    ASSERT_TRUE(direct.has_value());
    ASSERT_TRUE(overlaid.has_value());
    EXPECT_NEAR(direct->residual, 0.0, 1.0e-12);
    EXPECT_DOUBLE_EQ(direct->adjusted_response, overlaid->adjusted_response);
    EXPECT_DOUBLE_EQ(direct->residual, overlaid->residual);
    const auto samples{ residual_detail::BuildSecondStageAdjustedSamples(context, 0, candidate_snapshot) };
    EXPECT_DOUBLE_EQ(samples.at(1).response, direct->adjusted_response);
    EXPECT_EQ(candidate_view.size(), 1U);

    // A later refresh must not change an already captured snapshot or overlay.
    const residual_detail::FitState next_state{ MakeGaussianResult(candidate_model) };
    context.frozen_background = residual_detail::BuildFrozenBackground(context, next_state, { { 0 } });
    ASSERT_TRUE(context.frozen_background);
    EXPECT_NE(context.frozen_background->response_by_atom, candidate_snapshot.frozen_background->response_by_atom);
    const auto old_response{ residual_detail::EvaluateResidualSample(context, sample_ref, candidate_snapshot) };
    ASSERT_TRUE(old_response.has_value());
    EXPECT_DOUBLE_EQ(old_response->residual, direct->residual);
    EXPECT_DOUBLE_EQ(overlay(sample_ref)->residual, direct->residual);
    const auto refreshed{ residual_detail::BuildResidualBaseline(context, next_state) };
    ASSERT_TRUE(refreshed.sample_list.front().back().has_value());
    EXPECT_NE(refreshed.sample_list.front().back()->residual, direct->residual);

}

TEST(EstimatorSecondStageDefenseTest, AuditObjectiveSourcesAgreeAcrossTailPartitions)
{
    audit_detail::SecondStageContext context;
    context.atom_list.resize(1);
    context.at(0).neighbor_atom_sample_offset_list = { 0, 0, 0 };

    const rg::GaussianModel3D model{ 8.0, 0.50, -0.10 };
    for (const auto distance : { 0.15, 0.45 })
    {
        context.at(0).raw_sampling_entries.emplace_back(
            LocalPotentialSample{
                model.ResponseAtDistance(distance),
                SamplingPoint{ distance }
            });
    }

    audit_detail::FitState state;
    state.emplace_back(MakeGaussianResult(model));
    const auto baseline{ audit_detail::BuildResidualBaseline(context, state) };
    const auto cluster_key_list{
        std::vector<audit_detail::ClusterKey>{ audit_detail::ClusterKey{ 0 } }
    };

    for (const auto distance_max : { 0.30, 1.0 })
    {
        const auto domain{
            audit_detail::BuildObjectiveDomain(
                context,
                baseline.model_snapshot,
                cluster_key_list,
                0.0,
                distance_max)
        };
        const auto snapshot_objective{
            audit_detail::EvaluateAuditObjective(
                domain,
                audit_detail::SnapshotResidualEvaluator{
                    context,
                    baseline.model_snapshot
                })
        };
        const auto baseline_objective{
            audit_detail::EvaluateAuditObjective(domain, baseline)
        };
        ASSERT_TRUE(snapshot_objective.has_value());
        ASSERT_TRUE(baseline_objective.has_value());
        EXPECT_DOUBLE_EQ(
            snapshot_objective->fit_range_residual_objective,
            baseline_objective->fit_range_residual_objective);
        EXPECT_DOUBLE_EQ(
            snapshot_objective->tail_validation_loss,
            baseline_objective->tail_validation_loss);
        EXPECT_DOUBLE_EQ(
            snapshot_objective->offset_plausibility_penalty,
            baseline_objective->offset_plausibility_penalty);
        EXPECT_DOUBLE_EQ(
            snapshot_objective->GetTotalObjective(),
            baseline_objective->GetTotalObjective());
    }
}

TEST(EstimatorSecondStageDefenseTest, TransformedBacktrackingIncludesOffset)
{
    const rg::GaussianModel3D previous{ 8.0, 0.50, -0.10 };
    const rg::GaussianModel3D endpoint{ 12.0, 0.75, 0.40 };
    const auto previous_coordinates{ previous.ToTransformedCoordinates() };
    const auto endpoint_coordinates{ endpoint.ToTransformedCoordinates() };
    ASSERT_TRUE(previous_coordinates.has_value());
    ASSERT_TRUE(endpoint_coordinates.has_value());

    double previous_offset_distance{
        std::abs(endpoint.GetOffset() - previous.GetOffset())
    };
    for (const auto factor : { 0.5, 0.25, 0.125 })
    {
        const auto candidate{
            rg::GaussianModel3D::FromTransformedCoordinates(
                *previous_coordinates +
                factor * (*endpoint_coordinates - *previous_coordinates))
        };
        ASSERT_TRUE(candidate.has_value());
        const auto offset_distance{
            std::abs(candidate->GetOffset() - previous.GetOffset())
        };
        EXPECT_LT(offset_distance, previous_offset_distance);
        previous_offset_distance = offset_distance;
    }
}

TEST(EstimatorSecondStageDefenseTest, TransformedExtrapolationKeepsPositiveShape)
{
    const auto left{
        rg::GaussianModel3D{ 8.0, 0.50, -0.10 }.ToTransformedCoordinates()
    };
    const auto right{
        rg::GaussianModel3D{ 9.0, 0.60, -0.15 }.ToTransformedCoordinates()
    };
    ASSERT_TRUE(left.has_value());
    ASSERT_TRUE(right.has_value());

    const auto extrapolated{
        rg::GaussianModel3D::FromTransformedCoordinates(
            2.0 * *right - *left)
    };
    ASSERT_TRUE(extrapolated.has_value());
    EXPECT_GT(extrapolated->GetAmplitude(), 0.0);
    EXPECT_GT(extrapolated->GetWidth(), 0.0);
    EXPECT_TRUE(std::isfinite(extrapolated->GetOffset()));
}

TEST(EstimatorSecondStageDefenseTest, TransformedChangeSeparatesPeakHeightAndWidth)
{
    const auto change{
        change_detail::CalculateTransformedChange(
            rg::GaussianModel3D{ 8.0, 1.0, 0.0 },
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 })
    };

    EXPECT_NEAR(
        0.0,
        change.at(rg::GaussianModel3D::LogPeakHeightCoordinateIndex()),
        1.0e-12);
    EXPECT_NEAR(
        std::log(2.0),
        change.at(rg::GaussianModel3D::LogWidthCoordinateIndex()),
        1.0e-12);
    EXPECT_DOUBLE_EQ(
        0.0,
        change.at(rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()));
}

TEST(EstimatorSecondStageDefenseTest, TransformedConvergenceIgnoresHiddenMaximumTail)
{
    std::vector<change_detail::TransformedChange> change_list(
        1000,
        change_detail::TransformedChange{});
    change_list.back().at(
        rg::GaussianModel3D::LogPeakHeightCoordinateIndex()) =
        2.0e-3;
    const auto summary{ change_detail::SummarizeTransformedChanges(change_list) };
    EXPECT_LT(
        summary.percentile_list.at(
            rg::GaussianModel3D::LogPeakHeightCoordinateIndex()),
        1.0e-4);
    EXPECT_GT(
        summary.maximum_list.at(
            rg::GaussianModel3D::LogPeakHeightCoordinateIndex()),
        1.0e-3);
    EXPECT_TRUE(change_detail::IsTransformedPercentileConverged(summary));
}

TEST(EstimatorSecondStageDefenseTest, ConvergenceCertificateKeepsAcceptedResidualIndependent)
{
    const auto make_summary = [](double percentile, double maximum)
    {
        change_detail::TransformedChangeSummary summary;
        summary.percentile_list.fill(percentile);
        summary.maximum_list.fill(maximum);
        summary.population_size_list.fill(100);
        return summary;
    };
    const auto small{ make_summary(5.0e-5, 2.0e-3) };
    const auto large{ make_summary(2.0e-4, 2.0e-3) };
    audit_detail::ConvergenceCertificate certificate;
    certificate.accepted_active_movement = small;
    certificate.operator_nominal_residual = large;
    certificate.solver_qualified = true;
    EXPECT_FALSE(certificate.ProductionConverged());
}

TEST(EstimatorSecondStageDefenseTest,
    ActiveCoordinatePopulationKeepsMovingSelectedNodesWithoutFixedDilution)
{
    constexpr std::size_t atom_count{ 1000 };
    constexpr std::size_t fixed_atom_count{ 990 };
    std::vector<change_detail::TransformedChange> change_list(
        atom_count,
        change_detail::TransformedChange{});
    for (std::size_t i = fixed_atom_count; i < atom_count; i++)
    {
        change_list.at(i).fill(5.0e-4);
    }
    std::vector<std::size_t> atom_index_list(atom_count);
    for (std::size_t i = 0; i < atom_count; i++) atom_index_list.at(i) = i;

    audit_detail::SuspiciousBlockActivity block_activity{
        audit_detail::SuspiciousUpdateMask(atom_count, 0),
        audit_detail::SuspiciousUpdateMask(atom_count, 0),
        audit_detail::SuspiciousUpdateMask(atom_count, 0)
    };
    for (std::size_t i = 0; i < fixed_atom_count; i++)
    {
        block_activity.shape_fixed_atom_mask.at(i) = 1;
        block_activity.offset_fixed_atom_mask.at(i) = 1;
    }
    const auto all_selected{
        change_detail::SummarizeTransformedChanges(change_list)
    };

    EXPECT_TRUE(change_detail::IsTransformedPercentileConverged(all_selected));

    const auto population{
        audit_detail::BuildActiveCoordinatePopulation(
            atom_index_list,
            block_activity)
    };
    const auto dual_shadow{
        audit_detail::SummarizeActiveDofChanges(change_list, population)
    };
    EXPECT_FALSE(change_detail::IsTransformedPercentileConverged(dual_shadow));
    EXPECT_EQ(
        dual_shadow.population_size_list.at(
            rg::GaussianModel3D::LogPeakHeightCoordinateIndex()),
        10U);
    EXPECT_EQ(
        dual_shadow.population_size_list.at(
            rg::GaussianModel3D::LogWidthCoordinateIndex()),
        10U);
    EXPECT_EQ(
        dual_shadow.population_size_list.at(
            rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()),
        10U);

    std::vector<audit_detail::SuspiciousGaussianAssessment> assessment_by_atom(
        atom_count);
    const auto benign_fixed_mask{
        audit_detail::BuildSuspiciousFailureAtomMask(
            block_activity,
            assessment_by_atom)
    };
    EXPECT_EQ(
        std::ranges::count(benign_fixed_mask, 1),
        0);
    assessment_by_atom.front().reason =
        audit_detail::SuspiciousGaussianReason::WidthGrowth;
    const auto suspicious_fixed_mask{
        audit_detail::BuildSuspiciousFailureAtomMask(
            block_activity,
            assessment_by_atom)
    };
    EXPECT_EQ(
        std::ranges::count(suspicious_fixed_mask, 1),
        1);
}

TEST(EstimatorSecondStageDefenseTest,
    ActiveCoordinatePopulationCountsEverySelectedOffset)
{
    constexpr std::size_t stable_atom_count{ 100 };
    constexpr std::size_t atom_count{ stable_atom_count + 1 };
    std::vector<change_detail::TransformedChange> change_list(
        atom_count,
        change_detail::TransformedChange{});
    change_list.back().at(
        rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) =
        5.0e-4;
    std::vector<std::size_t> atom_index_list(atom_count);
    std::iota(atom_index_list.begin(), atom_index_list.end(), 0);
    audit_detail::SuspiciousBlockActivity activity{
        audit_detail::SuspiciousUpdateMask(atom_count, 1),
        audit_detail::SuspiciousUpdateMask(atom_count, 0),
        audit_detail::SuspiciousUpdateMask(atom_count, 0)
    };
    const auto population{
        audit_detail::BuildActiveCoordinatePopulation(
            atom_index_list,
            activity)
    };
    const auto audit{
        audit_detail::SummarizeActiveDofChanges(change_list, population)
    };

    EXPECT_EQ(
        audit.population_size_list.at(
            rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()),
        atom_count);
    EXPECT_TRUE(change_detail::IsTransformedPercentileConverged(audit));
    EXPECT_TRUE(change_detail::IsTransformedChangeMaterial(change_list.back(), 1.0e-4));
}

TEST(EstimatorSecondStageDefenseTest, ActiveCoordinatePopulationPreservesExtremeAndNonFiniteMembers)
{
    const std::vector<std::size_t> atom_index_list{ 0, 1, 2 };
    audit_detail::SuspiciousBlockActivity activity{
        audit_detail::SuspiciousUpdateMask(3, 1),
        audit_detail::SuspiciousUpdateMask(3, 0),
        audit_detail::SuspiciousUpdateMask(3, 0)
    };
    const auto population{
        audit_detail::BuildActiveCoordinatePopulation(
            atom_index_list,
            activity)
    };
    std::vector<change_detail::TransformedChange> change_list(
        3,
        change_detail::TransformedChange{});
    change_list.at(1).at(
        rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) =
        2.0e-3;
    const auto extreme{
        audit_detail::SummarizeActiveDofChanges(change_list, population)
    };
    EXPECT_DOUBLE_EQ(
        extreme.maximum_list.at(
            rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()),
        2.0e-3);
    EXPECT_FALSE(change_detail::IsTransformedPercentileConverged(extreme));

    change_list.at(1).at(
        rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) =
        std::numeric_limits<double>::quiet_NaN();
    const auto nonfinite{
        audit_detail::SummarizeActiveDofChanges(change_list, population)
    };
    EXPECT_TRUE(std::isinf(nonfinite.maximum_list.at(
        rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex())));
    EXPECT_FALSE(change_detail::IsTransformedPercentileConverged(nonfinite));

    activity.shape_fixed_atom_mask = { 1, 0, 1 };
    activity.offset_fixed_atom_mask.assign(3, 1);
    const auto shape_population{ audit_detail::BuildActiveCoordinatePopulation(
        atom_index_list, activity) };
    EXPECT_EQ(shape_population.active_shape_atom_index_list, (std::vector<std::size_t>{ 1 }));
    for (const auto parameter_index : std::array<std::size_t, 2>{
        rg::GaussianModel3D::LogPeakHeightCoordinateIndex(),
        rg::GaussianModel3D::LogWidthCoordinateIndex() })
    {
        change_list.assign(3, {});
        change_list.at(1).at(parameter_index) = std::numeric_limits<double>::quiet_NaN();
        const auto shape_nonfinite{ audit_detail::SummarizeActiveDofChanges(change_list, shape_population) };
        EXPECT_EQ(shape_nonfinite.population_size_list.at(parameter_index), 1U);
        EXPECT_FALSE(std::isfinite(shape_nonfinite.percentile_list.at(parameter_index)));
        EXPECT_FALSE(change_detail::IsTransformedPercentileConverged(shape_nonfinite));
        EXPECT_EQ(shape_nonfinite.population_size_list.at(
            rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()), 0U);
    }
}

TEST(EstimatorSecondStageDefenseTest, IndependentOffsetActivityRequiresItsOwnClusterQualification)
{
    const std::vector<std::size_t> atoms{ 0, 1 };
    audit_detail::SuspiciousBlockActivity activity{
        { 1, 1 }, { 0, 1 }, { 0, 0 } };
    const auto population{ audit_detail::BuildActiveCoordinatePopulation(atoms, activity) };
    EXPECT_EQ(population.active_offset_atom_index_list, (std::vector<std::size_t>{ 0 }));
    const std::vector<change_detail::TransformedChange> changes(2);
    const auto summary{ audit_detail::SummarizeActiveDofChanges(changes, population) };
    EXPECT_TRUE(change_detail::IsTransformedPercentileConverged(summary));
    const std::vector<std::optional<rg::RHBMEstimationStatus>> status(2);
    EXPECT_FALSE(audit_detail::AreActiveCoordinatesSolverQualified(
        atoms, { atoms }, activity, status, {}));
    audit_detail::ClusterHealthMap health;
    health.emplace(atoms, audit_detail::ClusterHealth{
        audit_detail::JointOffsetSolveStatus::Converged });
    EXPECT_TRUE(audit_detail::AreActiveCoordinatesSolverQualified(
        atoms, { atoms }, activity, status, health));
    EXPECT_FALSE(audit_detail::AreActiveCoordinatesSolverQualified(
        atoms, {}, activity, status, health));
}

TEST(EstimatorSecondStageDefenseTest, ConvergenceCertificateSeparatesAcceptedAndNominalPopulations)
{
    const std::vector<std::size_t> atom_index_list{ 0, 1, 2 };
    audit_detail::SuspiciousBlockActivity activity{
        audit_detail::SuspiciousUpdateMask{ 1, 1, 0 },
        audit_detail::SuspiciousUpdateMask{ 1, 1, 0 },
        audit_detail::SuspiciousUpdateMask(3, 0)
    };
    audit_detail::SuspiciousBlockActivity nominal_activity{
        audit_detail::SuspiciousUpdateMask(3, 0),
        audit_detail::SuspiciousUpdateMask(3, 0),
        audit_detail::SuspiciousUpdateMask(3, 0)
    };
    const auto accepted_population{ audit_detail::BuildActiveCoordinatePopulation(
        atom_index_list,
        activity) };
    const auto nominal_population{ audit_detail::BuildActiveCoordinatePopulation(
        atom_index_list,
        nominal_activity) };
    std::vector<change_detail::TransformedChange> changes(
        3, change_detail::TransformedChange{});
    changes.at(0).fill(2.0e-3);
    changes.at(1).fill(2.0e-3);
    changes.at(2).fill(5.0e-5);

    audit_detail::ConvergenceCertificate certificate;
    certificate.accepted_active_movement =
        audit_detail::SummarizeActiveDofChanges(changes, accepted_population);
    certificate.operator_nominal_residual =
        audit_detail::SummarizeActiveDofChanges(changes, nominal_population);
    certificate.solver_qualified = true;

    EXPECT_EQ(
        certificate.accepted_active_movement.population_size_list.at(
            rg::GaussianModel3D::LogPeakHeightCoordinateIndex()),
        1U);
    EXPECT_EQ(
        certificate.accepted_active_movement.population_size_list.at(
            rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()),
        1U);
    EXPECT_EQ(
        certificate.operator_nominal_residual.population_size_list.at(
            rg::GaussianModel3D::LogPeakHeightCoordinateIndex()),
        3U);
    EXPECT_EQ(
        certificate.operator_nominal_residual.population_size_list.at(
            rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()),
        3U);
    EXPECT_TRUE(change_detail::IsTransformedPercentileConverged(
        certificate.accepted_active_movement));
    EXPECT_FALSE(change_detail::IsTransformedPercentileConverged(
        certificate.operator_nominal_residual));
    EXPECT_FALSE(certificate.ProductionConverged());

    changes.at(0).fill(5.0e-5);
    changes.at(1).fill(5.0e-5);
    certificate.operator_nominal_residual =
        audit_detail::SummarizeActiveDofChanges(changes, nominal_population);
    EXPECT_TRUE(certificate.ProductionConverged());
}

TEST(EstimatorSecondStageDefenseTest, ConvergenceCertificateQualifiesIndependentOffsetActivity)
{
    const std::vector<std::size_t> atom_index_list{ 0, 1 };
    audit_detail::SuspiciousBlockActivity activity{
        audit_detail::SuspiciousUpdateMask(2, 0),
        audit_detail::SuspiciousUpdateMask{ 0, 1 },
        audit_detail::SuspiciousUpdateMask(2, 0)
    };
    const std::vector<audit_detail::ClusterKey> cluster_key_list{
        atom_index_list
    };
    const auto population{ audit_detail::BuildActiveCoordinatePopulation(
        atom_index_list,
        activity) };
    std::vector<change_detail::TransformedChange> changes(
        2, change_detail::TransformedChange{});

    audit_detail::ConvergenceCertificate certificate;
    certificate.accepted_active_movement =
        audit_detail::SummarizeActiveDofChanges(changes, population);
    const audit_detail::SuspiciousBlockActivity nominal{ { 0, 0 }, { 0, 0 }, { 0, 0 } };
    const auto nominal_population{ audit_detail::BuildActiveCoordinatePopulation(atom_index_list, nominal) };
    certificate.operator_nominal_residual =
        audit_detail::SummarizeActiveDofChanges(changes, nominal_population);
    const std::vector<std::optional<rg::RHBMEstimationStatus>>
        local_refit_status_by_atom(2, rg::RHBMEstimationStatus::SUCCESS);
    audit_detail::ClusterHealthMap health_by_key;
    health_by_key.emplace(
        atom_index_list,
        audit_detail::ClusterHealth{
            audit_detail::JointOffsetSolveStatus::Converged
        });
    certificate.solver_qualified =
        audit_detail::AreActiveCoordinatesSolverQualified(
            atom_index_list,
            cluster_key_list,
            activity,
            local_refit_status_by_atom,
            health_by_key);

    EXPECT_TRUE(certificate.solver_qualified);
    EXPECT_TRUE(certificate.StrictOperatorPassed());
    EXPECT_TRUE(certificate.ProductionConverged());
    changes.at(1).at(rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) = 2.0e-3;
    certificate.operator_nominal_residual =
        audit_detail::SummarizeActiveDofChanges(changes, nominal_population);
    EXPECT_FALSE(certificate.StrictOperatorPassed());
    EXPECT_FALSE(certificate.ProductionConverged());
}

TEST(EstimatorSecondStageDefenseTest, ConvergenceCertificateAllFixedStillRequiresCompleteOperator)
{
    const auto make_summary = [](double value, std::size_t population)
    {
        change_detail::TransformedChangeSummary summary;
        summary.percentile_list.fill(value);
        summary.maximum_list.fill(value);
        summary.population_size_list.fill(population);
        return summary;
    };

    audit_detail::ConvergenceCertificate certificate;
    certificate.accepted_active_movement = make_summary(0.0, 0);
    certificate.operator_nominal_residual = make_summary(5.0e-5, 4);
    certificate.solver_qualified = true;

    EXPECT_TRUE(change_detail::IsTransformedPercentileConverged(
        certificate.accepted_active_movement));
    EXPECT_TRUE(certificate.StrictOperatorPassed());
    EXPECT_TRUE(certificate.ProductionConverged());

    certificate.operator_complete = false;
    EXPECT_FALSE(certificate.StrictOperatorPassed());
    EXPECT_FALSE(certificate.ProductionConverged());
}

TEST(EstimatorSecondStageDefenseTest, NonFiniteChangeFailsPercentilePredicate)
{
    change_detail::TransformedChangeSummary summary;
    summary.percentile_list.fill(0.0);
    summary.maximum_list.fill(0.0);
    summary.population_size_list.fill(1);
    summary.percentile_list.at(
        rg::GaussianModel3D::LogPeakHeightCoordinateIndex()) =
        std::numeric_limits<double>::infinity();

    EXPECT_FALSE(change_detail::IsTransformedPercentileConverged(summary));
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

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    for (std::size_t i = 0; i < selected_atoms.size(); i++)
    {
        const auto fitted_model{ GetEstimateModel(*selected_atoms.at(i)) };
        ExpectGaussianModelsNear(fitted_model, previous_model_list.at(i), 1.0e-12);
    }
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingFallsBackWhenJointOffsetSamplesAreNonFinite)
{
    auto model{ BuildNonFiniteJointOffsetDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    const auto previous_model{ GetEstimateModel(*atom) };
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    ExpectGaussianModelsNear(GetEstimateModel(*atom), previous_model, 1.0e-12);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, SystemBuildFailureDoesNotBlockRemoteCluster)
{
    auto model{ BuildSeparatedSystemBuildFailureDefenseModel() };
    const auto initial_remote_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 2, 4)
    };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

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

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

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

TEST(EstimatorSecondStageDefenseTest, PersistentQuarantineReasonRequiresStableReasonAndReleasesOnProbation)
{
    const audit_detail::QuarantineTarget target{
        audit_detail::QuarantineTargetKind::ShapeAtom,
        { 0 }
    };
    audit_detail::QuarantineFailureStateMap state_by_target;
    const auto observe = [&](
        audit_detail::SuspiciousGaussianReason reason,
        std::size_t accepted_iteration_count)
    {
        return audit_detail::UpdateQuarantineFailureState(
            {
                {
                    target,
                    audit_detail::StabilizationTerminalFailure{
                        audit_detail::StabilizationTerminalReason::GuardInfeasible,
                        reason
                    }
                }
            },
            {},
            accepted_iteration_count,
            state_by_target);
    };

    EXPECT_TRUE(observe(audit_detail::SuspiciousGaussianReason::WidthGrowth, 1)
        .entered_target_list.empty());
    EXPECT_NE(
        state_by_target.at(target).lifecycle,
        audit_detail::QuarantineLifecycle::Exhausted);
    EXPECT_TRUE(observe(audit_detail::SuspiciousGaussianReason::WidthGrowth, 2)
        .entered_target_list.empty());
    EXPECT_TRUE(observe(
        audit_detail::SuspiciousGaussianReason::AmplitudeOffsetCompensation,
        3).entered_target_list.empty());
    EXPECT_EQ(state_by_target.at(target).stable_iteration_count, 1U);
    for (std::size_t accepted_iteration = 4;
        accepted_iteration < 7;
        accepted_iteration++)
    {
        EXPECT_TRUE(observe(
            audit_detail::SuspiciousGaussianReason::AmplitudeOffsetCompensation,
            accepted_iteration).entered_target_list.empty());
    }
    const auto entered{
        observe(
            audit_detail::SuspiciousGaussianReason::AmplitudeOffsetCompensation,
            7)
    };
    ASSERT_EQ(entered.entered_target_list, (std::vector{ target }));
    ASSERT_EQ(
        state_by_target.at(target).lifecycle,
        audit_detail::QuarantineLifecycle::Quarantined);
    EXPECT_EQ(
        state_by_target.at(target).next_probation_iteration,
        7U + audit_detail::kQuarantineProbationCooldown);

    state_by_target.at(target).lifecycle =
        audit_detail::QuarantineLifecycle::Probation;
    const auto released{
        audit_detail::UpdateQuarantineFailureState(
            {},
            { target },
            9,
            state_by_target)
    };
    EXPECT_EQ(released.released_target_list, (std::vector{ target }));
    EXPECT_TRUE(state_by_target.empty());

    state_by_target.emplace(
        target,
        audit_detail::QuarantineFailureState{
            audit_detail::StabilizationTerminalFailure{
                audit_detail::StabilizationTerminalReason::GuardInfeasible,
                audit_detail::SuspiciousGaussianReason::WidthGrowth
            },
            audit_detail::kPersistentQuarantineFailureIterationLimit,
            0,
            0,
            audit_detail::QuarantineLifecycle::Quarantined
        });
    for (std::size_t probation = 1;
        probation <= audit_detail::kQuarantineMaximumProbationCount;
        probation++)
    {
        state_by_target.at(target).lifecycle =
            audit_detail::QuarantineLifecycle::Probation;
        const auto failed{
            audit_detail::UpdateQuarantineFailureState(
                {
                    {
                        target,
                        audit_detail::StabilizationTerminalFailure{
                            audit_detail::StabilizationTerminalReason::GuardInfeasible,
                            audit_detail::SuspiciousGaussianReason::WidthGrowth
                        }
                    }
                },
                {},
                9 + probation,
                state_by_target)
        };
        EXPECT_EQ(failed.failed_probation_target_list, (std::vector{ target }));
        EXPECT_EQ(state_by_target.at(target).probation_count, probation);
    }
    EXPECT_EQ(
        state_by_target.at(target).lifecycle,
        audit_detail::QuarantineLifecycle::Exhausted);

    const audit_detail::QuarantineTarget offset_target{
        audit_detail::QuarantineTargetKind::OffsetAtom, { 0 }
    };
    const audit_detail::QuarantineTarget other_offset_target{
        audit_detail::QuarantineTargetKind::OffsetAtom, { 1 }
    };
    auto quarantined{ state_by_target.at(target) };
    quarantined.lifecycle = audit_detail::QuarantineLifecycle::Quarantined;
    for (const auto & released_target : { target, offset_target })
    {
        state_by_target = {
            { target, quarantined },
            { offset_target, quarantined },
            { other_offset_target, quarantined }
        };
        state_by_target.at(released_target).lifecycle =
            audit_detail::QuarantineLifecycle::Probation;
        const auto independent_release{
            audit_detail::UpdateQuarantineFailureState(
                {}, { released_target }, 20, state_by_target)
        };
        EXPECT_EQ(independent_release.released_target_list,
            (std::vector{ released_target }));
        EXPECT_FALSE(state_by_target.contains(released_target));
        EXPECT_EQ(state_by_target.size(), 2U);
        for (const auto & [remaining_target, state] : state_by_target)
        {
            EXPECT_NE(remaining_target, released_target);
            EXPECT_EQ(state.lifecycle, audit_detail::QuarantineLifecycle::Quarantined);
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

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    ExpectGaussianModelsNear(
        GetEstimateModel(*selected_atoms.front()),
        previous_empty_model,
        1.0e-12);
    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 1, 3),
        initial_remote_error);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingRejectsTerminalFiniteNonphysicalProfile)
{
    auto model{ BuildFiniteNonphysicalProfileDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    const auto previous_model{ GetEstimateModel(*atom) };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    const auto fitted_model{ GetEstimateModel(*atom) };
    EXPECT_GT(fitted_model.GetAmplitude(), 0.0);
    EXPECT_GT(fitted_model.GetWidth(), 0.0);
    ExpectGaussianModelsNear(fitted_model, previous_model, 1.0e-12);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingAppliesCollinearRidgeGuard)
{
    auto model{ BuildNearCollinearDefenseModel() };
    const auto initial_error{ CalculateSelectedAtomResponseMeanSquaredError(*model) };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    const auto fitted_error{ CalculateSelectedAtomResponseMeanSquaredError(*model) };
    const auto tolerance{ 1.0e-3 * std::max(initial_error, 1.0) };
    EXPECT_LE(fitted_error, initial_error + tolerance);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(
    EstimatorSecondStageDefenseTest,
    RunSecondStageLocalFittingPersistsFinalModelAndPeelingWithoutGroupFitting)
{
    auto model{ BuildIndependentOffsetDefenseModel() };
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    const auto options{ MakeSecondStageOptions() };
    const auto initial_analysis_view{ model->GetAnalysisView() };
    const auto group_key_list{ initial_analysis_view.CollectAtomGroupKeys(
        rg::FittingStage::Second) };
    std::vector<rg::GaussianModel3D> initial_group_mean_list;
    std::vector<rg::GaussianModel3D> initial_group_mdpde_list;
    std::vector<rg::GaussianModel3D> initial_group_prior_list;
    initial_group_mean_list.reserve(group_key_list.size());
    initial_group_mdpde_list.reserve(group_key_list.size());
    initial_group_prior_list.reserve(group_key_list.size());
    for (const auto group_key : group_key_list)
    {
        initial_group_mean_list.emplace_back(
            initial_analysis_view.GetAtomGroupMean(
                rg::FittingStage::Second, group_key));
        initial_group_mdpde_list.emplace_back(
            initial_analysis_view.GetAtomGroupMDPDE(
                rg::FittingStage::Second, group_key));
        initial_group_prior_list.emplace_back(
            initial_analysis_view.GetAtomGroupPrior(
                rg::FittingStage::Second, group_key));
    }

    rt::RunSecondStageLocalFitting(*model, options);

    ExpectPeelingSamplingEntriesMatchFinalModels(*model);
    const auto analysis_view{ model->GetAnalysisView() };
    ASSERT_EQ(
        analysis_view.CollectAtomGroupKeys(rg::FittingStage::Second),
        group_key_list);
    for (std::size_t group_index = 0;
        group_index < group_key_list.size();
        group_index++)
    {
        const auto group_key{ group_key_list.at(group_index) };
        const auto & atom_list{
            analysis_view.GetAtomObjectList(
                rg::FittingStage::Second, group_key)
        };
        EXPECT_DOUBLE_EQ(
            analysis_view.GetAtomAlphaG(
                rg::FittingStage::Second, group_key),
            0.0);
        ExpectGaussianModelsNear(
            analysis_view.GetAtomGroupMean(
                rg::FittingStage::Second, group_key),
            initial_group_mean_list.at(group_index),
            0.0);
        ExpectGaussianModelsNear(
            analysis_view.GetAtomGroupMDPDE(
                rg::FittingStage::Second, group_key),
            initial_group_mdpde_list.at(group_index),
            0.0);
        ExpectGaussianModelsNear(
            analysis_view.GetAtomGroupPrior(
                rg::FittingStage::Second, group_key),
            initial_group_prior_list.at(group_index),
            0.0);
        for (const auto * atom : atom_list)
        {
            EXPECT_FALSE(
                rg::AtomLocalPotentialView::For(*atom)
                    .GetGaussianResult(FittingStage::Second)
                    .posterior.has_value());
        }
    }
}

TEST(
    EstimatorSecondStageDefenseTest,
    RunSecondStageLocalFittingUsesFrozenClusterBackgroundWithoutGroupOrResidueKeys)
{
    const std::array seeds{ rg::GaussianModel3D{ 5.0, 0.50, 0.05 },
        rg::GaussianModel3D{ 7.0, 0.60, 0.15 } };
    const std::array truth{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 }, rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 5.0, 0.48, 0.10 }, rg::GaussianModel3D{ 7.0, 0.63, 0.10 },
        rg::GaussianModel3D{ 2.5, 0.45, 0.05 }, rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 } };
    auto options{ MakeSecondStageOptions() };
    options.exclude_hydrogen = true;
    for (const auto & [shared_cluster, shared_contributor] : {
        std::pair{ false, false }, std::pair{ true, false }, std::pair{ false, true } })
    {
        std::array<rg::GaussianModel3D, 2> reference_selected;
        std::array<rg::GaussianModel3D, 2> reference_background;
        for (const double scale : { 1.0, 100.0 })
        {
            auto scaled_seeds{ seeds };
            auto scaled_truth{ truth };
            for (auto & model : scaled_seeds)
                model = { scale * model.GetAmplitude(), model.GetWidth(), scale * model.GetOffset() };
            for (auto & model : scaled_truth)
                model = { scale * model.GetAmplitude(), model.GetWidth(), scale * model.GetOffset() };
            auto serial{ BuildUnselectedContributorDefenseModel(
                scaled_seeds, scaled_truth, false, shared_cluster, shared_contributor) };
            auto parallel{ BuildUnselectedContributorDefenseModel(
                scaled_seeds, scaled_truth, true, shared_cluster, shared_contributor) };
            auto logged{ BuildUnselectedContributorDefenseModel(
                scaled_seeds, scaled_truth, false, shared_cluster, shared_contributor) };
            EXPECT_EQ(rg::data_internal::GetGroupKey(serial->FindAtomPtr(3)),
                rg::data_internal::GetGroupKey(serial->FindAtomPtr(4)));
            EXPECT_NE(rg::data_internal::GetGroupKey(serial->FindAtomPtr(3)),
                rg::data_internal::GetGroupKey(parallel->FindAtomPtr(3)));
            options.thread_size = 1;
            ASSERT_TRUE(rt::RunSecondStageLocalFitting(*serial, options));
            options.thread_size = 2;
            ASSERT_TRUE(rt::RunSecondStageLocalFitting(*parallel, options));
            const auto previous_level{ Logger::GetLogLevel() };
            Logger::SetLogLevel(LogLevel::Debug);
            testing::internal::CaptureStdout();
            options.thread_size = 1;
            options.quiet_mode = false;
            const bool completed{ rt::RunSecondStageLocalFitting(*logged, options) };
            const auto output{ testing::internal::GetCapturedStdout() };
            auto alternate_logged{ BuildUnselectedContributorDefenseModel(
                scaled_seeds, scaled_truth, true, shared_cluster, shared_contributor) };
            testing::internal::CaptureStdout();
            const bool alternate_completed{
                rt::RunSecondStageLocalFitting(*alternate_logged, options) };
            const auto alternate_output{ testing::internal::GetCapturedStdout() };
            auto relabeled_logged{ BuildUnselectedContributorDefenseModel(
                scaled_seeds, scaled_truth, false, shared_cluster, shared_contributor, true) };
            testing::internal::CaptureStdout();
            const bool relabeled_completed{
                rt::RunSecondStageLocalFitting(*relabeled_logged, options) };
            const auto relabeled_output{ testing::internal::GetCapturedStdout() };
            Logger::SetLogLevel(previous_level);
            options.quiet_mode = true;
            ASSERT_TRUE(completed);
            ASSERT_TRUE(alternate_completed);
            ASSERT_TRUE(relabeled_completed);
            const auto audit_records = [](const std::string & log)
            {
                std::vector<std::string> records;
                std::istringstream lines{ log };
                std::string line;
                while (std::getline(lines, line))
                {
                    for (const std::string marker : {
                        "Convergence safeguard audit:", "Second-stage audit terminal:",
                        "Second-stage audit terminal atom:", "Second-stage local fitting summary:",
                        "Local-fitting atom cutoff:", "Adaptive local-fitting topology rebuild:" })
                    {
                        const auto position{ line.find(marker) };
                        if (position != std::string::npos)
                            records.emplace_back(line.substr(position));
                    }
                }
                return records;
            };
            EXPECT_EQ(audit_records(output), audit_records(alternate_output));
            EXPECT_EQ(audit_records(output), audit_records(relabeled_output));
            EXPECT_NE(output.find(shared_cluster ?
                "initial components/max atoms/ratio = 1/2/1.00" :
                "initial components/max atoms/ratio = 2/1/0.50"), std::string::npos);
            if (shared_contributor)
            {
                EXPECT_NE(output.find("candidate/retained/cut edges = 0/0/0"), std::string::npos);
            }
            std::array<std::optional<rg::GaussianModel3D>, 2> first_background;
            std::array<std::optional<rg::GaussianModel3D>, 2> last_background;
            const auto terminal_position{ output.find("Second-stage audit terminal:") };
            ASSERT_NE(terminal_position, std::string::npos);
            const auto attempt_count{
                std::stoull(output.substr(output.find(", try=", terminal_position) + 6)) };
            EXPECT_GT(attempt_count, 1U);
            const auto first_iteration_end{ output.find("Convergence safeguard audit:") };
            ASSERT_NE(first_iteration_end, std::string::npos);
            std::array<std::size_t, 2> background_count{};
            std::array<std::size_t, 2> first_iteration_background_count{};
            const std::string marker{ "Second-stage frozen background: target=" };
            std::size_t position{ 0 };
            while ((position = output.find(marker, position)) != std::string::npos)
            {
                position += marker.size();
                const auto target{ static_cast<std::size_t>(std::stoi(output.substr(position)) - 1) };
                ASSERT_LT(target, 2U);
                background_count.at(target)++;
                if (position < first_iteration_end) first_iteration_background_count.at(target)++;
                const auto parameter_begin{ output.find("A/W/C=", position) + 6 };
                auto parameters{ output.substr(parameter_begin, output.find(",", parameter_begin) - parameter_begin) };
                std::replace(parameters.begin(), parameters.end(), '/', ' ');
                std::istringstream values{ parameters };
                double amplitude{ 0.0 }, width{ 0.0 }, offset{ 0.0 };
                ASSERT_TRUE(static_cast<bool>(values >> amplitude >> width >> offset));
                last_background.at(target) = rg::GaussianModel3D{ amplitude, width, offset };
                if (!first_background.at(target).has_value()) first_background.at(target) = last_background.at(target);
            }
            for (std::size_t target = 0; target < 2; target++)
            {
                ASSERT_TRUE(first_background.at(target).has_value());
                ASSERT_TRUE(last_background.at(target).has_value());
                EXPECT_EQ(first_iteration_background_count.at(target), 1U);
                EXPECT_EQ(background_count.at(target), attempt_count);
                const auto expected_initial{ shared_cluster ?
                    *median_detail::BuildGaussianParameterMedian({ scaled_seeds[0], scaled_seeds[1] }) : scaled_seeds[target] };
                ExpectGaussianModelsNear(*first_background.at(target), expected_initial, 1.0e-12 * scale);
                const int serial_id{ static_cast<int>(target + 1) };
                const auto selected{ GetEstimateModel(*serial->FindAtomPtr(serial_id)) };
                ExpectGaussianModelsNear(selected, GetEstimateModel(*parallel->FindAtomPtr(serial_id)), 1.0e-10 * scale);
                ExpectGaussianModelsNear(selected, GetEstimateModel(*logged->FindAtomPtr(serial_id)), 1.0e-10 * scale);
                ExpectGaussianModelsNear(selected, GetEstimateModel(*alternate_logged->FindAtomPtr(serial_id)), 0.0);
                ExpectGaussianModelsNear(selected, GetEstimateModel(*relabeled_logged->FindAtomPtr(serial_id)), 0.0);
                const auto view{ rg::AtomLocalPotentialView::For(*serial->FindAtomPtr(serial_id)) };
                const auto raw{ view.GetRawSamplingEntries(false) };
                const auto peeled{ view.GetPeelingSamplingEntries(false) };
                const auto alternate{ rg::AtomLocalPotentialView::For(*parallel->FindAtomPtr(serial_id)).GetPeelingSamplingEntries(false) };
                const auto relabeled{ rg::AtomLocalPotentialView::For(*relabeled_logged->FindAtomPtr(serial_id)).GetPeelingSamplingEntries(false) };
                ASSERT_EQ(raw.size(), peeled.size());
                ASSERT_EQ(alternate.size(), peeled.size());
                ASSERT_EQ(relabeled.size(), peeled.size());
                for (std::size_t row = 0; row < raw.size(); row++)
                {
                    double background{ 0.0 };
                    for (const int contributor_id : { 3, 4, 6, 7 })
                    {
                        const auto distance{ Distance(raw.at(row).point.position,
                            serial->FindAtomPtr(contributor_id)->GetPosition()) };
                        if (distance <= 2.5)
                            background += last_background.at(target)->ResponseAtDistance(distance);
                    }
                    const auto * selected_neighbor{ serial->FindAtomPtr(target == 0 ? 2 : 1) };
                    const auto selected_distance{ Distance(raw.at(row).point.position, selected_neighbor->GetPosition()) };
                    if (selected_distance <= 2.5)
                        background += GetEstimateModel(*selected_neighbor).ResponseAtDistance(selected_distance);
                    EXPECT_NEAR(peeled.at(row).response, raw.at(row).response - background, 1.0e-10 * scale);
                    EXPECT_NEAR(peeled.at(row).response, alternate.at(row).response, 1.0e-10 * scale);
                    EXPECT_DOUBLE_EQ(peeled.at(row).response, relabeled.at(row).response);
                }
                EXPECT_EQ(view.GetNeighborCountForPeeling(), shared_cluster ? 3 :
                    shared_contributor && target == 1 ? 2 : 1);
                const rg::GaussianModel3D normalized_selected{
                    selected.GetAmplitude() / scale, selected.GetWidth(), selected.GetOffset() / scale };
                const rg::GaussianModel3D normalized_background{
                    last_background.at(target)->GetAmplitude() / scale, last_background.at(target)->GetWidth(),
                    last_background.at(target)->GetOffset() / scale };
                if (scale == 1.0)
                {
                    reference_selected.at(target) = normalized_selected;
                    reference_background.at(target) = normalized_background;
                }
                else
                {
                    ExpectGaussianModelsNear(normalized_selected, reference_selected.at(target), 1.0e-6);
                    ExpectGaussianModelsNear(normalized_background, reference_background.at(target), 1.0e-6);
                }
            }
            if (shared_cluster)
                ExpectGaussianModelsNear(*last_background[0], *last_background[1], 0.0);
            EXPECT_EQ(serial->GetSelectedAtomCount(), 2U);
            for (int id = 3; id <= 7; id++)
            {
                EXPECT_FALSE(rg::AtomLocalPotentialView::For(*serial->FindAtomPtr(id)).IsAvailable());
                EXPECT_FALSE(rg::AtomLocalPotentialView::For(*parallel->FindAtomPtr(id)).IsAvailable());
            }
        }
    }
    auto excluded{ BuildUnselectedContributorDefenseModel(seeds, truth) };
    auto included{ BuildUnselectedContributorDefenseModel(seeds, truth) };
    ASSERT_TRUE(rt::RunSecondStageLocalFitting(*excluded, options));
    options.exclude_hydrogen = false;
    ASSERT_TRUE(rt::RunSecondStageLocalFitting(*included, options));
    EXPECT_EQ(rg::AtomLocalPotentialView::For(*included->FindAtomPtr(1)).GetNeighborCountForPeeling(), 2);
    EXPECT_NE(rg::AtomLocalPotentialView::For(*excluded->FindAtomPtr(1)).GetPeelingSamplingEntries(false).front().response,
        rg::AtomLocalPotentialView::For(*included->FindAtomPtr(1)).GetPeelingSamplingEntries(false).front().response);

}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingJointlyPolishesClusterParameters)
{
    auto model{ BuildJointPolishDefenseModel() };
    const auto & atom_list{ model->GetSelectedAtoms() };
    std::vector<rg::GaussianModel3D> initial_model_list;
    for (const auto * atom : atom_list)
    {
        initial_model_list.emplace_back(GetEstimateModel(*atom));
    }
    const auto initial_error{ CalculateSelectedAtomResponseMeanSquaredError(*model) };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    const auto fitted_error{ CalculateSelectedAtomResponseMeanSquaredError(*model) };
    EXPECT_LT(fitted_error, initial_error);
    bool amplitude_changed{ false };
    bool width_changed{ false };
    bool offset_changed{ false };
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        const auto fitted{ GetEstimateModel(*atom_list.at(i)) };
        const auto & initial{ initial_model_list.at(i) };
        amplitude_changed = amplitude_changed ||
            std::abs(fitted.GetAmplitude() - initial.GetAmplitude()) > 1.0e-6;
        width_changed = width_changed ||
            std::abs(fitted.GetWidth() - initial.GetWidth()) > 1.0e-6;
        offset_changed = offset_changed ||
            std::abs(fitted.GetOffset() - initial.GetOffset()) > 1.0e-6;
    }
    EXPECT_TRUE(amplitude_changed);
    EXPECT_TRUE(width_changed);
    EXPECT_TRUE(offset_changed);
    EXPECT_GT(
        std::abs(
            GetEstimateModel(*atom_list.at(0)).GetOffset() -
            GetEstimateModel(*atom_list.at(1)).GetOffset()),
        1.0e-6);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, SameChemicalKeyAtomsKeepIndependentOffsetsAndGroupKeyInvariantResults)
{
    auto original{ BuildIndependentOffsetDefenseModel() };
    auto relabeled{ BuildIndependentOffsetDefenseModel(1.0, true) };
    const auto initial_error{ CalculateSelectedAtomResponseMeanSquaredError(*original) };
    const auto run_logged = [](rg::ModelObject & model)
    {
        auto options{ MakeSecondStageOptions() };
        options.quiet_mode = false;
        const auto previous_level{ Logger::GetLogLevel() };
        Logger::SetLogLevel(LogLevel::Debug);
        testing::internal::CaptureStdout();
        const bool completed{ rt::RunSecondStageLocalFitting(model, options) };
        const auto output{ testing::internal::GetCapturedStdout() };
        Logger::SetLogLevel(previous_level);
        EXPECT_TRUE(completed);
        std::vector<std::string> evidence;
        std::istringstream lines{ output };
        std::string line;
        while (std::getline(lines, line))
        {
            for (const std::string marker : {
                "Convergence safeguard audit:", "Second-stage audit terminal:",
                "Second-stage audit terminal atom:" })
            {
                const auto position{ line.find(marker) };
                if (position != std::string::npos) evidence.emplace_back(line.substr(position));
            }
        }
        EXPECT_FALSE(evidence.empty());
        return evidence;
    };
    const auto original_trace{ run_logged(*original) };
    EXPECT_EQ(original_trace, run_logged(*relabeled));
    const auto & atoms{ original->GetSelectedAtoms() };
    const auto & changed{ relabeled->GetSelectedAtoms() };
    ASSERT_EQ(atoms.size(), 2U);
    EXPECT_EQ(rg::data_internal::GetGroupKey(atoms.at(0)), rg::data_internal::GetGroupKey(atoms.at(1)));
    EXPECT_NE(rg::data_internal::GetGroupKey(changed.at(0)), rg::data_internal::GetGroupKey(changed.at(1)));
    EXPECT_GT(std::abs(GetEstimateModel(*atoms.at(0)).GetOffset() -
        GetEstimateModel(*atoms.at(1)).GetOffset()), 1.0e-5);
    for (std::size_t atom = 0; atom < atoms.size(); atom++)
    {
        ExpectGaussianModelsNear(GetEstimateModel(*atoms.at(atom)), GetEstimateModel(*changed.at(atom)), 1.0e-12);
        const auto peeled{ rg::AtomLocalPotentialView::For(*atoms.at(atom)).GetPeelingSamplingEntries(false) };
        const auto other{ rg::AtomLocalPotentialView::For(*changed.at(atom)).GetPeelingSamplingEntries(false) };
        ASSERT_EQ(peeled.size(), other.size());
        for (std::size_t row = 0; row < peeled.size(); row++)
            EXPECT_DOUBLE_EQ(peeled.at(row).response, other.at(row).response);
    }
    ExpectPeelingSamplingEntriesMatchFinalModels(*original);
    EXPECT_LT(CalculateSelectedAtomResponseMeanSquaredError(*original), initial_error);
    ExpectSelectedAtomEstimatesAreFinite(*original);
}

TEST(EstimatorSecondStageDefenseTest, IndependentOffsetJointPolishIsIntensityScaleInvariant)
{
    constexpr double scale{ 100.0 };
    auto base_model{ BuildIndependentOffsetDefenseModel() };
    auto scaled_model{ BuildIndependentOffsetDefenseModel(scale) };

    rt::RunSecondStageLocalFitting(*base_model, MakeSecondStageOptions());
    rt::RunSecondStageLocalFitting(*scaled_model, MakeSecondStageOptions());

    const auto & base_atoms{ base_model->GetSelectedAtoms() };
    const auto & scaled_atoms{ scaled_model->GetSelectedAtoms() };
    ASSERT_EQ(base_atoms.size(), scaled_atoms.size());
    ASSERT_EQ(base_atoms.size(), 2U);
    EXPECT_NE(
        GetEstimateModel(*base_atoms.at(0)).GetOffset(),
        GetEstimateModel(*base_atoms.at(1)).GetOffset());
    EXPECT_NE(
        GetEstimateModel(*scaled_atoms.at(0)).GetOffset(),
        GetEstimateModel(*scaled_atoms.at(1)).GetOffset());
    for (std::size_t i = 0; i < base_atoms.size(); i++)
    {
        const auto base{ GetEstimateModel(*base_atoms.at(i)) };
        const auto scaled{ GetEstimateModel(*scaled_atoms.at(i)) };
        EXPECT_NEAR(
            base.GetAmplitude() * scale,
            scaled.GetAmplitude(),
            std::max(1.0e-8, std::abs(scaled.GetAmplitude()) * 5.0e-5));
        EXPECT_NEAR(
            base.GetWidth(),
            scaled.GetWidth(),
            std::max(1.0e-8, std::abs(scaled.GetWidth()) * 5.0e-5));
        EXPECT_NEAR(
            base.GetOffset() * scale,
            scaled.GetOffset(),
            std::max(1.0e-8, std::abs(scaled.GetOffset()) * 5.0e-5));
    }
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingUpdatesHealthyVariablesAcrossClusters)
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

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

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

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingMatchesSerialAndParallelSelection)
{
    auto serial_model{ BuildSeparatedRollbackDefenseModel() };
    auto parallel_model{ BuildSeparatedRollbackDefenseModel() };
    auto serial_options{ MakeSecondStageOptions() };
    auto parallel_options{ MakeSecondStageOptions() };
    serial_options.thread_size = 1;
    parallel_options.thread_size = 2;

    EXPECT_EQ(
        rt::RunSecondStageLocalFitting(*serial_model, serial_options),
        rt::RunSecondStageLocalFitting(*parallel_model, parallel_options));

    const auto & serial_atoms{ serial_model->GetSelectedAtoms() };
    const auto & parallel_atoms{ parallel_model->GetSelectedAtoms() };
    ASSERT_EQ(serial_atoms.size(), parallel_atoms.size());
    for (std::size_t i = 0; i < serial_atoms.size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*serial_atoms.at(i)),
            GetEstimateModel(*parallel_atoms.at(i)),
            1.0e-12);
    }
}

TEST(
    EstimatorSecondStageDefenseTest,
    BoundaryComponentReconciliationMatchesSerialAndParallelSelection)
{
    auto serial_model{ BuildBoundaryComponentConflictDefenseModel() };
    auto parallel_model{ BuildBoundaryComponentConflictDefenseModel() };
    auto serial_options{ MakeSecondStageOptions() };
    auto parallel_options{ MakeSecondStageOptions() };
    serial_options.thread_size = 1;
    parallel_options.thread_size = 2;

    const auto previous_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    serial_options.quiet_mode = false;
    testing::internal::CaptureStdout();
    const auto serial_completed{ rt::RunSecondStageLocalFitting(*serial_model, serial_options) };
    const auto output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_level);
    EXPECT_EQ(serial_completed, rt::RunSecondStageLocalFitting(*parallel_model, parallel_options));
    const auto cutoff_position{ output.find("Local-fitting atom cutoff: atoms=103, limit=100, clusters=") };
    ASSERT_NE(cutoff_position, std::string::npos);
    const auto maximum_position{ output.find(", max-atoms=", cutoff_position) };
    const auto cuts_position{ output.find(", cutoff-edges=", cutoff_position) };
    ASSERT_NE(maximum_position, std::string::npos);
    ASSERT_NE(cuts_position, std::string::npos);
    EXPECT_LE(std::stoull(output.substr(maximum_position + 12)), 100U);
    EXPECT_GT(std::stoull(output.substr(cuts_position + 15)), 0U);
    EXPECT_NE(output.find("Boundary-component reconciliation: clusters/atoms/boundary-samples = 2/101/"),
        std::string::npos);
    const auto & serial_atoms{ serial_model->GetSelectedAtoms() };
    const auto & parallel_atoms{ parallel_model->GetSelectedAtoms() };
    ASSERT_EQ(serial_atoms.size(), parallel_atoms.size());
    for (std::size_t i = 0; i < serial_atoms.size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*serial_atoms.at(i)),
            GetEstimateModel(*parallel_atoms.at(i)),
            1.0e-12);
    }
}

TEST(
    EstimatorSecondStageDefenseTest,
    BoundaryComponentReconciliationIsIntensityScaleInvariant)
{
    constexpr double intensity_scale{ 100.0 };
    auto base_model{ BuildBoundaryComponentConflictDefenseModel() };
    auto scaled_model{
        BuildBoundaryComponentConflictDefenseModel(intensity_scale)
    };

    rt::RunSecondStageLocalFitting(*base_model, MakeSecondStageOptions());
    rt::RunSecondStageLocalFitting(*scaled_model, MakeSecondStageOptions());
    const auto & base_atoms{ base_model->GetSelectedAtoms() };
    const auto & scaled_atoms{ scaled_model->GetSelectedAtoms() };
    ASSERT_EQ(base_atoms.size(), scaled_atoms.size());
    for (std::size_t i = 0; i < base_atoms.size(); i++)
    {
        const auto base{ GetEstimateModel(*base_atoms.at(i)) };
        const auto scaled{ GetEstimateModel(*scaled_atoms.at(i)) };
        EXPECT_NEAR(
            base.GetAmplitude() * intensity_scale,
            scaled.GetAmplitude(),
            std::max(1.0e-8, std::abs(scaled.GetAmplitude()) * 5.0e-5));
        EXPECT_NEAR(base.GetWidth(), scaled.GetWidth(), 1.0e-5);
        EXPECT_NEAR(
            base.GetOffset() * intensity_scale,
            scaled.GetOffset(),
            std::max(1.0e-8, std::abs(scaled.GetOffset()) * 1.0e-4));
    }
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingIsIntensityScaleInvariant)
{
    constexpr double scale{ 100.0 };
    auto base_model{ BuildNearCollinearDefenseModel() };
    auto scaled_model{ BuildNearCollinearDefenseModel(scale) };

    rt::RunSecondStageLocalFitting(*base_model, MakeSecondStageOptions());
    rt::RunSecondStageLocalFitting(*scaled_model, MakeSecondStageOptions());

    const auto & base_atoms{ base_model->GetSelectedAtoms() };
    const auto & scaled_atoms{ scaled_model->GetSelectedAtoms() };
    ASSERT_EQ(base_atoms.size(), scaled_atoms.size());
    for (std::size_t i = 0; i < base_atoms.size(); i++)
    {
        const auto base{ GetEstimateModel(*base_atoms.at(i)) };
        const auto scaled{ GetEstimateModel(*scaled_atoms.at(i)) };
        EXPECT_NEAR(
            base.GetAmplitude() * scale,
            scaled.GetAmplitude(),
            std::max(1.0e-8, std::abs(scaled.GetAmplitude()) * 1.0e-5));
        EXPECT_NEAR(base.GetWidth(), scaled.GetWidth(), 1.0e-6);
        EXPECT_NEAR(
            base.GetOffset() * scale,
            scaled.GetOffset(),
            std::max(1.0e-8, std::abs(scaled.GetOffset()) * 5.0e-5));
    }
}

TEST(
    EstimatorSecondStageDefenseTest,
    BoundaryJointCorrectionMatchesSerialParallelAndIntensityScaling)
{
    constexpr double intensity_scale{ 100.0 };
    auto serial_model{ BuildBoundaryJointCorrectionDefenseModel() };
    auto parallel_model{ BuildBoundaryJointCorrectionDefenseModel() };
    auto scaled_model{
        BuildBoundaryJointCorrectionDefenseModel(intensity_scale)
    };
    auto serial_options{ MakeSecondStageOptions() };
    auto parallel_options{ MakeSecondStageOptions() };
    serial_options.thread_size = 1;
    parallel_options.thread_size = 2;

    EXPECT_EQ(
        rt::RunSecondStageLocalFitting(*serial_model, serial_options),
        rt::RunSecondStageLocalFitting(*parallel_model, parallel_options));
    rt::RunSecondStageLocalFitting(*scaled_model, MakeSecondStageOptions());

    const auto & serial_atoms{ serial_model->GetSelectedAtoms() };
    const auto & parallel_atoms{ parallel_model->GetSelectedAtoms() };
    const auto & scaled_atoms{ scaled_model->GetSelectedAtoms() };
    ASSERT_EQ(serial_atoms.size(), parallel_atoms.size());
    ASSERT_EQ(serial_atoms.size(), scaled_atoms.size());
    for (std::size_t i = 0; i < serial_atoms.size(); i++)
    {
        const auto serial{ GetEstimateModel(*serial_atoms.at(i)) };
        const auto parallel{ GetEstimateModel(*parallel_atoms.at(i)) };
        const auto scaled{ GetEstimateModel(*scaled_atoms.at(i)) };
        ExpectGaussianModelsNear(serial, parallel, 1.0e-12);
        EXPECT_NEAR(
            serial.GetAmplitude() * intensity_scale,
            scaled.GetAmplitude(),
            std::max(1.0e-8, std::abs(scaled.GetAmplitude()) * 5.0e-5));
        EXPECT_NEAR(serial.GetWidth(), scaled.GetWidth(), 5.0e-6);
        EXPECT_NEAR(
            serial.GetOffset() * intensity_scale,
            scaled.GetOffset(),
            std::max(1.0e-8, std::abs(scaled.GetOffset()) * 1.0e-4));
    }
}

TEST(
    EstimatorSecondStageDefenseTest,
    BoundaryComponentReconciliationBacktracksAndPreservesRemoteCluster)
{
    auto model{ BuildBoundaryComponentConflictDefenseModel() };
    const auto initial_remote_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 101, 103)
    };
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());
    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 101, 103),
        initial_remote_error);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(
    EstimatorSecondStageDefenseTest,
    RunSecondStageLocalFittingDampsOffsetStepIntoInitialTrustRadius)
{
    const rg::GaussianModel3D initial_model{ 6.0, 0.55, 0.0 };
    auto truth_coordinates{ initial_model.ToTransformedCoordinates() };
    ASSERT_TRUE(truth_coordinates.has_value());
    (*truth_coordinates)(static_cast<Eigen::Index>(
        rg::GaussianModel3D::OffsetToPeakRatioCoordinateIndex())) = 1.25;
    const auto truth_model{
        rg::GaussianModel3D::FromTransformedCoordinates(
            *truth_coordinates)
    };
    ASSERT_TRUE(truth_model.has_value());

    auto model{
        BuildDefenseModel(
            { std::array<double, 3>{ 0.0, 0.0, 0.0 } },
            { Spot::O },
            { Element::OXYGEN },
            { *truth_model },
            initial_model)
    };
    const auto previous_model{
        GetEstimateModel(*model->GetSelectedAtoms().front())
    };
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    const auto fitted_model{
        GetEstimateModel(*model->GetSelectedAtoms().front())
    };
    EXPECT_NE(fitted_model.GetOffset(), previous_model.GetOffset());
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(
    EstimatorSecondStageDefenseTest,
    MissingGroupSeedsUsesLocalMdpdeAndGlobalFallbackWithoutGroupFitting)
{
    auto model{ BuildNearCollinearDefenseModel() };
    auto options{ MakeSecondStageOptions() };
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    for (auto * atom : model->GetSelectedAtoms())
    {
        auto result{
            rg::AtomLocalPotentialView::For(*atom).GetGaussianResult(
                FittingStage::Second)
        };
        ASSERT_TRUE(seed_detail::IsValidSecondStageGaussianModel(
            result.mdpde.GetModel()));
        ASSERT_TRUE(seed_detail::IsValidSecondStageGaussianModel(
            result.ols.GetModel()));
        if (atom == model->GetSelectedAtoms().front())
        {
            result.mdpde = rg::GaussianModel3DWithUncertainty{
                rg::GaussianModel3D{ 0.0, 0.0, 0.0 },
                rg::GaussianModel3DUncertainty{}
            };
        }
        result.posterior.reset();
        analysis.SetAtomLocalGaussianResult(
            FittingStage::Second,
            *atom,
            std::move(result));
        LocalPotentialSampleList sentinel_peeling_sampling_entries{
            LocalPotentialSample{
                100.0 + static_cast<double>(atom->GetSerialID()),
                SamplingPoint{ 0.5, atom->GetPosition(), true }
            }
        };
        analysis.SetAtomLocalPeelingSamplingEntries(
            *atom, sentinel_peeling_sampling_entries);
        analysis.SetAtomLocalNeighborCountForPeeling(*atom, 99);
    }
    const auto previous_analysis_view{ model->GetAnalysisView() };
    std::vector<rg::GaussianModel3D> previous_group_prior_list;
    for (const auto group_key : previous_analysis_view.CollectAtomGroupKeys(
        rg::FittingStage::Second))
    {
        previous_group_prior_list.emplace_back(
            previous_analysis_view.GetAtomGroupPrior(
                rg::FittingStage::Second, group_key));
    }

    const auto peeling_applied{ rt::RunSecondStageLocalFitting(*model, options) };

    EXPECT_TRUE(peeling_applied);

    for (std::size_t i = 0; i < model->GetSelectedAtoms().size(); i++)
    {
        const auto * atom{ model->GetSelectedAtoms().at(i) };
        EXPECT_TRUE(seed_detail::IsValidSecondStageGaussianModel(
            GetEstimateModel(*atom)));
        const auto peeling_sampling_entries{
            rg::AtomLocalPotentialView::For(
                *atom)
                .GetPeelingSamplingEntries(false)
        };
        EXPECT_GT(peeling_sampling_entries.size(), 1U);
        EXPECT_NE(
            rg::AtomLocalPotentialView::For(
                *atom)
                .GetNeighborCountForPeeling(),
            99);
        EXPECT_FALSE(
            rg::AtomLocalPotentialView::For(*atom)
                .GetGaussianResult(FittingStage::Second)
                .posterior.has_value());
    }
    const auto final_analysis_view{ model->GetAnalysisView() };
    const auto group_key_list{ final_analysis_view.CollectAtomGroupKeys(
        rg::FittingStage::Second) };
    ASSERT_EQ(group_key_list.size(), previous_group_prior_list.size());
    for (std::size_t i = 0; i < group_key_list.size(); i++)
    {
        ExpectGaussianModelsNear(
            final_analysis_view.GetAtomGroupPrior(
                rg::FittingStage::Second, group_key_list.at(i)),
            previous_group_prior_list.at(i),
            0.0);
    }
}
