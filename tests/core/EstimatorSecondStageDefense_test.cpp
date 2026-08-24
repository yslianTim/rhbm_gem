#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "core/detail/FittingModel.hpp"
#include "core/detail/CouplingGraph.hpp"
#include "core/detail/JointFitting.hpp"
#include "core/detail/CandidateSelection.hpp"
#include "core/detail/IterationProcess.hpp"
#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

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
    const std::array<float, 3> & lhs,
    const std::array<float, 3> & rhs)
{
    const auto dx{ static_cast<double>(lhs.at(0) - rhs.at(0)) };
    const auto dy{ static_cast<double>(lhs.at(1) - rhs.at(1)) };
    const auto dz{ static_cast<double>(lhs.at(2) - rhs.at(2)) };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::vector<coupling_detail::ResidueKey>
MakeUniqueResidueKeys(std::size_t count)
{
    std::vector<coupling_detail::ResidueKey> key_list;
    key_list.reserve(count);
    for (std::size_t index = 0; index < count; index++)
    {
        key_list.emplace_back("A", static_cast<int>(index));
    }
    return key_list;
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

std::optional<std::array<std::size_t, 4>> ParsePolishProgressCounts(
    std::string_view row)
{
    std::array<std::size_t, 4> separator_position{};
    std::size_t search_position{ 0 };
    for (auto & position : separator_position)
    {
        position = row.find('|', search_position);
        if (position == std::string_view::npos) return std::nullopt;
        search_position = position + 1;
    }

    auto cell{
        row.substr(
            separator_position.at(2) + 1,
            separator_position.at(3) - separator_position.at(2) - 1)
    };
    while (!cell.empty() && cell.front() == ' ') cell.remove_prefix(1);
    while (!cell.empty() && cell.back() == ' ') cell.remove_suffix(1);

    std::array<std::size_t, 4> count{};
    std::size_t value_start{ 0 };
    for (std::size_t value_index = 0; value_index < count.size(); value_index++)
    {
        const auto value_end{
            value_index + 1 == count.size() ?
                cell.size() : cell.find('/', value_start)
        };
        if (value_end == std::string_view::npos) return std::nullopt;
        count.at(value_index) = static_cast<std::size_t>(
            std::stoull(std::string{
                cell.substr(value_start, value_end - value_start) }));
        value_start = value_end + 1;
    }
    return count;
}

std::unique_ptr<rg::AtomObject> MakeAtom(
    int serial_id,
    Spot spot,
    Element element,
    const std::array<float, 3> & position)
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
    result.posterior = result.mdpde;
    return result;
}

LocalPotentialSampleList BuildSamples(
    const rg::AtomObject & target_atom,
    const std::vector<rg::AtomObject *> & atom_list,
    const std::vector<rg::GaussianModel3D> & truth_model_list)
{
    const std::array<std::array<float, 3>, 6> direction_list{
        std::array<float, 3>{ 1.0F, 0.0F, 0.0F },
        std::array<float, 3>{ -1.0F, 0.0F, 0.0F },
        std::array<float, 3>{ 0.0F, 1.0F, 0.0F },
        std::array<float, 3>{ 0.0F, -1.0F, 0.0F },
        std::array<float, 3>{ 0.0F, 0.0F, 1.0F },
        std::array<float, 3>{ 0.0F, 0.0F, -1.0F }
    };
    const std::array<float, 4> radius_list{ 0.15F, 0.35F, 0.65F, 0.95F };

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
            sample_list.emplace_back(LocalPotentialSample{
                static_cast<float>(response),
                point
            });
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
    const std::vector<std::size_t> & group_id_list,
    const std::vector<rg::GaussianModel3D> & model_list,
    const std::vector<double> & target_offset_list)
{
    if (group_id_list.size() != model_list.size() ||
        model_list.size() != target_offset_list.size())
    {
        throw std::invalid_argument(
            "Joint offset fixture input sizes are inconsistent.");
    }

    offset_detail::SecondStageContext context;
    context.selected_atom_list.resize(model_list.size());
    std::unordered_map<std::size_t, std::size_t> dense_group_id_by_id;
    for (std::size_t atom_index = 0;
        atom_index < model_list.size();
        atom_index++)
    {
        auto & atom_context{ context.at(atom_index) };
        const auto group_iter{
            dense_group_id_by_id.try_emplace(
                group_id_list.at(atom_index),
                dense_group_id_by_id.size()).first
        };
        atom_context.group_id = group_iter->second;
        SamplingPoint point;
        point.distance = 0.35F;
        atom_context.raw_sampling_entries.emplace_back(LocalPotentialSample{
            static_cast<float>(
                model_list.at(atom_index).SignalAtDistance(point.distance) +
                target_offset_list.at(atom_index) *
                    model_list.at(atom_index).OffsetBasisAtDistance(
                        point.distance)),
            point
        });
        atom_context.neighbor_atom_sample_offset_list = { 0, 0 };
    }
    return {
        std::move(context),
        offset_detail::SecondStageModelSnapshot{
            model_list,
            {}
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
    const std::vector<std::size_t> & group_id_list,
    const std::vector<rg::GaussianModel3D> & base_model_list,
    const std::vector<rg::GaussianModel3D> & target_model_list)
{
    if (group_id_list.size() != base_model_list.size() ||
        base_model_list.size() != target_model_list.size() ||
        base_model_list.empty())
    {
        throw std::invalid_argument(
            "Joint polish fixture input sizes are inconsistent.");
    }

    constexpr std::array<float, 5> distance_list{
        0.0F,
        0.15F,
        0.30F,
        0.45F,
        0.60F
    };
    JointPolishFixture fixture;
    fixture.context.selected_atom_list.resize(base_model_list.size());
    std::unordered_map<std::size_t, std::size_t> selected_group_id_by_id;
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
                    static_cast<float>(
                        target_model_list.at(atom_index).ResponseAtDistance(
                            distance)),
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

        const auto group_position_iter{
            selected_group_id_by_id.find(
                group_id_list.at(atom_index))
        };
        const auto group_position{
            group_position_iter ==
                selected_group_id_by_id.end() ?
                fixture.context.selected_atom_index_list_by_group.size() :
                group_position_iter->second
        };
        if (group_position_iter ==
            selected_group_id_by_id.end())
        {
            selected_group_id_by_id.emplace(
                group_id_list.at(atom_index),
                group_position);
            fixture.context.selected_atom_index_list_by_group.emplace_back();
        }
        atom_context.group_id = group_position;
        fixture.context.selected_atom_index_list_by_group.at(group_position)
            .emplace_back(atom_index);
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
            point.distance = static_cast<float>(radius);
            point.position = {
                static_cast<float>(radius),
                0.0F,
                0.0F
            };
            sample_list.emplace_back(LocalPotentialSample{
                static_cast<float>(
                    zero_offset_response + previous_offset_response),
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
    return audit_detail::EvaluateSuspiciousGaussianUpdate(
        sample_entries,
        candidate_model,
        options,
        previous_baseline,
        audit_detail::SuspiciousUpdateMode::PostRefit);
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

    auto non_finite_sample_list{ sample_list };
    non_finite_sample_list.front().response =
        std::numeric_limits<float>::quiet_NaN();
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
        audit_detail::EvaluateSuspiciousOffsetUpdate(
            sample_list,
            previous_model,
            large_offset_model,
            options),
        audit_detail::SuspiciousGaussianReason::OffsetMagnitude);

    const auto wide_model{ MakeGaussianWithCenterSignal(0.1, 2.0) };
    EXPECT_EQ(
        audit_detail::EvaluateSuspiciousOffsetUpdate(
            sample_list,
            previous_model,
            wide_model,
            options),
        audit_detail::SuspiciousGaussianReason::None);
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
        audit_detail::EvaluateSuspiciousOffsetUpdate(
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
        audit_detail::EvaluateSuspiciousOffsetUpdate(
            noisy_sample_list,
            previous_model,
            candidate_with_center_offset(1.3),
            options),
        audit_detail::SuspiciousGaussianReason::None);
    EXPECT_EQ(
        audit_detail::EvaluateSuspiciousOffsetUpdate(
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
        audit_detail::EvaluateSuspiciousOffsetUpdate(
            zero_mad_sample_list,
            previous_model,
            candidate_with_center_offset(1.2),
            options),
        audit_detail::SuspiciousGaussianReason::None);
    EXPECT_EQ(
        audit_detail::EvaluateSuspiciousOffsetUpdate(
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
        audit_detail::EvaluateSuspiciousOffsetUpdate(
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
        audit_detail::EvaluateSuspiciousOffsetUpdate(
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
        audit_detail::EvaluateSuspiciousOffsetUpdate(
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
        audit_detail::EvaluateSuspiciousOffsetUpdate(
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
        audit_detail::EvaluateSuspiciousOffsetUpdate(
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
        audit_detail::EvaluateSuspiciousOffsetUpdate(
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

TEST(EstimatorSecondStageDefenseTest, SuspiciousRollbackExpandsOnlyWithinSharedOffsetGroup)
{
    EXPECT_EQ(
        audit_detail::ExpandSuspiciousSharedOffsetGroups(
            std::vector<std::size_t>{ 10, 10, 20, 30, 20 },
            std::vector<char>{ 0, 1, 0, 0, 0 }),
        (std::vector<char>{ 1, 1, 0, 0, 0 }));
    EXPECT_EQ(
        audit_detail::ExpandSuspiciousSharedOffsetGroups(
            std::vector<std::size_t>{ 10, 20, 10, 20 },
            std::vector<char>{ 0, 1, 0, 0 }),
        (std::vector<char>{ 0, 1, 0, 1 }));
}

std::unique_ptr<rg::ModelObject> BuildDefenseModel(
    const std::vector<std::array<float, 3>> & position_list,
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
        auto local_editor{ analysis.EnsureAtomLocalPotential(*atom) };
        local_editor.SetAlphaR(FittingStage::Second, 0.0);
        local_editor.SetGaussianResult(
            FittingStage::Second,
            MakeGaussianResult(initial_model));
        local_editor.SetRawSamplingEntries(
            BuildSamples(*atom, selected_atoms, truth_model_list));
    }
    return model;
}

std::unique_ptr<rg::ModelObject> BuildUnselectedContributorDefenseModel(
    const rg::GaussianModel3D & model_value)
{
    const std::vector<std::array<float, 3>> position_list{
        { 0.0F, 0.0F, 0.0F },
        { 0.9F, 0.0F, 0.0F },
        { 0.45F, 0.0F, 0.0F },
        { 0.60F, 0.0F, 0.0F },
        { 0.75F, 0.0F, 0.0F },
        { 4.50F, 0.0F, 0.0F },
        { 10.0F, 0.0F, 0.0F }
    };
    const std::vector<Spot> spot_list{
        Spot::C,
        Spot::C,
        Spot::C,
        Spot::O,
        Spot::C,
        Spot::O,
        Spot::O
    };
    const std::vector<Element> element_list{
        Element::CARBON,
        Element::CARBON,
        Element::CARBON,
        Element::OXYGEN,
        Element::HYDROGEN,
        Element::OXYGEN,
        Element::OXYGEN
    };
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
    for (int serial_id = 3; serial_id <= 7; serial_id++)
    {
        model->SetAtomSelected(serial_id, false);
    }

    std::vector<rg::AtomObject *> all_atoms;
    std::vector<rg::GaussianModel3D> truth_models;
    for (const auto & atom : model->GetAtomList())
    {
        all_atoms.emplace_back(atom.get());
        truth_models.emplace_back(model_value);
    }

    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    analysis.InitializeGroupAlpha(rg::FittingStage::Second, 0.0);
    for (auto * atom : model->GetSelectedAtoms())
    {
        auto local_editor{ analysis.EnsureAtomLocalPotential(*atom) };
        local_editor.SetAlphaR(FittingStage::Second, 0.0);
        local_editor.SetGaussianResult(
            FittingStage::Second,
            MakeGaussianResult(model_value));
        local_editor.SetRawSamplingEntries(
            BuildSamples(*atom, all_atoms, truth_models));
    }
    return model;
}

void ExpectUnselectedContributorPeeling(
    const rg::ModelObject & model,
    const rg::GaussianModel3D & global_fallback_model,
    bool exclude_hydrogen)
{
    const auto & selected_atoms{ model.GetSelectedAtoms() };
    std::vector<rg::GaussianModel3D> selected_models;
    std::vector<GroupKey> selected_group_keys;
    selected_models.reserve(selected_atoms.size());
    selected_group_keys.reserve(selected_atoms.size());
    for (const auto * atom : selected_atoms)
    {
        selected_models.emplace_back(
            rg::AtomLocalPotentialView::RequireFor(*atom).GetEstimateMDPDE(
                FittingStage::Second));
        selected_group_keys.emplace_back(
            rg::data_internal::GetGroupKey(atom));
    }

    for (const auto * target_atom : selected_atoms)
    {
        const auto view{ rg::AtomLocalPotentialView::RequireFor(*target_atom) };
        const auto raw_entries{ view.GetRawSamplingEntries(false) };
        const auto peeling_entries{ view.GetPeelingSamplingEntries(false) };
        ASSERT_EQ(raw_entries.size(), peeling_entries.size());
        for (std::size_t sample_index = 0;
            sample_index < raw_entries.size();
            sample_index++)
        {
            const auto & raw_sample{ raw_entries.at(sample_index) };
            auto expected_response{ static_cast<double>(raw_sample.response) };
            for (const auto * neighbor : target_atom->FindNeighborAtoms(5.0))
            {
                if (exclude_hydrogen &&
                    neighbor->GetElement() == Element::HYDROGEN)
                {
                    continue;
                }
                const auto sample_distance{
                    Distance(raw_sample.point.position, neighbor->GetPosition())
                };
                if (sample_distance > 2.5) continue;

                const auto selected_iter{
                    std::find(selected_atoms.begin(), selected_atoms.end(), neighbor)
                };
                rg::GaussianModel3D contributor_model{ global_fallback_model };
                if (selected_iter != selected_atoms.end())
                {
                    contributor_model =
                        rg::AtomLocalPotentialView::RequireFor(*neighbor)
                            .GetEstimateMDPDE(FittingStage::Second);
                }
                else
                {
                    std::vector<rg::GaussianModel3D> group_models;
                    const auto neighbor_group_key{
                        rg::data_internal::GetGroupKey(neighbor)
                    };
                    for (std::size_t i = 0; i < selected_atoms.size(); i++)
                    {
                        if (selected_group_keys.at(i) == neighbor_group_key)
                        {
                            group_models.emplace_back(selected_models.at(i));
                        }
                    }
                    const auto group_median{
                        median_detail::BuildGaussianParameterMedian(
                            group_models)
                    };
                    if (group_median.has_value())
                    {
                        contributor_model = *group_median;
                    }
                }
                expected_response -= contributor_model.ResponseAtDistance(
                    sample_distance);
            }
            EXPECT_FLOAT_EQ(
                peeling_entries.at(sample_index).response,
                static_cast<float>(expected_response));
        }
    }
}

std::unique_ptr<rg::ModelObject> BuildNearCollinearDefenseModel(
    double intensity_scale = 1.0)
{
    return BuildDefenseModel(
        {
            std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
            std::array<float, 3>{ 1.0e-4F, 0.0F, 0.0F }
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
            std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
            std::array<float, 3>{ 0.8F, 0.0F, 0.0F }
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
    analysis.EnsureAtomLocalPotential(*atom_list.at(0))
        .SetGaussianResult(FittingStage::Second, MakeGaussianResult(
            rg::GaussianModel3D{ 5.8, 0.47, 0.15 }));
    analysis.EnsureAtomLocalPotential(*atom_list.at(1))
        .SetGaussianResult(FittingStage::Second, MakeGaussianResult(
            rg::GaussianModel3D{ 4.7, 0.68, -0.08 }));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildSharedOffsetJointPolishDefenseModel(
    double intensity_scale = 1.0)
{
    auto model{ BuildDefenseModel(
        {
            std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
            std::array<float, 3>{ 0.8F, 0.0F, 0.0F }
        },
        { Spot::C, Spot::C },
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
                0.05 * intensity_scale
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
            rg::AtomLocalPotentialView::RequireFor(*atom_list.at(i))
                .GetGaussianResult(FittingStage::Second).mdpde.GetModel()
        };
        analysis.EnsureAtomLocalPotential(*atom_list.at(i))
            .SetGaussianResult(FittingStage::Second, MakeGaussianResult(
                current.WithOffset(initial_offset_list.at(i))));
    }
    return model;
}

void MakeAtomSamplesSuspicious(rg::ModelObject & model, std::size_t atom_index)
{
    const auto & selected_atoms{ model.GetSelectedAtoms() };
    auto * target_atom{ selected_atoms.at(atom_index) };
    const auto target_position{ target_atom->GetPosition() };

    auto raw_sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*target_atom)
            .GetRawSamplingEntries(false)
    };
    raw_sampling_entries.resize(256);
    raw_sampling_entries.front().response = 0.0F;
    raw_sampling_entries.front().point.distance = 0.0F;
    raw_sampling_entries.front().point.position = target_position;
    for (std::size_t i = 1; i < raw_sampling_entries.size(); i++)
    {
        auto & sample{ raw_sampling_entries.at(i) };
        const auto response_scale{
            0.5F + 0.5F * static_cast<float>(i) /
                static_cast<float>(raw_sampling_entries.size())
        };
        sample.response = std::numeric_limits<float>::max() * response_scale;
        sample.point.position = target_position;
        sample.point.position.at(0) += 100.0F;
        sample.point.distance = 100.0F;
    }

    auto analysis{ model.EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*target_atom).SetRawSamplingEntries(
        std::move(raw_sampling_entries));
}

std::unique_ptr<rg::ModelObject> BuildSeparatedRollbackDefenseModel()
{
    auto model{
        BuildDefenseModel(
            {
                std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 1.0e-4F, 0.0F, 0.0F },
                std::array<float, 3>{ 10.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 10.0001F, 0.0F, 0.0F }
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
    std::vector<std::array<float, 3>> position_list;
    std::vector<Spot> spot_list;
    std::vector<Element> element_list;
    std::vector<rg::GaussianModel3D> truth_model_list;
    for (std::size_t i = 0; i < 13; i++)
    {
        const auto x_position{
            i < 11 ? 0.45F * static_cast<float>(i) :
                20.0F + 0.45F * static_cast<float>(i - 11)
        };
        position_list.push_back({ x_position, 0.0F, 0.0F });
        spot_list.push_back(i % 2 == 0 ? Spot::C : Spot::O);
        element_list.push_back(
            i % 2 == 0 ? Element::CARBON : Element::OXYGEN);
        const auto truth_model{
            i >= 11 ? rg::GaussianModel3D{ 7.5, 0.70, 0.10 } :
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

std::unique_ptr<rg::ModelObject> BuildSeparatedSystemBuildFailureDefenseModel()
{
    auto model{ BuildSeparatedRollbackDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    auto raw_sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*atom)
            .GetRawSamplingEntries(false)
    };
    raw_sampling_entries.front().response =
        std::numeric_limits<float>::quiet_NaN();
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*atom).SetRawSamplingEntries(
        std::move(raw_sampling_entries));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildTerminalWithPersistentLocalRefitFallbackDefenseModel()
{
    auto model{ BuildSeparatedRollbackDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().at(2) };
    auto raw_sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*atom)
            .GetRawSamplingEntries(false)
    };
    raw_sampling_entries.resize(1);
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*atom).SetRawSamplingEntries(
        std::move(raw_sampling_entries));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildSeparatedLocalRefitFallbackDefenseModel()
{
    auto model{
        BuildDefenseModel(
            {
                std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 1.0e-4F, 0.0F, 0.0F },
                std::array<float, 3>{ 10.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 10.0001F, 0.0F, 0.0F }
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
    auto * atom{ model->GetSelectedAtoms().front() };
    auto raw_sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*atom)
            .GetRawSamplingEntries(false)
    };
    LocalPotentialSampleList fallback_sampling_entries{
        raw_sampling_entries.at(0),
        raw_sampling_entries.at(6)
    };
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*atom).SetRawSamplingEntries(
        std::move(fallback_sampling_entries));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildSeparatedEmptyJointOffsetDefenseModel()
{
    auto model{
        BuildDefenseModel(
            {
                std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 10.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 10.0001F, 0.0F, 0.0F }
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
    analysis.EnsureAtomLocalPotential(*model->GetSelectedAtoms().front())
        .SetRawSamplingEntries({});
    return model;
}

std::unique_ptr<rg::ModelObject> BuildPostRefitRollbackChainDefenseModel()
{
    auto model{
        BuildDefenseModel(
            {
                std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 2.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 4.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 6.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 8.0F, 0.0F, 0.0F }
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
    analysis.EnsureAtomLocalPotential(*model->GetSelectedAtoms().front())
        .SetGaussianResult(FittingStage::Second, MakeGaussianResult(
            rg::GaussianModel3D{ 1.0e300, 0.55, 0.0 }));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildNonFiniteJointOffsetDefenseModel()
{
    auto model{
        BuildDefenseModel(
            { std::array<float, 3>{ 0.0F, 0.0F, 0.0F } },
            { Spot::O },
            { Element::OXYGEN },
            { rg::GaussianModel3D{ 8.0, 0.5, -0.1 } },
            rg::GaussianModel3D{ 7.0, 0.5, 0.0 })
    };
    auto * atom{ model->GetSelectedAtoms().front() };
    auto raw_sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*atom)
            .GetRawSamplingEntries(false)
    };
    raw_sampling_entries.front().response =
        std::numeric_limits<float>::quiet_NaN();
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*atom).SetRawSamplingEntries(
        std::move(raw_sampling_entries));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildFiniteNonphysicalProfileDefenseModel()
{
    const rg::GaussianModel3D initial_model{ 6.0, 0.55, 0.0 };
    auto model{
        BuildDefenseModel(
            { std::array<float, 3>{ 0.0F, 0.0F, 0.0F } },
            { Spot::O },
            { Element::OXYGEN },
            { initial_model },
            initial_model)
    };
    auto * atom{ model->GetSelectedAtoms().front() };
    auto raw_sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*atom)
            .GetRawSamplingEntries(false)
    };
    for (auto & sample : raw_sampling_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto outer_bias{ distance > 0.2 ? 12.0 : 8.0 };
        sample.response = static_cast<float>(
            initial_model.SignalAtDistance(distance) + outer_bias);
    }
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*atom).SetRawSamplingEntries(
        std::move(raw_sampling_entries));
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
        const auto local_view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
        for (const auto & sample : local_view.GetRawSamplingEntries(false))
        {
            double fitted_response{ 0.0 };
            for (const auto * fitted_atom : selected_atoms)
            {
                const auto fitted_view{
                    rg::AtomLocalPotentialView::RequireFor(*fitted_atom)
                };
                fitted_response += fitted_view.GetEstimateMDPDE(
                    FittingStage::Second).ResponseAtDistance(
                    Distance(sample.point.position, fitted_atom->GetPosition()));
            }
            const auto residual{ static_cast<double>(sample.response) - fitted_response };
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
    return rg::AtomLocalPotentialView::RequireFor(atom).GetEstimateMDPDE(
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
            rg::AtomLocalPotentialView::RequireFor(*atom).GetEstimateMDPDE(
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
            rg::AtomLocalPotentialView::RequireFor(*target_atom)
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
            auto expected_response{ static_cast<double>(raw_sample.response) };
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

            EXPECT_FLOAT_EQ(
                peeling_sample.response,
                static_cast<float>(expected_response));
            EXPECT_FLOAT_EQ(
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

TEST(EstimatorSecondStageDefenseTest, SeedSelectionUsesConfiguredFallbackPriority)
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
    seed_detail::SecondStageSeedCandidates candidates;
    candidates.group_posterior = make_candidate(1.0);
    candidates.group_prior = make_candidate(2.0);
    candidates.group_median = make_candidate(3.0);
    candidates.global_median = make_candidate(4.0);

    const auto expect_source = [&](seed_detail::SecondStageSeedSource source)
    {
        const auto selection{ seed_detail::SelectSecondStageSeed(candidates) };
        ASSERT_TRUE(selection.has_value());
        EXPECT_EQ(selection->source, source);
    };
    expect_source(seed_detail::SecondStageSeedSource::GroupPosterior);
    candidates.group_posterior = invalid_candidate;
    expect_source(seed_detail::SecondStageSeedSource::GroupPrior);
    candidates.group_prior = invalid_candidate;
    expect_source(seed_detail::SecondStageSeedSource::GroupMedian);
    candidates.group_median = invalid_candidate;
    expect_source(seed_detail::SecondStageSeedSource::GlobalMedian);
    candidates.global_median = invalid_candidate;
    EXPECT_FALSE(seed_detail::SelectSecondStageSeed(candidates).has_value());
}

TEST(EstimatorSecondStageDefenseTest, SeedSelectionReturnsCompleteSourceModelAndUncertainty)
{
    seed_detail::SecondStageSeedCandidates candidates;
    candidates.group_posterior =
        rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 6.0, 0.55, -0.2 },
            rg::GaussianModel3DUncertainty{ 0.1, 0.02, 0.03 }
    };
    candidates.group_prior = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 4.0, 0.65, 0.4 },
        rg::GaussianModel3DUncertainty{}
    };

    const auto selection{ seed_detail::SelectSecondStageSeed(candidates) };
    ASSERT_TRUE(selection.has_value());
    EXPECT_EQ(selection->source, seed_detail::SecondStageSeedSource::GroupPosterior);
    EXPECT_DOUBLE_EQ(selection->model.GetModel().GetAmplitude(), 6.0);
    EXPECT_DOUBLE_EQ(selection->model.GetModel().GetWidth(), 0.55);
    EXPECT_DOUBLE_EQ(selection->model.GetModel().GetOffset(), -0.2);
    EXPECT_DOUBLE_EQ(
        selection->model.GetStandardDeviationModel().GetAmplitude(),
        0.1);
    EXPECT_DOUBLE_EQ(
        selection->model.GetStandardDeviationModel().GetWidth(),
        0.02);
    EXPECT_DOUBLE_EQ(
        selection->model.GetStandardDeviationModel().GetOffset(),
        0.03);
}

TEST(EstimatorSecondStageDefenseTest, SeedInitializationBuildsMediansOnlyFromDirectSources)
{
    auto model{
        BuildDefenseModel(
            {
                std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 6.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 12.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 18.0F, 0.0F, 0.0F }
            },
            { Spot::C, Spot::C, Spot::O, Spot::N },
            {
                Element::CARBON,
                Element::CARBON,
                Element::OXYGEN,
                Element::NITROGEN
            },
            {
                rg::GaussianModel3D{ 2.0, 0.5, 0.2 },
                rg::GaussianModel3D{ 3.0, 0.5, 0.2 },
                rg::GaussianModel3D{ 10.0, 0.7, -0.2 },
                rg::GaussianModel3D{ 6.0, 0.6, 0.0 }
            },
            rg::GaussianModel3D{ 20.0, 0.8, 0.3 })
    };

    auto analysis{ model->EditAnalysis() };
    const auto & atom_list{ model->GetSelectedAtoms() };
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        auto result{ MakeGaussianResult(
            rg::GaussianModel3D{
                20.0 + static_cast<double>(i),
                0.8,
                0.3
            }) };
        result.posterior.reset();
        if (i == 0)
        {
            result.posterior = rg::GaussianModel3DWithUncertainty{
                rg::GaussianModel3D{ 2.0, 0.5, 0.2 },
                rg::GaussianModel3DUncertainty{ 0.1, 0.02, 0.03 }
            };
        }
        else if (i == 2)
        {
            result.posterior = rg::GaussianModel3DWithUncertainty{
                rg::GaussianModel3D{ 10.0, 0.7, -0.2 },
                rg::GaussianModel3DUncertainty{ 0.2, 0.03, 0.04 }
            };
        }
        analysis.EnsureAtomLocalPotential(*atom_list.at(i))
            .SetGaussianResult(
                FittingStage::Second,
                std::move(result));
    }

    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(
        out.find(
            "Selected second-stage initial seeds = 4, sources = "
            "group-posterior:2, group-prior:0, group-median:1, global-median:1."),
        std::string::npos);
    EXPECT_NE(
        out.find(
            "atom index = 1, source = group-median, "
            "original MDPDE A/B/C = 2.10e+01/8.00e-01/3.00e-01, "
            "selected A/B/C = 2.00e+00/5.00e-01/2.00e-01."),
        std::string::npos);
    EXPECT_NE(
        out.find(
            "atom index = 3, source = global-median, "
            "original MDPDE A/B/C = 2.30e+01/8.00e-01/3.00e-01, "
            "selected A/B/C = 6.00e+00/6.00e-01/0.00e+00."),
        std::string::npos);
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
    EXPECT_FALSE(audit_detail::IsBetterAuditObjective(
        std::numeric_limits<double>::infinity(), 1.0, tolerance));
    EXPECT_TRUE(audit_detail::IsBetterAuditObjective(
        1.0, std::numeric_limits<double>::infinity(), tolerance));
    EXPECT_THROW(
        audit_detail::IsBetterAuditObjective(
            0.8,
            1.0,
            audit_detail::ObjectiveTolerance{ -1.0, 0.0 }),
        std::invalid_argument);
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
}

TEST(EstimatorSecondStageDefenseTest, AuditObjectiveProgressGuardChecksPreviousAndBest)
{
    const audit_detail::ObjectiveTolerance tolerance{
        1.0e-8,
        1.0e-3
    };
    const audit_detail::ObjectiveBreakdown best_one{ 1.0, 0.0, 0.0 };
    const audit_detail::ObjectiveBreakdown best_below{ 0.99, 0.0, 0.0 };
    const audit_detail::ObjectiveBreakdown best_infinite{
        std::numeric_limits<double>::infinity(),
        0.0,
        0.0
    };
    EXPECT_TRUE(audit_detail::IsAuditObjectiveAcceptableForProgress(
        1.0005, 1.0, &best_one, tolerance));
    EXPECT_FALSE(audit_detail::IsAuditObjectiveAcceptableForProgress(
        1.002, 1.0, &best_one, tolerance));
    EXPECT_FALSE(audit_detail::IsAuditObjectiveAcceptableForProgress(
        1.0, 1.0, &best_below, tolerance));
    EXPECT_FALSE(audit_detail::IsAuditObjectiveAcceptableForProgress(
        std::numeric_limits<double>::infinity(),
        1.0,
        nullptr,
        tolerance));
    EXPECT_FALSE(audit_detail::IsAuditObjectiveAcceptableForProgress(
        1.0,
        std::numeric_limits<double>::infinity(),
        nullptr,
        tolerance));
    EXPECT_FALSE(audit_detail::IsAuditObjectiveAcceptableForProgress(
        1.0,
        1.0,
        &best_infinite,
        tolerance));
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
    EXPECT_THROW(
        audit_detail::IsAuditObjectiveAcceptableForProgress(
            1.0,
            1.0,
            nullptr,
            audit_detail::ObjectiveTolerance{ 0.0, -1.0 }),
        std::invalid_argument);
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
    EXPECT_THROW(
        audit_detail::CalculateClusterAtomWeight(0, 4),
        std::invalid_argument);
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

TEST(EstimatorSecondStageDefenseTest, TrustRegionDampingCapsLargeTransformedStep)
{
    std::vector<Eigen::Vector3d> previous{
        Eigen::Vector3d::Zero()
    };
    auto candidate{ previous };
    candidate.at(0) << 1.0, 0.35, 0.5;

    const auto capped{
        trust_detail::LimitTrustRegionDamping(
            previous, candidate, 1.0, 1.0)
    };
    EXPECT_DOUBLE_EQ(capped.effective_damping, 0.5);
    EXPECT_DOUBLE_EQ(capped.step_norm, 1.0);

    const auto inside{
        trust_detail::LimitTrustRegionDamping(
            previous, candidate, 0.25, 1.0)
    };
    EXPECT_DOUBLE_EQ(inside.effective_damping, 0.25);
    EXPECT_DOUBLE_EQ(inside.step_norm, 0.5);
}

TEST(EstimatorSecondStageDefenseTest, TrustRegionDampingIsIntensityScaleInvariant)
{
    const auto encode = [](const rg::GaussianModel3D & model)
    {
        const auto estimation{
            change_detail::EncodeTransformedCoordinates(model)
        };
        EXPECT_TRUE(estimation.has_value());
        return estimation.value_or(Eigen::Vector3d::Zero());
    };
    const std::vector<Eigen::Vector3d> base_previous{
        encode(rg::GaussianModel3D{ 2.0, 0.8, 0.2 })
    };
    const std::vector<Eigen::Vector3d> base_candidate{
        encode(rg::GaussianModel3D{ 4.0, 1.0, 0.4 })
    };
    constexpr double intensity_scale{ 1.0e5 };
    const std::vector<Eigen::Vector3d> scaled_previous{
        encode(rg::GaussianModel3D{
            2.0 * intensity_scale,
            0.8,
            0.2 * intensity_scale })
    };
    const std::vector<Eigen::Vector3d> scaled_candidate{
        encode(rg::GaussianModel3D{
            4.0 * intensity_scale,
            1.0,
            0.4 * intensity_scale })
    };

    const auto base{
        trust_detail::LimitTrustRegionDamping(
            base_previous, base_candidate, 1.0, 0.5)
    };
    const auto scaled{
        trust_detail::LimitTrustRegionDamping(
            scaled_previous, scaled_candidate, 1.0, 0.5)
    };
    EXPECT_NEAR(base.effective_damping, scaled.effective_damping, 1.0e-12);
    EXPECT_NEAR(base.step_norm, scaled.step_norm, 1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, TrustRegionStateReconcilesShrinksGrowsAndSaturates)
{
    trust_detail::TrustRegionStateSet state;
    const trust_detail::ClusterKey key{ 0 };
    state.Reconcile({ key });
    EXPECT_DOUBLE_EQ(state.GetRadius(key), 1.0);

    for (const auto expected : { 0.5, 0.25, 0.125, 0.0625 })
    {
        const auto update{ state.Shrink({ key }) };
        EXPECT_EQ(
            update.changed_key_list,
            std::vector<trust_detail::ClusterKey>{ key });
        EXPECT_TRUE(update.saturated_key_list.empty());
        EXPECT_DOUBLE_EQ(state.GetRadius(key), expected);
    }
    const auto saturated{ state.Shrink({ key }) };
    EXPECT_TRUE(saturated.changed_key_list.empty());
    EXPECT_EQ(
        saturated.saturated_key_list,
        std::vector<trust_detail::ClusterKey>{ key });

    state.Grow({ key });
    EXPECT_DOUBLE_EQ(state.GetRadius(key), 0.125);

    const trust_detail::ClusterKey replacement_key{ 1 };
    state.Reconcile({ replacement_key });
    EXPECT_THROW(state.GetRadius(key), std::invalid_argument);
    EXPECT_DOUBLE_EQ(state.GetRadius(replacement_key), 1.0);
}

TEST(EstimatorSecondStageDefenseTest, ExhaustedRejectionsAreExcludedFromRadiusShrink)
{
    using Key = trust_detail::ClusterKey;
    const Key exhausted_key{ 0 };
    const Key retryable_key{ 1 };
    const auto partition{
        trust_detail::PartitionRejectedClusters(
            { exhausted_key, retryable_key },
            { exhausted_key })
    };
    ASSERT_EQ(partition.exhausted_key_list, std::vector<Key>{ exhausted_key });
    ASSERT_EQ(partition.retryable_key_list, std::vector<Key>{ retryable_key });

    trust_detail::TrustRegionStateSet state;
    state.Reconcile({ exhausted_key, retryable_key });
    const auto update{ state.Shrink(partition.retryable_key_list) };

    EXPECT_DOUBLE_EQ(state.GetRadius(exhausted_key), 1.0);
    EXPECT_DOUBLE_EQ(state.GetRadius(retryable_key), 0.5);
    EXPECT_EQ(update.changed_key_list, std::vector<Key>{ retryable_key });

    const auto all_exhausted{
        trust_detail::PartitionRejectedClusters(
            { exhausted_key },
            { exhausted_key })
    };
    const auto exhausted_update{
        state.Shrink(all_exhausted.retryable_key_list)
    };
    EXPECT_TRUE(exhausted_update.changed_key_list.empty());
    EXPECT_TRUE(exhausted_update.saturated_key_list.empty());
    EXPECT_DOUBLE_EQ(state.GetRadius(exhausted_key), 1.0);
}

TEST(EstimatorSecondStageDefenseTest, AllRejectedResolutionDistinguishesTerminalCases)
{
    using Key = trust_detail::ClusterKey;
    using Resolution = trust_detail::AllRejectedResolution;
    const Key first_key{ 0 };
    const Key second_key{ 1 };

    const auto exhausted{
        trust_detail::PartitionRejectedClusters(
            { first_key, second_key },
            { first_key, second_key })
    };
    EXPECT_EQ(
        trust_detail::ResolveAllRejected(false, exhausted, {}),
        Resolution::BacktrackingExhausted);

    const auto retryable{
        trust_detail::PartitionRejectedClusters(
            { first_key, second_key },
            {})
    };
    trust_detail::TrustRegionRadiusUpdate changed;
    changed.changed_key_list = { first_key };
    EXPECT_EQ(
        trust_detail::ResolveAllRejected(
            false,
            retryable,
            changed),
        Resolution::Retry);

    trust_detail::TrustRegionRadiusUpdate saturated;
    saturated.saturated_key_list = { first_key, second_key };
    EXPECT_EQ(
        trust_detail::ResolveAllRejected(
            false,
            retryable,
            saturated),
        Resolution::MinimumRadius);

    const auto mixed{
        trust_detail::PartitionRejectedClusters(
            { first_key, second_key },
            { first_key })
    };
    saturated.saturated_key_list = { second_key };
    EXPECT_EQ(
        trust_detail::ResolveAllRejected(
            false,
            mixed,
            saturated),
        Resolution::NoRetryProgress);
    EXPECT_EQ(
        trust_detail::ResolveAllRejected(
            true,
            mixed,
            changed),
        Resolution::MaximumIterations);
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

TEST(EstimatorSecondStageDefenseTest,
    JointOffsetParameterizationSharesColumnsAndUsesGroupMedianSeeds)
{
    Eigen::VectorXd initial_atom_offset(5);
    initial_atom_offset << 1.0, 4.0, 3.0, 2.0, 6.0;
    const auto parameterization{
        offset_detail::BuildJointOffsetParameterization(
            std::vector<std::size_t>{ 20, 10, 20, 20, 10 },
            initial_atom_offset)
    };

    ASSERT_TRUE(parameterization.has_value());
    EXPECT_EQ(parameterization->group_position_by_atom.size(), 5U);
    EXPECT_EQ(parameterization->seed_offset.size(), 2);
    EXPECT_EQ(parameterization->OffsetColumn(0), 1);
    EXPECT_EQ(parameterization->OffsetColumn(1), 0);
    EXPECT_EQ(
        parameterization->OffsetColumn(0),
        parameterization->OffsetColumn(2));
    EXPECT_EQ(
        parameterization->OffsetColumn(0),
        parameterization->OffsetColumn(3));
    EXPECT_EQ(
        parameterization->OffsetColumn(1),
        parameterization->OffsetColumn(4));
    EXPECT_DOUBLE_EQ(
        parameterization->seed_offset(parameterization->OffsetColumn(0)),
        2.0);
    EXPECT_DOUBLE_EQ(
        parameterization->seed_offset(parameterization->OffsetColumn(1)),
        5.0);

    Eigen::VectorXd group_offset{ Eigen::VectorXd::Zero(2) };
    group_offset(parameterization->OffsetColumn(0)) = 2.5;
    group_offset(parameterization->OffsetColumn(1)) = 5.5;
    const auto atom_offset{ parameterization->ExpandOffsets(group_offset) };
    EXPECT_EQ(atom_offset.size(), 5);
    EXPECT_DOUBLE_EQ(atom_offset(0), 2.5);
    EXPECT_DOUBLE_EQ(atom_offset(1), 5.5);
    EXPECT_DOUBLE_EQ(atom_offset(2), 2.5);
    EXPECT_DOUBLE_EQ(atom_offset(3), 2.5);
    EXPECT_DOUBLE_EQ(atom_offset(4), 5.5);

    EXPECT_FALSE(
        offset_detail::BuildJointOffsetParameterization(
            std::vector<std::size_t>{ 20 },
            atom_offset).has_value());
    auto non_finite_atom_offset{ initial_atom_offset };
    non_finite_atom_offset(0) = std::numeric_limits<double>::infinity();
    EXPECT_FALSE(
        offset_detail::BuildJointOffsetParameterization(
            std::vector<std::size_t>{ 20, 10, 20, 20, 10 },
            non_finite_atom_offset).has_value());
}

TEST(EstimatorSecondStageDefenseTest,
    JointOffsetParameterizationUsesDeterministicGroupColumn)
{
    Eigen::VectorXd atom_offset(3);
    atom_offset << 1.0, 4.0, 3.0;
    const auto parameterization{
        offset_detail::BuildJointOffsetParameterization(
            std::vector<std::size_t>{ 20, 10, 20 },
            atom_offset)
    };
    ASSERT_TRUE(parameterization.has_value());

    Eigen::VectorXd reordered_atom_offset(3);
    reordered_atom_offset << 4.0, 1.0, 3.0;
    const auto reordered{
        offset_detail::BuildJointOffsetParameterization(
            std::vector<std::size_t>{ 10, 20, 20 },
            reordered_atom_offset)
    };
    ASSERT_TRUE(reordered.has_value());
    EXPECT_EQ(reordered->OffsetColumn(0), 0);
    EXPECT_EQ(reordered->OffsetColumn(1), 1);
    EXPECT_EQ(reordered->OffsetColumn(2), 1);
}

TEST(EstimatorSecondStageDefenseTest, LocalRefitHealthTracksStationarity)
{
    EXPECT_TRUE(health_detail::IsLocalRefitStatusStationarityEligible(
        rg::RHBMEstimationStatus::SUCCESS));
    EXPECT_FALSE(health_detail::IsLocalRefitStatusStationarityEligible(
        rg::RHBMEstimationStatus::MAX_ITERATIONS_REACHED));
    for (const auto status : {
        rg::RHBMEstimationStatus::NUMERICAL_FALLBACK,
        rg::RHBMEstimationStatus::INSUFFICIENT_DATA,
        rg::RHBMEstimationStatus::SINGLE_MEMBER })
    {
        EXPECT_FALSE(health_detail::IsLocalRefitStatusStationarityEligible(status));
    }
    EXPECT_THROW(
        health_detail::IsLocalRefitStatusStationarityEligible(
            static_cast<rg::RHBMEstimationStatus>(-1)),
        std::logic_error);
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetHealthSeparatesHardFailureFromStationarity)
{
    using Status = offset_detail::JointOffsetSolveStatus;

    EXPECT_TRUE(health_detail::ClusterHealth{ Status::Converged }.IsStationarityEligible());
    EXPECT_FALSE(offset_detail::IsJointOffsetSolveHardFailure(Status::Converged));

    for (const auto status : {
        Status::IrlsObjectiveDeteriorated,
        Status::IrlsMaximumIterationsReached })
    {
        EXPECT_FALSE(health_detail::ClusterHealth{ status }.IsStationarityEligible());
        EXPECT_FALSE(offset_detail::IsJointOffsetSolveHardFailure(status));
    }

    for (const auto status : {
        Status::SystemBuildFailed,
        Status::EmptySystem,
        Status::InitialSolveFailed,
        Status::IrlsSolveFailed })
    {
        EXPECT_FALSE(health_detail::ClusterHealth{ status }.IsStationarityEligible());
        EXPECT_TRUE(offset_detail::IsJointOffsetSolveHardFailure(status));
    }

    const auto invalid_status{ static_cast<Status>(-1) };
    EXPECT_FALSE(health_detail::ClusterHealth{ invalid_status }.IsStationarityEligible());
    EXPECT_THROW(
        offset_detail::IsJointOffsetSolveHardFailure(invalid_status),
        std::logic_error);
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetEstimatorSharesGroupOffsets)
{
    auto fixture{
        BuildJointOffsetEstimationFixture(
            { 20, 20 },
            {
                rg::GaussianModel3D{ 6.0, 0.55, 0.0 },
                rg::GaussianModel3D{ 7.0, 0.60, 4.0 }
            },
            { 2.0, 2.0 })
    };
    offset_detail::ReusableWeightedRidgeSolver solver;
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
    EXPECT_NEAR(result.offset(0), 2.0, 1.0e-5);
    EXPECT_DOUBLE_EQ(result.offset(0), result.offset(1));
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetEstimatorKeepsIndependentGroups)
{
    auto fixture{
        BuildJointOffsetEstimationFixture(
            { 20, 10 },
            {
                rg::GaussianModel3D{ 6.0, 0.55, 1.0 },
                rg::GaussianModel3D{ 7.0, 0.60, 3.0 }
            },
            { 1.0, 3.0 })
    };
    offset_detail::ReusableWeightedRidgeSolver solver;
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
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetEstimatorReportsBuildAndEmptyFailures)
{
    auto empty_fixture{
        BuildJointOffsetEstimationFixture(
            { 20 },
            { rg::GaussianModel3D{ 6.0, 0.55, 2.0 } },
            { 2.0 })
    };
    empty_fixture.first.at(0).raw_sampling_entries.clear();
    empty_fixture.first.at(0).neighbor_atom_sample_offset_list.clear();
    offset_detail::ReusableWeightedRidgeSolver empty_solver;
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
            { 20 },
            { rg::GaussianModel3D{ 6.0, 0.55, 2.0 } },
            { 2.0 })
    };
    invalid_fixture.first.at(0).raw_sampling_entries.at(0).response =
        std::numeric_limits<float>::infinity();
    offset_detail::ReusableWeightedRidgeSolver invalid_solver;
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
}

TEST(EstimatorSecondStageDefenseTest, JointOffsetStatusTextCoversAllStatuses)
{
    using Status = offset_detail::JointOffsetSolveStatus;
    EXPECT_STREQ(
        offset_detail::GetJointOffsetSolveStatusText(Status::Converged),
        "converged");
    EXPECT_STREQ(
        offset_detail::GetJointOffsetSolveStatusText(Status::SystemBuildFailed),
        "system-build-failed");
    EXPECT_STREQ(
        offset_detail::GetJointOffsetSolveStatusText(Status::EmptySystem),
        "empty-system");
    EXPECT_STREQ(
        offset_detail::GetJointOffsetSolveStatusText(Status::InitialSolveFailed),
        "initial-solve-failed");
    EXPECT_STREQ(
        offset_detail::GetJointOffsetSolveStatusText(Status::IrlsSolveFailed),
        "irls-solve-failed");
    EXPECT_STREQ(
        offset_detail::GetJointOffsetSolveStatusText(
            Status::IrlsObjectiveDeteriorated),
        "irls-objective-deteriorated");
    EXPECT_STREQ(
        offset_detail::GetJointOffsetSolveStatusText(
            Status::IrlsMaximumIterationsReached),
        "irls-maximum-iterations-reached");
    EXPECT_THROW(
        offset_detail::GetJointOffsetSolveStatusText(
            static_cast<Status>(-1)),
        std::logic_error);
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
        ASSERT_EQ(base_change.value_list.size(), scaled_change.value_list.size());
        for (std::size_t i = 0; i < base_change.value_list.size(); i++)
        {
            EXPECT_NEAR(
                base_change.value_list.at(i),
                scaled_change.value_list.at(i),
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
        change_detail::EncodeTransformedCoordinates(
            reference_state.at(0).mdpde.GetModel())
    };
    ASSERT_TRUE(reference_coordinates.has_value());
    auto small_coordinates{ *reference_coordinates };
    auto large_coordinates{ *reference_coordinates };
    small_coordinates(change_detail::kLogWidthChangeIndex) += 0.099;
    large_coordinates(change_detail::kLogWidthChangeIndex) += 0.101;
    const auto small_model{
        change_detail::DecodeTransformedCoordinates(small_coordinates)
    };
    const auto large_model{
        change_detail::DecodeTransformedCoordinates(large_coordinates)
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
        const auto transformed{
            change_detail::EncodeTransformedCoordinates(model)
        };
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
                    change_detail::DecodeTransformedCoordinates(lower)
                };
                const auto upper_model{
                    change_detail::DecodeTransformedCoordinates(upper)
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

TEST(EstimatorSecondStageDefenseTest,
    GroupMedianModelListUsesComponentMediansAndKeepsGroupsSeparate)
{
    const std::vector<std::size_t> group_id_list{ 10, 10, 10, 20, 20, 30 };
    const std::vector<rg::GaussianModel3D> model_list{
        rg::GaussianModel3D{ 1.0, 0.20, 1.0 },
        rg::GaussianModel3D{ 9.0, 0.60, 3.0 },
        rg::GaussianModel3D{ 5.0, 0.40, 2.0 },
        rg::GaussianModel3D{ 4.0, 0.50, -2.0 },
        rg::GaussianModel3D{ 8.0, 0.90, 2.0 },
        rg::GaussianModel3D{ 7.0, 0.80, 0.4 }
    };

    const auto median_model_list{
        median_detail::BuildGroupMedianModelList(
            group_id_list,
            model_list)
    };

    ASSERT_EQ(median_model_list.size(), model_list.size());
    for (const auto atom_position : std::array<std::size_t, 3>{ 0, 1, 2 })
    {
        ExpectGaussianModelsNear(
            median_model_list.at(atom_position),
            rg::GaussianModel3D{ 5.0, 0.40, 2.0 },
            1.0e-12);
    }
    for (const auto atom_position : std::array<std::size_t, 2>{ 3, 4 })
    {
        ExpectGaussianModelsNear(
            median_model_list.at(atom_position),
            rg::GaussianModel3D{ 6.0, 0.70, 0.0 },
            1.0e-12);
    }
    ExpectGaussianModelsNear(
        median_model_list.at(5),
        model_list.at(5),
        1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest,
    GroupMedianModelListIgnoresInvalidMembersAndFallsBackToSnapshot)
{
    const rg::GaussianModel3D valid_model{ 6.0, 0.55, 0.2 };
    const rg::GaussianModel3D invalid_model{ -1.0, 0.60, 9.0 };
    const std::vector<rg::GaussianModel3D> model_list{
        valid_model,
        invalid_model,
        invalid_model
    };
    const auto median_model_list{
        median_detail::BuildGroupMedianModelList(
            std::vector<std::size_t>{ 10, 10, 20 },
            model_list)
    };

    ASSERT_EQ(median_model_list.size(), model_list.size());
    ExpectGaussianModelsNear(median_model_list.at(0), valid_model, 1.0e-12);
    ExpectGaussianModelsNear(median_model_list.at(1), valid_model, 1.0e-12);
    ExpectGaussianModelsNear(median_model_list.at(2), invalid_model, 1.0e-12);
    EXPECT_THROW(
        median_detail::BuildGroupMedianModelList(
            std::vector<std::size_t>{ 10 },
            model_list),
        std::invalid_argument);
}

TEST(EstimatorSecondStageDefenseTest,
    SharedOffsetDampedModelsInterpolatePhysicalGroupOffset)
{
    const std::vector<std::size_t> group_id_list{ 10, 10, 20 };
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
    const auto previous_shared_offset_list{
        median_detail::BuildGroupMedianOffsetList(
            group_id_list,
            previous_model_list)
    };
    const auto raw_shared_offset_list{
        median_detail::BuildGroupMedianOffsetList(
            group_id_list,
            raw_model_list)
    };

    std::vector<rg::GaussianModel3D> candidate_model_list;
    for (const auto damping : std::array<double, 3>{ 0.0, 0.25, 1.0 })
    {
        ASSERT_TRUE(median_detail::TryBuildSharedOffsetDampedModelList(
            previous_model_list,
            raw_model_list,
            previous_shared_offset_list,
            raw_shared_offset_list,
            damping,
            candidate_model_list));
        const auto expected_first_group_offset{ 0.5 + damping * 0.1 };
        const auto expected_second_group_offset{ -1.0 + damping * 3.0 };
        EXPECT_NEAR(
            candidate_model_list.at(0).GetOffset(),
            expected_first_group_offset,
            1.0e-12);
        EXPECT_NEAR(
            candidate_model_list.at(1).GetOffset(),
            expected_first_group_offset,
            1.0e-12);
        EXPECT_NEAR(
            candidate_model_list.at(2).GetOffset(),
            expected_second_group_offset,
            1.0e-12);

        for (std::size_t atom_position = 0;
            atom_position < candidate_model_list.size();
            atom_position++)
        {
            const auto previous_coordinates{
                change_detail::EncodeTransformedCoordinates(
                    previous_model_list.at(atom_position))
            };
            const auto raw_coordinates{
                change_detail::EncodeTransformedCoordinates(
                    raw_model_list.at(atom_position))
            };
            const auto candidate_coordinates{
                change_detail::EncodeTransformedCoordinates(
                    candidate_model_list.at(atom_position))
            };
            ASSERT_TRUE(previous_coordinates.has_value());
            ASSERT_TRUE(raw_coordinates.has_value());
            ASSERT_TRUE(candidate_coordinates.has_value());
            for (const auto parameter_index : std::array<std::size_t, 2>{
                change_detail::kLogPeakHeightChangeIndex,
                change_detail::kLogWidthChangeIndex })
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

TEST(EstimatorSecondStageDefenseTest, GroupMedianModelsAreIntensityScaleInvariant)
{
    constexpr double scale{ 100.0 };
    const std::vector<std::size_t> group_id_list{ 10, 10, 10 };
    const std::vector<rg::GaussianModel3D> model_list{
        rg::GaussianModel3D{ 3.0, 0.40, -0.2 },
        rg::GaussianModel3D{ 5.0, 0.60, 0.1 },
        rg::GaussianModel3D{ 7.0, 0.80, 0.4 }
    };
    std::vector<rg::GaussianModel3D> scaled_model_list;
    scaled_model_list.reserve(model_list.size());
    for (const auto & model : model_list)
    {
        scaled_model_list.emplace_back(rg::GaussianModel3D{
            scale * model.GetAmplitude(),
            model.GetWidth(),
            scale * model.GetOffset()
        });
    }

    const auto median_model_list{
        median_detail::BuildGroupMedianModelList(
            group_id_list,
            model_list)
    };
    const auto scaled_median_model_list{
        median_detail::BuildGroupMedianModelList(
            group_id_list,
            scaled_model_list)
    };
    ASSERT_EQ(median_model_list.size(), scaled_median_model_list.size());
    for (std::size_t atom_position = 0;
        atom_position < median_model_list.size();
        atom_position++)
    {
        EXPECT_DOUBLE_EQ(
            scale * median_model_list.at(atom_position).GetAmplitude(),
            scaled_median_model_list.at(atom_position).GetAmplitude());
        EXPECT_DOUBLE_EQ(
            median_model_list.at(atom_position).GetWidth(),
            scaled_median_model_list.at(atom_position).GetWidth());
        EXPECT_DOUBLE_EQ(
            scale * median_model_list.at(atom_position).GetOffset(),
            scaled_median_model_list.at(atom_position).GetOffset());
    }
}

TEST(EstimatorSecondStageDefenseTest,
    IndividualRefitUsesGroupMedianModelForTargetOffsetResponse)
{
    const auto median_model_list{
        median_detail::BuildGroupMedianModelList(
            std::vector<std::size_t>{ 10, 10, 10 },
            std::vector<rg::GaussianModel3D>{
                rg::GaussianModel3D{ 3.0, 0.70, 0.2 },
                rg::GaussianModel3D{ 5.0, 0.80, 0.3 },
                rg::GaussianModel3D{ 7.0, 0.90, 0.4 }
            })
    };
    ASSERT_EQ(median_model_list.size(), 3U);
    const auto & offset_model{ median_model_list.front() };
    const rg::GaussianModel3D truth_shape{ 6.0, 0.55, 0.0 };
    LocalPotentialSampleList sample_list;
    for (const auto distance : std::array<float, 5>{
        0.0F, 0.15F, 0.30F, 0.45F, 0.60F })
    {
        sample_list.emplace_back(LocalPotentialSample{
            static_cast<float>(
                truth_shape.SignalAtDistance(distance) +
                offset_model.GetOffset() *
                    offset_model.OffsetBasisAtDistance(distance)),
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
    JointPolishParameterizationSharesOneOffsetColumnPerGroup)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 1.0 },
        rg::GaussianModel3D{ 7.0, 0.60, 4.0 },
        rg::GaussianModel3D{ 8.0, 0.65, 3.0 }
    };
    const auto parameterization{
        polish_detail::BuildJointPolishParameterization(
            std::vector<std::size_t>{ 20, 10, 20 },
            base_model_list)
    };

    ASSERT_TRUE(parameterization.has_value());
    EXPECT_EQ(parameterization->group_position_by_atom.size(), 3U);
    EXPECT_EQ(parameterization->seed_parameter.size(), 8);
    EXPECT_EQ(parameterization->ShapeColumn(0, 0), 0);
    EXPECT_EQ(parameterization->ShapeColumn(1, 0), 2);
    EXPECT_EQ(parameterization->ShapeColumn(2, 0), 4);
    EXPECT_EQ(
        parameterization->OffsetColumn(0),
        parameterization->OffsetColumn(2));
    EXPECT_NE(
        parameterization->OffsetColumn(0),
        parameterization->OffsetColumn(1));
}

TEST(EstimatorSecondStageDefenseTest, JointPolishSharedOffsetSeedUsesGroupMedians)
{
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 1.0 },
        rg::GaussianModel3D{ 7.0, 0.60, 4.0 },
        rg::GaussianModel3D{ 8.0, 0.65, 3.0 },
        rg::GaussianModel3D{ 9.0, 0.70, 2.0 },
        rg::GaussianModel3D{ 10.0, 0.75, 6.0 }
    };
    const auto parameterization{
        polish_detail::BuildJointPolishParameterization(
            std::vector<std::size_t>{ 20, 10, 20, 20, 10 },
            base_model_list)
    };

    ASSERT_TRUE(parameterization.has_value());
    EXPECT_EQ(parameterization->seed_parameter.size(), 12);
    EXPECT_DOUBLE_EQ(
        parameterization->seed_parameter(
            parameterization->OffsetColumn(0)),
        2.0);
    EXPECT_DOUBLE_EQ(
        parameterization->seed_parameter(
            parameterization->OffsetColumn(1)),
        5.0);

    const auto seed_model_list{ parameterization->DecodeSeedModels() };
    const auto zero_direction{
        Eigen::VectorXd::Zero(parameterization->seed_parameter.size())
    };
    const auto zero_model_list{
        parameterization->DecodeModels(zero_direction, 0.0)
    };
    ASSERT_TRUE(seed_model_list.has_value());
    ASSERT_TRUE(zero_model_list.has_value());
    ASSERT_EQ(seed_model_list->size(), zero_model_list->size());
    for (std::size_t atom_position = 0;
        atom_position < seed_model_list->size();
        atom_position++)
    {
        ExpectGaussianModelsNear(
            seed_model_list->at(atom_position),
            zero_model_list->at(atom_position),
            1.0e-12);
    }

    Eigen::VectorXd direction{
        Eigen::VectorXd::Zero(parameterization->seed_parameter.size())
    };
    direction(parameterization->OffsetColumn(1)) = 2.0;
    direction(parameterization->OffsetColumn(0)) = -1.0;
    const auto candidate_model_list{
        parameterization->DecodeModels(direction, 0.5)
    };
    ASSERT_TRUE(candidate_model_list.has_value());
    EXPECT_DOUBLE_EQ(candidate_model_list->at(0).GetOffset(), 1.5);
    EXPECT_DOUBLE_EQ(candidate_model_list->at(1).GetOffset(), 6.0);
    EXPECT_DOUBLE_EQ(candidate_model_list->at(2).GetOffset(), 1.5);
    EXPECT_DOUBLE_EQ(candidate_model_list->at(3).GetOffset(), 1.5);
    EXPECT_DOUBLE_EQ(candidate_model_list->at(4).GetOffset(), 6.0);
    EXPECT_FALSE(
        polish_detail::BuildJointPolishParameterization(
            std::vector<std::size_t>{ 20 },
            base_model_list).has_value());
}

TEST(
    EstimatorSecondStageDefenseTest,
    JointPolishDirectionAndProposalShareGroupOffset)
{
    const std::vector<std::size_t> group_id_list{ 20, 20 };
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.10 }
    };
    const std::vector<rg::GaussianModel3D> target_model_list{
        rg::GaussianModel3D{ 6.5, 0.60, 0.30 },
        rg::GaussianModel3D{ 4.0, 0.65, 0.30 }
    };
    auto fixture{
        BuildJointPolishFixture(
            group_id_list,
            base_model_list,
            target_model_list)
    };
    const auto parameterization{
        polish_detail::BuildJointPolishParameterization(
            group_id_list,
            base_model_list)
    };
    ASSERT_TRUE(parameterization.has_value());
    polish_detail::ReusableWeightedRidgeSolver direction_solver;
    const polish_detail::FitStatePatch base_patch;
    const polish_detail::FitStateView base_state_view{ fixture.state, base_patch };
    const auto direction{
        polish_detail::BuildJointPolishDirection(
            fixture.context,
            base_state_view,
            polish_detail::ClusterKey{ 0, 1 },
            fixture.sample_ref_list,
            { 1.0, 1.0 },
            *parameterization,
            direction_solver)
    };
    ASSERT_TRUE(direction.has_value());
    EXPECT_TRUE(direction->allFinite());
    EXPECT_GT(direction->norm(), 1.0e-8);

    polish_detail::ReusableWeightedRidgeSolver proposal_solver;
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
    EXPECT_DOUBLE_EQ(
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
    JointPolishProposalRejectsEmptyInvalidUnchangedAndOutOfRegionInputs)
{
    const std::vector<std::size_t> group_id_list{ 20, 20 };
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.10 }
    };
    const std::vector<rg::GaussianModel3D> target_model_list{
        rg::GaussianModel3D{ 6.5, 0.60, 0.30 },
        rg::GaussianModel3D{ 4.0, 0.65, 0.30 }
    };
    auto fixture{
        BuildJointPolishFixture(
            group_id_list,
            base_model_list,
            target_model_list)
    };
    const polish_detail::FitStatePatch base_patch;
    const polish_detail::FitStateView base_state_view{ fixture.state, base_patch };
    const auto key{ polish_detail::ClusterKey{ 0, 1 } };

    polish_detail::ReusableWeightedRidgeSolver empty_solver;
    EXPECT_FALSE(
        polish_detail::BuildJointPolishProposal(
            fixture.context,
            base_state_view,
            key,
            {},
            { 1.0, 1.0 },
            empty_solver,
            4.0).has_value());

    auto invalid_fixture{ fixture };
    invalid_fixture.state.at(0).mdpde =
        rg::GaussianModel3DWithUncertainty{
            rg::GaussianModel3D{ 0.0, 0.55, 0.10 },
            rg::GaussianModel3DUncertainty{}
        };
    const polish_detail::FitStatePatch invalid_patch;
    const polish_detail::FitStateView invalid_state_view{
        invalid_fixture.state,
        invalid_patch
    };
    polish_detail::ReusableWeightedRidgeSolver invalid_solver;
    EXPECT_FALSE(
        polish_detail::BuildJointPolishProposal(
            invalid_fixture.context,
            invalid_state_view,
            key,
            invalid_fixture.sample_ref_list,
            { 1.0, 1.0 },
            invalid_solver,
            4.0).has_value());

    const std::vector<rg::GaussianModel3D> unchanged_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, 0.10 }
    };
    auto unchanged_fixture{
        BuildJointPolishFixture(
            group_id_list,
            unchanged_model_list,
            unchanged_model_list)
    };
    const polish_detail::FitStatePatch unchanged_patch;
    const polish_detail::FitStateView unchanged_state_view{
        unchanged_fixture.state,
        unchanged_patch
    };
    polish_detail::ReusableWeightedRidgeSolver unchanged_solver;
    EXPECT_FALSE(
        polish_detail::BuildJointPolishProposal(
            unchanged_fixture.context,
            unchanged_state_view,
            key,
            unchanged_fixture.sample_ref_list,
            { 1.0, 1.0 },
            unchanged_solver,
            4.0).has_value());

    polish_detail::ReusableWeightedRidgeSolver trust_region_solver;
    EXPECT_FALSE(
        polish_detail::BuildJointPolishProposal(
            fixture.context,
            base_state_view,
            key,
            fixture.sample_ref_list,
            { 1.0, 1.0 },
            trust_region_solver,
            0.01).has_value());
}

TEST(
    EstimatorSecondStageDefenseTest,
    JointPolishProposalUsesFitStateViewBaseForTrustRegionOrigin)
{
    const std::vector<std::size_t> group_id_list{ 20, 20 };
    const std::vector<rg::GaussianModel3D> base_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.10 },
        rg::GaussianModel3D{ 4.5, 0.70, -0.10 }
    };
    const std::vector<rg::GaussianModel3D> target_model_list{
        rg::GaussianModel3D{ 6.5, 0.60, 0.30 },
        rg::GaussianModel3D{ 4.0, 0.65, 0.30 }
    };
    auto fixture{
        BuildJointPolishFixture(
            group_id_list,
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
        polish_detail::BuildJointPolishParameterization(
            group_id_list,
            std::vector<rg::GaussianModel3D>{
                base_state_view.GetModel(0),
                base_state_view.GetModel(1) })
    };
    ASSERT_TRUE(patched_parameterization.has_value());
    const auto patched_seed{ patched_parameterization->DecodeSeedModels() };
    ASSERT_TRUE(patched_seed.has_value());
    EXPECT_DOUBLE_EQ(patched_seed->at(0).GetAmplitude(), 10.0);

    polish_detail::ReusableWeightedRidgeSolver solver;
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

TEST(EstimatorSecondStageDefenseTest, SharedOffsetJacobianMatchesFiniteDifference)
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
        const auto transformed{
            change_detail::EncodeTransformedCoordinates(model)
        };
        ASSERT_TRUE(transformed.has_value());
        for (const auto distance : distance_list)
        {
            const auto evaluation{
        polish_detail::EvaluateSharedOffsetResponse(
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
                    change_detail::kOffsetToPeakRatioChangeIndex)) = 0.0;
                upper(static_cast<Eigen::Index>(
                    change_detail::kOffsetToPeakRatioChangeIndex)) = 0.0;
                const auto lower_shape{
                    change_detail::DecodeTransformedCoordinates(lower)
                };
                const auto upper_shape{
                    change_detail::DecodeTransformedCoordinates(upper)
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
    const auto topology{ builder.BuildTopology(MakeUniqueResidueKeys(2)) };
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
        scaled_builder.BuildTopology(MakeUniqueResidueKeys(2))
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
        tiny_builder.BuildTopology(MakeUniqueResidueKeys(2))
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

    const auto topology{ builder.BuildTopology(MakeUniqueResidueKeys(2)) };
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
    EXPECT_EQ(topology.residue_cutoff_summary.maximum_residue_count_limit, 10U);
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

    const auto topology{ builder.BuildTopology(MakeUniqueResidueKeys(2)) };
    EXPECT_FALSE(topology.summary.uses_weighted_graph);
    EXPECT_TRUE(topology.summary.threshold_sensitivity_list.empty());
    ASSERT_EQ(topology.sample_dependency_list.size(), 1U);
    EXPECT_EQ(
        topology.sample_dependency_list.front().contributor_atom_index_list,
        (std::vector<std::size_t>{ 0, 1 }));
    EXPECT_TRUE(HasCouplingNeighbor(topology, 0, 1));
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphValidatesSampleAndBuildOptions)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    coupling_detail::CouplingGraphBuilder builder{ 2 };
    EXPECT_THROW(
        AddCouplingGraphSample(builder, { 0, 0 }, { { 2, unit } }),
        std::invalid_argument);

    AddCouplingGraphSample(builder, { 0, 0 }, { { 0, unit } });
    auto invalid_weight_options{
        coupling_detail::CouplingGraphOptions{}
    };
    invalid_weight_options.minimum_weight = -0.01;
    EXPECT_THROW(
        builder.BuildTopology(MakeUniqueResidueKeys(2), invalid_weight_options),
        std::invalid_argument);

    auto invalid_residue_limit_options{
        coupling_detail::CouplingGraphOptions{}
    };
    invalid_residue_limit_options.maximum_residue_count = 0;
    EXPECT_THROW(
        builder.BuildTopology(
            MakeUniqueResidueKeys(2),
            invalid_residue_limit_options),
        std::invalid_argument);
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphSummaryIncludesResidueComponents)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    coupling_detail::CouplingGraphBuilder builder{ 2 };
    AddCouplingGraphSample(builder, { 0, 0 }, { { 0, unit } });
    AddCouplingGraphSample(builder, { 1, 0 }, { { 1, unit } });

    const auto topology{
        builder.BuildTopology({ { "A", 1 }, { "A", 1 } })
    };
    EXPECT_TRUE(topology.adjacency_list.at(0).empty());
    EXPECT_TRUE(topology.adjacency_list.at(1).empty());
    EXPECT_EQ(topology.summary.component_count, 1U);
    EXPECT_EQ(topology.summary.maximum_component_size, 2U);
    EXPECT_DOUBLE_EQ(topology.summary.maximum_component_ratio, 1.0);
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphCutsWeakAndCancelledEdges)
{
    const Eigen::Vector3d unit{ 1.0, 0.0, 0.0 };
    coupling_detail::CouplingGraphBuilder weak_builder{ 2 };
    AddCouplingGraphSample(weak_builder, { 0, 0 }, { { 0, unit }, { 1, unit } });
    AddCouplingGraphSample(weak_builder, { 0, 1 }, { { 0, 10.0 * unit } });
    AddCouplingGraphSample(weak_builder, { 1, 0 }, { { 1, 10.0 * unit } });
    const auto weak_topology{
        weak_builder.BuildTopology(MakeUniqueResidueKeys(2))
    };
    EXPECT_FALSE(HasCouplingNeighbor(weak_topology, 0, 1));
    EXPECT_EQ(weak_topology.summary.candidate_edge_count, 1U);
    EXPECT_EQ(weak_topology.summary.cut_edge_count, 1U);

    coupling_detail::CouplingGraphBuilder cancelled_builder{ 2 };
    AddCouplingGraphSample(cancelled_builder, { 0, 0 }, { { 0, unit }, { 1, unit } });
    AddCouplingGraphSample(cancelled_builder, { 0, 1 }, { { 0, unit }, { 1, -unit } });
    const auto cancelled_topology{
        cancelled_builder.BuildTopology(MakeUniqueResidueKeys(2))
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
            MakeUniqueResidueKeys(2),
            options,
            previous_topology);
    };

    coupling_detail::GraphTopology absent_previous;
    absent_previous.adjacency_list.resize(2);
    const auto absent_midpoint{ build_topology(0.05, &absent_previous) };
    EXPECT_FALSE(HasCouplingNeighbor(absent_midpoint, 0, 1));
    const auto added{ build_topology(0.061, &absent_previous) };
    EXPECT_TRUE(HasCouplingNeighbor(added, 0, 1));
    const auto retained_midpoint{ build_topology(0.05, &added) };
    EXPECT_TRUE(HasCouplingNeighbor(retained_midpoint, 0, 1));
    const auto removed{ build_topology(0.039, &retained_midpoint) };
    EXPECT_FALSE(HasCouplingNeighbor(removed, 0, 1));
}

TEST(EstimatorSecondStageDefenseTest, CouplingGraphValidatesAdaptiveHysteresisInputs)
{
    coupling_detail::CouplingGraphBuilder builder{ 1 };
    coupling_detail::CouplingGraphOptions options;
    options.minimum_weight = 0.05;
    options.retained_edge_minimum_weight = 0.051;
    EXPECT_THROW(
        builder.BuildTopology(MakeUniqueResidueKeys(1), options),
        std::invalid_argument);

    coupling_detail::CouplingGraphBuilder previous_size_builder{ 1 };
    options.retained_edge_minimum_weight = 0.04;
    coupling_detail::GraphTopology wrong_size_previous;
    wrong_size_previous.adjacency_list.resize(2);
    EXPECT_THROW(
        previous_size_builder.BuildTopology(
            MakeUniqueResidueKeys(1),
            options,
            &wrong_size_previous),
        std::invalid_argument);
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
    const auto topology{
        builder.BuildTopology(MakeUniqueResidueKeys(7), options)
    };
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
    EXPECT_EQ(formal_threshold.component_count, formal_partition.sample_id_list_by_key.size());
    EXPECT_EQ(formal_threshold.maximum_component_size, 2U);
    EXPECT_EQ(topology.summary.component_count, 4U);
    EXPECT_EQ(topology.summary.maximum_component_size, 2U);
    EXPECT_NEAR(topology.summary.maximum_component_ratio, 2.0 / 7.0, 1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, CouplingPartitionCutsWeakBridgeAndDuplicatesBoundarySample)
{
    coupling_detail::GraphTopology topology;
    topology.adjacency_list.resize(3);
    topology.adjacency_list.at(0).push_back(1);
    topology.adjacency_list.at(1).push_back(0);
    topology.sample_dependency_list = {
        { { 0, 0 }, { 0, 1 } },
        { { 1, 0 }, { 1, 2 } }
    };

    const auto partition{
        coupling_detail::BuildGraphPartition(topology, { 0, 1, 2 })
    };
    ASSERT_EQ(partition.sample_id_list_by_key.size(), 2U);
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 0, 1 }), 1U);
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 2 }), 1U);
    EXPECT_EQ(partition.boundary_sample_count, 1U);
    ASSERT_EQ(partition.boundary_sample_dependency_list.size(), 1U);
    EXPECT_EQ(
        partition.boundary_sample_dependency_list.front().sample_id,
        (coupling_detail::SampleRef{ 1, 0 }));
    EXPECT_EQ(
        partition.boundary_sample_dependency_list.front().cluster_key_list,
        (std::vector<audit_detail::ClusterKey>{ { 0, 1 }, { 2 } }));
    EXPECT_EQ(partition.sample_id_list_by_key.at({ 0, 1 }).size(), 2U);
    EXPECT_EQ(partition.sample_id_list_by_key.at({ 2 }).size(), 1U);

    const auto key_list{
        coupling_detail::BuildGraphClusterKeyList(partition)
    };
    EXPECT_EQ(key_list, (std::vector<std::vector<std::size_t>>{ { 0, 1 }, { 2 } }));
    const auto affected_sample_list{
        coupling_detail::BuildGraphAffectedSampleUnion(
            partition,
            key_list)
    };
    EXPECT_EQ(affected_sample_list.size(), 2U);
    EXPECT_EQ(affected_sample_list.at(0).atom_index, 0U);
    EXPECT_EQ(affected_sample_list.at(0).sample_index, 0U);
    EXPECT_EQ(affected_sample_list.at(1).atom_index, 1U);
    EXPECT_EQ(affected_sample_list.at(1).sample_index, 0U);

    const auto inactive_partition{
        coupling_detail::BuildGraphPartition(topology, { 2, 0 })
    };
    EXPECT_EQ(inactive_partition.sample_id_list_by_key.count({ 0 }), 1U);
    EXPECT_EQ(inactive_partition.sample_id_list_by_key.count({ 2 }), 1U);
    EXPECT_EQ(inactive_partition.boundary_sample_count, 0U);
}

TEST(EstimatorSecondStageDefenseTest, BoundaryReconciliationComponentsUseAcceptedSharedSamples)
{
    const audit_detail::ClusterKey key_a{ 0 };
    const audit_detail::ClusterKey key_b{ 1 };
    const audit_detail::ClusterKey key_c{ 2 };
    const audit_detail::ClusterKey key_d{ 3 };
    const audit_detail::ClusterKey key_e{ 4 };
    const coupling_detail::SampleRef sample_ab{ 0, 0 };
    const coupling_detail::SampleRef sample_bc{ 1, 0 };
    const coupling_detail::SampleRef sample_de{ 3, 0 };
    coupling_detail::CouplingGraphPartition partition;
    partition.sample_id_list_by_key = {
        { key_a, { sample_ab } },
        { key_b, { sample_ab, sample_bc } },
        { key_c, { sample_bc } },
        { key_d, { sample_de } },
        { key_e, { sample_de } }
    };
    partition.boundary_sample_dependency_list = {
        { sample_ab, { key_a, key_b } },
        { sample_bc, { key_b, key_c } },
        { sample_de, { key_d, key_e } }
    };
    partition.boundary_sample_count =
        partition.boundary_sample_dependency_list.size();

    const auto component_list{
        coupling_detail::BuildBoundaryReconciliationComponents(
            partition,
            { key_e, key_c, key_a, key_d, key_b })
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
        component_list.at(1).key_list,
        (std::vector<audit_detail::ClusterKey>{ key_d, key_e }));
    EXPECT_EQ(
        component_list.at(1).affected_sample_ref_list,
        (std::vector<coupling_detail::SampleRef>{ sample_de }));
    EXPECT_EQ(component_list.at(1).boundary_sample_count, 1U);

    EXPECT_TRUE(
        coupling_detail::BuildBoundaryReconciliationComponents(
            partition,
            { key_c, key_a }).empty());
    EXPECT_EQ(
        coupling_detail::BuildBoundaryReconciliationComponents(
            partition,
            { key_c, key_b, key_a }),
        std::vector<coupling_detail::BoundaryReconciliationComponent>{
            component_list.front()
        });
}

TEST(EstimatorSecondStageDefenseTest, CouplingResidueCutoffKeepsResiduesWholeAndBounded)
{
    coupling_detail::GraphTopology topology;
    topology.adjacency_list.resize(13);
    for (std::size_t atom_index = 1; atom_index < 12; atom_index++)
    {
        topology.retained_edge_list.emplace_back(
            coupling_detail::GraphWeightedEdge{
                atom_index,
                atom_index + 1,
                1.0 - 0.01 * static_cast<double>(atom_index)
            });
    }
    topology.sample_dependency_list = {
        { { 0, 0 }, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 } }
    };
    std::vector<coupling_detail::ResidueKey> residue_key_list{
        { "A", 1 },
        { "A", 1 },
        { "B", 1 }
    };
    for (int residue_id = 2; residue_id <= 11; residue_id++)
    {
        residue_key_list.emplace_back("A", residue_id);
    }
    topology.residue_key_by_atom_index = residue_key_list;

    const auto capped_topology{
        coupling_detail::ApplyGraphResidueCutoff(
            std::move(topology),
            10)
    };
    EXPECT_EQ(capped_topology.residue_cutoff_summary.residue_count, 12U);
    EXPECT_GE(capped_topology.residue_cutoff_summary.cluster_count, 2U);
    EXPECT_LE(capped_topology.residue_cutoff_summary.maximum_residue_count, 10U);
    EXPECT_GT(capped_topology.residue_cutoff_summary.cut_edge_count, 0U);

    std::vector<std::size_t> active_index_list(13);
    for (std::size_t atom_index = 0; atom_index < active_index_list.size(); atom_index++)
    {
        active_index_list.at(atom_index) = atom_index;
    }
    const auto partition{
        coupling_detail::BuildGraphPartition(
            capped_topology,
            active_index_list)
    };
    EXPECT_EQ(partition.boundary_sample_count, 1U);
    EXPECT_GE(partition.sample_id_list_by_key.size(), 2U);

    bool found_first_residue{ false };
    for (const auto & [key, sample_id_list] : partition.sample_id_list_by_key)
    {
        EXPECT_EQ(sample_id_list.size(), 1U);
        std::set<coupling_detail::ResidueKey> residue_key_set;
        for (const auto atom_index : key)
        {
            residue_key_set.emplace(residue_key_list.at(atom_index));
        }
        EXPECT_LE(residue_key_set.size(), 10U);
        if (std::find(key.begin(), key.end(), 0U) == key.end()) continue;
        found_first_residue = true;
        EXPECT_NE(std::find(key.begin(), key.end(), 1U), key.end());
    }
    EXPECT_TRUE(found_first_residue);

    std::reverse(active_index_list.begin(), active_index_list.end());
    const auto reversed_partition{
        coupling_detail::BuildGraphPartition(
            capped_topology,
            active_index_list)
    };
    EXPECT_EQ(
        partition.sample_id_list_by_key.size(),
        reversed_partition.sample_id_list_by_key.size());
    auto expected_iter{ partition.sample_id_list_by_key.begin() };
    auto actual_iter{ reversed_partition.sample_id_list_by_key.begin() };
    for (; expected_iter != partition.sample_id_list_by_key.end();
        expected_iter++, actual_iter++)
    {
        EXPECT_EQ(expected_iter->first, actual_iter->first);
    }
}

TEST(EstimatorSecondStageDefenseTest, CouplingResidueCutoffPrioritizesStrongEdges)
{
    coupling_detail::GraphTopology topology;
    topology.adjacency_list.resize(3);
    topology.retained_edge_list = {
        { 1, 2, 0.80 },
        { 0, 1, 0.90 }
    };
    topology.residue_key_by_atom_index = {
        { "A", 1 }, { "A", 2 }, { "A", 3 }
    };
    const auto capped_topology{
        coupling_detail::ApplyGraphResidueCutoff(
            std::move(topology),
            2)
    };
    const auto partition{
        coupling_detail::BuildGraphPartition(
            capped_topology,
            { 2, 1, 0 })
    };
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 0, 1 }), 1U);
    EXPECT_EQ(partition.sample_id_list_by_key.count({ 2 }), 1U);
    EXPECT_TRUE(HasCouplingNeighbor(capped_topology, 0, 1));
    EXPECT_FALSE(HasCouplingNeighbor(capped_topology, 1, 2));
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
    auto binary_topology{
        builder.BuildTopology(MakeUniqueResidueKeys(2), fallback_options)
    };
    EXPECT_FALSE(binary_topology.summary.uses_weighted_graph);
    EXPECT_TRUE(binary_topology.summary.threshold_sensitivity_list.empty());
    const auto binary_partition{
        coupling_detail::BuildGraphPartition(
            binary_topology,
            { 0, 1 })
    };
    EXPECT_EQ(binary_partition.sample_id_list_by_key.count({ 0, 1 }), 1U);
    binary_topology.residue_key_by_atom_index = { { "A", 1 }, { "A", 2 } };
    const auto capped_binary_topology{
        coupling_detail::ApplyGraphResidueCutoff(
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
    const auto overflow_topology{
        overflow_builder.BuildTopology(MakeUniqueResidueKeys(2))
    };
    EXPECT_FALSE(overflow_topology.summary.uses_weighted_graph);
}

TEST(EstimatorSecondStageDefenseTest, TransformedCoordinatesRoundTrip)
{
    const rg::GaussianModel3D model{ 8.5, 0.65, -0.12 };
    const auto encoded{
        change_detail::EncodeTransformedCoordinates(model)
    };
    ASSERT_TRUE(encoded.has_value());

    const auto decoded{
        change_detail::DecodeTransformedCoordinates(*encoded)
    };
    ASSERT_TRUE(decoded.has_value());
    ExpectGaussianModelsNear(model, *decoded, 1.0e-12);
}

TEST(EstimatorSecondStageDefenseTest, TransformedDampingIsIntensityScaleInvariant)
{
    const rg::GaussianModel3D previous{ 8.0, 0.50, -0.10 };
    const rg::GaussianModel3D current{ 9.0, 0.60, -0.15 };
    constexpr double damping{ 0.25 };

    const auto damp = [&](const rg::GaussianModel3D & lhs,
                          const rg::GaussianModel3D & rhs)
    {
        const auto lhs_coordinates{
            change_detail::EncodeTransformedCoordinates(lhs)
        };
        const auto rhs_coordinates{
            change_detail::EncodeTransformedCoordinates(rhs)
        };
        EXPECT_TRUE(lhs_coordinates.has_value());
        EXPECT_TRUE(rhs_coordinates.has_value());
        return change_detail::DecodeTransformedCoordinates(
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

    backtracking_detail::SecondStageContext context;
    context.selected_atom_list.resize(previous_model_list.size());
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
        context,
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

    const auto candidate_state{
        backtracking_detail::FitStateView{ previous_state, candidate_patch }.Materialize()
    };
    for (const auto atom_index : candidate_patch.atom_index_list)
    {
        ExpectGaussianModelsNear(
            candidate_state.at(atom_index).mdpde.GetModel(),
            candidate_patch.mdpde_list.at(atom_index).GetModel(),
            1.0e-12);
    }
    EXPECT_DOUBLE_EQ(
        candidate_state.at(0).mdpde.GetModel().GetOffset(),
        candidate_state.at(1).mdpde.GetModel().GetOffset());

    const auto merged_provenance{
        workspace.BuildCandidatePolishProvenance(
            std::vector<char>{ 0, 1 },
            std::vector<char>{ 1, 0 })
    };
    EXPECT_EQ(merged_provenance, (std::vector<char>{ 1, 0 }));
}

TEST(EstimatorSecondStageDefenseTest,
    BacktrackingWorkspaceDistinguishesInvalidAndExhaustedSteps)
{
    backtracking_detail::SecondStageContext context;
    context.selected_atom_list.resize(1);
    backtracking_detail::FitState previous_state(1);
    previous_state.at(0).mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 8.0, 0.50, -0.10 },
        rg::GaussianModel3DUncertainty{ 0.1, 0.02, 0.03 }
    };

    backtracking_detail::FitState endpoint_state(1);
    endpoint_state.at(0).mdpde = rg::GaussianModel3DWithUncertainty{
        rg::GaussianModel3D{ 0.0, 0.0, 0.0 },
        rg::GaussianModel3DUncertainty{ 0.2, 0.04, 0.05 }
    };
    const auto invalid_endpoint_patch{
        backtracking_detail::FitStatePatch::FromState(
            endpoint_state,
            std::vector<std::size_t>{ 0 })
    };
    backtracking_detail::BacktrackingWorkspace invalid_workspace{
        context,
        previous_state,
        invalid_endpoint_patch,
        1.0e-4
    };
    const auto invalid_step{ invalid_workspace.BuildNextCandidate() };
    EXPECT_EQ(
        invalid_step.status,
        backtracking_detail::BacktrackingStepStatus::InvalidCandidate);
    EXPECT_EQ(invalid_step.trial_number, 1U);

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
        context,
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

    backtracking_detail::SecondStageContext empty_context;
    backtracking_detail::FitState empty_state;
    backtracking_detail::BacktrackingWorkspace factor_workspace{
        empty_context,
        empty_state,
        backtracking_detail::FitStatePatch{},
        0.0
    };
    std::size_t ready_count{ 0 };
    double expected_factor{ 0.5 };
    std::size_t expected_trial_number{ 2 };
    backtracking_detail::BacktrackingStep factor_step;
    do
    {
        factor_step = factor_workspace.BuildNextCandidate();
        if (factor_step.status == backtracking_detail::BacktrackingStepStatus::CandidateReady)
        {
            ready_count++;
            EXPECT_DOUBLE_EQ(factor_step.factor, expected_factor);
            EXPECT_EQ(factor_step.trial_number, expected_trial_number);
            expected_factor *= 0.5;
            expected_trial_number++;
        }
    }
    while (factor_step.status == backtracking_detail::BacktrackingStepStatus::CandidateReady);
    EXPECT_EQ(
        factor_step.status,
        backtracking_detail::BacktrackingStepStatus::Exhausted);
    EXPECT_GT(ready_count, 40U);
    EXPECT_LT(
        factor_step.factor,
        std::numeric_limits<double>::epsilon());
}

TEST(EstimatorSecondStageDefenseTest, ResidualBaselineAndOverlayAgreeForCandidate)
{
    residual_detail::SecondStageContext context;
    context.selected_atom_list.resize(1);
    context.selected_atom_index_list_by_group.emplace_back(
        std::vector<std::size_t>{ 0 });
    context.at(0).neighbor_atom_sample_offset_list = { 0, 0, 0 };

    const rg::GaussianModel3D previous_model{ 8.0, 0.50, -0.10 };
    const rg::GaussianModel3D candidate_model{ 10.0, 0.60, 0.20 };
    for (const auto distance : { 0.15F, 0.45F })
    {
        context.at(0).raw_sampling_entries.emplace_back(
            LocalPotentialSample{
                static_cast<float>(candidate_model.ResponseAtDistance(distance)),
                SamplingPoint{ distance }
            });
    }

    residual_detail::FitState previous_state;
    previous_state.emplace_back(MakeGaussianResult(previous_model));
    const auto baseline{
        residual_detail::BuildResidualBaseline(
            context,
            previous_state)
    };
    ASSERT_TRUE(baseline.sample_list.at(0).at(1).has_value());
    EXPECT_DOUBLE_EQ(
        baseline.sample_list.at(0).at(1)->adjusted_response,
        static_cast<double>(context.at(0).raw_sampling_entries.at(1).response));
    EXPECT_NEAR(
        baseline.sample_list.at(0).at(1)->residual,
        baseline.sample_list.at(0).at(1)->adjusted_response -
            previous_model.ResponseAtDistance(0.45),
        1.0e-6);

    residual_detail::FitStatePatch patch;
    patch.atom_index_list = { 0 };
    auto candidate_result{ MakeGaussianResult(candidate_model) };
    patch.mdpde_list.emplace_back(candidate_result.mdpde);
    const residual_detail::FitStateView candidate_view{
        previous_state,
        patch
    };
    const residual_detail::CandidateEvaluationOverlay overlay{
        context,
        baseline,
        candidate_view
    };
    const residual_detail::SampleRef sample_ref{ 0, 1 };
    const auto direct{
        residual_detail::EvaluateResidualSample(
            context,
            candidate_view,
            sample_ref,
            baseline.model_snapshot)
    };
    const auto overlaid{ overlay(sample_ref) };
    ASSERT_TRUE(direct.has_value());
    ASSERT_TRUE(overlaid.has_value());
    EXPECT_DOUBLE_EQ(direct->adjusted_response, overlaid->adjusted_response);
    EXPECT_DOUBLE_EQ(direct->residual, overlaid->residual);
}

TEST(EstimatorSecondStageDefenseTest, AuditObjectiveSourcesAgreeAcrossTailPartitions)
{
    audit_detail::SecondStageContext context;
    context.selected_atom_list.resize(1);
    context.selected_atom_index_list_by_group.emplace_back(
        std::vector<std::size_t>{ 0 });
    context.at(0).neighbor_atom_sample_offset_list = { 0, 0, 0 };

    const rg::GaussianModel3D model{ 8.0, 0.50, -0.10 };
    for (const auto distance : { 0.15F, 0.45F })
    {
        context.at(0).raw_sampling_entries.emplace_back(
            LocalPotentialSample{
                static_cast<float>(model.ResponseAtDistance(distance)),
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
    const auto previous_coordinates{
        change_detail::EncodeTransformedCoordinates(previous)
    };
    const auto endpoint_coordinates{
        change_detail::EncodeTransformedCoordinates(endpoint)
    };
    ASSERT_TRUE(previous_coordinates.has_value());
    ASSERT_TRUE(endpoint_coordinates.has_value());

    double previous_offset_distance{
        std::abs(endpoint.GetOffset() - previous.GetOffset())
    };
    for (const auto factor : { 0.5, 0.25, 0.125 })
    {
        const auto candidate{
            change_detail::DecodeTransformedCoordinates(
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
        change_detail::EncodeTransformedCoordinates(
            rg::GaussianModel3D{ 8.0, 0.50, -0.10 })
    };
    const auto right{
        change_detail::EncodeTransformedCoordinates(
            rg::GaussianModel3D{ 9.0, 0.60, -0.15 })
    };
    ASSERT_TRUE(left.has_value());
    ASSERT_TRUE(right.has_value());

    const auto extrapolated{
        change_detail::DecodeTransformedCoordinates(
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
        change.value_list.at(change_detail::kLogPeakHeightChangeIndex),
        1.0e-12);
    EXPECT_NEAR(
        std::log(2.0),
        change.value_list.at(change_detail::kLogWidthChangeIndex),
        1.0e-12);
    EXPECT_DOUBLE_EQ(
        0.0,
        change.value_list.at(change_detail::kOffsetToPeakRatioChangeIndex));
}

TEST(EstimatorSecondStageDefenseTest, InvalidTransformedCoordinatesProduceInfiniteChange)
{
    const rg::GaussianModel3D valid_model{ 1.0, 0.5, 0.0 };
    const std::array<rg::GaussianModel3D, 4> invalid_model_list{
        rg::GaussianModel3D{ 0.0, 0.5, 0.0 },
        rg::GaussianModel3D{ 1.0, 0.0, 0.0 },
        rg::GaussianModel3D{
            std::numeric_limits<double>::quiet_NaN(),
            0.5,
            0.0
        },
        rg::GaussianModel3D{
            std::numeric_limits<double>::denorm_min(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max()
        }
    };

    for (const auto & invalid_model : invalid_model_list)
    {
        EXPECT_FALSE(seed_detail::IsValidSecondStageGaussianModel(invalid_model));
        const auto change{
            change_detail::CalculateTransformedChange(
                invalid_model,
                valid_model)
        };
        ASSERT_EQ(change_detail::kTransformedChangeSize, change.value_list.size());
        for (const auto value : change.value_list)
        {
            EXPECT_TRUE(std::isinf(value));
        }
    }
}

TEST(EstimatorSecondStageDefenseTest, TransformedConvergenceRejectsHiddenMaximumTail)
{
    std::vector<alg::ParameterChange> change_list(
        1000,
        alg::ParameterChange{ std::vector<double>(3, 0.0) });
    change_list.back().value_list.at(change_detail::kLogPeakHeightChangeIndex) =
        2.0e-3;
    std::vector<std::size_t> index_list(change_list.size());
    for (std::size_t i = 0; i < index_list.size(); i++)
    {
        index_list.at(i) = i;
    }

    const auto percentile_stats{
        alg::SummarizeParameterChangeStats(change_list, index_list, 0.99)
    };
    const auto maximum_list{
        change_detail::SummarizeMaximumTransformedChanges(
            change_list,
            index_list)
    };
    EXPECT_LT(
        percentile_stats.percentile_list.at(
            change_detail::kLogPeakHeightChangeIndex),
        1.0e-4);
    EXPECT_FALSE(change_detail::IsTransformedChangeConverged(
        percentile_stats,
        maximum_list));
}

TEST(EstimatorSecondStageDefenseTest, PostRefitRollbackRestoresCompleteLongChain)
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
        ExpectGaussianModelsNear(
            GetEstimateModel(*selected_atoms.at(i)),
            previous_model_list.at(i),
            1.0e-12);
    }
}

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingFallsBackWhenJointOffsetSamplesAreNonFinite)
{
    auto model{ BuildNonFiniteJointOffsetDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    const auto previous_model{ GetEstimateModel(*atom) };
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    ExpectGaussianModelsNear(GetEstimateModel(*atom), previous_model, 1.0e-12);
    ExpectSelectedAtomEstimatesAreFinite(*model);
    EXPECT_NE(out.find("objective-unavailable"), std::string::npos);
    EXPECT_EQ(out.find("objective = not-evaluated"), std::string::npos);
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

TEST(EstimatorSecondStageDefenseTest, LocalRefitFallbackDoesNotBlockRemoteCluster)
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

TEST(EstimatorSecondStageDefenseTest, TerminalFallbackPreservesAffectedCluster)
{
    auto model{ BuildTerminalWithPersistentLocalRefitFallbackDefenseModel() };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    const std::array<rg::GaussianModel3D, 2> previous_terminal_model_list{
        GetEstimateModel(*selected_atoms.at(0)),
        GetEstimateModel(*selected_atoms.at(1))
    };
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };
    const std::string err{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    for (std::size_t i = 0; i < previous_terminal_model_list.size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*selected_atoms.at(i)),
            previous_terminal_model_list.at(i),
            1.0e-12);
    }
    EXPECT_NE(
        out.find("Reset second-stage objective domain"),
        std::string::npos);
    EXPECT_NE(
        err.find("offsets finite = "),
        std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, PersistentTerminalReasonRequiresStableReason)
{
    const audit_detail::ClusterKey key{ 0 };
    const std::vector<audit_detail::ClusterKey> accepted_key_list{
        key
    };
    const rg::GaussianModel3D model{ 8.0, 0.50, -0.10 };
    audit_detail::FitState previous_state;
    previous_state.emplace_back(MakeGaussianResult(model));
    const auto assembled_state{ previous_state };
    audit_detail::ClusterHealthMap health_by_key;
    health_by_key.emplace(
        key,
        audit_detail::ClusterHealth{
            audit_detail::JointOffsetSolveStatus::Converged });
    std::vector<char> suspicious_atom_mask{ 1 };
    audit_detail::PersistentTerminalFailureStateMap state_by_key;

    for (std::size_t stable_count = 1;
        stable_count < audit_detail::kPersistentTerminalFailureIterationLimit;
        stable_count++)
    {
        const auto terminal_failure_by_key{
            audit_detail::UpdatePersistentTerminalFailureState(
                accepted_key_list,
                suspicious_atom_mask,
                health_by_key,
                assembled_state,
                previous_state,
                state_by_key)
        };
        EXPECT_TRUE(terminal_failure_by_key.empty());
        ASSERT_EQ(state_by_key.size(), 1U);
        EXPECT_EQ(
            state_by_key.at(key).stable_iteration_count,
            stable_count);
    }

    const auto terminal_failure_by_key{
        audit_detail::UpdatePersistentTerminalFailureState(
            accepted_key_list,
            suspicious_atom_mask,
            health_by_key,
            assembled_state,
            previous_state,
            state_by_key)
    };
    ASSERT_EQ(terminal_failure_by_key.size(), 1U);
    ASSERT_TRUE(state_by_key.empty());
    ASSERT_TRUE(std::holds_alternative<
        audit_detail::PersistentSuspiciousRollbackReason>(
            terminal_failure_by_key.at(key)));
    EXPECT_EQ(
        std::get<audit_detail::PersistentSuspiciousRollbackReason>(
            terminal_failure_by_key.at(key)),
        (audit_detail::PersistentSuspiciousRollbackReason{ 0 }));

    suspicious_atom_mask.at(0) = 0;
    health_by_key.at(key).joint_offset_status =
        audit_detail::JointOffsetSolveStatus::SystemBuildFailed;
    const auto changed_reason_terminal_failure_by_key{
        audit_detail::UpdatePersistentTerminalFailureState(
            accepted_key_list,
            suspicious_atom_mask,
            health_by_key,
            assembled_state,
            previous_state,
            state_by_key)
    };
    EXPECT_TRUE(changed_reason_terminal_failure_by_key.empty());
    ASSERT_EQ(state_by_key.size(), 1U);
    EXPECT_EQ(state_by_key.at(key).stable_iteration_count, 1U);
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

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingRollsBackFiniteNonphysicalProfile)
{
    auto model{ BuildFiniteNonphysicalProfileDefenseModel() };
    auto * atom{ model->GetSelectedAtoms().front() };
    const auto previous_model{ GetEstimateModel(*atom) };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    ExpectGaussianModelsNear(GetEstimateModel(*atom), previous_model, 1.0e-12);
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
    RunSecondStageLocalFittingPersistsFinalModelPeelingBeforeGroupFitting)
{
    auto model{ BuildSharedOffsetJointPolishDefenseModel() };
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    analysis.InitializeGroupAlpha(rg::FittingStage::Second, 0.0);
    const auto options{ MakeSecondStageOptions() };

    rt::RunSecondStageLocalFitting(*model, options);

    ExpectPeelingSamplingEntriesMatchFinalModels(*model);
    const auto analysis_view{ model->GetAnalysisView() };
    for (const auto group_key : analysis_view.CollectAtomGroupKeys(
        rg::FittingStage::Second))
    {
        const auto & atom_list{
            analysis_view.GetAtomObjectList(
                rg::FittingStage::Second, group_key)
        };
        std::vector<LocalPotentialSampleList> sample_entries_list;
        std::vector<rg::LocalGaussianResult> member_result_list;
        sample_entries_list.reserve(atom_list.size());
        member_result_list.reserve(atom_list.size());
        for (const auto * atom : atom_list)
        {
            const auto local_view{
                rg::AtomLocalPotentialView::RequireFor(*atom)
            };
            sample_entries_list.emplace_back(
                local_view.GetPeelingSamplingEntries(false));
            member_result_list.emplace_back(
                local_view.GetGaussianResult(FittingStage::Second));
        }
        const auto expected_group_result{
            rt::EstimateGroupGaussian(
                sample_entries_list,
                member_result_list,
                analysis_view.GetAtomAlphaG(
                    rg::FittingStage::Second, group_key),
                options)
        };

        ExpectGaussianModelsNear(
            analysis_view.GetAtomGroupMean(
                rg::FittingStage::Second, group_key),
            expected_group_result.mean,
            1.0e-10);
        ExpectGaussianModelsNear(
            analysis_view.GetAtomGroupMDPDE(
                rg::FittingStage::Second, group_key),
            expected_group_result.mdpde,
            1.0e-10);
        ExpectGaussianModelsNear(
            analysis_view.GetAtomGroupPrior(
                rg::FittingStage::Second, group_key),
            expected_group_result.prior.GetModel(),
            1.0e-10);
        ASSERT_EQ(
            expected_group_result.member_results.size(),
            atom_list.size());
        for (std::size_t i = 0; i < atom_list.size(); i++)
        {
            const auto actual_result{
                rg::AtomLocalPotentialView::RequireFor(*atom_list.at(i))
                    .GetGaussianResult(FittingStage::Second)
            };
            ASSERT_TRUE(actual_result.posterior.has_value());
            ExpectGaussianModelsNear(
                actual_result.posterior->GetModel(),
                expected_group_result.member_results.at(i).mdpde.GetModel(),
                1.0e-10);
        }
    }
}

TEST(
    EstimatorSecondStageDefenseTest,
    RunSecondStageLocalFittingIncludesEffectiveUnselectedContributors)
{
    const rg::GaussianModel3D seed_model{ 6.0, 0.55, 0.10 };
    std::optional<float> include_hydrogen_response;
    for (const bool exclude_hydrogen : { false, true })
    {
        auto model{ BuildUnselectedContributorDefenseModel(seed_model) };
        auto options{ MakeSecondStageOptions() };
        options.quiet_mode = false;
        options.exclude_hydrogen = exclude_hydrogen;
        const auto previous_log_level{ Logger::GetLogLevel() };
        Logger::SetLogLevel(LogLevel::Debug);
        testing::internal::CaptureStdout();
        const auto peeling_applied{ rt::RunSecondStageLocalFitting(*model, options) };
        const std::string out{ testing::internal::GetCapturedStdout() };
        Logger::SetLogLevel(previous_log_level);

        EXPECT_EQ(model->GetSelectedAtomCount(), 2U);
        EXPECT_TRUE(peeling_applied);
        const auto expected_neighbor_count{ exclude_hydrogen ? 3 : 4 };
        for (const int serial_id : { 1, 2 })
        {
            EXPECT_EQ(
                rg::AtomLocalPotentialView::RequireFor(
                    *model->FindAtomPtr(serial_id))
                    .GetNeighborCountForPeeling(),
                expected_neighbor_count);
        }
        for (int serial_id = 3; serial_id <= 7; serial_id++)
        {
            EXPECT_FALSE(rg::AtomLocalPotentialView::For(
                *model->FindAtomPtr(serial_id)).IsAvailable());
        }
        ExpectUnselectedContributorPeeling(
            *model,
            seed_model,
            exclude_hydrogen);

        EXPECT_NE(
            out.find(
                exclude_hydrogen ?
                    "Unselected second-stage neighbor seeds = 2, "
                    "sources = group-median:1, global-median:1." :
                    "Unselected second-stage neighbor seeds = 3, "
                    "sources = group-median:2, global-median:1."),
            std::string::npos);
        EXPECT_NE(
            out.find(
                "Unselected second-stage neighbor seed selection: serial ID = 3"),
            std::string::npos);
        EXPECT_NE(
            out.find(
                "Unselected second-stage neighbor seed selection: serial ID = 4"),
            std::string::npos);
        EXPECT_EQ(
            out.find(
                "Unselected second-stage neighbor seed selection: serial ID = 6"),
            std::string::npos);
        EXPECT_EQ(
            out.find(
                "Unselected second-stage neighbor seed selection: serial ID = 7"),
            std::string::npos);
        const auto hydrogen_log_position{
            out.find(
                "Unselected second-stage neighbor seed selection: serial ID = 5")
        };
        if (exclude_hydrogen)
        {
            EXPECT_EQ(hydrogen_log_position, std::string::npos);
        }
        else
        {
            EXPECT_NE(hydrogen_log_position, std::string::npos);
        }

        const auto first_response{
            rg::AtomLocalPotentialView::RequireFor(
                *model->GetSelectedAtoms().front())
                .GetPeelingSamplingEntries(false).front().response
        };
        if (exclude_hydrogen)
        {
            ASSERT_TRUE(include_hydrogen_response.has_value());
            EXPECT_NE(first_response, *include_hydrogen_response);
        }
        else
        {
            include_hydrogen_response = first_response;
        }
    }

    auto quiet_model{ BuildUnselectedContributorDefenseModel(seed_model) };
    auto quiet_options{ MakeSecondStageOptions() };
    quiet_options.quiet_mode = true;
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*quiet_model, quiet_options);
    const std::string quiet_out{ testing::internal::GetCapturedStdout() };
    EXPECT_EQ(
        quiet_out.find("Unselected second-stage neighbor seeds"),
        std::string::npos);
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

TEST(EstimatorSecondStageDefenseTest, SharedOffsetGroupParticipatesInAcceptedJointPolish)
{
    auto model{ BuildSharedOffsetJointPolishDefenseModel() };
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    const auto initial_error{ CalculateSelectedAtomResponseMeanSquaredError(*model) };

    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };

    const auto & atom_list{ model->GetSelectedAtoms() };
    ASSERT_EQ(atom_list.size(), 2U);
    EXPECT_EQ(atom_list.at(0)->GetAtomKey(), atom_list.at(1)->GetAtomKey());
    EXPECT_DOUBLE_EQ(
        GetEstimateModel(*atom_list.at(0)).GetOffset(),
        GetEstimateModel(*atom_list.at(1)).GetOffset());
    EXPECT_LT(CalculateSelectedAtomResponseMeanSquaredError(*model), initial_error);
    EXPECT_NE(out.find("1/1/0/0"), std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorSecondStageDefenseTest, SharedOffsetJointPolishIsIntensityScaleInvariant)
{
    constexpr double scale{ 100.0 };
    auto base_model{ BuildSharedOffsetJointPolishDefenseModel() };
    auto scaled_model{ BuildSharedOffsetJointPolishDefenseModel(scale) };

    rt::RunSecondStageLocalFitting(*base_model, MakeSecondStageOptions());
    rt::RunSecondStageLocalFitting(*scaled_model, MakeSecondStageOptions());

    const auto & base_atoms{ base_model->GetSelectedAtoms() };
    const auto & scaled_atoms{ scaled_model->GetSelectedAtoms() };
    ASSERT_EQ(base_atoms.size(), scaled_atoms.size());
    ASSERT_EQ(base_atoms.size(), 2U);
    EXPECT_DOUBLE_EQ(
        GetEstimateModel(*base_atoms.at(0)).GetOffset(),
        GetEstimateModel(*base_atoms.at(1)).GetOffset());
    EXPECT_DOUBLE_EQ(
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

TEST(EstimatorSecondStageDefenseTest, RunSecondStageLocalFittingAcceptsRemoteClusterWhenLocalClusterRollsBack)
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
        EXPECT_NEAR(base.GetWidth(), scaled.GetWidth(), 1.0e-6);
        EXPECT_NEAR(
            base.GetOffset() * intensity_scale,
            scaled.GetOffset(),
            std::max(1.0e-8, std::abs(scaled.GetOffset()) * 5.0e-5));
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

TEST(EstimatorSecondStageDefenseTest, ObjectiveSeparatesFitAndTailDomains)
{
    auto model{
        BuildDefenseModel(
            { std::array<float, 3>{ 0.0F, 0.0F, 0.0F } },
            { Spot::O },
            { Element::OXYGEN },
            { rg::GaussianModel3D{ 8.0, 0.45, -0.1 } },
            rg::GaussianModel3D{ 6.0, 0.60, 0.0 })
    };
    auto options{ MakeSecondStageOptions() };
    options.distance_max = 0.5;
    options.quiet_mode = false;
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(out.find("unique fit/tail samples = 12/12"), std::string::npos);
    EXPECT_NE(
        out.find("fixed tail scale median/p99/max = "),
        std::string::npos);
    EXPECT_EQ(
        out.find("fixed tail scale median/p99/max = unavailable"),
        std::string::npos);
}

TEST(EstimatorSecondStageDefenseTest, ObjectiveDomainCountsCutBoundarySamplesOnce)
{
    auto model{
        BuildDefenseModel(
            {
                std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 3.4F, 0.0F, 0.0F }
            },
            { Spot::C, Spot::O },
            { Element::CARBON, Element::OXYGEN },
            {
                rg::GaussianModel3D{ 8.0, 0.90, 0.20 },
                rg::GaussianModel3D{ 3.0, 0.80, -0.10 }
            },
            rg::GaussianModel3D{ 5.5, 0.30, 0.0 })
    };
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(out.find("candidate/retained/cut edges = 1/0/1"), std::string::npos);
    EXPECT_NE(out.find("clusters = 2"), std::string::npos);
    EXPECT_NE(out.find("unique fit/tail samples = 48/0"), std::string::npos);
    EXPECT_NE(
        out.find("Boundary-component reconciliation:"),
        std::string::npos);
}

TEST(
    EstimatorSecondStageDefenseTest,
    BoundaryComponentReconciliationBacktracksAndPreservesRemoteCluster)
{
    auto model{ BuildBoundaryComponentConflictDefenseModel() };
    const auto initial_remote_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model, 11, 13)
    };
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);
    EXPECT_NE(
        out.find(
            "Boundary-component reconciliation: clusters/atoms/boundary-samples = 2/11/"),
        std::string::npos);
    EXPECT_NE(
        out.find("trials/factor/accepted/exhausted = 2/0.5/yes/no"),
        std::string::npos);
    EXPECT_NE(
        out.find("trials/factor/accepted/exhausted = 1/-/no/yes"),
        std::string::npos);
    EXPECT_NE(out.find("| 1/2"), std::string::npos);
    EXPECT_LT(
        CalculateSelectedAtomResponseMeanSquaredError(*model, 11, 13),
        initial_remote_error);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(
    EstimatorSecondStageDefenseTest,
    RunSecondStageLocalFittingDampsOffsetStepIntoInitialTrustRadius)
{
    const rg::GaussianModel3D initial_model{ 6.0, 0.55, 0.0 };
    auto truth_coordinates{
        change_detail::EncodeTransformedCoordinates(initial_model)
    };
    ASSERT_TRUE(truth_coordinates.has_value());
    (*truth_coordinates)(static_cast<Eigen::Index>(
        change_detail::kOffsetToPeakRatioChangeIndex)) = 1.25;
    const auto truth_model{
        change_detail::DecodeTransformedCoordinates(
            *truth_coordinates)
    };
    ASSERT_TRUE(truth_model.has_value());

    auto model{
        BuildDefenseModel(
            { std::array<float, 3>{ 0.0F, 0.0F, 0.0F } },
            { Spot::O },
            { Element::OXYGEN },
            { *truth_model },
            initial_model)
    };
    const auto previous_model{
        GetEstimateModel(*model->GetSelectedAtoms().front())
    };
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    const auto fitted_model{
        GetEstimateModel(*model->GetSelectedAtoms().front())
    };
    EXPECT_NE(fitted_model.GetOffset(), previous_model.GetOffset());
    EXPECT_NE(out.find("accepted_iterations="), std::string::npos);
    EXPECT_EQ(out.find("accepted_iterations=0"), std::string::npos);
    EXPECT_EQ(
        out.find("previous-shared-offset-projection-outside-trust-region"),
        std::string::npos);
    EXPECT_EQ(out.find("objective = not-evaluated"), std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(
    EstimatorSecondStageDefenseTest,
    PreObjectiveTrustFailureDoesNotReportUnavailableObjective)
{
    const rg::GaussianModel3D initial_model{ 6.0, 0.55, 0.0 };
    auto model{
        BuildDefenseModel(
            {
                std::array<float, 3>{ 0.0F, 0.0F, 0.0F },
                std::array<float, 3>{ 1.0e-4F, 0.0F, 0.0F }
            },
            { Spot::C, Spot::C },
            { Element::CARBON, Element::CARBON },
            { initial_model, initial_model },
            initial_model)
    };
    const auto & atom_list{ model->GetSelectedAtoms() };
    auto analysis{ model->EditAnalysis() };
    analysis.EnsureAtomLocalPotential(*atom_list.at(0))
        .SetGaussianResult(FittingStage::Second, MakeGaussianResult(
            initial_model.WithOffset(-2.0)));
    analysis.EnsureAtomLocalPotential(*atom_list.at(1))
        .SetGaussianResult(FittingStage::Second, MakeGaussianResult(
            initial_model.WithOffset(2.0)));
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(
        out.find("previous-shared-offset-projection-outside-trust-region"),
        std::string::npos);
    EXPECT_NE(out.find("objective = not-evaluated"), std::string::npos);
    EXPECT_EQ(out.find("objective-unavailable"), std::string::npos);
    EXPECT_EQ(out.find("fit/tail samples = 0/0"), std::string::npos);
}

TEST(
    EstimatorSecondStageDefenseTest,
    MissingPosteriorAndPriorSkipsSecondStageDespiteValidLocalSeeds)
{
    auto model{ BuildNearCollinearDefenseModel() };
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    std::vector<rg::GaussianModel3D> previous_model_list;
    std::vector<LocalPotentialSampleList> previous_peeling_sampling_entries_list;
    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    analysis.InitializeGroupAlpha(rg::FittingStage::Second, 0.0);
    for (auto * atom : model->GetSelectedAtoms())
    {
        auto result{
            rg::AtomLocalPotentialView::RequireFor(*atom).GetGaussianResult(
                FittingStage::Second)
        };
        ASSERT_TRUE(seed_detail::IsValidSecondStageGaussianModel(
            result.mdpde.GetModel()));
        ASSERT_TRUE(seed_detail::IsValidSecondStageGaussianModel(
            result.ols.GetModel()));
        result.posterior.reset();
        auto local_editor{ analysis.EnsureAtomLocalPotential(*atom) };
        local_editor.SetGaussianResult(
            FittingStage::Second,
            std::move(result));
        LocalPotentialSampleList sentinel_peeling_sampling_entries{
            LocalPotentialSample{
                static_cast<float>(100.0 + previous_model_list.size()),
                SamplingPoint{ 0.5F, atom->GetPosition(), true }
            }
        };
        local_editor.SetPeelingSamplingEntries(
            sentinel_peeling_sampling_entries);
        local_editor.SetNeighborCountForPeeling(99);
        previous_model_list.emplace_back(GetEstimateModel(*atom));
        previous_peeling_sampling_entries_list.emplace_back(
            std::move(sentinel_peeling_sampling_entries));
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

    testing::internal::CaptureStdout();
    const auto peeling_applied{ rt::RunSecondStageLocalFitting(*model, options) };
    const std::string out{ testing::internal::GetCapturedStdout() };

    EXPECT_FALSE(peeling_applied);

    for (std::size_t i = 0; i < model->GetSelectedAtoms().size(); i++)
    {
        ExpectGaussianModelsNear(
            GetEstimateModel(*model->GetSelectedAtoms().at(i)),
            previous_model_list.at(i),
            0.0);
        const auto peeling_sampling_entries{
            rg::AtomLocalPotentialView::RequireFor(
                *model->GetSelectedAtoms().at(i))
                .GetPeelingSamplingEntries(false)
        };
        ASSERT_EQ(peeling_sampling_entries.size(), 1U);
        EXPECT_FLOAT_EQ(
            peeling_sampling_entries.front().response,
            previous_peeling_sampling_entries_list.at(i).front().response);
        EXPECT_NE(
            rg::AtomLocalPotentialView::RequireFor(
                *model->GetSelectedAtoms().at(i))
                .GetNeighborCountForPeeling(),
            99);
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
    EXPECT_NE(out.find("stop_reason=no-valid-seed"), std::string::npos);
    EXPECT_NE(out.find("final_uses_polish=unavailable"), std::string::npos);
    EXPECT_EQ(
        out.find("Second-stage best-iteration application:"),
        std::string::npos);
    EXPECT_NE(
        out.find("final_state_source=unavailable"),
        std::string::npos);
}

TEST(EstimatorSecondStageDefenseTest, NonQuietSecondStageReportsAcceptedJointPolish)
{
    auto model{ BuildJointPolishDefenseModel() };
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(out.find("final_uses_polish=yes"), std::string::npos);
    bool found_accepted_polish{ false };
    bool found_skipped_polish{ false };
    for (std::size_t row_start = out.find('\r');
        row_start != std::string::npos;
        row_start = out.find('\r', row_start))
    {
        row_start++;
        const auto row_end{ out.find_first_of("\r\n", row_start) };
        ASSERT_NE(row_end, std::string::npos);
        const std::string_view row{
            out.data() + row_start,
            row_end - row_start
        };
        const auto polish_count{ ParsePolishProgressCounts(row) };
        if (!polish_count.has_value()) continue;
        found_accepted_polish = found_accepted_polish || polish_count->at(1) > 0;
        found_skipped_polish = found_skipped_polish || polish_count->at(3) > 0;
        row_start = row_end;
    }
    EXPECT_TRUE(found_accepted_polish);
    EXPECT_TRUE(found_skipped_polish);
}

TEST(EstimatorSecondStageDefenseTest, NonQuietSecondStageReportsNoPolishForEmptySelection)
{
    auto model{ BuildJointPolishDefenseModel() };
    model->SelectAllAtoms(false);
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;

    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };

    EXPECT_NE(out.find("final_uses_polish=no"), std::string::npos);
    EXPECT_NE(
        out.find("final_state_source=latest-validated"),
        std::string::npos);
}

TEST(EstimatorSecondStageDefenseTest, NonQuietSecondStageLogsEveryOuterAttempt)
{
    auto model{ BuildSeparatedRollbackDefenseModel() };
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    const auto count_occurrences = [&](std::string_view text)
    {
        std::size_t count{ 0 };
        for (std::size_t position = 0;
            (position = out.find(text, position)) != std::string::npos;
            position += text.size())
        {
            count++;
        }
        return count;
    };
    EXPECT_EQ(count_occurrences("Try/Acc"), 1U);
    EXPECT_NE(out.find("Atom A/T"), std::string::npos);
    EXPECT_NE(out.find("Cluster A/R"), std::string::npos);
    EXPECT_NE(out.find("Polish E/A/R/S"), std::string::npos);
    EXPECT_NE(out.find("Suspicious"), std::string::npos);
    EXPECT_NE(out.find("dMax A/R"), std::string::npos);
    EXPECT_EQ(count_occurrences("Local-fitting residue cutoff:"), 1U);
    EXPECT_NE(
        out.find("Adaptive local-fitting topology rebuild:"),
        std::string::npos);
    EXPECT_NE(
        out.find("% Rebuild local-fitting coupling topology"),
        std::string::npos);
    EXPECT_NE(
        out.find("topology_rebuilds/partition_changes="),
        std::string::npos);
    EXPECT_NE(
        out.find("boundary_reconciliations/backtracked/rejected="),
        std::string::npos);
    EXPECT_NE(
        out.find("boundary_reconciliation_ms="),
        std::string::npos);
    EXPECT_NE(
        out.find("iteration/candidate/topology/total_ms="),
        std::string::npos);

    const auto header_start{ out.find("Try/Acc") };
    ASSERT_NE(header_start, std::string::npos);
    const auto header_end{ out.find('\n', header_start) };
    ASSERT_NE(header_end, std::string::npos);
    const std::string_view header{
        out.data() + header_start,
        header_end - header_start
    };
    std::vector<std::string_view> progress_row_list;
    for (std::size_t row_start = out.find('\r');
        row_start != std::string::npos;
        row_start = out.find('\r', row_start))
    {
        row_start++;
        const auto row_end{ out.find_first_of("\r\n", row_start) };
        ASSERT_NE(row_end, std::string::npos);
        const std::string_view row{
            out.data() + row_start,
            row_end - row_start
        };
        if (ParsePolishProgressCounts(row).has_value())
        {
            progress_row_list.emplace_back(row);
        }
        row_start = row_end;
    }
    const auto separator_position_list = [](std::string_view row)
    {
        std::vector<std::size_t> position_list;
        for (std::size_t position = row.find('|');
            position != std::string::npos;
            position = row.find('|', position + 1))
        {
            position_list.emplace_back(position);
        }
        return position_list;
    };
    const auto header_separator_position_list{
        separator_position_list(header)
    };
    ASSERT_EQ(header_separator_position_list.size(), 5U);
    ASSERT_EQ(progress_row_list.size(), 4U);
    bool found_skipped_polish{ false };
    for (const auto row : progress_row_list)
    {
        EXPECT_EQ(row.size(), header.size());
        EXPECT_EQ(
            separator_position_list(row),
            header_separator_position_list);

        const auto polish_count{ ParsePolishProgressCounts(row) };
        ASSERT_TRUE(polish_count.has_value());
        EXPECT_EQ(
            polish_count->at(0),
            polish_count->at(1) + polish_count->at(2) + polish_count->at(3));
        found_skipped_polish = found_skipped_polish || polish_count->at(3) > 0;
    }
    EXPECT_TRUE(found_skipped_polish);
    EXPECT_NE(
        progress_row_list.front().find("3.55e-02/4.14e-02"),
        std::string::npos);

    const std::string summary_prefix{
        "Second-stage local fitting summary: accepted_iterations="
    };
    const auto summary_position{ out.find(summary_prefix) };
    ASSERT_NE(summary_position, std::string::npos);
    ASSERT_GT(summary_position, 0U);
    EXPECT_EQ(out.at(summary_position - 1), '\n');
    const auto accepted_start{ summary_position + summary_prefix.size() };
    const auto accepted_end{ out.find(',', accepted_start) };
    ASSERT_NE(accepted_end, std::string::npos);
    const auto accepted_iteration_count{
        static_cast<std::size_t>(std::stoull(
            out.substr(accepted_start, accepted_end - accepted_start)))
    };
    EXPECT_EQ(accepted_iteration_count, 4U);
    EXPECT_NE(
        out.find(
            "best_iteration=1, stop_reason=audit-patience"),
        std::string::npos);
    EXPECT_TRUE(
        out.find("final_uses_polish=yes") != std::string::npos ||
        out.find("final_uses_polish=no") != std::string::npos);
    EXPECT_NE(
        out.find("final_state_source=best-audit"),
        std::string::npos);
    ExpectPeelingSamplingEntriesMatchFinalModels(*model);
}

TEST(EstimatorSecondStageDefenseTest, QuietSecondStageSuppressesIterationTable)
{
    auto model{ BuildSeparatedRollbackDefenseModel() };
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());
    const std::string out{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_EQ(out.find("Try/Acc"), std::string::npos);
    EXPECT_EQ(std::count(out.begin(), out.end(), '\r'), 0);
    EXPECT_EQ(
        out.find("Rebuild local-fitting coupling topology"),
        std::string::npos);
}
