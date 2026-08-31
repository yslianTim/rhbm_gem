#include "core/detail/FittingModel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

namespace rhbm_gem::core::detail {

namespace {

constexpr double kTransformedChangePercentile{ 0.99 };
constexpr std::array<double, kTransformedChangeSize>
    kTrustRegionParameterScale{ 0.50, 0.35, 1.0 };
constexpr double kTrustRegionBoundaryTolerance{ 1.0e-12 };

std::vector<double> SummarizeMaximumTransformedChanges(
    const std::vector<algorithm::ParameterChange> & change_list,
    const std::vector<std::size_t> & index_list)
{
    std::vector<double> maximum_list(kTransformedChangeSize, 0.0);
    for (const auto index : index_list)
    {
        if (index >= change_list.size() ||
            change_list.at(index).value_list.size() != kTransformedChangeSize)
        {
            throw std::invalid_argument(
                "Local fitting maximum transformed change input is inconsistent.");
        }
        for (std::size_t parameter_index = 0;
            parameter_index < kTransformedChangeSize;
            parameter_index++)
        {
            maximum_list.at(parameter_index) = std::max(
                maximum_list.at(parameter_index),
                change_list.at(index).value_list.at(parameter_index));
        }
    }
    return maximum_list;
}

double CalculateSecondStageAdjustedResponse(
    const AtomContext & atom_context,
    std::size_t sample_index,
    const SecondStageModelSnapshot & model_snapshot)
{
    auto response_value{
        static_cast<double>(
            atom_context.raw_sampling_entries.at(sample_index).response)
    };
    for (const auto & neighbor_atom_sample : atom_context.Neighbors(sample_index))
    {
        response_value -= ResolveNeighborAtomModel(
            neighbor_atom_sample,
            model_snapshot).ResponseAtDistance(neighbor_atom_sample.distance);
    }
    return response_value;
}

template<typename State>
FittedGaussianSnapshot BuildFittedGaussianSnapshotImpl(const State & state)
{
    FittedGaussianSnapshot snapshot;
    snapshot.reserve(state.size());
    for (std::size_t i = 0; i < state.size(); i++)
    {
        snapshot.emplace_back(GetFitModel(state, i));
    }
    return snapshot;
}

template<typename State>
std::optional<ResidualSample> EvaluateResidualSampleImpl(
    const SecondStageContext & context,
    const State & state,
    const SampleRef & sample_ref,
    const SecondStageModelSnapshot & model_snapshot)
{
    const auto & atom_context{ context.at(sample_ref.atom_index) };
    const auto & sample{
        atom_context.raw_sampling_entries.at(sample_ref.sample_index)
    };
    const auto adjusted_response{
        CalculateSecondStageAdjustedResponse(
            atom_context,
            sample_ref.sample_index,
            model_snapshot)
    };
    const auto expected_response{
        GetFitModel(state, sample_ref.atom_index).ResponseAtDistance(
            static_cast<double>(sample.point.distance))
    };
    const auto residual{ adjusted_response - expected_response };
    if (!std::isfinite(adjusted_response) || !std::isfinite(residual))
    {
        return std::nullopt;
    }
    return ResidualSample{ adjusted_response, residual };
}

template<typename CurrentState, typename PreviousState>
TransformedChangeSummary SummarizeTransformedChangesImpl(
    const CurrentState & current_state,
    const PreviousState & previous_state,
    const std::vector<std::size_t> & index_list)
{
    std::vector<algorithm::ParameterChange> change_list;
    change_list.reserve(index_list.size());
    for (const auto i : index_list)
    {
        change_list.emplace_back(CalculateTransformedChange(
            GetFitModel(current_state, i),
            GetFitModel(previous_state, i)));
    }

    std::vector<std::size_t> local_index_list(change_list.size());
    for (std::size_t i = 0; i < local_index_list.size(); i++)
    {
        local_index_list.at(i) = i;
    }
    TransformedChangeSummary summary{
        algorithm::SummarizeParameterChangeStats(
            change_list,
            local_index_list,
            kTransformedChangePercentile),
        SummarizeMaximumTransformedChanges(change_list, local_index_list)
    };
    summary.population_size_list.fill(index_list.size());
    return summary;
}

} // namespace

FitStatePatch FitStatePatch::FromState(
    const FitState & state,
    ClusterKey atom_index_list)
{
    std::ranges::sort(atom_index_list);
    atom_index_list.erase(
        std::ranges::unique(atom_index_list).begin(),
        atom_index_list.end());

    FitStatePatch patch;
    patch.atom_index_list = std::move(atom_index_list);
    patch.mdpde_list.reserve(patch.atom_index_list.size());
    for (const auto atom_index : patch.atom_index_list)
    {
        patch.mdpde_list.emplace_back(state.at(atom_index).mdpde);
    }
    return patch;
}

const GaussianModel3DWithUncertainty * FitStatePatch::Find(
    std::size_t atom_index) const
{
    const auto iter{
        std::ranges::lower_bound(atom_index_list, atom_index)
    };
    if (iter == atom_index_list.end() || *iter != atom_index)
    {
        return nullptr;
    }
    return &mdpde_list.at(static_cast<std::size_t>(
        std::distance(atom_index_list.begin(), iter)));
}

void FitStatePatch::ApplyTo(FitState & state) const
{
    if (atom_index_list.size() != mdpde_list.size())
    {
        throw std::invalid_argument(
            "Local fitting state patch sizes are inconsistent.");
    }
    for (std::size_t i = 0; i < atom_index_list.size(); i++)
    {
        state.at(atom_index_list.at(i)).mdpde = mdpde_list.at(i);
    }
}

std::optional<ResidualSample> SnapshotResidualEvaluator::operator()(
    const SampleRef & sample_ref) const
{
    return EvaluateResidualSample(
        context,
        sample_ref,
        model_snapshot);
}

const GaussianModel3D & GetFitModel(const FitState & state, std::size_t atom_index)
{
    return state.at(atom_index).mdpde.GetModel();
}

const GaussianModel3D & GetFitModel(const FitStateView & state, std::size_t atom_index)
{
    return state.GetModel(atom_index);
}

const GaussianModel3D & GetFitModel(
    const FittedGaussianSnapshot & state,
    std::size_t atom_index)
{
    return state.at(atom_index);
}

FittedGaussianSnapshot BuildFittedGaussianSnapshot(const FitState & state)
{
    return BuildFittedGaussianSnapshotImpl(state);
}

FittedGaussianSnapshot BuildFittedGaussianSnapshot(const FitStateView & state)
{
    return BuildFittedGaussianSnapshotImpl(state);
}

std::optional<ResidualSample> EvaluateResidualSample(
    const SecondStageContext & context,
    const FitStateView & state,
    const SampleRef & sample_ref,
    const SecondStageModelSnapshot & model_snapshot)
{
    return EvaluateResidualSampleImpl(context, state, sample_ref, model_snapshot);
}

std::optional<ResidualSample> EvaluateResidualSample(
    const SecondStageContext & context,
    const SampleRef & sample_ref,
    const SecondStageModelSnapshot & model_snapshot)
{
    return EvaluateResidualSampleImpl(
        context,
        model_snapshot.selected,
        sample_ref,
        model_snapshot);
}

TransformedChangeSummary SummarizeTransformedChanges(
    const FitState & current_state,
    const FitState & previous_state,
    const std::vector<std::size_t> & index_list)
{
    return SummarizeTransformedChangesImpl(
        current_state,
        previous_state,
        index_list);
}

TransformedChangeSummary SummarizeTransformedChanges(
    const FitStateView & current_state,
    const FittedGaussianSnapshot & previous_state,
    const std::vector<std::size_t> & index_list)
{
    return SummarizeTransformedChangesImpl(
        current_state,
        previous_state,
        index_list);
}

TransformedChangeSummary SummarizeTransformedChangesByParameter(
    const std::vector<algorithm::ParameterChange> & change_list,
    const TransformedChangeIndexListByParameter & index_list_by_parameter)
{
    TransformedChangeSummary summary;
    summary.percentile_stats.percentile_list.resize(kTransformedChangeSize, 0.0);
    summary.maximum_list.resize(kTransformedChangeSize, 0.0);
    for (std::size_t parameter_index = 0;
        parameter_index < kTransformedChangeSize;
        parameter_index++)
    {
        const auto & index_list{ index_list_by_parameter.at(parameter_index) };
        summary.population_size_list.at(parameter_index) = index_list.size();
        std::vector<double> parameter_change_list;
        parameter_change_list.reserve(index_list.size());
        for (const auto index : index_list)
        {
            if (index >= change_list.size() ||
                change_list.at(index).value_list.size() != kTransformedChangeSize)
            {
                throw std::invalid_argument(
                    "Local fitting masked transformed change input is inconsistent.");
            }
            parameter_change_list.emplace_back(
                change_list.at(index).value_list.at(parameter_index));
        }
        summary.percentile_stats.percentile_list.at(parameter_index) =
            array_helper::ComputePercentile(
                parameter_change_list,
                kTransformedChangePercentile);
        summary.maximum_list.at(parameter_index) =
            parameter_change_list.empty() ? 0.0 :
                *std::ranges::max_element(parameter_change_list);
    }
    return summary;
}

TransformedChangeSummary SummarizeTransformedChangesByParameter(
    const FitState & current_state,
    const FitState & previous_state,
    const TransformedChangeIndexListByParameter & index_list_by_parameter)
{
    if (current_state.size() != previous_state.size())
    {
        throw std::invalid_argument(
            "Local fitting masked transformed change state sizes are inconsistent.");
    }
    std::vector<algorithm::ParameterChange> change_list;
    change_list.reserve(current_state.size());
    for (std::size_t i = 0; i < current_state.size(); i++)
    {
        change_list.emplace_back(CalculateTransformedChange(
            GetFitModel(current_state, i),
            GetFitModel(previous_state, i)));
    }
    return SummarizeTransformedChangesByParameter(
        change_list,
        index_list_by_parameter);
}

std::optional<Eigen::Vector3d> EncodeTransformedCoordinates(const GaussianModel3D & model)
{
    const auto amplitude{ model.GetAmplitude() };
    const auto width{ model.GetWidth() };
    const auto offset{ model.GetOffset() };
    if (!std::isfinite(amplitude) || amplitude <= 0.0 ||
        !std::isfinite(width) || width <= 0.0 ||
        !std::isfinite(offset))
    {
        return std::nullopt;
    }

    const auto log_width{ std::log(width) };
    const auto log_peak_height{
        std::log(amplitude) - 1.5 * std::log(Constants::two_pi) - 3.0 * log_width
    };
    double offset_to_peak_ratio{ 0.0 };
    if (offset != 0.0)
    {
        const auto log_abs_offset_to_peak_ratio{
            std::log(std::abs(offset)) +
            0.5 * std::log(4.0 / Constants::two_pi) - log_width - log_peak_height
        };
        if (log_abs_offset_to_peak_ratio > std::log(std::numeric_limits<double>::max()))
        {
            return std::nullopt;
        }
        offset_to_peak_ratio = std::copysign(std::exp(log_abs_offset_to_peak_ratio), offset);
    }

    if (!std::isfinite(log_peak_height) ||
        !std::isfinite(log_width) ||
        !std::isfinite(offset_to_peak_ratio))
    {
        return std::nullopt;
    }

    return Eigen::Vector3d{ log_peak_height, log_width, offset_to_peak_ratio };
}

std::optional<GaussianModel3D> DecodeTransformedCoordinates(const Eigen::Vector3d & coordinates)
{
    if (!coordinates.allFinite()) return std::nullopt;

    const auto log_peak_height{
        coordinates(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex))
    };
    const auto log_width{
        coordinates(static_cast<Eigen::Index>(kLogWidthChangeIndex))
    };
    const auto offset_to_peak_ratio{
        coordinates(static_cast<Eigen::Index>(kOffsetToPeakRatioChangeIndex))
    };
    const auto log_amplitude{
        log_peak_height + 1.5 * std::log(Constants::two_pi) + 3.0 * log_width
    };
    const auto amplitude{ std::exp(log_amplitude) };
    const auto width{ std::exp(log_width) };
    double offset{ 0.0 };
    if (offset_to_peak_ratio != 0.0)
    {
        const auto log_abs_offset{
            std::log(std::abs(offset_to_peak_ratio)) +
            log_peak_height + log_width - 0.5 * std::log(4.0 / Constants::two_pi)
        };
        offset = std::copysign(std::exp(log_abs_offset), offset_to_peak_ratio);
    }

    if (!std::isfinite(amplitude) || amplitude <= 0.0 ||
        !std::isfinite(width) || width <= 0.0 ||
        !std::isfinite(offset))
    {
        return std::nullopt;
    }
    return GaussianModel3D{ amplitude, width, offset };
}

bool IsValidSecondStageGaussianModel(const GaussianModel3D & model)
{
    return EncodeTransformedCoordinates(model).has_value();
}

GaussianModel3DWithUncertainty WithPreservedUncertaintyOffset(
    const GaussianModel3DWithUncertainty & gaussian,
    double offset)
{
    return {
        gaussian.GetModel().WithOffset(offset),
        gaussian.GetStandardDeviationModel()
    };
}

std::optional<GaussianModel3D> BuildGaussianParameterMedian(
    const std::vector<GaussianModel3D> & model_list)
{
    std::vector<double> amplitude_list;
    std::vector<double> width_list;
    std::vector<double> offset_list;
    amplitude_list.reserve(model_list.size());
    width_list.reserve(model_list.size());
    offset_list.reserve(model_list.size());
    for (const auto & model : model_list)
    {
        if (!IsValidSecondStageGaussianModel(model)) continue;
        amplitude_list.emplace_back(model.GetAmplitude());
        width_list.emplace_back(model.GetWidth());
        offset_list.emplace_back(model.GetOffset());
    }
    if (amplitude_list.empty()) return std::nullopt;

    const GaussianModel3D median_model{
        array_helper::ComputeMedian(amplitude_list),
        array_helper::ComputeMedian(width_list),
        array_helper::ComputeMedian(offset_list)
    };
    return IsValidSecondStageGaussianModel(median_model) ?
        std::optional<GaussianModel3D>{ median_model } : std::nullopt;
}

std::vector<GaussianModel3D> BuildGroupMedianModelList(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & model_list)
{
    if (group_id_by_atom_position.size() != model_list.size())
    {
        throw std::invalid_argument("Local fitting group median inputs are inconsistent.");
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> atom_position_list_by_group;
    atom_position_list_by_group.reserve(model_list.size());
    for (std::size_t atom_position = 0; atom_position < model_list.size(); atom_position++)
    {
        atom_position_list_by_group[group_id_by_atom_position.at(atom_position)]
            .emplace_back(atom_position);
    }

    auto group_median_model_list{ model_list };
    std::vector<GaussianModel3D> group_model_list;
    group_model_list.reserve(model_list.size());
    for (const auto & entry : atom_position_list_by_group)
    {
        const auto & atom_position_list{ entry.second };
        group_model_list.clear();
        for (const auto atom_position : atom_position_list)
        {
            group_model_list.emplace_back(model_list.at(atom_position));
        }
        const auto median_model{ BuildGaussianParameterMedian(group_model_list) };
        if (!median_model.has_value()) continue;
        for (const auto atom_position : atom_position_list)
        {
            group_median_model_list.at(atom_position) = *median_model;
        }
    }
    return group_median_model_list;
}

std::vector<double> BuildGroupMedianOffsetList(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & model_list)
{
    if (group_id_by_atom_position.size() != model_list.size())
    {
        throw std::invalid_argument("Local fitting group median inputs are inconsistent.");
    }

    std::unordered_map<std::size_t, std::vector<std::size_t>> atom_position_list_by_group;
    atom_position_list_by_group.reserve(model_list.size());
    std::vector<double> offset_list;
    offset_list.reserve(model_list.size());
    for (std::size_t atom_position = 0;
        atom_position < model_list.size();
        atom_position++)
    {
        atom_position_list_by_group[group_id_by_atom_position.at(atom_position)]
            .emplace_back(atom_position);
        offset_list.emplace_back(model_list.at(atom_position).GetOffset());
    }

    std::vector<double> valid_offset_list;
    valid_offset_list.reserve(model_list.size());
    for (const auto & entry : atom_position_list_by_group)
    {
        valid_offset_list.clear();
        for (const auto atom_position : entry.second)
        {
            const auto & model{ model_list.at(atom_position) };
            if (IsValidSecondStageGaussianModel(model))
            {
                valid_offset_list.emplace_back(model.GetOffset());
            }
        }
        if (valid_offset_list.empty()) continue;
        const auto median{ array_helper::ComputeMedian(valid_offset_list) };
        if (!std::isfinite(median)) continue;
        for (const auto atom_position : entry.second)
        {
            offset_list.at(atom_position) = median;
        }
    }
    return offset_list;
}

std::optional<std::vector<GaussianModel3D>> BuildSharedOffsetDampedModelList(
    const std::vector<GaussianModel3D> & previous_model_list,
    const std::vector<GaussianModel3D> & raw_model_list,
    const std::vector<double> & previous_shared_offset_list,
    const std::vector<double> & raw_shared_offset_list,
    double damping)
{
    if (raw_model_list.size() != previous_model_list.size() ||
        previous_shared_offset_list.size() != previous_model_list.size() ||
        raw_shared_offset_list.size() != previous_model_list.size() ||
        !std::isfinite(damping) || damping < 0.0 || damping > 1.0)
    {
        return std::nullopt;
    }

    std::vector<GaussianModel3D> candidate_model_list;
    candidate_model_list.reserve(previous_model_list.size());
    for (std::size_t atom_position = 0; atom_position < previous_model_list.size(); atom_position++)
    {
        const auto previous_coordinates{
            EncodeTransformedCoordinates(previous_model_list.at(atom_position))
        };
        const auto raw_coordinates{
            EncodeTransformedCoordinates(raw_model_list.at(atom_position))
        };
        if (!previous_coordinates.has_value() || !raw_coordinates.has_value())
        {
            return std::nullopt;
        }

        Eigen::Vector3d shape_coordinates{
            (*previous_coordinates)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)) +
                damping * (
                    (*raw_coordinates)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)) -
                    (*previous_coordinates)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex))),
            (*previous_coordinates)(static_cast<Eigen::Index>(kLogWidthChangeIndex)) +
                damping * (
                    (*raw_coordinates)(static_cast<Eigen::Index>(kLogWidthChangeIndex)) -
                    (*previous_coordinates)(static_cast<Eigen::Index>(kLogWidthChangeIndex))),
            0.0
        };
        const auto shape_model{ DecodeTransformedCoordinates(shape_coordinates) };
        if (!shape_model.has_value())
        {
            return std::nullopt;
        }

        const auto previous_shared_offset{ previous_shared_offset_list.at(atom_position) };
        const auto raw_shared_offset{ raw_shared_offset_list.at(atom_position) };
        const auto candidate_model{
            shape_model->WithOffset(previous_shared_offset + damping * (raw_shared_offset - previous_shared_offset))
        };
        if (!IsValidSecondStageGaussianModel(candidate_model))
        {
            return std::nullopt;
        }
        candidate_model_list.emplace_back(candidate_model);
    }
    return candidate_model_list;
}

static std::optional<SharedOffsetResponse> EvaluateValidSharedOffsetResponse(
    const GaussianModel3D & model,
    double distance)
{
    if (!std::isfinite(distance) || distance < 0.0) return std::nullopt;

    const auto width{ model.GetWidth() };
    const auto evaluation{ model.EvaluateAtDistance(distance) };
    if (!std::isfinite(evaluation.signal) ||
        !std::isfinite(evaluation.offset_basis) ||
        !std::isfinite(evaluation.response))
    {
        return std::nullopt;
    }

    const double center_offset_basis_scale{ std::sqrt(2.0 / M_PI) };
    const auto normalized_distance{ distance / width };
    auto log_width_derivative{
        evaluation.signal * normalized_distance * normalized_distance
    };
    if (distance < 1.0e-5)
    {
        log_width_derivative -= model.GetOffset() * center_offset_basis_scale / width;
    }
    else
    {
        const auto exponent{ -0.5 * normalized_distance * normalized_distance };
        log_width_derivative -= model.GetOffset() * center_offset_basis_scale / width * std::exp(exponent);
    }

    const Eigen::Vector2d shape_jacobian{ evaluation.signal, log_width_derivative };
    if (!shape_jacobian.allFinite()) return std::nullopt;

    return SharedOffsetResponse{
        evaluation.response,
        shape_jacobian,
        evaluation.offset_basis
    };
}

std::optional<SharedOffsetResponse> EvaluateSharedOffsetResponse(
    const GaussianModel3D & model,
    double distance)
{
    if (!EncodeTransformedCoordinates(model).has_value()) return std::nullopt;
    return EvaluateValidSharedOffsetResponse(model, distance);
}

std::optional<TransformedModelInvariants> BuildTransformedModelInvariants(const GaussianModel3D & model)
{
    const auto transformed{ EncodeTransformedCoordinates(model) };
    if (!transformed.has_value()) return std::nullopt;

    const auto peak_height{
        std::exp((*transformed)(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)))
    };
    if (!std::isfinite(peak_height)) return std::nullopt;

    return TransformedModelInvariants{ model, peak_height };
}

std::optional<Eigen::Vector3d> EvaluateTransformedJacobian(
    const TransformedModelInvariants & invariants,
    double distance)
{
    const auto shared_offset_evaluation{
        EvaluateValidSharedOffsetResponse(invariants.model, distance)
    };
    if (!shared_offset_evaluation.has_value()) return std::nullopt;

    const auto & model{ invariants.model };
    const auto width{ model.GetWidth() };
    const double center_offset_basis_scale{ std::sqrt(2.0 / M_PI) };
    Eigen::Vector3d jacobian{ Eigen::Vector3d::Zero() };
    jacobian(static_cast<Eigen::Index>(kLogPeakHeightChangeIndex)) =
        shared_offset_evaluation->shape_jacobian(0) +
        model.GetOffset() * shared_offset_evaluation->offset_jacobian;
    jacobian(static_cast<Eigen::Index>(kLogWidthChangeIndex)) =
        shared_offset_evaluation->shape_jacobian(1) +
        model.GetOffset() * shared_offset_evaluation->offset_jacobian;
    jacobian(static_cast<Eigen::Index>(kOffsetToPeakRatioChangeIndex)) =
        invariants.peak_height * width *
        shared_offset_evaluation->offset_jacobian / center_offset_basis_scale;
    if (!jacobian.allFinite()) return std::nullopt;

    return jacobian;
}

RHBMExecutionOptions MakeExecutionOptions(int thread_size)
{
    return RHBMExecutionOptions{
        .quiet_mode = false,
        .thread_size = thread_size
    };
}

static float CalculateAdjustedResponse(
    double sample_response,
    double distance,
    const GaussianModel3D & offset_model)
{
    const auto evaluation{ offset_model.EvaluateAtDistance(distance) };
    return static_cast<float>(sample_response - (evaluation.response - evaluation.signal));
}

LocalPotentialSampleList BuildSamplesForZeroOffsetGaussianFit(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model)
{
    LocalPotentialSampleList adjusted_sampling_entries;
    adjusted_sampling_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto response{
            CalculateAdjustedResponse(
                static_cast<double>(sample.response),
                distance,
                model)
        };
        adjusted_sampling_entries.emplace_back(LocalPotentialSample{ response, sample.point });
    }
    return adjusted_sampling_entries;
}

static LocalGaussianResult DecodeLocalGaussianResult(
    double alpha_r,
    const RHBMBetaEstimateResult & fit_result,
    double offset)
{
    const auto ols_model{
        linearization_service::DecodeParameterVector(fit_result.beta_ols).WithOffset(offset)
    };
    const auto mdpde_model{
        linearization_service::DecodeParameterVector(fit_result.beta_mdpde).WithOffset(offset)
    };
    return LocalGaussianResult{
        alpha_r,
        GaussianModel3DWithUncertainty{ ols_model, GaussianModel3DUncertainty{} },
        GaussianModel3DWithUncertainty{ mdpde_model, GaussianModel3DUncertainty{} },
        std::nullopt,
        false,
        0.0,
        fit_result
    };
}

LocalGaussianDesignTemplate BuildLocalGaussianDesignTemplate(
    const LocalPotentialSampleList & sample_entries,
    double range_min,
    double range_max)
{
    numeric_validation::RequireFiniteNonNegativeRange(range_min, range_max, "fit range");

    LocalGaussianDesignTemplate design_template;
    design_template.source_sample_count = sample_entries.size();
    design_template.source_sample_index_list.reserve(sample_entries.size());
    design_template.distance_list.reserve(sample_entries.size());
    for (std::size_t sample_index = 0; sample_index < sample_entries.size(); sample_index++)
    {
        const auto distance{
            static_cast<double>(sample_entries.at(sample_index).point.distance)
        };
        if (distance < range_min || distance > range_max) continue;
        design_template.source_sample_index_list.emplace_back(sample_index);
        design_template.distance_list.emplace_back(distance);
    }

    const auto row_count{
        static_cast<Eigen::Index>(design_template.distance_list.size())
    };
    design_template.design_matrix = RHBMDesignMatrix::Zero(row_count, 2);
    for (Eigen::Index row = 0; row < row_count; row++)
    {
        const auto distance{
            design_template.distance_list.at(static_cast<std::size_t>(row))
        };
        design_template.design_matrix(row, 0) = 1.0;
        design_template.design_matrix(row, 1) = -0.5 * distance * distance;
    }
    return design_template;
}

RHBMMemberDataset BuildLocalGaussianPreparedDataset(
    const LocalGaussianDesignTemplate & design_template,
    const std::vector<double> & sample_response_list,
    const GaussianModel3D & offset_model)
{
    numeric_validation::RequireFinite(offset_model.GetOffset(), "offset");
    if (sample_response_list.size() != design_template.source_sample_count ||
        design_template.source_sample_index_list.size() != design_template.distance_list.size() ||
        design_template.design_matrix.rows() != static_cast<Eigen::Index>(design_template.distance_list.size()) ||
        design_template.design_matrix.cols() != 2)
    {
        throw std::invalid_argument("Prepared local Gaussian design inputs are inconsistent.");
    }

    const auto candidate_count{
        static_cast<Eigen::Index>(design_template.distance_list.size())
    };
    RHBMMemberDataset dataset;
    dataset.X = RHBMDesignMatrix::Zero(candidate_count, 2);
    dataset.y = RHBMResponseVector::Zero(candidate_count);
    Eigen::Index retained_count{ 0 };
    for (std::size_t row = 0; row < design_template.distance_list.size(); row++)
    {
        const auto sample_index{ design_template.source_sample_index_list.at(row) };
        if (sample_index >= sample_response_list.size())
        {
            throw std::invalid_argument("Prepared local Gaussian sample index is out of range.");
        }
        const auto adjusted_response{
            static_cast<double>(CalculateAdjustedResponse(
                sample_response_list.at(sample_index),
                design_template.distance_list.at(row),
                offset_model))
        };
        if (adjusted_response <= 0.0) continue;
        numeric_validation::RequireFinite(
            adjusted_response,
            "response",
            "Member dataset contains non-finite value.");
        dataset.X.row(retained_count) = design_template.design_matrix.row(
            static_cast<Eigen::Index>(row));
        dataset.y(retained_count) = std::log(adjusted_response);
        retained_count++;
    }

    if (retained_count == 0)
    {
        dataset.X = RHBMDesignMatrix::Zero(1, 2);
        dataset.y = RHBMResponseVector::Zero(1);
    }
    else
    {
        dataset.X.conservativeResize(retained_count, 2);
        dataset.y.conservativeResize(retained_count);
    }

    return dataset;
}

LocalGaussianResult EstimateLocalGaussianPrepared(
    const LocalGaussianDesignTemplate & design_template,
    const std::vector<double> & sample_response_list,
    double alpha_r,
    int thread_size,
    const GaussianModel3D & offset_model)
{
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");
    auto dataset{
        BuildLocalGaussianPreparedDataset(design_template, sample_response_list, offset_model)
    };
    const auto result{
        rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, MakeExecutionOptions(thread_size))
    };
    return DecodeLocalGaussianResult(alpha_r, result, offset_model.GetOffset());
}

const GaussianModel3D & ResolveNeighborAtomModel(
    const NeighborAtomSample & neighbor_atom_sample,
    const SecondStageModelSnapshot & model_snapshot)
{
    return neighbor_atom_sample.is_selected ?
        GetFitModel(model_snapshot.selected, neighbor_atom_sample.atom_index) :
        GetFitModel(model_snapshot.unselected, neighbor_atom_sample.atom_index);
}

static FittedGaussianSnapshot BuildUnselectedAtomContributorSnapshot(
    const SecondStageContext & context,
    const FittedGaussianSnapshot & selected_snapshot)
{
    if (selected_snapshot.size() != context.size())
    {
        throw std::invalid_argument(
            "Second-stage selected contributor snapshot size is inconsistent.");
    }

    std::vector<std::optional<GaussianModel3D>> median_model_by_group(
        context.selected_atom_index_list_by_group.size());
    std::vector<GaussianModel3D> model_list;
    for (std::size_t group_id = 0;
        group_id < context.selected_atom_index_list_by_group.size();
        group_id++)
    {
        const auto & atom_index_list{
            context.selected_atom_index_list_by_group.at(group_id)
        };
        model_list.clear();
        model_list.reserve(atom_index_list.size());
        for (const auto atom_index : atom_index_list)
        {
            model_list.emplace_back(GetFitModel(selected_snapshot, atom_index));
        }
        median_model_by_group.at(group_id) = BuildGaussianParameterMedian(model_list);
    }

    FittedGaussianSnapshot snapshot;
    snapshot.reserve(context.unselected_atom_list.size());
    for (const auto & unselected_atom_contributor : context.unselected_atom_list)
    {
        if (unselected_atom_contributor.selected_group_id.has_value() &&
            median_model_by_group.at(*unselected_atom_contributor.selected_group_id).has_value())
        {
            snapshot.emplace_back(
                *median_model_by_group.at(*unselected_atom_contributor.selected_group_id));
            continue;
        }
        if (!IsValidSecondStageGaussianModel(unselected_atom_contributor.initial_seed))
        {
            throw std::logic_error(
                "Second-stage unselected contributor seed is unavailable.");
        }
        snapshot.emplace_back(unselected_atom_contributor.initial_seed);
    }
    return snapshot;
}

SecondStageModelSnapshot BuildSecondStageModelSnapshot(
    const SecondStageContext & context,
    FittedGaussianSnapshot selected_snapshot)
{
    auto unselected_snapshot{
        BuildUnselectedAtomContributorSnapshot(context, selected_snapshot)
    };
    return SecondStageModelSnapshot{
        std::move(selected_snapshot),
        std::move(unselected_snapshot)
    };
}

SecondStageModelSnapshot BuildSecondStageModelSnapshot(
    const SecondStageContext & context,
    const FitState & state)
{
    if (context.size() != state.size())
    {
        throw std::invalid_argument("Local fitting context and state sizes are inconsistent.");
    }
    return BuildSecondStageModelSnapshot(context, BuildFittedGaussianSnapshot(state));
}

SecondStageAdjustedResponseCache BuildSecondStageAdjustedResponseCache(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot)
{
    SecondStageAdjustedResponseCache cache(context.size());
    for (std::size_t i = 0; i < context.size(); i++)
    {
        const auto & atom_context{ context.at(i) };
        const auto sample_count{ atom_context.raw_sampling_entries.size() };
        auto & response_list{ cache.at(i) };
        response_list.reserve(sample_count);
        for (std::size_t j = 0; j < sample_count; j++)
        {
            response_list.emplace_back(CalculateSecondStageAdjustedResponse(atom_context, j, model_snapshot));
        }
    }
    return cache;
}

LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const AtomContext & atom_context,
    const std::vector<double> & adjusted_response_list)
{
    if (adjusted_response_list.size() != atom_context.raw_sampling_entries.size())
    {
        throw std::invalid_argument("Second-stage adjusted response count is inconsistent.");
    }
    LocalPotentialSampleList adjusted_sampling_entries;
    adjusted_sampling_entries.reserve(adjusted_response_list.size());
    for (std::size_t i = 0; i < adjusted_response_list.size(); i++)
    {
        auto sample{ atom_context.raw_sampling_entries.at(i) };
        sample.response = static_cast<float>(adjusted_response_list.at(i));
        adjusted_sampling_entries.emplace_back(sample);
    }
    return adjusted_sampling_entries;
}

LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const AtomContext & atom_context,
    const SecondStageModelSnapshot & model_snapshot)
{
    LocalPotentialSampleList adjusted_sampling_entries;
    adjusted_sampling_entries.reserve(atom_context.raw_sampling_entries.size());
    for (std::size_t i = 0; i < atom_context.raw_sampling_entries.size(); i++)
    {
        auto sample{ atom_context.raw_sampling_entries.at(i) };
        sample.response = static_cast<float>(CalculateSecondStageAdjustedResponse(atom_context, i, model_snapshot));
        adjusted_sampling_entries.emplace_back(sample);
    }
    return adjusted_sampling_entries;
}

ResidualBaseline BuildResidualBaseline(const SecondStageContext & context, const FitState & state)
{
    ResidualBaseline baseline{
        BuildSecondStageModelSnapshot(context, state),
        std::vector<std::vector<std::optional<ResidualSample>>>(context.size())
    };
    for (std::size_t i = 0; i < context.size(); i++)
    {
        const auto & atom_context{ context.at(i) };
        const auto sample_count{ atom_context.raw_sampling_entries.size() };
        baseline.sample_list.at(i).reserve(sample_count);
        for (std::size_t j = 0; j < sample_count; j++)
        {
            baseline.sample_list.at(i).emplace_back(
                EvaluateResidualSample(
                    context,
                    SampleRef{ i, j },
                    baseline.model_snapshot));
        }
    }
    return baseline;
}

static algorithm::ParameterChange MakeInfiniteTransformedChange()
{
    return algorithm::ParameterChange{
        std::vector<double>(
            kTransformedChangeSize,
            std::numeric_limits<double>::infinity())
    };
}

algorithm::ParameterChange CalculateTransformedChange(
    const GaussianModel3D & current,
    const GaussianModel3D & previous)
{
    const auto current_coordinates{
        EncodeTransformedCoordinates(current)
    };
    const auto previous_coordinates{
        EncodeTransformedCoordinates(previous)
    };
    if (!current_coordinates.has_value() || !previous_coordinates.has_value())
    {
        return MakeInfiniteTransformedChange();
    }

    algorithm::ParameterChange change{
        std::vector<double>(kTransformedChangeSize, 0.0)
    };
    for (std::size_t i = 0; i < kTransformedChangeSize; i++)
    {
        const auto eigen_index{ static_cast<Eigen::Index>(i) };
        const auto value{
            std::abs(
                (*current_coordinates)(eigen_index) -
                (*previous_coordinates)(eigen_index))
        };
        if (!std::isfinite(value))
        {
            return MakeInfiniteTransformedChange();
        }
        change.value_list.at(i) = value;
    }
    return change;
}

double GetMaximumTransformedChange(const std::vector<double> & value_list)
{
    if (value_list.empty()) return 0.0;
    if (value_list.size() != kTransformedChangeSize)
    {
        throw std::invalid_argument(
            "Local fitting transformed change input is inconsistent.");
    }
    return std::ranges::max(value_list);
}

bool IsTransformedChangeMaterial(
    const algorithm::ParameterChange & change,
    double minimum_change)
{
    if (!std::isfinite(minimum_change) || minimum_change < 0.0)
    {
        throw std::invalid_argument(
            "Local fitting transformed change threshold is invalid.");
    }
    if (change.value_list.size() != kTransformedChangeSize)
    {
        throw std::invalid_argument(
            "Local fitting transformed change input is inconsistent.");
    }
    for (const auto value : change.value_list)
    {
        if (std::isfinite(value) && value >= minimum_change) return true;
    }
    return false;
}

double GetMaximumTransformedChange(const TransformedChangeSummary & summary)
{
    return GetMaximumTransformedChange(summary.maximum_list);
}

bool IsTransformedPercentileConverged(const TransformedChangeSummary & summary)
{
    for (const auto value : summary.percentile_stats.percentile_list)
    {
        if (!std::isfinite(value) || value >= kTransformedChangeTolerance) return false;
    }
    return true;
}

bool IsTrustRegionStepWithinRadius(double step_norm, double radius)
{
    return std::isfinite(step_norm) &&
        std::isfinite(radius) &&
        radius > 0.0 &&
        step_norm <= radius + kTrustRegionBoundaryTolerance;
}

static double CalculateScaledTransformedStepNorm(
    const Eigen::Vector3d & previous_estimation,
    const Eigen::Vector3d & candidate_estimation)
{
    double step_norm{ 0.0 };
    for (std::size_t parameter_index = 0;
        parameter_index < kTransformedChangeSize;
        parameter_index++)
    {
        const auto eigen_index{ static_cast<Eigen::Index>(parameter_index) };
        step_norm = std::max(
            step_norm,
            std::abs(
                candidate_estimation(eigen_index) -
                previous_estimation(eigen_index)) /
                kTrustRegionParameterScale.at(parameter_index));
    }
    return step_norm;
}

std::optional<double> CalculateModelTrustRegionStepNorm(
    const std::vector<GaussianModel3D> & previous_model_list,
    const std::vector<GaussianModel3D> & candidate_model_list)
{
    if (candidate_model_list.size() != previous_model_list.size())
    {
        return std::nullopt;
    }
    double step_norm{ 0.0 };
    for (std::size_t atom_position = 0; atom_position < previous_model_list.size(); atom_position++)
    {
        const auto previous{
            EncodeTransformedCoordinates(previous_model_list.at(atom_position))
        };
        const auto candidate{
            EncodeTransformedCoordinates(candidate_model_list.at(atom_position))
        };
        if (!previous.has_value() || !candidate.has_value())
        {
            return std::nullopt;
        }
        step_norm = std::max(
            step_norm,
            CalculateScaledTransformedStepNorm(*previous, *candidate));
    }
    return std::isfinite(step_norm) ? std::optional<double>{ step_norm } : std::nullopt;
}

} // namespace rhbm_gem::core::detail
