#include "core/detail/GaussianModelOperations.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rhbm_gem::core::detail {

namespace {

constexpr double kTransformedChangePercentile{ 0.99 };
constexpr double kTrustRegionBoundaryTolerance{ 1.0e-12 };
constexpr TransformedChange kTrustRegionParameterScale{ 0.50, 0.35, 1.0 };

using AtomPositionListByGroup = std::unordered_map<std::size_t, std::vector<std::size_t>>;

AtomPositionListByGroup BuildAtomPositionListByGroup(
    const std::vector<std::size_t> & group_id_by_atom_position,
    std::size_t model_count)
{
    if (group_id_by_atom_position.size() != model_count)
    {
        throw std::invalid_argument("Local fitting group median inputs are inconsistent.");
    }

    AtomPositionListByGroup atom_position_list_by_group;
    atom_position_list_by_group.reserve(group_id_by_atom_position.size());
    for (std::size_t atom_position = 0; atom_position < group_id_by_atom_position.size(); atom_position++)
    {
        atom_position_list_by_group[group_id_by_atom_position.at(atom_position)].emplace_back(atom_position);
    }
    return atom_position_list_by_group;
}

std::optional<SharedOffsetResponse> EvaluateValidSharedOffsetResponse(const GaussianModel3D & model, double distance)
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
    auto log_width_derivative{ evaluation.signal * normalized_distance * normalized_distance };
    if (distance < 1.0e-5)
    {
        log_width_derivative -= model.GetOffset() * center_offset_basis_scale / width;
    }
    else
    {
        const auto exponent{ -0.5 * normalized_distance * normalized_distance };
        log_width_derivative -= model.GetOffset() * center_offset_basis_scale / width * std::exp(exponent);
    }

    if (!std::isfinite(log_width_derivative)) return std::nullopt;

    return SharedOffsetResponse{
        evaluation.response,
        Eigen::Vector2d{ evaluation.signal, log_width_derivative },
        evaluation.offset_basis
    };
}

TransformedChange MakeInfiniteTransformedChange()
{
    TransformedChange change;
    change.fill(std::numeric_limits<double>::infinity());
    return change;
}

} // namespace

bool IsValidSecondStageGaussianModel(const GaussianModel3D & model)
{
    return model.ToTransformedCoordinates().has_value();
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

std::optional<GaussianModel3D> BuildGaussianParameterMedian(const std::vector<GaussianModel3D> & model_list)
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
    if (!IsValidSecondStageGaussianModel(median_model)) return std::nullopt;
    return median_model;
}

std::vector<GaussianModel3D> BuildGroupMedianModelList(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & model_list)
{
    const auto atom_position_list_by_group{
        BuildAtomPositionListByGroup(group_id_by_atom_position, model_list.size())
    };
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
    const auto atom_position_list_by_group{
        BuildAtomPositionListByGroup(group_id_by_atom_position, model_list.size())
    };
    std::vector<double> offset_list;
    offset_list.reserve(model_list.size());
    for (const auto & model : model_list)
    {
        offset_list.emplace_back(model.GetOffset());
    }

    std::vector<double> valid_offset_list;
    valid_offset_list.reserve(model_list.size());
    for (const auto & entry : atom_position_list_by_group)
    {
        const auto & atom_position_list{ entry.second };
        valid_offset_list.clear();
        for (const auto atom_position : atom_position_list)
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
        for (const auto atom_position : atom_position_list)
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
            previous_model_list.at(atom_position).ToTransformedCoordinates()
        };
        const auto raw_coordinates{
            raw_model_list.at(atom_position).ToTransformedCoordinates()
        };
        if (!previous_coordinates.has_value() || !raw_coordinates.has_value()) return std::nullopt;

        GaussianModel3D::TransformedCoordinates shape_coordinates{
            std::lerp(
                (*previous_coordinates)(static_cast<Eigen::Index>(
                    GaussianModel3D::LogPeakHeightCoordinateIndex())),
                (*raw_coordinates)(static_cast<Eigen::Index>(
                    GaussianModel3D::LogPeakHeightCoordinateIndex())),
                damping),
            std::lerp(
                (*previous_coordinates)(static_cast<Eigen::Index>(
                    GaussianModel3D::LogWidthCoordinateIndex())),
                (*raw_coordinates)(static_cast<Eigen::Index>(
                    GaussianModel3D::LogWidthCoordinateIndex())),
                damping),
            0.0
        };
        const auto shape_model{
            GaussianModel3D::FromTransformedCoordinates(shape_coordinates)
        };
        if (!shape_model.has_value()) return std::nullopt;

        const auto candidate_model{
            shape_model->WithOffset(
                std::lerp(
                    previous_shared_offset_list.at(atom_position),
                    raw_shared_offset_list.at(atom_position),
                    damping))
        };
        if (!IsValidSecondStageGaussianModel(candidate_model)) return std::nullopt;
        candidate_model_list.emplace_back(candidate_model);
    }
    return candidate_model_list;
}

std::optional<SharedOffsetResponse> EvaluateSharedOffsetResponse(const GaussianModel3D & model, double distance)
{
    if (!IsValidSecondStageGaussianModel(model)) return std::nullopt;
    return EvaluateValidSharedOffsetResponse(model, distance);
}

std::optional<TransformedModelInvariants> BuildTransformedModelInvariants(const GaussianModel3D & model)
{
    const auto transformed{ model.ToTransformedCoordinates() };
    if (!transformed.has_value()) return std::nullopt;

    const auto peak_height{
        std::exp((*transformed)(static_cast<Eigen::Index>(GaussianModel3D::LogPeakHeightCoordinateIndex())))
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
    jacobian(static_cast<Eigen::Index>(
        GaussianModel3D::LogPeakHeightCoordinateIndex())) =
        shared_offset_evaluation->shape_jacobian(0) +
        model.GetOffset() * shared_offset_evaluation->offset_jacobian;
    jacobian(static_cast<Eigen::Index>(
        GaussianModel3D::LogWidthCoordinateIndex())) =
        shared_offset_evaluation->shape_jacobian(1) +
        model.GetOffset() * shared_offset_evaluation->offset_jacobian;
    jacobian(static_cast<Eigen::Index>(
        GaussianModel3D::OffsetToPeakRatioCoordinateIndex())) =
        invariants.peak_height * width *
        shared_offset_evaluation->offset_jacobian / center_offset_basis_scale;
    if (!jacobian.allFinite()) return std::nullopt;

    return jacobian;
}

TransformedChange CalculateTransformedChange(
    const GaussianModel3D & current,
    const GaussianModel3D & previous)
{
    const auto current_coordinates{ current.ToTransformedCoordinates() };
    const auto previous_coordinates{ previous.ToTransformedCoordinates() };
    if (!current_coordinates.has_value() || !previous_coordinates.has_value())
    {
        return MakeInfiniteTransformedChange();
    }

    TransformedChange change{};
    for (std::size_t i = 0; i < change.size(); i++)
    {
        const auto eigen_index{ static_cast<Eigen::Index>(i) };
        const auto value{
            std::abs((*current_coordinates)(eigen_index) - (*previous_coordinates)(eigen_index))
        };
        if (!std::isfinite(value))
        {
            return MakeInfiniteTransformedChange();
        }
        change.at(i) = value;
    }
    return change;
}

bool IsTransformedChangeMaterial(const TransformedChange & change, double minimum_change)
{
    if (!std::isfinite(minimum_change) || minimum_change < 0.0)
    {
        throw std::invalid_argument("Local fitting transformed change threshold is invalid.");
    }
    for (const auto value : change)
    {
        if (std::isfinite(value) && value >= minimum_change) return true;
    }
    return false;
}

TransformedChangeSummary SummarizeTransformedChanges(const std::vector<TransformedChange> & change_list)
{
    TransformedChangeSummary summary;
    summary.population_size_list.fill(change_list.size());
    for (std::size_t parameter_index = 0;
        parameter_index < summary.percentile_list.size(); parameter_index++)
    {
        std::vector<double> parameter_change_list;
        parameter_change_list.reserve(change_list.size());
        for (const auto & change : change_list)
        {
            parameter_change_list.emplace_back(change.at(parameter_index));
        }
        summary.percentile_list.at(parameter_index) =
            array_helper::ComputePercentile(parameter_change_list, kTransformedChangePercentile);
        summary.maximum_list.at(parameter_index) =
            parameter_change_list.empty() ? 0.0 : *std::ranges::max_element(parameter_change_list);
    }
    return summary;
}

bool IsTransformedPercentileConverged(const TransformedChangeSummary & summary)
{
    for (const auto value : summary.percentile_list)
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

std::optional<double> CalculateModelTrustRegionStepNorm(
    const std::vector<GaussianModel3D> & previous_model_list,
    const std::vector<GaussianModel3D> & candidate_model_list)
{
    if (candidate_model_list.size() != previous_model_list.size()) return std::nullopt;
    double step_norm{ 0.0 };
    for (std::size_t atom_position = 0; atom_position < previous_model_list.size(); atom_position++)
    {
        const auto change{
            CalculateTransformedChange(
                candidate_model_list.at(atom_position),
                previous_model_list.at(atom_position))
        };
        for (std::size_t parameter_index = 0;
            parameter_index < kTrustRegionParameterScale.size(); parameter_index++)
        {
            step_norm = std::max(
                step_norm,
                change.at(parameter_index) / kTrustRegionParameterScale.at(parameter_index));
        }
    }
    return std::isfinite(step_norm) ? std::optional<double>{ step_norm } : std::nullopt;
}

} // namespace rhbm_gem::core::detail
