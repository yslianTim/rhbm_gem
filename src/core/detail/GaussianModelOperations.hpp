#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem::core::detail {

constexpr double kTransformedChangeTolerance{ 1.0e-4 };

using TransformedChange =
    std::array<double, GaussianModel3D::TransformedCoordinateSize()>;
bool IsValidSecondStageGaussianModel(const GaussianModel3D & model);

GaussianModel3DWithUncertainty WithPreservedUncertaintyOffset(
    const GaussianModel3DWithUncertainty & gaussian,
    double offset);

std::optional<GaussianModel3D> BuildGaussianParameterMedian(
    const std::vector<GaussianModel3D> & model_list);

std::vector<GaussianModel3D> BuildGroupMedianModelList(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & model_list);

std::vector<double> BuildGroupMedianOffsetList(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & model_list);

std::optional<std::vector<GaussianModel3D>> BuildSharedOffsetDampedModelList(
    const std::vector<GaussianModel3D> & previous_model_list,
    const std::vector<GaussianModel3D> & raw_model_list,
    const std::vector<double> & previous_shared_offset_list,
    const std::vector<double> & raw_shared_offset_list,
    double damping);

struct SharedOffsetResponse
{
    double response{ 0.0 };
    Eigen::Vector2d shape_jacobian{ Eigen::Vector2d::Zero() };
    double offset_jacobian{ 0.0 };
};

std::optional<SharedOffsetResponse> EvaluateSharedOffsetResponse(
    const GaussianModel3D & model,
    double distance);

struct TransformedModelInvariants
{
    GaussianModel3D model{};
    double peak_height{ 0.0 };
};

std::optional<TransformedModelInvariants> BuildTransformedModelInvariants(const GaussianModel3D & model);

std::optional<Eigen::Vector3d> EvaluateTransformedJacobian(
    const TransformedModelInvariants & invariants,
    double distance);

struct TransformedChangeSummary
{
    TransformedChange percentile_list{};
    TransformedChange maximum_list{};
    std::array<
        std::size_t,
        GaussianModel3D::TransformedCoordinateSize()> population_size_list{};
};

TransformedChange CalculateTransformedChange(
    const GaussianModel3D & current,
    const GaussianModel3D & previous);

bool IsTransformedChangeMaterial(const TransformedChange & change, double minimum_change);
TransformedChangeSummary SummarizeTransformedChanges(const std::vector<TransformedChange> & change_list);

bool IsTransformedPercentileConverged(const TransformedChangeSummary & summary);
bool IsTrustRegionStepWithinRadius(double step_norm, double radius);

std::optional<double> CalculateModelTrustRegionStepNorm(
    const std::vector<GaussianModel3D> & previous_model_list,
    const std::vector<GaussianModel3D> & candidate_model_list);

} // namespace rhbm_gem::core::detail
