#pragma once

#include "WidthBounds.hpp"
#include <rhbm_gem/utils/domain/JointEstimationTypes.hpp>
#include <optional>
#include <string>
#include <vector>

namespace rhbm_gem { class MapObject; }

namespace rhbm_gem::core::joint_component {

struct WidthBoundCalibration
{
    std::string status{"CalibrationUnavailable"};
    std::string unavailable_reason;
    std::size_t target_count{}, qualified_donor_count{};
    double qualified_donor_coverage{};
    std::optional<double> median_width;
    std::optional<double> tau_max4, eta_max4;
};

struct BoundedInitialization
{
    bool available{};
    std::string reason;
    double tau_max{};
    WidthBounds bounds;
    std::vector<double> tau_values;
    std::size_t projection_count{}, lower_projection_count{};
    double maximum_projection_magnitude{};
};

// Calibration consumes only the captured InitialStage/initialization record. It
// deliberately has no simulation-truth input.
WidthBoundCalibration CalibrateGlobalWidth(const JointInitialization &);
std::optional<double> IsotropicVoxelSpacing(const MapObject &);

// Builds a feasible working vector without changing the captured InitialStage
// estimates or the caller's JointInitialization.
BoundedInitialization MakeBoundedInitialization(
    const JointInitialization &, const WidthBoundCalibration &, std::optional<double> tau_min = {});

} // namespace rhbm_gem::core::joint_component
