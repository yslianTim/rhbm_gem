#pragma once

#include <rhbm_gem/utils/domain/JointEstimationTypes.hpp>
#include <optional>
#include <string>
#include <vector>

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
    std::vector<double> tau_values;
    std::size_t projection_count{};
    double maximum_projection_magnitude{};
};

// Calibration consumes only the captured InitialStage/initialization record. It
// deliberately has no simulation-truth input.
WidthBoundCalibration CalibrateGlobalWidth(const JointInitialization &);

// Builds a feasible working vector without changing the captured InitialStage
// estimates or the caller's JointInitialization.
BoundedInitialization MakeBoundedInitialization(
    const JointInitialization &, const WidthBoundCalibration &);

} // namespace rhbm_gem::core::joint_component
