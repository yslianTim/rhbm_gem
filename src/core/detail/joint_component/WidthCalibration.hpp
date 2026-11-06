#pragma once

#include <rhbm_gem/utils/domain/JointEstimationTypes.hpp>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace rhbm_gem::core::joint_component {

struct WidthBoundCalibration
{
    std::string dataset_identity;
    std::string status{"CalibrationUnavailable"};
    std::string unavailable_reason;
    std::size_t target_count{}, qualified_donor_count{};
    double qualified_donor_coverage{};
    std::size_t native_success_count{}, refined_success_count{}, unqualified_count{};
    std::size_t invalid_width_count{}, fallback_count{}, initialization_failure_count{};
    std::map<std::string,std::size_t> excluded_donor_reasons;
    std::optional<double> minimum_width, median_width, p90, p95, p99, maximum_width;
    std::optional<double> tau_max4, eta_max4;
    std::size_t initial_violations4{};
};

struct InitialTauProjection
{
    std::size_t joint_index{};
    std::string atom_id;
    std::optional<double> original_first_stage_tau;
    double original_initial_tau{}, bounded_initial_tau{};
    bool projection_applied{};
};

struct BoundedInitialization
{
    bool available{};
    std::string reason;
    double tau_max{};
    std::vector<double> tau_values;
    std::vector<InitialTauProjection> projections;
    std::size_t projection_count{};
    double maximum_projection_magnitude{};
};

// Calibration consumes only the captured FirstStage/initialization record. It
// deliberately has no simulation-truth input.
WidthBoundCalibration CalibrateGlobalWidth(const JointInitialization &, std::string dataset_identity);

// Shared math after caller-specific eligibility: donors must be finite positive
// qualified widths. Preserves the donor minimum, coverage, Type-7 and 4*p95 rule.
WidthBoundCalibration CalibrateQualifiedWidthDonors(WidthBoundCalibration, std::vector<double> donors);

// Builds a feasible working vector without changing the captured FirstStage
// estimates or the caller's JointInitialization.
BoundedInitialization MakeBoundedInitialization(
    const JointInitialization &, const WidthBoundCalibration &);

} // namespace rhbm_gem::core::joint_component
