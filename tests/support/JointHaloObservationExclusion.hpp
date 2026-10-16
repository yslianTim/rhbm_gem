#pragma once

#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <Eigen/Core>
#include <array>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace joint_test::halo_exclusion {

inline constexpr std::array<double, 6> DefaultRadii{0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
inline constexpr std::array<double, 6> RadialShellEdges{0.0, 0.5, 1.0, 1.5, 2.0, 2.5};

using Position = std::array<double, 3>;
using RowMask = std::vector<std::uint8_t>;

struct RadialCoverage
{
    std::size_t original{}, retained{};
};

struct TargetCoverage
{
    std::size_t contributor_index{};
    std::string atom_id;
    double nearest_halo_distance{};
    std::size_t original_rows{}, retained_rows{};
    std::array<RadialCoverage, 5> radial{};
};

struct CoverageSummary
{
    std::size_t original_rows{}, excluded_rows{}, retained_rows{};
    std::size_t mixed_target_halo_rows{}, excluded_mixed_target_halo_rows{};
    std::size_t target_components{}, targets_without_retained_rows{};
    std::size_t rows_shared_by_targets{}, retained_rows_shared_by_targets{}, max_targets_per_row{};
    std::array<std::size_t, 6> target_radial_diversity{};
    std::array<std::size_t, 5> halo_proximity_bins{};
    std::array<std::size_t, 6> observation_halo_distance_bins{};
    std::array<RadialCoverage, 5> radial{};
    std::vector<TargetCoverage> targets;
};

struct DesignSpectrum
{
    std::size_t rows{}, columns{}, rank{};
    double minimum_singular_value{}, maximum_singular_value{};
    double condition{};
    bool available{}, full_rank{}, ill_conditioned{};
};

struct LocalDesignMetrics
{
    DesignSpectrum ac, abc, projected_width;
    double minimum_projected_width_fraction{};
};

struct MaskedInput
{
    rhbm_gem::core::JointProblemInput input;
    // Each compact fitting row maps back to its stable row in the original domain.
    std::vector<std::size_t> source_rows;
    std::vector<std::size_t> source_contributors;
};

struct NoiseSummary
{
    double mean{}, standard_deviation{};
};

std::vector<double> NearestHaloDistances(std::span<const Position> observations,
    std::span<const Position> halos);
RowMask RetainedRows(std::span<const double> nearest_halo_distances, double radius);
// A row is excluded exactly when a non-target contributor has support on it.
RowMask HaloSupportExclusionRows(const rhbm_gem::core::JointProblemInput &);
NoiseSummary FillGaussianNoise(std::span<const double> source, std::span<double> destination,
    double sigma, std::uint64_t seed);
MaskedInput BuildMaskedInput(const rhbm_gem::core::JointProblemInput &,
    std::span<const std::uint8_t> retained_rows, bool retain_halos);
CoverageSummary AnalyzeCoverage(const rhbm_gem::core::JointProblemInput &,
    std::span<const std::uint8_t> retained_rows, std::span<const double> target_halo_distances,
    std::span<const double> observation_halo_distances);
LocalDesignMetrics AnalyzeLocalDesign(const Eigen::MatrixXd & ac, const Eigen::MatrixXd & log_width);

}
