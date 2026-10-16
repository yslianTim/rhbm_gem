#include "JointHaloObservationExclusion.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace joint_test::halo_exclusion {
namespace {
std::size_t Shell(double square)
{
    const double distance = std::sqrt(square);
    for (std::size_t shell = 0; shell < 4; ++shell)
        if (distance < RadialShellEdges[shell + 1]) return shell;
    return 4;
}

class DisjointSets
{
    std::vector<std::size_t> parent_, rank_;
public:
    explicit DisjointSets(std::size_t count) : parent_(count), rank_(count)
    { std::iota(parent_.begin(), parent_.end(), 0); }
    std::size_t Find(std::size_t x)
    { return parent_[x] == x ? x : parent_[x] = Find(parent_[x]); }
    void Join(std::size_t a, std::size_t b)
    {
        a = Find(a); b = Find(b);
        if (a == b) return;
        if (rank_[a] < rank_[b]) std::swap(a, b);
        parent_[b] = a;
        if (rank_[a] == rank_[b]) ++rank_[a];
    }
};

}

std::vector<double> NearestHaloDistances(std::span<const Position> observations,
    std::span<const Position> halos)
{
    std::vector<double> result;
    result.reserve(observations.size());
    for (const auto & observation : observations)
    {
        double nearest_square = std::numeric_limits<double>::infinity();
        for (const auto & halo : halos)
        {
            double square{};
            for (std::size_t axis = 0; axis < 3; ++axis)
            {
                const double delta = observation[axis] - halo[axis];
                square += delta * delta;
            }
            nearest_square = std::min(nearest_square, square);
        }
        result.push_back(std::sqrt(nearest_square));
    }
    return result;
}

RowMask RetainedRows(std::span<const double> nearest_halo_distances, double radius)
{
    if (!std::isfinite(radius) || radius < 0.0)
        throw std::invalid_argument("Halo exclusion radius must be finite and nonnegative.");
    RowMask retained;
    retained.reserve(nearest_halo_distances.size());
    for (const double distance : nearest_halo_distances)
    {
        if (std::isnan(distance) || distance < 0.0)
            throw std::invalid_argument("Invalid nearest-halo distance.");
        retained.push_back(distance >= radius);
    }
    return retained;
}

CoverageSummary AnalyzeCoverage(const rhbm_gem::core::JointProblemInput & input,
    std::span<const std::uint8_t> retained_rows, std::span<const double> target_halo_distances,
    std::span<const double> observation_halo_distances)
{
    if (!input.selection_domain || input.selection_domain->target_indices.empty() ||
        retained_rows.size() != input.observations.size() ||
        observation_halo_distances.size() != input.observations.size() ||
        target_halo_distances.size() != input.selection_domain->target_indices.size())
        throw std::invalid_argument("Invalid halo exclusion coverage dimensions.");

    CoverageSummary result;
    result.original_rows = input.observations.size();
    result.targets.resize(target_halo_distances.size());
    std::vector<std::uint8_t> target_row(result.original_rows), halo_row(result.original_rows);
    std::vector<std::size_t> target_incidence(result.original_rows);
    std::vector<std::size_t> row_owner(result.original_rows, std::numeric_limits<std::size_t>::max());
    DisjointSets components(result.targets.size());
    for (std::size_t row = 0; row < retained_rows.size(); ++row)
    {
        if (retained_rows[row]) ++result.retained_rows;
        else ++result.excluded_rows;
        const double distance = observation_halo_distances[row];
        const std::size_t bin = distance < 0.5 ? 0 : distance < 1.0 ? 1 : distance < 1.5 ? 2 :
            distance < 2.0 ? 3 : distance < 2.5 ? 4 : 5;
        ++result.observation_halo_distance_bins[bin];
    }

    const auto & targets = input.selection_domain->target_indices;
    for (std::size_t target = 0; target < targets.size(); ++target)
    {
        auto & coverage = result.targets[target];
        coverage.contributor_index = targets[target];
        coverage.atom_id = input.atom_ids.at(targets[target]);
        coverage.nearest_halo_distance = target_halo_distances[target];
        const auto & support = input.support.at(targets[target]);
        std::array<bool, 5> has_shell{};
        for (const auto & entry : support)
        {
            const auto row = entry.row;
            ++coverage.original_rows;
            ++coverage.radial[Shell(entry.squared_distance)].original;
            target_row.at(row) = 1;
            ++target_incidence.at(row);
            if (retained_rows[row])
            {
                ++coverage.retained_rows;
                ++coverage.radial[Shell(entry.squared_distance)].retained;
                has_shell[Shell(entry.squared_distance)] = true;
                if (row_owner[row] == std::numeric_limits<std::size_t>::max()) row_owner[row] = target;
                else components.Join(target, row_owner[row]);
            }
        }
        if (coverage.retained_rows == 0) ++result.targets_without_retained_rows;
        const auto occupied_shells = static_cast<std::size_t>(std::count(has_shell.begin(), has_shell.end(), true));
        result.target_radial_diversity[occupied_shells]++;
        for (std::size_t shell = 0; shell < coverage.radial.size(); ++shell)
        {
            result.radial[shell].original += coverage.radial[shell].original;
            result.radial[shell].retained += coverage.radial[shell].retained;
        }
        const double distance = coverage.nearest_halo_distance;
        const std::size_t proximity = distance < 1.0 ? 0 : distance < 1.5 ? 1 : distance < 2.0 ? 2 : distance < 2.5 ? 3 : 4;
        ++result.halo_proximity_bins[proximity];
    }

    for (std::size_t atom = 0; atom < input.support.size(); ++atom)
    {
        if (std::binary_search(targets.begin(), targets.end(), atom)) continue;
        for (const auto & entry : input.support[atom]) halo_row.at(entry.row) = 1;
    }
    for (std::size_t row = 0; row < result.original_rows; ++row)
    {
        result.max_targets_per_row = std::max(result.max_targets_per_row, target_incidence[row]);
        if (target_incidence[row] > 1)
        {
            ++result.rows_shared_by_targets;
            if (retained_rows[row]) ++result.retained_rows_shared_by_targets;
        }
        if (target_row[row] && halo_row[row])
        {
            ++result.mixed_target_halo_rows;
            if (!retained_rows[row]) ++result.excluded_mixed_target_halo_rows;
        }
    }
    for (std::size_t target = 0; target < result.targets.size(); ++target)
        if (components.Find(target) == target) ++result.target_components;
    return result;
}

}
