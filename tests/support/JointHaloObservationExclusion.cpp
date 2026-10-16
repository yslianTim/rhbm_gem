#include "JointHaloObservationExclusion.hpp"
#include "core/detail/joint_component/LocalConditioning.hpp"
#include <algorithm>
#include <cmath>
#include <Eigen/SVD>
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

Eigen::MatrixXd NormalizeColumns(const Eigen::MatrixXd & matrix)
{
    Eigen::MatrixXd normalized = matrix;
    for (Eigen::Index column = 0; column < matrix.cols(); ++column)
    {
        const double norm = normalized.col(column).norm();
        if (norm > 0.0) normalized.col(column) /= norm;
    }
    return normalized;
}

DesignSpectrum Spectrum(const Eigen::MatrixXd & matrix)
{
    DesignSpectrum result;
    result.rows = static_cast<std::size_t>(matrix.rows());
    result.columns = static_cast<std::size_t>(matrix.cols());
    if (matrix.rows() == 0 || matrix.cols() == 0 || !matrix.allFinite()) return result;
    const Eigen::MatrixXd normalized = NormalizeColumns(matrix);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(normalized);
    if (svd.info() != Eigen::Success || svd.singularValues().size() == 0 || !svd.singularValues().allFinite())
        return result;
    const auto singular = svd.singularValues();
    const double maximum = singular(0);
    const double threshold = std::numeric_limits<double>::epsilon() *
        static_cast<double>(std::max(normalized.rows(), normalized.cols())) * maximum;
    result.available = true;
    result.maximum_singular_value = maximum;
    result.minimum_singular_value = normalized.rows() < normalized.cols() ? 0.0 : singular(singular.size() - 1);
    result.rank = static_cast<std::size_t>((singular.array() > threshold).count());
    result.full_rank = normalized.rows() >= normalized.cols() && result.rank == result.columns &&
        result.minimum_singular_value > 0.0;
    if (result.full_rank)
    {
        result.condition = maximum / result.minimum_singular_value;
        result.ill_conditioned = result.condition * std::numeric_limits<double>::epsilon() >=
            rhbm_gem::core::joint_component::LocalConditioningRoundoffAmplificationWarning;
    }
    return result;
}
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

MaskedInput BuildMaskedInput(const rhbm_gem::core::JointProblemInput & source,
    std::span<const std::uint8_t> retained_rows, bool retain_halos)
{
    if (!source.selection_domain || retained_rows.size() != source.observations.size() ||
        source.support.size() != source.atom_ids.size())
        throw std::invalid_argument("Invalid masked joint problem dimensions.");

    MaskedInput result;
    auto & input = result.input;
    std::vector<std::size_t> row_map(source.observations.size(), std::numeric_limits<std::size_t>::max());
    for (std::size_t row = 0; row < retained_rows.size(); ++row)
    {
        if (!retained_rows[row]) continue;
        row_map[row] = input.observations.size();
        input.observations.push_back(source.observations[row]);
        input.row_ids.push_back(source.row_ids[row]);
    }

    const auto & targets = source.selection_domain->target_indices;
    std::vector<bool> is_target(source.atom_ids.size());
    for (const auto target : targets) is_target.at(target) = true;
    std::vector<std::size_t> source_to_masked(source.atom_ids.size(), std::numeric_limits<std::size_t>::max());
    for (std::size_t contributor = 0; contributor < source.atom_ids.size(); ++contributor)
    {
        std::vector<rhbm_gem::core::JointSupport> support;
        support.reserve(source.support[contributor].size());
        for (const auto & entry : source.support[contributor])
            if (retained_rows[entry.row]) support.push_back({row_map[entry.row], entry.squared_distance});
        if (!is_target[contributor] && (!retain_halos || support.empty())) continue;

        source_to_masked[contributor] = input.atom_ids.size();
        result.source_contributors.push_back(contributor);
        input.atom_ids.push_back(source.atom_ids[contributor]);
        input.support.push_back(std::move(support));
    }

    input.selection_domain.emplace();
    input.selection_domain->target_indices.reserve(targets.size());
    for (const auto target : targets)
        input.selection_domain->target_indices.push_back(source_to_masked.at(target));
    if (retain_halos)
    {
        const auto & prior = source.selection_domain->fixed_halo_selection;
        if (!prior) throw std::invalid_argument("Original problem lacks fixed-halo selection provenance.");
        const auto retained_halo_count = input.atom_ids.size() - targets.size();
        input.selection_domain->contributor_policy = std::string(rhbm_gem::FixedHaloContributorPolicy);
        input.selection_domain->fixed_halo_selection = rhbm_gem::JointFixedHaloSelectionProvenance{
            1, prior->candidate_halos, retained_halo_count, prior->candidate_halos - retained_halo_count};
    }
    return result;
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

LocalDesignMetrics AnalyzeLocalDesign(const Eigen::MatrixXd & ac, const Eigen::MatrixXd & log_width)
{
    if (ac.rows() != log_width.rows() || ac.cols() != 2 * log_width.cols())
        throw std::invalid_argument("Invalid local A/C and width design dimensions.");
    LocalDesignMetrics result;
    result.ac = Spectrum(ac);
    const Eigen::MatrixXd normalized_ac = NormalizeColumns(ac);
    const Eigen::MatrixXd normalized_width = NormalizeColumns(log_width);
    Eigen::MatrixXd full(ac.rows(), ac.cols() + log_width.cols());
    full << normalized_ac, normalized_width;
    result.abc = Spectrum(full);
    if (ac.rows() == 0 || log_width.cols() == 0 || !normalized_ac.allFinite() || !normalized_width.allFinite())
        return result;
    Eigen::JacobiSVD<Eigen::MatrixXd> ac_svd(normalized_ac, Eigen::ComputeThinU);
    if (ac_svd.info() != Eigen::Success || !ac_svd.singularValues().allFinite()) return result;
    const double maximum = ac_svd.singularValues().size() ? ac_svd.singularValues()(0) : 0.0;
    const double threshold = std::numeric_limits<double>::epsilon() *
        static_cast<double>(std::max(normalized_ac.rows(), normalized_ac.cols())) * maximum;
    const auto ac_rank = static_cast<Eigen::Index>((ac_svd.singularValues().array() > threshold).count());
    Eigen::MatrixXd projected = normalized_width;
    if (ac_rank > 0)
    {
        const Eigen::MatrixXd q = ac_svd.matrixU().leftCols(ac_rank);
        projected.noalias() -= q * (q.transpose() * normalized_width);
    }
    result.minimum_projected_width_fraction = std::numeric_limits<double>::infinity();
    for (Eigen::Index column = 0; column < projected.cols(); ++column)
        result.minimum_projected_width_fraction = std::min(result.minimum_projected_width_fraction,
            projected.col(column).norm());
    result.projected_width = Spectrum(projected);
    return result;
}
}
