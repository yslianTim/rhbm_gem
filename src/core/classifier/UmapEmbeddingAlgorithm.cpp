#include "core/classifier/UmapEmbeddingAlgorithm.hpp"

#include <rhbm_gem/utils/algorithm/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <knncolle/knncolle.hpp>
#include <umappp/umappp.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

namespace rhbm_gem::core::classifier {

namespace {

constexpr double kRadiusRatio{ 0.1 };
constexpr double kSubnormalCoordinateScale{ 0x1p54 };

struct UmapPoint
{
    std::array<double, kOutputDimensionCount> position;
    std::size_t index;

    const auto & GetPosition() const { return position; }
};

using Tree = KDTreeAlgorithm<UmapPoint>;

std::vector<UmapPoint> MakeUmapPoints(std::span<const double> embedding)
{
    std::vector<UmapPoint> points;
    points.reserve(embedding.size() / kOutputDimensionCount);
    for (std::size_t offset = 0; offset < embedding.size(); offset += kOutputDimensionCount)
    {
        points.push_back({ { embedding[offset], embedding[offset + 1] }, points.size() });
    }
    return points;
}

enum class MeanShiftStep { Updated, EmptyNeighborhood, ZeroWeight, NonFiniteUpdate };

MeanShiftStep ShiftUmapPoint(
    const KDNode<UmapPoint> * root,
    UmapPoint & query,
    double radius,
    std::vector<UmapPoint *> & neighbors)
{
    constexpr double alpha{ 2.0 };
    constexpr double bandwidth_ratio{ 0.5 };
    constexpr double exponent_factor{ alpha / (2.0 * bandwidth_ratio * bandwidth_ratio) };
    Tree::RangeSearch(root, &query, radius, neighbors);
    if (neighbors.empty()) return MeanShiftStep::EmptyNeighborhood;
    // At zero radius all neighbors coincide with the query, so their mean is its position.
    if (radius == 0.0) return MeanShiftStep::Updated;

    // A subnormal radius can have an unrepresentable reciprocal.
    const double inverse_radius{ radius >= std::numeric_limits<double>::min() ? 1.0 / radius : 0.0 };
    // Each source coordinate is at most radius away from the query. Leave headroom for summation.
    const double safe_query_magnitude{
        (std::numeric_limits<double>::max() / 4.0) / static_cast<double>(neighbors.size()) - radius
    };
    const bool use_running_mean{ inverse_radius == 0.0
        || std::abs(query.position[0]) > safe_query_magnitude
        || std::abs(query.position[1]) > safe_query_magnitude };
    std::sort(neighbors.begin(), neighbors.end(), [](const auto * a, const auto * b) {
        return a->index < b->index;
    });
    std::array<double, kOutputDimensionCount> mean{};
    double weight_sum{};
    for (const auto * neighbor : neighbors)
    {
        const double dx{ query.position[0] - neighbor->position[0] };
        const double dy{ query.position[1] - neighbor->position[1] };
        const double ux{ inverse_radius == 0.0 ? dx / radius : dx * inverse_radius };
        const double uy{ inverse_radius == 0.0 ? dy / radius : dy * inverse_radius };
        if (!std::isfinite(ux) || !std::isfinite(uy)) return MeanShiftStep::NonFiniteUpdate;
        const double weight{ std::exp(-exponent_factor * (ux * ux + uy * uy)) };
        weight_sum += weight;
        if (use_running_mean)
        {
            if (weight == 0.0) continue;
            const double fraction{ weight / weight_sum };
            for (std::size_t axis = 0; axis < kOutputDimensionCount; ++axis)
                mean[axis] = std::lerp(mean[axis], neighbor->position[axis], fraction);
        }
        else
        {
            for (std::size_t axis = 0; axis < kOutputDimensionCount; ++axis)
                mean[axis] += weight * neighbor->position[axis];
        }
    }
    if (weight_sum == 0.0) return MeanShiftStep::ZeroWeight;
    if (!use_running_mean)
        for (auto & coordinate : mean) coordinate /= weight_sum;
    if (!std::isfinite(weight_sum) || !HasFiniteCoordinates(mean)) return MeanShiftStep::NonFiniteUpdate;

    query.position = mean;
    return MeanShiftStep::Updated;
}

} // namespace

std::vector<double> RunUmap(
    const UmapFeatureMatrix & features,
    int effective_neighbors,
    double min_dist,
    int num_epochs,
    std::uint64_t random_seed,
    int job_count)
{
    const auto point_count{ static_cast<int>(features.values.size() / features.column_count) };
    std::vector<double> embedding(static_cast<std::size_t>(point_count) * kOutputDimensionCount);
    knncolle::VptreeBuilder<int, double, double> neighbor_builder(
        std::make_shared<knncolle::EuclideanDistance<double, double>>());
    umappp::Options options;
    options.num_neighbors = effective_neighbors;
    options.min_dist = min_dist;
    if (num_epochs > 0) options.num_epochs = num_epochs;
    const auto seed{ static_cast<umappp::RngEngine::result_type>(random_seed) };
    options.initialize_seed = seed;
    options.initialize_spectral_irlba_options.seed = seed;
    options.optimize_seed = seed;
    options.num_threads = job_count;
    options.num_threads_spectral = job_count;
    options.num_threads_optimize = 1;

    auto status{ umappp::initialize(features.column_count, point_count, features.values.data(), neighbor_builder,
        kOutputDimensionCount, embedding.data(), options) };
    status.run(embedding.data());
    return embedding;
}

UmapSpatialScale GetUmapSpatialScale(std::span<const double> embedding)
{
    if (embedding.size() % kOutputDimensionCount != 0)
        throw std::invalid_argument("UMAP spatial scale requires two coordinates per point.");
    if (!HasFiniteCoordinates(embedding))
        throw std::runtime_error("UMAP mean shift received a non-finite embedding coordinate.");
    if (embedding.empty()) return {};

    std::array lower{ embedding[0], embedding[1] };
    auto upper{ lower };
    for (std::size_t offset = kOutputDimensionCount; offset < embedding.size(); offset += kOutputDimensionCount)
    {
        for (std::size_t axis = 0; axis < kOutputDimensionCount; ++axis)
        {
            lower[axis] = std::min(lower[axis], embedding[offset + axis]);
            upper[axis] = std::max(upper[axis], embedding[offset + axis]);
        }
    }
    const double dx{ upper[0] - lower[0] };
    const double dy{ upper[1] - lower[1] };
    if (!std::isfinite(dx) || !std::isfinite(dy))
        throw std::runtime_error("UMAP mean shift bounding box span exceeds the finite double range.");

    const double diagonal{ std::hypot(dx, dy) };
    const double radius{ std::isfinite(diagonal) ? kRadiusRatio * diagonal
        : std::hypot(kRadiusRatio * dx, kRadiusRatio * dy) };
    // Keep the diagonal precise even when both spans are subnormal.
    const double scaled_diagonal{ radius < std::numeric_limits<double>::min()
        ? std::hypot(dx * kSubnormalCoordinateScale, dy * kSubnormalCoordinateScale) : 0.0 };
    return { lower, diagonal, radius, scaled_diagonal };
}

UmapMeanShiftResult RunUmapMeanShift(
    const std::vector<double> & embedding,
    const std::vector<LocalFittingFeatureRow> & feature_rows,
    [[maybe_unused]] int thread_size)
{
    constexpr int iteration_count{ 100 };
    const auto point_count{ embedding.size() / kOutputDimensionCount };
    if (embedding.size() % kOutputDimensionCount != 0 || point_count != feature_rows.size())
    {
        throw std::invalid_argument("UMAP mean shift requires two coordinates per feature row.");
    }
    const auto scale{ GetUmapSpatialScale(embedding) };
    if (embedding.empty()) return {};

    auto points{ MakeUmapPoints(embedding) };
    const double radius{ scale.radius };
    const auto root{ Tree::BuildKDTree(points) };
    auto shifted{ embedding };
    std::size_t empty_neighborhoods{};
    std::size_t zero_weights{};
    bool nonfinite_update{};
    std::exception_ptr worker_error;
#ifdef USE_OPENMP
    #pragma omp parallel num_threads(thread_size) reduction(+:empty_neighborhoods,zero_weights) reduction(||:nonfinite_update)
#endif
    {
        std::vector<UmapPoint *> neighbors;
#ifdef USE_OPENMP
        #pragma omp for schedule(static)
#endif
        for (std::size_t point = 0; point < point_count; ++point)
        {
            auto query{ points[point] };
            try
            {
                for (int iteration = 0; iteration < iteration_count; ++iteration)
                {
                    const auto step{ ShiftUmapPoint(root.get(), query, radius, neighbors) };
                    empty_neighborhoods += step == MeanShiftStep::EmptyNeighborhood;
                    zero_weights += step == MeanShiftStep::ZeroWeight;
                    if (step == MeanShiftStep::NonFiniteUpdate)
                    {
                        nonfinite_update = true;
                        break;
                    }
                }
            }
            catch (...)
            {
#ifdef USE_OPENMP
                #pragma omp critical(rhbm_umap_mean_shift_exception)
#endif
                {
                    if (!worker_error) worker_error = std::current_exception();
                }
            }
            const auto offset{ point * kOutputDimensionCount };
            shifted[offset] = query.position[0];
            shifted[offset + 1] = query.position[1];
        }
    }
    if (worker_error) std::rethrow_exception(worker_error);
    if (nonfinite_update)
    {
        throw std::runtime_error("Weighted mean shift produced a non-finite update.");
    }
    if (empty_neighborhoods != 0 || zero_weights != 0)
    {
        Logger::Log(LogLevel::Warning, "UMAP mean shift kept the current position for "
            + std::to_string(empty_neighborhoods) + " updates with no neighbors and "
            + std::to_string(zero_weights) + " updates with zero total neighbor weight.");
    }
    return { std::move(shifted), scale };
}

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
std::array<double, kOutputDimensionCount> RunUmapMeanShiftStepForTesting(
    std::span<const double> embedding,
    std::array<double, kOutputDimensionCount> position,
    double radius)
{
    auto points{ MakeUmapPoints(embedding) };
    const auto root{ Tree::BuildKDTree(points) };
    UmapPoint query{ position, 0 };
    std::vector<UmapPoint *> neighbors;
    ShiftUmapPoint(root.get(), query, radius, neighbors);
    return query.position;
}
#endif

UmapBinning BinMeanShiftEmbedding(
    const std::vector<LocalFittingFeatureRow> & feature_rows,
    const std::vector<double> & mean_shift_embedding,
    const UmapSpatialScale & scale)
{
    if (mean_shift_embedding.size() % kOutputDimensionCount != 0
        || mean_shift_embedding.size() / kOutputDimensionCount != feature_rows.size())
        throw std::invalid_argument("UMAP binning requires two coordinates per feature row.");
    if (!HasFiniteCoordinates(mean_shift_embedding))
        throw std::runtime_error("UMAP binning received a non-finite mean shift coordinate.");

    constexpr double min_bin_index{ static_cast<double>(std::numeric_limits<std::int64_t>::min()) };
    UmapBinning result;
    result.scale = scale;
    result.point_bin_numbers.resize(feature_rows.size());
    for (std::size_t point = 0; point < feature_rows.size(); ++point)
    {
        std::array<std::int64_t, kOutputDimensionCount> bin_index{};
        for (std::size_t axis = 0; scale.diagonal != 0.0 && axis < kOutputDimensionCount; ++axis)
        {
            const double offset{ mean_shift_embedding[point * kOutputDimensionCount + axis] - scale.origin[axis] };
            // Normalize by D first when R is subnormal or has rounded to zero.
            const double relative{ scale.scaled_diagonal > 0.0
                ? (offset * kSubnormalCoordinateScale / scale.scaled_diagonal) / kRadiusRatio
                : offset / scale.radius };
            const double index{ std::floor(relative / kUmapBinWidthRatio) };
            if (!(index >= min_bin_index && index < -min_bin_index))
            {
                throw std::runtime_error("Mean shift bin index is out of range for atom serial "
                    + std::to_string(feature_rows[point].serial_id) + ".");
            }
            bin_index[axis] = static_cast<std::int64_t>(index);
        }
        result.bins[bin_index].points.push_back(point);
    }
    std::size_t bin_number{};
    for (auto & [index, bin] : result.bins)
    {
        bin.number = ++bin_number;
        for (const auto point : bin.points) result.point_bin_numbers[point] = bin.number;
    }
    return result;
}

bool HasFiniteCoordinates(std::span<const double> coordinates)
{
    return std::all_of(coordinates.begin(), coordinates.end(),
        [](double coordinate) { return std::isfinite(coordinate); });
}

} // namespace rhbm_gem::core::classifier
