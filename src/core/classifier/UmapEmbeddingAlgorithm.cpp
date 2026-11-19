#include "core/classifier/UmapEmbeddingAlgorithm.hpp"

#include <rhbm_gem/utils/domain/Logger.hpp>

#include <knncolle/knncolle.hpp>
#include <umappp/umappp.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

namespace rhbm_gem::core::classifier {

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

std::vector<double> RunUmapMeanShift(
    const std::vector<double> & embedding,
    const std::vector<LocalFittingFeatureRow> & feature_rows,
    [[maybe_unused]] int thread_size)
{
    constexpr int iteration_count{ 100 };
    constexpr int requested_neighbors{ 200 };
    constexpr double alpha{ 2.0 };
    constexpr double sigma{ 1.5 };
    const int point_count{ static_cast<int>(embedding.size() / kOutputDimensionCount) };
    const int neighbors{ std::min(requested_neighbors, point_count) };
    if (neighbors != requested_neighbors)
    {
        Logger::Log(LogLevel::Warning,
            "UMAP mean shift neighbor count was limited to "
                + std::to_string(neighbors) + " for " + std::to_string(point_count)
                + " data rows.");
    }

    knncolle::VptreeBuilder<int, double, double> builder(
        std::make_shared<knncolle::EuclideanDistance<double, double>>());
    const auto index{ builder.build_unique(knncolle::SimpleMatrix<int, double>(
        kOutputDimensionCount, point_count, embedding.data())) };
    auto shifted{ embedding };
    std::atomic<bool> nonfinite_update{ false };
#ifdef USE_OPENMP
    #pragma omp parallel num_threads(thread_size)
#endif
    {
        auto searcher{ index->initialize() };
        std::vector<int> nearest;
        for (int iteration = 0; iteration < iteration_count; ++iteration)
        {
#ifdef USE_OPENMP
            #pragma omp for schedule(static)
#endif
            for (int point = 0; point < point_count; ++point)
            {
                const auto offset{ static_cast<std::size_t>(point) * kOutputDimensionCount };
                // Coordinate queries include coincident source points, including self.
                searcher->search(shifted.data() + offset, neighbors, &nearest, nullptr);
                std::array<double, kOutputDimensionCount> mean{};
                double weight_sum{};
                for (const int neighbor : nearest)
                {
                    const auto source{ static_cast<std::size_t>(neighbor) * kOutputDimensionCount };
                    const double dx{ shifted[offset] - embedding[source] };
                    const double dy{ shifted[offset + 1] - embedding[source + 1] };
                    const double distance_squared{ dx * dx + dy * dy };
                    const double weight{ std::exp(-alpha * distance_squared / (2.0 * sigma * sigma)) };
                    weight_sum += weight;
                    for (std::size_t dimension = 0; dimension < kOutputDimensionCount; ++dimension)
                    {
                        mean[dimension] += weight * embedding[source + dimension];
                    }
                }
                if (weight_sum == 0.0)
                {
                    Logger::Log(LogLevel::Warning,
                        "UMAP mean shift has zero total neighbor weight for atom serial "
                            + std::to_string(feature_rows[static_cast<std::size_t>(point)].serial_id)
                            + "; keeping its current position.");
                    continue;
                }
                for (auto & coordinate : mean) coordinate /= weight_sum;
                if (!std::isfinite(weight_sum)
                    || !HasFiniteCoordinates(mean))
                {
                    nonfinite_update.store(true, std::memory_order_relaxed);
                    continue;
                }
                for (std::size_t dimension = 0; dimension < kOutputDimensionCount; ++dimension)
                {
                    shifted[offset + dimension] = mean[dimension];
                }
            }
        }
    }
    if (nonfinite_update.load(std::memory_order_relaxed))
    {
        throw std::runtime_error("Weighted mean shift produced a non-finite update.");
    }
    return shifted;
}

UmapBinning BinMeanShiftEmbedding(
    const std::vector<LocalFittingFeatureRow> & feature_rows,
    const std::vector<double> & mean_shift_embedding)
{
    constexpr double bin_width{ 0.5 };
    constexpr double min_bin_index{ static_cast<double>(std::numeric_limits<std::int64_t>::min()) };
    UmapBinning result;
    result.point_bin_numbers.resize(feature_rows.size());
    for (std::size_t point = 0; point < feature_rows.size(); ++point)
    {
        std::array<std::int64_t, kOutputDimensionCount> bin_index{};
        for (std::size_t axis = 0; axis < kOutputDimensionCount; ++axis)
        {
            const double index{ std::floor(mean_shift_embedding[point * kOutputDimensionCount + axis] / bin_width) };
            if (index < min_bin_index || index >= -min_bin_index)
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
