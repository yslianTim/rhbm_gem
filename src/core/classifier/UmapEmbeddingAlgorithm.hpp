#pragma once

#include "core/classifier/LocalFittingFeatures.hpp"
#include "core/classifier/UmapEmbeddingTypes.hpp"

#include <span>

namespace rhbm_gem::core::classifier {

std::vector<double> RunUmap(
    const UmapFeatureMatrix & features,
    int effective_neighbors,
    double min_dist,
    int num_epochs,
    std::uint64_t random_seed,
    int job_count);

std::vector<double> RunUmapMeanShift(
    const std::vector<double> & embedding,
    const std::vector<LocalFittingFeatureRow> & feature_rows,
    int thread_size);

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
std::array<double, kOutputDimensionCount> RunUmapMeanShiftStepForTesting(
    std::span<const double> embedding,
    std::array<double, kOutputDimensionCount> position,
    double radius);
#endif

UmapBinning BinMeanShiftEmbedding(
    const std::vector<LocalFittingFeatureRow> & feature_rows,
    const std::vector<double> & mean_shift_embedding);

bool HasFiniteCoordinates(std::span<const double> coordinates);

} // namespace rhbm_gem::core::classifier
