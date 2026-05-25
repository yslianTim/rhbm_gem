#pragma once

#include <tuple>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

namespace rhbm_gem::local_potential_series
{

enum class QScoreReference
{
    Fixed = 0,
    MDPDE = 1
};

std::tuple<float, float> ComputeDistanceRange(
    const LocalPotentialSampleList & sampling_entries,
    double margin_rate = 0.0);
std::tuple<float, float> ComputeResponseRange(
    const LocalPotentialSampleList & sampling_entries,
    double margin_rate = 0.0);
SeriesPointList BuildBinnedDistanceResponseSeries(
    const LocalPotentialSampleList & sampling_entries,
    int bin_size = 15,
    double x_min = 0.0,
    double x_max = 1.5);
double ComputeMapValueNearCenter(const LocalPotentialSampleList & sampling_entries);
double ComputeQScore(
    const LocalPotentialSampleList & sampling_entries,
    const LocalGaussianResult & gaussian_result,
    QScoreReference reference);

} // namespace rhbm_gem::local_potential_series
