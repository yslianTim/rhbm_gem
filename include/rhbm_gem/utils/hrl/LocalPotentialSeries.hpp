#pragma once

#include <tuple>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

namespace rhbm_gem::local_potential_series
{

std::tuple<double, double> ComputeDistanceRange(
    const LocalPotentialSampleList & sampling_entries,
    double margin_rate = 0.0);
std::tuple<double, double> ComputeResponseRange(
    const LocalPotentialSampleList & sampling_entries,
    double margin_rate = 0.0);
SeriesPointList BuildBinnedDistanceResponseSeries(
    const LocalPotentialSampleList & sampling_entries,
    int bin_size = 15,
    double x_min = 0.0,
    double x_max = 1.5);
double ComputeMapValueNearCenter(const LocalPotentialSampleList & sampling_entries);

} // namespace rhbm_gem::local_potential_series
