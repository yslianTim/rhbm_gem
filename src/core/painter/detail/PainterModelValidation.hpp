#pragma once

#include <stdexcept>
#include <cmath>
#include <limits>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <string>
#include <string_view>

#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem::painter_internal {

// Empty or constant subsets need a non-degenerate frame; no data points are added.
inline auto ComputePlotRange(const std::vector<double> & values, double margin, double minimum_range=0.1)
{
    auto range=array_helper::ComputeScalingRangeTuple(values,margin,minimum_range);
    auto & [low,high]=range;
    if (!std::isfinite(low) || !std::isfinite(high)) return decltype(range){0.0,1.0};
    if (!(high>low)) {const double pad=std::max(0.5,std::abs(low)*0.05); low-=pad; high+=pad;}
    return range;
}

// Candle distributions use available Original-grid responses in the histogram's
// [x_min, x_max) domain, not per-atom binned medians (or synthetic empty-bin zeros).
inline auto ComputePeeledResponsePlotRange(
    const ModelAnalysisView & view, const std::vector<GroupKey> & group_keys,
    double x_min, double x_max)
{
    double low{ std::numeric_limits<double>::infinity() };
    double high{ -std::numeric_limits<double>::infinity() };
    for (const auto key : group_keys)
    {
        for (const auto * atom : view.GetAtomObjectList(key))
        {
            const auto entry{ AtomLocalPotentialView::For(*atom) };
            if (!entry.IsAvailable()) continue;
            const auto & peeling{ entry.GetPostFitPeeling() };
            if (!peeling) continue;
            for (const auto & sample : peeling->samples)
            {
                if (!sample.response) continue;
                const auto distance{ std::sqrt(sample.observation.squared_distance) };
                if (distance < x_min || distance >= x_max) continue;
                low = std::min(low, *sample.response);
                high = std::max(high, *sample.response);
            }
        }
    }
    return ComputePlotRange(low <= high ? std::vector<double>{low, high} : std::vector<double>{}, 0.20);
}

inline void RequireLocalAnalyzedModel(
    const ModelObject & model_object,
    std::string_view painter_name)
{
    if (model_object.GetSelectedAtomCount() == 0)
    {
        throw std::runtime_error(
            std::string(painter_name) + " requires at least one selected atom.");
    }

    for (const auto * atom : model_object.GetSelectedAtoms())
    {
        if (!AtomLocalPotentialView::For(*atom).IsAvailable())
        {
            throw std::runtime_error(
                std::string(painter_name)
                + " requires analysis-ready selected atoms with local potential entries.");
        }
    }
}

inline void RequireGroupedAnalyzedModel(
    const ModelObject & model_object,
    std::string_view painter_name)
{
    RequireLocalAnalyzedModel(model_object, painter_name);
    if (model_object.GetAnalysisView().HasGroupedAnalysisData())
    {
        return;
    }

    throw std::runtime_error(
        std::string(painter_name)
        + " requires grouped analysis data. Run potential analysis before painting.");
}

} // namespace rhbm_gem::painter_internal
