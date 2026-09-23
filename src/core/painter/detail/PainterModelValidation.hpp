#pragma once

#include <stdexcept>
#include <cmath>
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
