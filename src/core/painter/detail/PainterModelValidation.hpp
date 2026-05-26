#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem::painter_internal {

inline void RequireSelectedAtomsHaveLocalEntries(
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

inline void RequireLocalAnalyzedModel(
    const ModelObject & model_object,
    std::string_view painter_name)
{
    RequireSelectedAtomsHaveLocalEntries(model_object, painter_name);
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
