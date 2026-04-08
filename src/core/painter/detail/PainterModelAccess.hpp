#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem {

class AtomPainter;
class ComparisonPainter;
class DemoPainter;
class GausPainter;
class ModelPainter;

namespace painter_internal {

struct LocalAnalyzedModelRef
{
    ModelObject & model_object;
};

struct GroupedAnalyzedModelRef
{
    ModelObject & model_object;
};

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
        if (!LocalPotentialView::For(*atom).IsAvailable())
        {
            throw std::runtime_error(
                std::string(painter_name)
                + " requires analysis-ready selected atoms with local potential entries.");
        }
    }
}

inline void RequireSelectedBondsHaveLocalEntries(
    const ModelObject & model_object,
    std::string_view painter_name)
{
    for (const auto * bond : model_object.GetSelectedBonds())
    {
        if (!LocalPotentialView::For(*bond).IsAvailable())
        {
            throw std::runtime_error(
                std::string(painter_name)
                + " requires analysis-ready selected bonds with local potential entries.");
        }
    }
}

inline void RequireGroupedAnalysisData(
    const ModelObject & model_object,
    std::string_view painter_name)
{
    if (model_object.GetAnalysisView().HasGroupedAnalysisData())
    {
        return;
    }

    throw std::runtime_error(
        std::string(painter_name)
        + " requires grouped analysis data. Run potential analysis before painting.");
}

inline LocalAnalyzedModelRef RequireLocalAnalyzedModel(
    ModelObject & model_object,
    std::string_view painter_name)
{
    RequireSelectedAtomsHaveLocalEntries(model_object, painter_name);
    RequireSelectedBondsHaveLocalEntries(model_object, painter_name);
    return LocalAnalyzedModelRef{ model_object };
}

inline GroupedAnalyzedModelRef RequireGroupedAnalyzedModel(
    ModelObject & model_object,
    std::string_view painter_name)
{
    RequireLocalAnalyzedModel(model_object, painter_name);
    RequireGroupedAnalysisData(model_object, painter_name);
    return GroupedAnalyzedModelRef{ model_object };
}

class PainterModelIngress
{
public:
    static void AddModel(AtomPainter & painter, LocalAnalyzedModelRef model_ref);
    static void AddModel(GausPainter & painter, GroupedAnalyzedModelRef model_ref);
    static void AddModel(ModelPainter & painter, GroupedAnalyzedModelRef model_ref);
    static void AddModel(ComparisonPainter & painter, GroupedAnalyzedModelRef model_ref);
    static void AddModel(DemoPainter & painter, GroupedAnalyzedModelRef model_ref);
    static void AddReferenceModel(
        ComparisonPainter & painter,
        GroupedAnalyzedModelRef model_ref,
        std::string_view label);
    static void AddReferenceModel(
        DemoPainter & painter,
        GroupedAnalyzedModelRef model_ref,
        std::string_view label);
};

} // namespace painter_internal
} // namespace rhbm_gem
