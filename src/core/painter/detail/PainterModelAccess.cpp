#include "core/painter/detail/PainterModelAccess.hpp"

#include <rhbm_gem/core/painter/AtomPainter.hpp>
#include <rhbm_gem/core/painter/ComparisonPainter.hpp>
#include <rhbm_gem/core/painter/DemoPainter.hpp>
#include <rhbm_gem/core/painter/GausPainter.hpp>
#include <rhbm_gem/core/painter/ModelPainter.hpp>

namespace rhbm_gem::painter_internal {

void PainterModelIngress::AddModel(AtomPainter & painter, LocalAnalyzedModelRef model_ref)
{
    if (painter.m_output_label.empty())
    {
        painter.m_output_label = model_ref.model_object.GetKeyTag();
    }
    for (auto * atom : model_ref.model_object.GetSelectedAtoms())
    {
        painter.AppendAtomObject(*atom);
    }
}

void PainterModelIngress::AddModel(GausPainter & painter, GroupedAnalyzedModelRef model_ref)
{
    painter.m_model_object_list.push_back(&model_ref.model_object);
}

void PainterModelIngress::AddModel(ModelPainter & painter, GroupedAnalyzedModelRef model_ref)
{
    painter.m_model_object_list.push_back(&model_ref.model_object);
}

void PainterModelIngress::AddModel(ComparisonPainter & painter, GroupedAnalyzedModelRef model_ref)
{
    painter.m_model_object_list.emplace_back(&model_ref.model_object);
    painter.m_resolution_list.emplace_back(model_ref.model_object.GetResolution());
}

void PainterModelIngress::AddModel(DemoPainter & painter, GroupedAnalyzedModelRef model_ref)
{
    painter.m_model_object_list.push_back(&model_ref.model_object);
}

void PainterModelIngress::AddReferenceModel(
    ComparisonPainter & painter,
    GroupedAnalyzedModelRef model_ref,
    std::string_view label)
{
    painter.m_ref_model_object_list_map[std::string(label)].push_back(&model_ref.model_object);
}

void PainterModelIngress::AddReferenceModel(
    DemoPainter & painter,
    GroupedAnalyzedModelRef model_ref,
    std::string_view label)
{
    painter.m_ref_model_object_list_map[std::string(label)].push_back(&model_ref.model_object);
}

} // namespace rhbm_gem::painter_internal
