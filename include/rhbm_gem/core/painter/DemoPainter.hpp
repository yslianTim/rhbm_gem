#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/core/painter/PainterBase.hpp>

#ifdef HAVE_ROOT
class TPad;
class TVirtualPad;
class TH2;
class TPaveText;
class TGraphErrors;
#endif

namespace rhbm_gem {

class ModelObject;
namespace painter_internal {
class PainterModelIngress;
}

class DemoPainter : public PainterBase
{
    std::vector<ModelObject *> m_model_object_list;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ref_model_object_list_map;

public:
    DemoPainter();
    ~DemoPainter();
    // Compatibility surface: internal call sites should pre-validate via painter detail helpers.
    void AddModel(ModelObject & data_object);
    // Compatibility surface: internal call sites should pre-validate via painter detail helpers.
    void AddReferenceModel(ModelObject & data_object, std::string_view label);
    void Painting() override;

private:
    friend class painter_internal::PainterModelIngress;

    void PainMapValueComparisonSingle(const std::string & name, ModelObject * model_object, ModelObject * ref_model_object);
    void PaintAtomMapValueExample(ModelObject * model_object, const std::string & name);
    void PaintGroupGausMainChainSummary(const std::vector<ModelObject *> & model_list, const std::string & name);
    void PaintGroupGausMainChainSingle(ModelObject * model_object, const std::string & name);
    void PaintAtomWidthScatterPlotSingle(ModelObject * model_object, const std::string & name, bool draw_box_plot=false);
    void PaintGroupGausToFSC(const std::vector<ModelObject *> & model_list, const std::string & name);
    void PaintGroupWidthScatterPlot(const std::vector<ModelObject *> & model_list, const std::string & name, int par_id=0, bool draw_box_plot=false);
    void PaintAtomGausMainChainDemo(ModelObject * model_object1, ModelObject * model_object2, const std::string & name, int par_id=0);
    void PaintAtomGausMainChainDemoSingle(ModelObject * model_object, const std::string & name, int par_id=0);
    void PaintGroupWidthAlphaCarbonDemo(const std::vector<ModelObject *> & model_list, const std::string & name);
    void PaintGroupGausMergeResidueDemo(const std::vector<ModelObject *> & model_list, const std::string & name);

#ifdef HAVE_ROOT
    void PrintGausResultPad(::TPad * pad, ::TH2 * hist, bool draw_x_axis, bool draw_title_label, bool is_right_side_pad);
    void PrintGausCorrelationPad(::TPad * pad, ::TH2 * hist, bool draw_x_axis, bool draw_title_label);
#endif
};

} // namespace rhbm_gem
