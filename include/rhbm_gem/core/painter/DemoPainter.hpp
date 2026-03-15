#pragma once

#include <memory>
#include <string>
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
class AtomClassifier;

class DemoPainter : public PainterBase
{
    std::string m_folder_path;
    std::vector<ModelObject *> m_model_object_list;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ref_model_object_list_map;
    std::unique_ptr<AtomClassifier> m_atom_classifier;

public:
    DemoPainter();
    ~DemoPainter();
    void SetFolder(const std::string & folder_path) override;
    void AddDataObject(DataObjectBase * data_object) override;
    void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label) override;
    void Painting() override;

private:
    void AppendModelObject(ModelObject & data_object);
    void AppendReferenceModelObject(ModelObject & data_object, const std::string & label);
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
    void PrintGausResultGlobalPad(::TPad * pad, ::TH2 * hist, double left_margin, double right_margin, double bottom_margin, double top_margin, bool is_right_side_pad);
    void PrintGausTitlePad(::TPad * pad, ::TPaveText * text, const std::string & title, float text_size);
    void PrintGausResultPad(::TPad * pad, ::TH2 * hist, bool draw_x_axis, bool draw_title_label, bool is_right_side_pad);
    void PrintGausCorrelationPad(::TPad * pad, ::TH2 * hist, bool draw_x_axis, bool draw_title_label);
    void BuildMapValueScatterGraph(GroupKey group_key, ::TGraphErrors * graph, ModelObject * model1, ModelObject * model2, int bin_size, double x_min, double x_max);
#endif
};

} // namespace rhbm_gem
