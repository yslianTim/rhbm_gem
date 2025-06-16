#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include "PainterBase.hpp"

class ModelObject;
class AtomClassifier;
enum class Residue : uint16_t;

#ifdef HAVE_ROOT
class TPad;
class TVirtualPad;
class TH2;
class TPaveText;
#endif

class DemoPainter : public PainterBase
{
    std::string m_folder_path;
    std::vector<ModelObject *> m_model_object_list;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ref_model_object_list_map;
    std::unique_ptr<AtomClassifier> m_atom_classifier;

public:
    DemoPainter(void);
    ~DemoPainter();
    void SetFolder(const std::string & folder_path) override;
    void AddDataObject(DataObjectBase * data_object) override;
    void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label) override;
    void Painting(void) override;

private:
    void PaintAtomMapValueExample(ModelObject * model_object, const std::string & name);
    void PaintGroupGausMainChainSummary(const std::vector<ModelObject *> & model_list, const std::string & name);
    void PaintGroupGausMainChainSingle(ModelObject * model_object, const std::string & name);
    void PaintAtomWidthScatterPlotSingle(ModelObject * model_object, const std::string & name, bool draw_box_plot=false);
    void PaintGroupGausToFSC(const std::string & name);
    void PaintGroupWidthScatterPlot(const std::vector<ModelObject *> & model_list, const std::string & name, int par_id=0, bool draw_box_plot=false);
    void PaintAtomGausMainChainDemo(ModelObject * model_object1, ModelObject * model_object2, const std::string & name, int par_id=0);
    void PaintAtomGausMainChainDemoSingle(ModelObject * model_object, const std::string & name, int par_id=0);

    #ifdef HAVE_ROOT
    void ModifyAxisLabelSideChain(TPad * pad, TH2 * hist, Residue residue, const std::vector<std::string> & label_list);
    void PrintIconMainChainPad(TPad * pad, TPaveText * text, double resolution, bool is_bottom_pad, bool is_top_pad);
    void PrintInfoMainChainPad(TPad * pad, TPaveText * text, const std::string & pdb_id, const std::string & emd_id, bool is_bottom_pad, bool is_top_pad);

    void PrintGausResultGlobalPad(TPad * pad, TH2 * hist, double left_margin, double right_margin, double bottom_margin, double top_margin, bool is_right_side_pad);

    void PrintGausTitlePad(TPad * pad, TPaveText * text, const std::string & title, float text_size);
    void PrintGausResultPad(TPad * pad, TH2 * hist, bool draw_x_axis, bool draw_title_label, bool is_right_side_pad);
    void PrintGausCorrelationPad(TPad * pad, TH2 * hist, bool draw_x_axis, bool draw_title_label);
    
    void PrintLeftSideChainPad(TPad * pad, TH2 * hist, Residue residue, const std::string & y_title, const std::vector<std::string> & label_list);
    void PrintRightSideChainPad(TPad * pad, TH2 * hist, Residue residue, const std::vector<std::string> & label_list);
    void PrintTitleSideChainPad(TPad * pad, TPaveText * text, const std::string & residue_name);

    void PrintGausResultInResidueIDPad(TPad * pad, TH2 * hist, int par_id=0);
    void PrintInfoInResidueIDPad(TVirtualPad * pad, TPaveText * text, const ModelObject * model_object, const std::string & chain_id, int residue_size);
    #endif

};