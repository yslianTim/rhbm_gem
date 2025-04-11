#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include "PainterBase.hpp"

class ModelObject;
class AtomClassifier;

#ifdef HAVE_ROOT
class TPad;
class TH2;
class TPaveText;
#endif

class ModelPainter : public PainterBase
{
    std::vector<ModelObject *> m_model_object_list;
    std::unordered_map<std::string, ModelObject *> m_ref_model_object_map;
    std::string m_folder_path;
    std::unique_ptr<AtomClassifier> m_atom_classifier;

public:
    ModelPainter(void);
    ~ModelPainter();
    void SetFolder(const std::string & folder_path) override;
    void AddDataObject(DataObjectBase * data_object) override;
    void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label) override;
    void Painting(void) override;

private:
    void PaintResidueClassGroupGausMainChain(ModelObject * model_object, const std::string & name);
    void PaintResidueClassGroupGausSideChain(ModelObject * model_object, const std::string & name);
    void PaintResidueClassGroupGausMainChain(const std::string & name);
    //void PaintResidueClassGroupGausSideChain(const std::string & name);
    void PaintResidueClassMapValue(ModelObject * model_object, const std::string & name);
    void PaintResidueClassKNN(ModelObject * model_object, const std::string & name);
    void PaintResidueClassXYPosition(ModelObject * model_object, const std::string & name);

    #ifdef HAVE_ROOT
    void PrintIconPad(TPad * pad, TPaveText * text);
    void PrintInfoPad(TPad * pad, TPaveText * text, const std::string & pdb_id, const std::string & emd_id);
    void PrintAmplitudePad(TPad * pad, TH2 * hist);
    void PrintWidthPad(TPad * pad, TH2 * hist);
    void PrintAmplitudeSummaryPad(TPad * pad, TH2 * hist);
    void PrintWidthSummaryPad(TPad * pad, TH2 * hist);
    void PrintGausSummaryPad(TPad * pad, TH2 * hist);

    void PrintInfoSideChainPad(TPad * pad, TPaveText * text, const std::string & residue_name);
    void PrintAmplitudeSideChainPad(TPad * pad, TH2 * hist, int residue, const std::vector<std::string> & label_list);
    void PrintWidthSideChainPad(TPad * pad, TH2 * hist, int residue, const std::vector<std::string> & label_list);
    void ModifyAxisLabelSideChain(TPad * pad, TH2 * hist, int residue, const std::vector<std::string> & label_list);

    void PrintIconMainChainPad(TPad * pad, TPaveText * text, const std::string & fsc, bool is_bottom_pad);
    void PrintInfoMainChainPad(TPad * pad, TPaveText * text, const std::string & pdb_id, const std::string & emd_id, bool is_bottom_pad);

    void PrintTitlePad(TPad * pad, TPaveText * text, const std::string & title);
    void PrintResultPad(TPad * pad, TH2 * hist, bool draw_x_axis);
    void PrintSummaryPad(TPad * pad, TH2 * hist, bool draw_x_axis);
    void PrintCorrelationPad(TPad * pad, TH2 * hist, bool draw_x_axis);
    void PrintLeftSideChainPad(TPad * pad, TH2 * hist, int residue, const std::string & y_title, const std::vector<std::string> & label_list);
    void PrintRightSideChainPad(TPad * pad, TH2 * hist, int residue, const std::vector<std::string> & label_list);
    void PrintTitleSideChainPad(TPad * pad, TPaveText * text, const std::string & residue_name);
    #endif

};