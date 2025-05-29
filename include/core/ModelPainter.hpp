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
class TH2;
class TPaveText;
#endif

class ModelPainter : public PainterBase
{
    std::vector<ModelObject *> m_model_object_list;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ref_model_object_map;
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
    void PaintGroupGausMainChain(ModelObject * model_object, const std::string & name, bool is_simulation=false);
    void PaintGroupChargeMainChain(ModelObject * model_object, const std::string & name, bool is_simulation=false);
    void PaintStructureClassGroupGausSideChain(ModelObject * model_object, const std::string & name);
    void PaintResidueClassGroupGausScatter(ModelObject * model_object, const std::string & name, int par_id=0);
    void PaintAtomMapValueMainChain(ModelObject * model_object, const std::string & name);
    void PaintResidueClassWidthScatterPlot(ModelObject * model_object, const std::string & name, int par_id=0, bool draw_box_plot=false);
    void PaintResidueClassXYPosition(ModelObject * model_object, const std::string & name);
    void PaintAtomXYPosition(ModelObject * model_object, const std::string & name);
    void PaintAtomGausScatter(ModelObject * model_object, const std::string & name, bool do_normalize=false);
    void PaintAtomGausMainChain(ModelObject * model_object, const std::string & name);
    void PaintAtomChargeMainChain(ModelObject * model_object, const std::string & name);
    void PaintAtomRankMainChain(ModelObject * model_object, const std::string & name);

    #ifdef HAVE_ROOT
    void PrintResolutionInfoPad(TPad * pad, TPaveText * text, double resolution);
    void PrintDataInfoPad(TPad * pad, TPaveText * text, const std::string & pdb_id, const std::string & emd_id);
    void PrintAmplitudePad(TPad * pad, TH2 * hist);
    void PrintWidthPad(TPad * pad, TH2 * hist);
    void PrintAmplitudeSummaryPad(TPad * pad, TH2 * hist);
    void PrintWidthSummaryPad(TPad * pad, TH2 * hist);
    void PrintGausSummaryPad(TPad * pad, TH2 * hist);
    void PrintCountSummaryPad(TPad * pad, TH2 * hist);

    void PrintInfoSideChainPad(TPad * pad, TPaveText * text, const std::string & residue_name);
    void PrintAmplitudeSideChainPad(TPad * pad, TH2 * hist, Residue residue, const std::vector<std::string> & label_list);
    void PrintWidthSideChainPad(TPad * pad, TH2 * hist, Residue residue, const std::vector<std::string> & label_list);
    void ModifyAxisLabelSideChain(TPad * pad, TH2 * hist, Residue residue, const std::vector<std::string> & label_list);
    #endif

};