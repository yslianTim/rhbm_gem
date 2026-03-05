#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "GlobalEnumClass.hpp"
#include "PainterBase.hpp"

#ifdef HAVE_ROOT
class TPad;
class TH2;
class TPaveText;
#endif

namespace rhbm_gem {

class ModelObject;
class AtomClassifier;
class BondClassifier;

class ModelPainter : public PainterBase
{
    std::string m_folder_path;
    std::vector<ModelObject *> m_model_object_list;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ref_model_object_list_map;
    std::unique_ptr<AtomClassifier> m_atom_classifier;
    std::unique_ptr<BondClassifier> m_bond_classifier;
    enum class IngestMode { Data, Reference };
    IngestMode m_ingest_mode{ IngestMode::Data };
    std::string m_ingest_label;

public:
    ModelPainter();
    ~ModelPainter();
    void SetFolder(const std::string & folder_path) override;
    void AddDataObject(DataObjectBase * data_object) override;
    void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label) override;
    void Painting() override;

private:
    void IngestModelObject(ModelObject & data_object);
    void PaintAtomGroupGausMainChain(ModelObject * model_object, const std::string & name);
    void PaintBondGroupGausMainChain(ModelObject * model_object, const std::string & name);
    void PaintAtomGroupGausNucleotideMainChain(ModelObject * model_object, const std::string & name);
    void PaintAtomMapValueMainChain(ModelObject * model_object, const std::string & name);
    void PaintBondMapValueMainChain(ModelObject * model_object, const std::string & name);
    void PaintGroupWidthScatterPlot(ModelObject * model_object, const std::string & name, int par_id=0, bool draw_box_plot=false);
    void PaintAtomXYPosition(ModelObject * model_object, const std::string & name);
    void PaintAtomGausScatterPlot(ModelObject * model_object, const std::string & name, bool do_normalize=false);
    void PaintBondGausScatterPlot(ModelObject * model_object, const std::string & name, bool do_normalize=false);
    void PaintAtomGausMainChain(ModelObject * model_object, const std::string & name);
    void PaintBondGausMainChain(ModelObject * model_object, const std::string & name);
    void PaintAtomRankMainChain(ModelObject * model_object, const std::string & name);

#ifdef HAVE_ROOT
    void PrintGausResultGlobalPad(::TPad * pad, ::TH2 * hist, double left_margin, double right_margin, double bottom_margin, double top_margin, bool is_right_side_pad);
    void PrintGausTitlePad(::TPad * pad, ::TPaveText * text, const std::string & title, float text_size);

    void PrintAmplitudePad(::TPad * pad, ::TH2 * hist);
    void PrintWidthPad(::TPad * pad, ::TH2 * hist);
    void PrintAmplitudeSummaryPad(::TPad * pad, ::TH2 * hist);
    void PrintAtomWidthSummaryPad(::TPad * pad, ::TH2 * hist);
    void PrintBondWidthSummaryPad(::TPad * pad, ::TH2 * hist);
    void PrintGausSummaryPad(::TPad * pad, ::TH2 * hist);

    void ModifyAxisLabelSideChain(::TPad * pad, ::TH2 * hist, Residue residue, const std::vector<std::string> & label_list);
#endif
};

} // namespace rhbm_gem
