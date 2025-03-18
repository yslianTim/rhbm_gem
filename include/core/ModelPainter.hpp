#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "PainterBase.hpp"

class ModelObject;

#ifdef HAVE_ROOT
class TPad;
class TH2;
class TPaveText;
#endif

class ModelPainter : public PainterBase
{
    ModelObject * m_model_object;
    std::string m_folder_path;

public:
    ModelPainter(void) = delete;
    explicit ModelPainter(ModelObject * model);
    ~ModelPainter();
    void SetFolder(const std::string & folder_path) override;
    void Painting(void) override;

private:
    void PaintResidueClassGroupGausMainChain(const std::string & name);
    void PaintResidueClassGroupGausSideChain(const std::string & name);

    #ifdef HAVE_ROOT
    void PrintIconPad(TPad * pad, TPaveText * text);
    void PrintInfoPad(TPad * pad, TPaveText * text, const std::string & pdb_id, const std::string & emd_id);
    void PrintAmplitudePad(TPad * pad, TH2 * hist);
    void PrintWidthPad(TPad * pad, TH2 * hist);
    void PrintAmplitudeSummaryPad(TPad * pad, TH2 * hist);
    void PrintWidthSummaryPad(TPad * pad, TH2 * hist);
    void PrintGausSummaryPad(TPad * pad, TH2 * hist);

    void PrintAmplitudeSideChainPad(TPad * pad, TH2 * hist, int residue, const std::vector<std::string> & label_list);
    void PrintWidthSideChainPad(TPad * pad, TH2 * hist, int residue, const std::vector<std::string> & label_list);
    void ModifyAxisLabelSideChain(TPad * pad, TH2 * hist, int residue, const std::vector<std::string> & label_list);

    #endif

};