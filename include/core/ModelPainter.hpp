#pragma once

#include <iostream>
#include <string>
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
    void PaintResidueClassGroupGaus(const std::string & name);

    #ifdef HAVE_ROOT
    void PrintTitlePad(TPad * pad, TPaveText * text, const std::string & pdb_id, const std::string & emd_id);
    void PrintAmplitudePad(TPad * pad, TH2 * hist);
    void PrintWidthPad(TPad * pad, TH2 * hist);
    void PrintAmplitudeSummaryPad(TPad * pad, TH2 * hist);
    void PrintWidthSummaryPad(TPad * pad, TH2 * hist);
    void PrintGausSummaryPad(TPad * pad, TH2 * hist);
    #endif

};