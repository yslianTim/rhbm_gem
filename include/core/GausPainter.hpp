#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>

#include "GlobalEnumClass.hpp"
#include "PainterBase.hpp"

class ModelObject;
class AtomClassifier;
class BondClassifier;

#ifdef HAVE_ROOT
class TPad;
class TH2;
class TPaveText;
#endif

class GausPainter : public PainterBase
{
    std::string m_folder_path;
    std::vector<ModelObject *> m_model_object_list;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ref_model_object_list_map;
    std::unique_ptr<AtomClassifier> m_atom_classifier;
    std::unique_ptr<BondClassifier> m_bond_classifier;

public:
    GausPainter(void);
    ~GausPainter();
    void SetFolder(const std::string & folder_path) override;
    void AddDataObject(DataObjectBase * data_object) override;
    void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label) override;
    void Painting(void) override;

private:
    void PaintAtomGroupGausSummary(ModelObject * model_object, const std::string & name);

    #ifdef HAVE_ROOT
    void RemodelPad(TPad * pad, TH2 * frame, const std::vector<std::string> & label_list);
    void RemodelElementPad(TPad * pad, TH2 * frame, const std::vector<std::string> & label_list);
    #endif

};
