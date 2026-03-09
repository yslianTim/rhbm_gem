#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/utils/GlobalEnumClass.hpp>
#include <rhbm_gem/core/PainterBase.hpp>

#ifdef HAVE_ROOT
class TPad;
class TH2;
class TAxis;
class TPaveText;
#endif

namespace rhbm_gem {

class ModelObject;
class AtomClassifier;
class BondClassifier;

class GausPainter : public PainterBase
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
    GausPainter();
    ~GausPainter();
    void SetFolder(const std::string & folder_path) override;
    void AddDataObject(DataObjectBase * data_object) override;
    void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label) override;
    void Painting() override;
#ifdef HAVE_ROOT
    static void RemodelAxisLabels(::TAxis * axis, const std::vector<std::string> & label_list, double angle, int align);
#endif

private:
    void IngestModelObject(ModelObject & data_object);
    void PaintAtomLocalGausSummary(ModelObject * model_object, const std::string & name);
    void PaintAtomGroupGausSummary(ModelObject * model_object, const std::string & name);
    void PaintAtomGroupMapValueAminoAcidMainChainComponent(ModelObject * model_object, const std::string & name);
    void PaintAtomGroupGausAminoAcidMainChainComponent(ModelObject * model_object, const std::string & name);
    void PaintAtomGroupGausAminoAcidMainChainComponentSimple(ModelObject * model_object, const std::string & name);
    void PaintAtomGroupGausAminoAcidMainChainStructure(ModelObject * model_object, const std::string & name);
    void PaintAtomLocalGausToSequenceAminoAcidMainChain(ModelObject * model_object, const std::string & name);

#ifdef HAVE_ROOT
    void RemodelFrameInPad(::TH2 * frame, ::TPad * pad, double x_tick_length, double y_tick_length);
    std::unique_ptr<::TPaveText> CreateDataInfoPaveText(ModelObject * model_object) const;
    std::unique_ptr<::TPaveText> CreateResolutionPaveText(ModelObject * model_object) const;
#endif
};

} // namespace rhbm_gem
