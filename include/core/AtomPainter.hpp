#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "DataObjectVisitor.hpp"
#include "PainterBase.hpp"

#ifdef HAVE_ROOT
class TPad;
class TH2;
class TPaveText;
#endif

namespace rhbm_gem {

class AtomObject;

class AtomPainter : public PainterBase, public DataObjectVisitor
{
    std::vector<AtomObject *> m_atom_object_list;
    std::unordered_map<std::string, AtomObject *> m_ref_atom_object_map;
    std::string m_folder_path;
    enum class IngestMode { Data, Reference };
    IngestMode m_ingest_mode{ IngestMode::Data };
    std::string m_ingest_label;

public:
    AtomPainter();
    ~AtomPainter();
    void SetFolder(const std::string & folder_path) override;
    void AddDataObject(DataObjectBase * data_object) override;
    void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label) override;
    void Painting() override;
    void VisitAtomObject(AtomObject & data_object) override;
    void VisitBondObject(BondObject & data_object) override;
    void VisitModelObject(ModelObject & data_object) override;
    void VisitMapObject(MapObject & data_object) override;

private:
    [[noreturn]] void ThrowInvalidType() const;
    void PaintDemoPlot(const std::string & name);
    void PaintAtomSamplingDataSummary(const std::string & name);
};

} // namespace rhbm_gem
