#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/core/painter/PainterBase.hpp>

#ifdef HAVE_ROOT
class TPad;
class TH2;
class TPaveText;
#endif

namespace rhbm_gem {

class AtomObject;

class AtomPainter : public PainterBase
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

private:
    void IngestAtomObject(AtomObject & data_object);
    void PaintDemoPlot(const std::string & name);
    void PaintAtomSamplingDataSummary(const std::string & name);
};

} // namespace rhbm_gem
