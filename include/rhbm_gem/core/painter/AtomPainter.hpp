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

public:
    AtomPainter();
    ~AtomPainter();
    void AddDataObject(DataObjectBase * data_object) override;
    void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label) override;
    void Painting() override;

private:
    void AppendAtomObject(AtomObject & data_object);
    void PaintDemoPlot(const std::string & name);
    void PaintAtomSamplingDataSummary(const std::string & name);
};

} // namespace rhbm_gem
