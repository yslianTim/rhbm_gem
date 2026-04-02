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
class ModelObject;
namespace painter_internal {
class PainterModelIngress;
}

class AtomPainter : public PainterBase
{
    std::vector<AtomObject *> m_atom_object_list;
    std::string m_output_label;

public:
    AtomPainter();
    ~AtomPainter();
    // Compatibility surface: internal call sites should pre-validate via painter detail helpers.
    void AddModel(ModelObject & data_object);
    void Painting() override;

private:
    friend class painter_internal::PainterModelIngress;

    void AppendAtomObject(AtomObject & data_object);
    void PaintDemoPlot(const std::string & name);
    void PaintAtomSamplingDataSummary(const std::string & name);
};

} // namespace rhbm_gem
