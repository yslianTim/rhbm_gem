#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/core/painter/PainterBase.hpp>

#ifdef HAVE_ROOT
class TPad;
class TH2;
class TPaveText;
class TGraphErrors;
#endif

namespace rhbm_gem {

class ModelObject;
namespace painter_internal {
class PainterModelIngress;
}

class ComparisonPainter : public PainterBase
{
    std::vector<ModelObject *> m_model_object_list;
    std::vector<double> m_resolution_list;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ref_model_object_list_map;

public:
    ComparisonPainter();
    ~ComparisonPainter();
    // Compatibility surface: internal call sites should pre-validate via painter detail helpers.
    void AddModel(ModelObject & data_object);
    // Compatibility surface: internal call sites should pre-validate via painter detail helpers.
    void AddReferenceModel(ModelObject & data_object, std::string_view label);
    void Painting() override;

private:
    friend class painter_internal::PainterModelIngress;

    void PaintGroupGausEstimateComparison(const std::string & name);
    void PaintGausEstimateResidueClassDenseComparison(const std::string & name);
    void PainMapValueComparison(const std::string & name, ModelObject * model_object, const std::vector<ModelObject *> & ref_model_object_list);

#ifdef HAVE_ROOT
    void BuildGausRatioToResolutionGraph(int par_id, size_t target_id, size_t reference_id, ::TGraphErrors * graph, const std::vector<ModelObject *> & model_list, const std::string & class_key, Residue residue=Residue::UNK);
    void BuildAmplitudeRatioToWidthGraph(size_t target_id, size_t reference_id, ::TGraphErrors * graph, const std::vector<ModelObject *> & model_list, const std::string & class_key, bool draw_index=false, Residue residue=Residue::UNK);
#endif
};

} // namespace rhbm_gem
