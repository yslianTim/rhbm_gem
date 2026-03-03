#pragma once

#include <cstddef>
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
class TGraphErrors;
#endif

namespace rhbm_gem {

class ModelObject;
class AtomClassifier;

class ComparisonPainter : public PainterBase
{
    std::string m_folder_path;
    std::vector<ModelObject *> m_model_object_list;
    std::vector<double> m_resolution_list;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ref_model_object_list_map;
    std::unique_ptr<AtomClassifier> m_atom_classifier;

public:
    ComparisonPainter(void);
    ~ComparisonPainter();
    void SetFolder(const std::string & folder_path) override;
    void AddDataObject(DataObjectBase * data_object) override;
    void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label) override;
    void Painting(void) override;

private:
    void PaintGroupGausEstimateComparison(const std::string & name);
    void PaintGausEstimateResidueClassDenseComparison(const std::string & name);
    void PainMapValueComparison(const std::string & name, ModelObject * model_object, const std::vector<ModelObject *> & ref_model_object_list);

#ifdef HAVE_ROOT
    void BuildGausRatioToResolutionGraph(int par_id, size_t target_id, size_t reference_id, ::TGraphErrors * graph, const std::vector<ModelObject *> & model_list, const std::string & class_key, Residue residue=Residue::UNK);
    void BuildAmplitudeRatioToWidthGraph(size_t target_id, size_t reference_id, ::TGraphErrors * graph, const std::vector<ModelObject *> & model_list, const std::string & class_key, bool draw_index=false, Residue residue=Residue::UNK);
    void BuildMapValueScatterGraph(GroupKey group_key, ::TGraphErrors * graph, ModelObject * model1, ModelObject * model2, int bin_size=15, double x_min=0.0, double x_max=1.5);
#endif
};

} // namespace rhbm_gem
