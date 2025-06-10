#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include "PainterBase.hpp"


class ModelObject;
class AtomClassifier;

#ifdef HAVE_ROOT
class TPad;
class TH2;
class TPaveText;
class TGraphErrors;
#endif

class ComparisonPainter : public PainterBase
{
    std::vector<ModelObject *> m_model_object_list;
    std::unordered_map<std::string, std::vector<ModelObject *>> m_ref_model_object_list_map;
    std::string m_folder_path;
    std::unique_ptr<AtomClassifier> m_atom_classifier;

public:
    ComparisonPainter(void);
    ~ComparisonPainter();
    void SetFolder(const std::string & folder_path) override;
    void AddDataObject(DataObjectBase * data_object) override;
    void AddReferenceDataObject(DataObjectBase * data_object, const std::string & label) override;
    void Painting(void) override;

private:
    void PaintSimulationGaus(const std::string & name);
    void PaintSimulationGausRatio(const std::string & name, const std::vector<ModelObject *> & model_list);
    void PaintSimulationGausEstimate(const std::string & name);
    void PaintGausEstimateComparison(const std::string & name);
    void PaintGausEstimateComparisonOld(const std::string & name);
    void PaintGausEstimateScatterComparison(const std::string & name);
    void PainMapValueComparison(const std::string & name, ModelObject * model_data, ModelObject * model_sim);
    void PainResidueClassGausComparison(const std::string & name, ModelObject * model_data, ModelObject * model_sim, int par_id);
    double CalculateErrorPropagation(double target_value, double reference_value, double target_error, double reference_error);

    #ifdef HAVE_ROOT
    void BuildRatioGraph(TGraphErrors * ratio_graph, const TGraphErrors * target_graph, std::vector<TGraphErrors *> reference_graph_list, bool draw_index=false);
    void BuildRatioGraph(TGraphErrors * ratio_graph, const TGraphErrors * target_graph, const TGraphErrors * reference_graph, bool draw_index=false);
    void BuildGausEstimateToBlurringWidthGraph(uint64_t group_key, TGraphErrors * graph, const std::vector<ModelObject *> & model_list, int par_id=0);
    void BuildAmplitudeToWidthGraph(uint64_t group_key, TGraphErrors * graph, const std::vector<ModelObject *> & model_list);
    void BuildAmplitudeRatioToWidthScatterGraph(size_t target_id, size_t reference_id, TGraphErrors * graph, ModelObject * model);
    void BuildMapValueScatterGraph(uint64_t group_key, TGraphErrors * graph, ModelObject * model1, ModelObject * model2, int bin_size=15, double x_min=0.0, double x_max=1.5);
    void BuildGausScatterGraph(const std::vector<uint64_t> & group_key_list, TGraphErrors * graph, ModelObject * model1, ModelObject * model2, const std::string & class_key, int par_id=0);
    #endif

};