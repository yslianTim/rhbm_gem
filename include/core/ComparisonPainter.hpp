#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include "PainterBase.hpp"
#include "AtomicInfoHelper.hpp"

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
    using ElementKeyType = GroupKeyMapping<ElementGroupClassifierTag>::type;
    using ResidueKeyType = GroupKeyMapping<ResidueGroupClassifierTag>::type;

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
    void PaintDataSimulationComparison(const std::string & name);
    double CalculateErrorPropagation(double target_value, double reference_value, double target_error, double reference_error);

    #ifdef HAVE_ROOT
    void BuildRatioGraph(TGraphErrors * ratio_graph, const TGraphErrors * target_graph, const TGraphErrors * reference_graph);
    void BuildGausEstimateToBlurringWidthGraph(ElementKeyType & group_key, TGraphErrors * graph, const std::vector<ModelObject *> & model_list, int par_id=0);
    void BuildAmplitudeRatioToWidthGraph(ElementKeyType & group_key, TGraphErrors * graph, const std::vector<ModelObject *> & model_list);
    #endif

};