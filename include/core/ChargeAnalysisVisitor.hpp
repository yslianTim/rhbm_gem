#pragma once

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include "DataObjectVisitorBase.hpp"

class AtomObject;

class ChargeAnalysisVisitor : public DataObjectVisitorBase
{
    unsigned int m_thread_size;
    double m_alpha_r, m_alpha_g;
    double m_x_min, m_x_max;
    std::string m_model_key_tag;
    std::string m_neutral_model_key_tag, m_positive_model_key_tag, m_negative_model_key_tag;
    std::vector<AtomObject *> m_selected_atom_list;

public:
    ChargeAnalysisVisitor(void);
    ~ChargeAnalysisVisitor();
    void VisitAtomObject(AtomObject * data_object) override;
    void VisitModelObject(ModelObject * data_object) override;
    void VisitMapObject(MapObject * data_object) override;
    void Analysis(DataObjectManager * data_manager) override;

    void SetThreadSize(unsigned int thread_size) { m_thread_size = thread_size; }
    void SetAlphaR(double alpha_r) { m_alpha_r = alpha_r; }
    void SetAlphaG(double alpha_g) { m_alpha_g = alpha_g; }
    void SetModelObjectKeyTag(const std::string & value) { m_model_key_tag = value; }
    void SetNeutralModelObjectKeyTag(const std::string & value) { m_neutral_model_key_tag = value; }
    void SetPositiveModelObjectKeyTag(const std::string & value) { m_positive_model_key_tag = value; }
    void SetNegativeModelObjectKeyTag(const std::string & value) { m_negative_model_key_tag = value; }
    void SetFitRange(double x_min, double x_max);

private:
    std::unordered_map<int, AtomObject *> BuildSerialIDAtomObjectMap(ModelObject * model_object);
    void RunChargeFitting(ModelObject * model_neutral, ModelObject * model_pos, ModelObject * model_neg);
    void PrintRegressionResult(const std::unordered_map<size_t, std::vector<std::tuple<float, float, float>>> & data);

};
