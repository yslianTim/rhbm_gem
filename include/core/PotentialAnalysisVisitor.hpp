#pragma once

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <set>
#include <tuple>
#include <any>
#include <unordered_map>
#include "DataObjectVisitorBase.hpp"
#include "AtomicInfoHelper.hpp"

class AtomSelector;
class SphereSampler;

class PotentialAnalysisVisitor : public DataObjectVisitorBase
{
    double m_alpha_r, m_alpha_g;
    double m_x_min, m_x_max;
    std::string m_map_key_tag, m_model_key_tag;
    std::shared_ptr<AtomSelector> m_atom_selector;
    std::shared_ptr<SphereSampler> m_sphere_sampler;
    std::vector<AtomObject *> m_selected_atom_list;

    using ElementKeyType = GroupKeyMapping<ElementGroupClassifierTag>::type;
    using ResidueKeyType = GroupKeyMapping<ResidueGroupClassifierTag>::type;
    std::set<ElementKeyType> element_class_group_set;
    std::set<ResidueKeyType> residue_class_group_set;

public:
    PotentialAnalysisVisitor(std::shared_ptr<AtomSelector> atom_selector, std::shared_ptr<SphereSampler> sphere_sampler);
    ~PotentialAnalysisVisitor();
    void VisitAtomObject(AtomObject * data_object) override;
    void VisitModelObject(ModelObject * data_object) override;
    void VisitMapObject(MapObject * data_object) override;
    void Analysis(DataObjectManager * data_manager) override;

    void SetFitRange(double x_min, double x_max);
    void SetAlphaR(double alpha_r) { m_alpha_r = alpha_r; }
    void SetAlphaG(double alpha_g) { m_alpha_g = alpha_g; }
    void SetMapObjectKeyTag(const std::string & value) { m_map_key_tag = value; }
    void SetModelObjectKeyTag(const std::string & value) { m_model_key_tag = value; }

private:
    void RunAtomClassification(const std::string & class_key, ModelObject * model_object);
    void RunPotentialFitting(const std::string & class_key, ModelObject * model_object);

};