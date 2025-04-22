#pragma once

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <set>
#include <tuple>
#include <unordered_map>
#include <functional>
#include "DataObjectVisitorBase.hpp"

class AtomSelector;
class SphereSampler;

class PotentialAnalysisVisitor : public DataObjectVisitorBase
{
    unsigned int m_thread_size;
    double m_alpha_r, m_alpha_g;
    double m_x_min, m_x_max;
    std::string m_map_key_tag, m_model_key_tag;
    std::shared_ptr<AtomSelector> m_atom_selector;
    std::shared_ptr<SphereSampler> m_sphere_sampler;
    std::vector<AtomObject *> m_selected_atom_list;

    std::set<uint64_t> element_class_group_set;
    std::set<uint64_t> residue_class_group_set;

public:
    PotentialAnalysisVisitor(std::shared_ptr<AtomSelector> atom_selector, std::shared_ptr<SphereSampler> sphere_sampler);
    ~PotentialAnalysisVisitor();
    void VisitAtomObject(AtomObject * data_object) override;
    void VisitModelObject(ModelObject * data_object) override;
    void VisitMapObject(MapObject * data_object) override;
    void Analysis(DataObjectManager * data_manager) override;

    void SetThreadSize(unsigned int thread_size) { m_thread_size = thread_size; }
    void SetFitRange(double x_min, double x_max);
    void SetAlphaR(double alpha_r) { m_alpha_r = alpha_r; }
    void SetAlphaG(double alpha_g) { m_alpha_g = alpha_g; }
    void SetMapObjectKeyTag(const std::string & value) { m_map_key_tag = value; }
    void SetModelObjectKeyTag(const std::string & value) { m_model_key_tag = value; }

private:
    const std::set<uint64_t> & GetGroupSet(const std::string & class_key);
    void RunAtomClassification(const std::string & class_key, ModelObject * model_object);
    void RunPotentialFitting(const std::string & class_key, ModelObject * model_object);

};