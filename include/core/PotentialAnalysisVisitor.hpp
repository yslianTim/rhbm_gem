#pragma once

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <tuple>
#include <unordered_map>
#include "DataObjectVisitorBase.hpp"
#include "TupleHashHelper.hpp"

class AtomSelector;
class SphereSampler;
class DataObjectBase;
struct AtomicPotentialEntry;
struct GroupPotentialEntry;

class PotentialAnalysisVisitor : public DataObjectVisitorBase
{
    double m_alpha_r, m_alpha_g;
    double m_x_min, m_x_max;
    std::string m_map_key_tag, m_model_key_tag;
    std::shared_ptr<AtomSelector> m_atom_selector;
    std::shared_ptr<SphereSampler> m_sphere_sampler;
    std::vector<AtomObject *> m_selected_atom_list;
    std::unordered_map<std::tuple<int, int, int, int>, std::unique_ptr<GroupPotentialEntry>, TupleHash<int, int, int, int>> m_grouping_gaus_entry_map;
    std::unordered_map<std::tuple<int, int, int, int>, std::vector<AtomObject *>, TupleHash<int, int, int, int>> m_grouping_atom_map;

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
    void RunAtomClassification(void);
    void RunPotentialFitting(void);

};