#pragma once

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include "DataObjectVisitorBase.hpp"
#include "PotentialAnalysisCommand.hpp"

class SphereSampler;

class PotentialAnalysisVisitor : public DataObjectVisitorBase
{
    bool m_is_asymmetry;
    unsigned int m_thread_size;
    double m_alpha_r, m_alpha_g;
    double m_x_min, m_x_max;
    std::string m_map_key_tag, m_model_key_tag;
    SphereSampler * m_sphere_sampler;
    std::vector<AtomObject *> m_selected_atom_list;

public:
    PotentialAnalysisVisitor(SphereSampler * sphere_sampler, const PotentialAnalysisCommand::Options & options);
    ~PotentialAnalysisVisitor() = default;
    void VisitAtomObject(AtomObject * data_object) override;
    void VisitModelObject(ModelObject * data_object) override;
    void VisitMapObject(MapObject * data_object) override;
    void VisitDataObjectManager(DataObjectManager * data_manager) override;

private:
    void RunPotentialFitting(ModelObject * model_object);

};
