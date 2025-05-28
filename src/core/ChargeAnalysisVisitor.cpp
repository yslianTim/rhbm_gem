#include "ChargeAnalysisVisitor.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "DataObjectManager.hpp"
#include "SphereSampler.hpp"
#include "HRLModelHelper.hpp"
#include "ScopeTimer.hpp"
#include "AtomicInfoHelper.hpp"
#include "KeyPacker.hpp"
#include "ArrayStats.hpp"
#include "AtomClassifier.hpp"

#include <iostream>
#include <tuple>
#include <fstream>
#include <unordered_set>

ChargeAnalysisVisitor::ChargeAnalysisVisitor(
    SphereSampler * sphere_sampler) :
    m_thread_size{ 1 },
    m_alpha_r{ 0.0 }, m_alpha_g{ 0.0 }, m_x_min{ 0.0 }, m_x_max{ 0.0 },
    m_sphere_sampler{ sphere_sampler }
{

}

ChargeAnalysisVisitor::~ChargeAnalysisVisitor()
{

}

void ChargeAnalysisVisitor::VisitAtomObject(AtomObject * data_object)
{
    if (data_object == nullptr) return;
    data_object->SetSelectedFlag(true);
}

void ChargeAnalysisVisitor::VisitModelObject(ModelObject * data_object)
{
    ScopeTimer timer("ChargeAnalysisVisitor::VisitModelObject");
    if (data_object == nullptr) return;
}

void ChargeAnalysisVisitor::VisitMapObject(MapObject * data_object)
{
    ScopeTimer timer("ChargeAnalysisVisitor::VisitMapObject");
    if (data_object == nullptr) return;
    
}

void ChargeAnalysisVisitor::Analysis(DataObjectManager * data_manager)
{
    std::cout <<"- Analysis..." << std::endl;
    try
    {
        const auto & model_object{ data_manager->GetDataObjectRef(m_model_key_tag) };
        model_object->Accept(this);

        const auto & map_object{ data_manager->GetDataObjectRef(m_map_key_tag) };
        map_object->Accept(this);


    }
    catch(const std::exception & e)
    {
        std::cerr << e.what() << std::endl;
    }
}

void ChargeAnalysisVisitor::SetAsymmetryFlag(bool value)
{
    m_is_asymmetry = value;
}

void ChargeAnalysisVisitor::SetThreadSize(unsigned int thread_size)
{
    m_thread_size = thread_size;
}

void ChargeAnalysisVisitor::SetFitRange(double x_min, double x_max)
{
    m_x_min = x_min;
    m_x_max = x_max;
}

void ChargeAnalysisVisitor::SetAlphaR(double alpha_r)
{
    m_alpha_r = alpha_r;
}

void ChargeAnalysisVisitor::SetAlphaG(double alpha_g)
{
    m_alpha_g = alpha_g;
}

void ChargeAnalysisVisitor::SetMapObjectKeyTag(const std::string & value)
{
    m_map_key_tag = value;
}

void ChargeAnalysisVisitor::SetModelObjectKeyTag(const std::string & value)
{
    m_model_key_tag = value;
}