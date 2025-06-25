#include "PotentialAnalysisCommand.hpp"
#include "DataObjectManager.hpp"
#include "PotentialAnalysisVisitor.hpp"
#include "SphereSampler.hpp"
#include "ModelObject.hpp"

#include <iostream>

PotentialAnalysisCommand::PotentialAnalysisCommand(void) :
    m_simulated_map_resolution{ 0.0 },
    m_sphere_sampler{ std::make_unique<SphereSampler>() }
{

}

PotentialAnalysisCommand::~PotentialAnalysisCommand()
{
    
}

void PotentialAnalysisCommand::Execute(void)
{
    std::cout << "PotentialAnalysisCommand::Execute() called." << std::endl;

    auto data_manager{ std::make_unique<DataObjectManager>(m_database_path) };
    try
    {
        data_manager->ProcessFile(m_model_file_path, "model");
        data_manager->ProcessFile(m_map_file_path, "map");
        if (m_is_simulation == true)
        {
            auto model_object{
                dynamic_cast<ModelObject *>(data_manager->GetDataObjectRef("model").get()) };
            if (model_object == nullptr)
            {
                throw std::runtime_error(
                    "PotentialAnalysisCommand::Execute(): invalid model object");
            }
            UpdateModelObjectForSimulation(model_object);
        }
    }
    catch (const std::exception & e)
    {
        std::cerr << e.what() << '\n';
        return;
    }

    auto analyzer{ std::make_unique<PotentialAnalysisVisitor>(m_sphere_sampler.get()) };
    analyzer->SetAsymmetryFlag(m_is_asymmetry);
    analyzer->SetThreadSize(static_cast<unsigned int>(m_thread_size));
    analyzer->SetModelObjectKeyTag("model");
    analyzer->SetMapObjectKeyTag("map");
    analyzer->SetFitRange(m_fit_range_min, m_fit_range_max);
    analyzer->SetAlphaR(m_alpha_r);
    analyzer->SetAlphaG(m_alpha_g);

    data_manager->Accept(analyzer.get());
    data_manager->SaveDataObject("model", m_saved_key_tag);
}

void PotentialAnalysisCommand::SetThreadSize(int value)
{
    m_thread_size = value;
    m_sphere_sampler->SetThreadSize(static_cast<unsigned int>(value));
}

void PotentialAnalysisCommand::SetSamplingSize(int value)
{
    m_sphere_sampler->SetSamplingSize(static_cast<unsigned int>(value));
}

void PotentialAnalysisCommand::SetSamplingRangeMinimum(double value)
{
    m_sphere_sampler->SetDistanceRangeMinimum(value);
}

void PotentialAnalysisCommand::SetSamplingRangeMaximum(double value)
{
    m_sphere_sampler->SetDistanceRangeMaximum(value);
}

void PotentialAnalysisCommand::UpdateModelObjectForSimulation(ModelObject * model_object)
{
    if (model_object == nullptr) return;
    if (m_simulated_map_resolution == 0.0)
    {
        std::cout <<"[Warning] The resolution of input simulated map hasn't been set.\n";
        std::cout <<"          Please give the corresponding resolution value for this map.\n";
        std::cout <<"          (-r, --sim-resolution)"<< std::endl;
    }
    model_object->SetEmdID("Simulation");
    model_object->SetResolution(m_simulated_map_resolution);
    model_object->SetResolutionMethod("Blurring Width");
}
