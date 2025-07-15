#include "PotentialAnalysisCommand.hpp"
#include "DataObjectManager.hpp"
#include "PotentialAnalysisVisitor.hpp"
#include "SphereSampler.hpp"
#include "ModelObject.hpp"
#include "Logger.hpp"

PotentialAnalysisCommand::PotentialAnalysisCommand(void) :
    m_thread_size{ 1 }, m_is_asymmetry{ false }, m_is_simulation{ false },
    m_fit_range_min{ 0.0 }, m_fit_range_max{ 1.0 }, m_alpha_r{ 0.1 }, m_alpha_g{ 0.2 },
    m_simulated_map_resolution{ 0.0 },
    m_database_path{ "database.sqlite" }, m_model_file_path{ "" }, m_map_file_path{ "" },
    m_saved_key_tag{ "model" },
    m_sphere_sampler{ std::make_unique<SphereSampler>() }
{
    m_sphere_sampler->SetThreadSize(static_cast<unsigned int>(m_thread_size));
    m_sphere_sampler->SetSamplingSize(1500);
    m_sphere_sampler->SetDistanceRangeMinimum(0.0);
    m_sphere_sampler->SetDistanceRangeMaximum(1.5);
}

PotentialAnalysisCommand::~PotentialAnalysisCommand()
{
    
}

void PotentialAnalysisCommand::Execute(void)
{
    Logger::Log(LogLevel::Info, "PotentialAnalysisCommand::Execute() called.");

    auto data_manager{ std::make_unique<DataObjectManager>(m_database_path) };
    try
    {
        data_manager->ProcessFile(m_model_file_path, "model");
        data_manager->ProcessFile(m_map_file_path, "map");
        if (m_is_simulation == true)
        {
            auto model_object{
                dynamic_cast<ModelObject *>(data_manager->GetDataObjectRef("model")) };
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
        Logger::Log(LogLevel::Error, e.what());
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
        Logger::Log(LogLevel::Warning,
            "[Warning] The resolution of input simulated map hasn't been set.\n"
            "          Please give the corresponding resolution value for this map.\n"
            "          (-r, --sim-resolution)");
    }
    model_object->SetEmdID("Simulation");
    model_object->SetResolution(m_simulated_map_resolution);
    model_object->SetResolutionMethod("Blurring Width");
}
