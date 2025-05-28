#include "ChargeAnalysisCommand.hpp"
#include "DataObjectManager.hpp"
#include "ChargeAnalysisVisitor.hpp"
#include "SphereSampler.hpp"

#include <iostream>

ChargeAnalysisCommand::ChargeAnalysisCommand(void) :
    m_sphere_sampler{ std::make_unique<SphereSampler>() }
{

}

ChargeAnalysisCommand::~ChargeAnalysisCommand()
{
    
}

void ChargeAnalysisCommand::Execute(void)
{
    std::cout << "ChargeAnalysisCommand::Execute() called." << std::endl;

    auto data_manager{ std::make_unique<DataObjectManager>(m_database_path) };
    data_manager->ProcessFile(m_model_file_path, "model");
    data_manager->ProcessFile(m_map_file_path, "map");
    data_manager->ProcessFile(m_sim_neutral_map_file_path, "map_neutral");
    data_manager->ProcessFile(m_sim_positive_map_file_path, "map_positive");
    data_manager->ProcessFile(m_sim_negative_map_file_path, "map_negative");

    m_sphere_sampler->Print();

    auto analyzer{ std::make_unique<ChargeAnalysisVisitor>(m_sphere_sampler.get()) };
    analyzer->SetAsymmetryFlag(m_is_asymmetry);
    analyzer->SetThreadSize(static_cast<unsigned int>(m_thread_size));
    analyzer->SetModelObjectKeyTag("model");
    analyzer->SetMapObjectKeyTag("map");
    analyzer->SetNeutralMapObjectKeyTag("map_neutral");
    analyzer->SetPositiveMapObjectKeyTag("map_positive");
    analyzer->SetNegativeMapObjectKeyTag("map_negative");
    analyzer->SetFitRange(m_fit_range_min, m_fit_range_max);
    analyzer->SetAlphaR(m_alpha_r);
    analyzer->SetAlphaG(m_alpha_g);

    data_manager->Accept(analyzer.get());

    data_manager->SaveDataObject("model", m_saved_key_tag);
}

void ChargeAnalysisCommand::SetThreadSize(int value)
{
    m_thread_size = value;
    m_sphere_sampler->SetThreadSize(static_cast<unsigned int>(value));
}

void ChargeAnalysisCommand::SetSamplingSize(int value)
{
    m_sphere_sampler->SetSamplingSize(static_cast<unsigned int>(value));
}

void ChargeAnalysisCommand::SetSamplingRangeMinimum(double value)
{
    m_sphere_sampler->SetDistanceRangeMinimum(value);
}

void ChargeAnalysisCommand::SetSamplingRangeMaximum(double value)
{
    m_sphere_sampler->SetDistanceRangeMaximum(value);
}