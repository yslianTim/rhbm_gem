#include "PotentialAnalysisCommand.hpp"
#include "DataObjectManager.hpp"
#include "PotentialAnalysisVisitor.hpp"
#include "AtomSelector.hpp"
#include "SphereSampler.hpp"

PotentialAnalysisCommand::PotentialAnalysisCommand(void) :
    m_atom_selector{ std::make_shared<AtomSelector>() },
    m_sphere_sampler{ std::make_shared<SphereSampler>() }
{

}

PotentialAnalysisCommand::~PotentialAnalysisCommand()
{
    
}

void PotentialAnalysisCommand::Execute(void)
{
    std::cout << "PotentialAnalysisCommand::Execute() called." << std::endl;

    auto data_manager{ std::make_unique<DataObjectManager>(m_database_path) };
    data_manager->ProcessFile(m_model_file_path, "model");
    data_manager->ProcessFile(m_map_file_path, "map");

    m_atom_selector->Print();
    m_sphere_sampler->Print();

    auto analyzer{ std::make_unique<PotentialAnalysisVisitor>(m_atom_selector, m_sphere_sampler) };
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

void PotentialAnalysisCommand::SetPickChainID(const std::string & value)
{
    m_atom_selector->PickChainID(value);
}

void PotentialAnalysisCommand::SetPickIndicator(const std::string & value)
{
    m_atom_selector->PickIndicator(value);
}

void PotentialAnalysisCommand::SetPickResidueType(const std::string & value)
{
    m_atom_selector->PickResidueType(value);
}

void PotentialAnalysisCommand::SetPickElementType(const std::string & value)
{
    m_atom_selector->PickElementType(value);
}

void PotentialAnalysisCommand::SetPickRemotenessType(const std::string & value)
{
    m_atom_selector->PickRemotenessType(value);
}

void PotentialAnalysisCommand::SetPickBranchType(const std::string & value)
{
    m_atom_selector->PickBranchType(value);
}

void PotentialAnalysisCommand::SetVetoChainID(const std::string & value)
{
    m_atom_selector->VetoChainID(value);
}

void PotentialAnalysisCommand::SetVetoIndicator(const std::string & value)
{
    m_atom_selector->VetoIndicator(value);
}

void PotentialAnalysisCommand::SetVetoResidueType(const std::string & value)
{
    m_atom_selector->VetoResidueType(value);
}

void PotentialAnalysisCommand::SetVetoElementType(const std::string & value)
{
    m_atom_selector->VetoElementType(value);
}

void PotentialAnalysisCommand::SetVetoRemotenessType(const std::string & value)
{
    m_atom_selector->VetoRemotenessType(value);
}

void PotentialAnalysisCommand::SetVetoBranchType(const std::string & value)
{
    m_atom_selector->VetoBranchType(value);
}