#include "PotentialAnalysisCommand.hpp"
#include "DataObjectManager.hpp"
#include "PotentialAnalysisVisitor.hpp"
#include "AtomSelector.hpp"
#include "SphereSampler.hpp"

PotentialAnalysisCommand::PotentialAnalysisCommand(void) :
    m_thread_size{ 1 }, m_sampling_size{ 1500 },
    m_sampling_range_min{ 0.0 }, m_sampling_range_max{ 1.5 },
    m_fit_range_min{ 0.0 }, m_fit_range_max{ 1.0 },
    m_alpha_r{ 0.1 }, m_alpha_g{ 0.2 },
    m_database_path{ "" }, m_model_file_path{ "" }, m_map_file_path{ "" },
    m_pick_chain_id{ "A,AA,BA" }, m_pick_indicator{ ".,A" },
    m_pick_residue{ "" }, m_pick_element{ "" },
    m_pick_remoteness{ "" }, m_pick_branch{ "" },
    m_veto_chain_id{ "" }, m_veto_indicator{ "" },
    m_veto_residue{ "" }, m_veto_element{ "H" },
    m_veto_remoteness{ "" }, m_veto_branch{ "" }
{

}

void PotentialAnalysisCommand::Execute(void)
{
    std::cout << "PotentialAnalysisCommand::Execute() called." << std::endl;

    auto data_manager{ std::make_unique<DataObjectManager>(m_database_path) };
    data_manager->ProcessFile(m_model_file_path, "model");
    data_manager->ProcessFile(m_map_file_path, "map");

    auto atom_selector{ std::make_shared<AtomSelector>() };
    atom_selector->PickChainID(m_pick_chain_id);
    atom_selector->PickIndicator(m_pick_indicator);
    atom_selector->PickResidueType(m_pick_residue);
    atom_selector->PickElementType(m_pick_element);
    atom_selector->PickRemotenessType(m_pick_remoteness);
    atom_selector->PickBranchType(m_pick_branch);
    atom_selector->VetoChainID(m_veto_chain_id);
    atom_selector->VetoIndicator(m_veto_indicator);
    atom_selector->VetoResidueType(m_veto_residue);
    atom_selector->VetoElementType(m_veto_element);
    atom_selector->VetoRemotenessType(m_veto_remoteness);
    atom_selector->VetoBranchType(m_veto_branch);

    auto sphere_sampler{ std::make_shared<SphereSampler>() };
    sphere_sampler->SetThreadSize(m_thread_size);
    sphere_sampler->SetSamplingSize(m_sampling_size);
    sphere_sampler->SetDistanceRangeMinimum(m_sampling_range_min);
    sphere_sampler->SetDistanceRangeMaximum(m_sampling_range_max);

    auto analyzer{ std::make_unique<PotentialAnalysisVisitor>(atom_selector, sphere_sampler) };
    analyzer->SetThreadSize(m_thread_size);
    analyzer->SetModelObjectKeyTag("model");
    analyzer->SetMapObjectKeyTag("map");
    analyzer->SetFitRange(m_fit_range_min, m_fit_range_max);
    analyzer->SetAlphaR(m_alpha_r);
    analyzer->SetAlphaG(m_alpha_g);

    data_manager->Accept(analyzer.get());

    data_manager->SaveDataObject("model");
}