#include "ChargeAnalysisCommand.hpp"
#include "DataObjectManager.hpp"
#include "ChargeAnalysisVisitor.hpp"

#include <iostream>
#include <sstream>

ChargeAnalysisCommand::ChargeAnalysisCommand(void)
{

}

ChargeAnalysisCommand::~ChargeAnalysisCommand()
{
    
}

void ChargeAnalysisCommand::Execute(void)
{
    std::cout << "ChargeAnalysisCommand::Execute() called." << std::endl;

    auto data_manager{ std::make_unique<DataObjectManager>(m_database_path) };
    data_manager->LoadDataObject(m_model_key_tag);
    data_manager->LoadDataObject(m_sim_neutral_model_key_tag);
    data_manager->LoadDataObject(m_sim_positive_model_key_tag);
    data_manager->LoadDataObject(m_sim_negative_model_key_tag);

    auto analyzer{ std::make_unique<ChargeAnalysisVisitor>() };
    analyzer->SetThreadSize(static_cast<unsigned int>(m_thread_size));
    analyzer->SetModelObjectKeyTag(m_model_key_tag);
    analyzer->SetNeutralModelObjectKeyTag(m_sim_neutral_model_key_tag);
    analyzer->SetPositiveModelObjectKeyTag(m_sim_positive_model_key_tag);
    analyzer->SetNegativeModelObjectKeyTag(m_sim_negative_model_key_tag);
    analyzer->SetFitRange(m_fit_range_min, m_fit_range_max);
    analyzer->SetAlphaR(m_alpha_r);
    analyzer->SetAlphaG(m_alpha_g);

    data_manager->Accept(analyzer.get());

    //data_manager->SaveDataObject(m_model_key_tag);
}