#include "ChargeAnalysisCommand.hpp"
#include "DataObjectManager.hpp"
#include "ChargeAnalysisVisitor.hpp"

#include <iostream>

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
    data_manager->ProcessFile(m_model_file_path, "model");
    data_manager->ProcessFile(m_map_file_path, "map");

    auto analyzer{ std::make_unique<ChargeAnalysisVisitor>() };
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