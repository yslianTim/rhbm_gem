#include "PotentialDisplayCommand.hpp"
#include "DataObjectManager.hpp"
#include "PotentialDisplayVisitor.hpp"

PotentialDisplayCommand::PotentialDisplayCommand(void) :
    m_database_path{ "" }, m_model_key_tag{ "" }
{

}

void PotentialDisplayCommand::Execute(void)
{
    std::cout << "PotentialDisplayCommand::Execute() called." << std::endl;

    auto data_manager{ std::make_unique<DataObjectManager>(m_database_path) };
    data_manager->LoadDataObject(m_model_key_tag);

    auto analyzer{ std::make_unique<PotentialDisplayVisitor>() };
    analyzer->SetModelObjectKeyTag(m_model_key_tag);

    data_manager->Accept(analyzer.get());

}