#include "PotentialDisplayCommand.hpp"
#include "DataObjectManager.hpp"
#include "PotentialDisplayVisitor.hpp"

PotentialDisplayCommand::PotentialDisplayCommand(void)
{

}

void PotentialDisplayCommand::Execute(void)
{
    std::cout << "PotentialDisplayCommand::Execute() called." << std::endl;

    auto data_manager{ std::make_unique<DataObjectManager>(m_database_path) };
    data_manager->LoadDataObject(m_model_key_tag);

    auto analyzer{ std::make_unique<PotentialDisplayVisitor>() };
    analyzer->SetModelObjectKeyTag(m_model_key_tag);
    analyzer->SetFolderPath(m_folder_path);

    data_manager->Accept(analyzer.get());

}