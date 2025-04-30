#include "PotentialDisplayCommand.hpp"
#include "DataObjectManager.hpp"
#include "PotentialDisplayVisitor.hpp"

#include <iostream>
#include <sstream>

PotentialDisplayCommand::PotentialDisplayCommand(void)
{

}

PotentialDisplayCommand::~PotentialDisplayCommand()
{

}

void PotentialDisplayCommand::Execute(void)
{
    std::cout << "PotentialDisplayCommand::Execute() called." << std::endl;
    std::cout << "Total number of model object sets to be display: "
              << m_model_key_tag_list.size() << std::endl;

    auto data_manager{ std::make_unique<DataObjectManager>(m_database_path) };
    LoadModelObjects(data_manager.get());
    LoadRefModelObjects(data_manager.get());

    auto model_displayer{ std::make_unique<PotentialDisplayVisitor>() };
    model_displayer->SetModelObjectKeyTagList(m_model_key_tag_list);
    model_displayer->SetRefModelObjectKeyTagListMap(m_ref_model_key_tag_list_map);
    model_displayer->SetFolderPath(m_folder_path);
    data_manager->Accept(model_displayer.get());
}

void PotentialDisplayCommand::SetModelKeyTagList(const std::string & value)
{
    std::stringstream ss(value);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        if (segment == "") continue;
        m_model_key_tag_list.emplace_back(segment);
    }
}

void PotentialDisplayCommand::SetRefModelKeyTagList(
    const std::string & map_key, const std::string & value)
{
    std::stringstream ss(value);
    std::string segment;
    while (std::getline(ss, segment, ','))
    {
        if (segment == "") continue;
        m_ref_model_key_tag_list_map[map_key].emplace_back(segment);
    }
}

void PotentialDisplayCommand::LoadModelObjects(DataObjectManager * data_manager)
{
    for (auto & key : m_model_key_tag_list)
    {
        data_manager->LoadDataObject(key);
    }
}

void PotentialDisplayCommand::LoadRefModelObjects(DataObjectManager * data_manager)
{
    for (auto & [map_key, key_tag_list] : m_ref_model_key_tag_list_map)
    {
        for (auto & key_tag : key_tag_list)
        {
            data_manager->LoadDataObject(key_tag);
        }
    }
}