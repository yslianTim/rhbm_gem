#include "PotentialDisplayCommand.hpp"
#include "DataObjectManager.hpp"
#include "PotentialDisplayVisitor.hpp"
#include "AtomSelector.hpp"

#include <iostream>
#include <sstream>

PotentialDisplayCommand::PotentialDisplayCommand(void) :
    m_atom_selector{ std::make_unique<AtomSelector>() }
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
    //LoadAdditionalReferenceModelObjects(data_manager.get());

    auto model_displayer{ std::make_unique<PotentialDisplayVisitor>(m_atom_selector.get()) };
    model_displayer->SetModelObjectKeyTagList(m_model_key_tag_list);
    model_displayer->SetRefModelObjectKeyTagListMap(m_ref_model_key_tag_list_map);
    model_displayer->SetFolderPath(m_folder_path);
    model_displayer->SetPainterChoice(m_painter_choice);
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

void PotentialDisplayCommand::SetPickChainID(const std::string & value)
{
    m_atom_selector->PickChainID(value);
}

void PotentialDisplayCommand::SetPickResidueType(const std::string & value)
{
    m_atom_selector->PickResidueType(value);
}

void PotentialDisplayCommand::SetPickElementType(const std::string & value)
{
    m_atom_selector->PickElementType(value);
}

void PotentialDisplayCommand::SetPickRemotenessType(const std::string & value)
{
    m_atom_selector->PickRemotenessType(value);
}

void PotentialDisplayCommand::SetVetoChainID(const std::string & value)
{
    m_atom_selector->VetoChainID(value);
}

void PotentialDisplayCommand::SetVetoResidueType(const std::string & value)
{
    m_atom_selector->VetoResidueType(value);
}

void PotentialDisplayCommand::SetVetoElementType(const std::string & value)
{
    m_atom_selector->VetoElementType(value);
}

void PotentialDisplayCommand::SetVetoRemotenessType(const std::string & value)
{
    m_atom_selector->VetoRemotenessType(value);
}

void PotentialDisplayCommand::LoadAdditionalReferenceModelObjects(DataObjectManager * data_manager)
{
    std::vector<std::string> key_tag_list
    {
        "amber_b15","amber_b25","amber_b35","amber_b45","amber_b55",
        "amber_b65","amber_b75","amber_b85","amber_b95"
    };
    m_ref_model_key_tag_list_map["amber95"] = key_tag_list;
    for (auto & key_tag : key_tag_list)
    {
        data_manager->LoadDataObject(key_tag);
    }
}
