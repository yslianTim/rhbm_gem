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
    LoadAdditionalReferenceModelObjects(data_manager.get());

    m_atom_selector->Print();

    auto model_displayer{ std::make_unique<PotentialDisplayVisitor>(m_atom_selector.get()) };
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
    /*
    std::vector<std::string> element_list{"Ap", "Cp", "Nn", "On"};
    std::vector<std::string> charge_list{"1", "3", "5", "7"};
    std::unordered_map<std::string, std::vector<std::string>> additional_list_map;
    for (auto & element : element_list)
    {
        for (auto & charge : charge_list)
        {
            auto model_class{ element + charge };
            std::vector<std::string> key_tag_list;
            key_tag_list.reserve(10);
            for (int i = 0; i < 10; i++)
            {
                key_tag_list.emplace_back(model_class + "_b" + std::to_string(i));
            }
            additional_list_map[model_class] = key_tag_list;
            m_ref_model_key_tag_list_map[model_class] = key_tag_list;
            
        }
    }
    for (auto & [model_class, key_tag_list] : additional_list_map)
    {
        std::cout <<"["<< model_class <<"]"<< std::endl;
        for (auto & key_tag : key_tag_list)
        {
            data_manager->LoadDataObject(key_tag);
        }
    }*/

    std::vector<std::string> key_tag_list
    {
        "amber_b0","amber_b1","amber_b2","amber_b3","amber_b4",
        "amber_b5","amber_b6","amber_b7","amber_b8","amber_b9"
    };
    m_ref_model_key_tag_list_map["amber95"] = key_tag_list;
    for (auto & key_tag : key_tag_list)
    {
        data_manager->LoadDataObject(key_tag);
    }

    std::vector<std::string> sim_key_tag_list
    {
        "Nn1_b0","Nn1_b1","Nn1_b2","Nn1_b3","Nn1_b4",
        "Nn1_b5","Nn1_b6","Nn1_b7","Nn1_b8","Nn1_b9"
    };
    m_ref_model_key_tag_list_map["sim_test"] = sim_key_tag_list;
    for (auto & key_tag : sim_key_tag_list)
    {
        data_manager->LoadDataObject(key_tag);
    }
}