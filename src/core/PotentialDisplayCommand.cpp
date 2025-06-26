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

void PotentialDisplayCommand::SetRefModelKeyTagListMap(const std::string & value)
{
    size_t pos{ 0 };
    size_t len{ value.size() };

    while (pos < len)
    {
        // Find '[' for start of group name
        if (value[pos] != '[') {
            throw std::runtime_error("Parser Error : expect '['");
        }
        // Find ']' for end of group name
        size_t end_name{ value.find(']', pos+1) };
        if (end_name == std::string::npos) {
            throw std::runtime_error("Parser Error : expect ']'");
        }
        std::string group_name{ value.substr(pos+1, end_name - (pos+1)) };

        // Find the start of members after ']'
        size_t start_members{ end_name + 1 };
        size_t end_block{ value.find(';', start_members) };
        if (end_block == std::string::npos) {
            end_block = len;
        }
        std::string members_str{ value.substr(start_members, end_block - start_members) };

        // Segment members_str by commas
        std::vector<std::string> members;
        std::istringstream iss(members_str);
        std::string member_token;
        while (std::getline(iss, member_token, ',')) {
            if (!member_token.empty()) {
                members.push_back(member_token);
            }
        }

        m_ref_model_key_tag_list_map.emplace(std::move(group_name), std::move(members));

        // Jump to the next block, which is after the semicolon
        pos = end_block + 1;
    }

    // Print the parsed model key tag list
    std::cout << "Parsed reference model key tag list: " << std::endl;
    for (const auto & [group_name, key_tags] : m_ref_model_key_tag_list_map)
    {
        std::cout << "Group: [" << group_name << "] -> ";
        for (const auto & key_tag : key_tags)
        {
            std::cout << key_tag << " ";
        }
        std::cout << std::endl;
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
