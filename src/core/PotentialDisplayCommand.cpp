#include "PotentialDisplayCommand.hpp"
#include "DataObjectManager.hpp"
#include "PotentialDisplayVisitor.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"

PotentialDisplayCommand::PotentialDisplayCommand(void) :
    m_atom_selector{ std::make_unique<AtomSelector>() }
{

}

void PotentialDisplayCommand::RegisterCLIOptions(CLI::App * cmd)
{
    cmd->add_option("-p,--painter", m_options.painter_choice,
        "Painter choice")->required();
    cmd->add_option("-k,--model-keylist", m_options.model_key_tag_list,
        "List of model key tag to be display")->required();
    cmd->add_option("-r,--ref-model-keylist", m_options.ref_model_key_tag_list,
        "List of reference model key tag to be display")->default_val(m_options.ref_model_key_tag_list);
    cmd->add_option("--pick-chain", m_options.pick_chain_id,
        "Pick chain ID")->default_val(m_options.pick_chain_id);
    cmd->add_option("--pick-residue", m_options.pick_residue,
        "Pick residue type")->default_val(m_options.pick_residue);
    cmd->add_option("--pick-element", m_options.pick_element,
        "Pick element type")->default_val(m_options.pick_element);
    cmd->add_option("--pick-remoteness", m_options.pick_remoteness,
        "Pick remoteness type")->default_val(m_options.pick_remoteness);
    cmd->add_option("--veto-chain", m_options.veto_chain_id,
        "Veto chain ID")->default_val(m_options.veto_chain_id);
    cmd->add_option("--veto-residue", m_options.veto_residue,
        "Veto residue type")->default_val(m_options.veto_residue);
    cmd->add_option("--veto-element", m_options.veto_element,
        "Veto element type")->default_val(m_options.veto_element);
    cmd->add_option("--veto-remoteness", m_options.veto_remoteness,
        "Veto remoteness type")->default_val(m_options.veto_remoteness);
    RegisterCommandOptions(cmd);
}

bool PotentialDisplayCommand::Execute(void)
{
    Logger::Log(LogLevel::Info, "PotentialDisplayCommand::Execute() called.");
    try
    {
        SetModelKeyTagList(m_options.model_key_tag_list);
        SetRefModelKeyTagListMap(m_options.ref_model_key_tag_list);
        SetPickChainID(m_options.pick_chain_id);
        SetPickResidueType(m_options.pick_residue);
        SetPickElementType(m_options.pick_element);
        SetPickRemotenessType(m_options.pick_remoteness);
        SetVetoChainID(m_options.veto_chain_id);
        SetVetoResidueType(m_options.veto_residue);
        SetVetoElementType(m_options.veto_element);
        SetVetoRemotenessType(m_options.veto_remoteness);
        Logger::Log(LogLevel::Info, "Total number of model object sets to be display: "
                    + std::to_string(m_model_key_tag_list.size()));

        auto data_manager{ std::make_unique<DataObjectManager>(m_options.database_path) };
        LoadModelObjects(data_manager.get());
        LoadRefModelObjects(data_manager.get());

        auto model_displayer{ std::make_unique<PotentialDisplayVisitor>(m_atom_selector.get()) };
        model_displayer->SetModelObjectKeyTagList(m_model_key_tag_list);
        model_displayer->SetRefModelObjectKeyTagListMap(m_ref_model_key_tag_list_map);
        model_displayer->SetFolderPath(m_options.folder_path);
        model_displayer->SetPainterChoice(m_options.painter_choice);

        data_manager->Accept(model_displayer.get());
    }
    catch(const std::exception & e)
    {
        Logger::Log(LogLevel::Error, e.what());
        return false;
    }
    return true;
}

void PotentialDisplayCommand::SetModelKeyTagList(const std::string & value)
{
    m_model_key_tag_list.clear();
    for (const auto & token : StringHelper::SplitStringLineFromDelimiter(value, ','))
    {
        m_model_key_tag_list.emplace_back(token);
    }
}

void PotentialDisplayCommand::SetRefModelKeyTagListMap(const std::string & value)
{
    m_ref_model_key_tag_list_map.clear();
    size_t pos{ 0 };
    size_t len{ value.size() };

    while (pos < len)
    {
        // Find '[' for start of group name
        if (value[pos] != '[')
        {
            throw std::runtime_error("Parser Error : expect '['");
        }
        // Find ']' for end of group name
        size_t end_name{ value.find(']', pos+1) };
        if (end_name == std::string::npos)
        {
            throw std::runtime_error("Parser Error : expect ']'");
        }
        std::string group_name{ value.substr(pos+1, end_name - (pos+1)) };

        // Find the start of members after ']'
        size_t start_members{ end_name + 1 };
        size_t end_block{ value.find(';', start_members) };
        if (end_block == std::string::npos)
        {
            end_block = len;
        }
        std::string members_string{ value.substr(start_members, end_block - start_members) };

        m_ref_model_key_tag_list_map.emplace(
            std::move(group_name),
            StringHelper::SplitStringLineFromDelimiter(members_string, ','));

        // Jump to the next block, which is after the semicolon
        pos = end_block + 1;
    }

    // Print the parsed model key tag list
    Logger::Log(LogLevel::Debug, "Parsed reference model key tag list:");
    for (const auto & [group_name, key_tags] : m_ref_model_key_tag_list_map)
    {
        std::string message{ "Group: [" + group_name + "] -> " };
        for (const auto & key_tag : key_tags)
        {
            message += key_tag + " ";
        }
        Logger::Log(LogLevel::Debug, message);
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
