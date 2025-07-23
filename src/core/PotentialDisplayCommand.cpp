#include "PotentialDisplayCommand.hpp"
#include "DataObjectManager.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "PainterBase.hpp"
#include "AtomPainter.hpp"
#include "ModelPainter.hpp"
#include "ComparisonPainter.hpp"
#include "DemoPainter.hpp"
#include "StringHelper.hpp"
#include "AtomSelector.hpp"
#include "Logger.hpp"
#include "CommandRegistry.hpp"

namespace {
CommandRegistrar<PotentialDisplayCommand> registrar_potential_display{
    "potential_display",
    "Run potential display"};
}

PotentialDisplayCommand::PotentialDisplayCommand(void) :
    m_atom_selector{ std::make_unique<AtomSelector>() }
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::PotentialDisplayCommand() called");
}

PotentialDisplayCommand::~PotentialDisplayCommand()
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::~PotentialDisplayCommand() called");
    m_atom_selector.reset();
}

void PotentialDisplayCommand::RegisterCLIOptions(CLI::App * cmd)
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::RegisterCLIOptions() called");
    cmd->add_option("-p,--painter", m_options.painter_choice,
        "Painter choice")->required();
    cmd->add_option("-k,--model-keylist", m_options.model_key_tag_list,
        "List of model key tag to be display")->required()->delimiter(',');
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
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::Execute() called");
    try
    {
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
                    + std::to_string(m_options.model_key_tag_list.size()));

        auto data_manager{ std::make_unique<DataObjectManager>(m_options.database_path) };
        LoadModelObjects(data_manager.get());
        LoadRefModelObjects(data_manager.get());

        for (auto & key : m_options.model_key_tag_list)
        {
            auto object{ data_manager->GetTypedDataObjectPtr<ModelObject>(key) };
            m_model_object_list.emplace_back(object);
        }

        for (auto & [class_key, tag_list] : m_ref_model_key_tag_list_map)
        {
            for (auto & tag : tag_list)
            {
                auto object{ data_manager->GetTypedDataObjectPtr<ModelObject>(tag) };
                m_ref_model_object_list_map[class_key].emplace_back(object);
            }
        }

        for (auto * model_object : m_model_object_list)
        {
            for (auto & atom : model_object->GetComponentsList())
            {
                bool selected_flag{
                    m_atom_selector->GetSelectionFlag(
                        atom->GetChainID(),
                        atom->GetResidue(),
                        atom->GetElement(),
                        atom->GetRemoteness())
                };
                atom->SetSelectedFlag(selected_flag);
            }
        }
        BuildOrderedModelObjectList();
        BuildOrderedRefModelObjectListMap();
        RunDisplay();
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
    m_options.model_key_tag_list = StringHelper::ParseListOption<std::string>(value, ',');
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
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::LoadModelObjects() called");
    for (auto & key : m_options.model_key_tag_list)
    {
        data_manager->LoadDataObject(key);
    }
}

void PotentialDisplayCommand::LoadRefModelObjects(DataObjectManager * data_manager)
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::LoadRefModelObjects() called");
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

void PotentialDisplayCommand::BuildOrderedModelObjectList(void)
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::BuildOrderedModelObjectList() called");
    m_ordered_model_object_list.clear();
    for (const auto & key : m_options.model_key_tag_list)
    {
        for (auto * model_object : m_model_object_list)
        {
            if (model_object->GetKeyTag() == key)
            {
                m_ordered_model_object_list.emplace_back(model_object);
                break;
            }
        }
    }
}

void PotentialDisplayCommand::BuildOrderedRefModelObjectListMap(void)
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::BuildOrderedRefModelObjectListMap() called");
    m_ordered_ref_model_object_list_map.clear();
    for (auto & [class_key, ordered_key_list] : m_ref_model_key_tag_list_map)
    {
        auto it{ m_ref_model_object_list_map.find(class_key) };
        if (it == m_ref_model_object_list_map.end()) continue;
        auto & model_object_list{ it->second };
        std::vector<ModelObject *> ordered_model_object_list;
        for (const auto & key : ordered_key_list)
        {
            for (auto * model_object : model_object_list)
            {
                if (model_object->GetKeyTag() == key)
                {
                    ordered_model_object_list.emplace_back(model_object);
                    break;
                }
            }
        }

        m_ordered_ref_model_object_list_map[class_key] = std::move(ordered_model_object_list);
    }
}

void PotentialDisplayCommand::RunDisplay(void)
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::RunDisplay() called");
    std::unique_ptr<PainterBase> painter{ nullptr };
    switch (m_options.painter_choice)
    {
        case 0:
            painter = std::make_unique<AtomPainter>();
            for (auto * model_object : m_ordered_model_object_list)
            {
                for (auto & atom : model_object->GetComponentsList())
                {
                    if (atom->GetSelectedFlag() == false) continue;
                    painter->AddDataObject(atom.get());
                }
            }
            m_atom_selector->Print();
            break;
        case 1:
            painter = std::make_unique<ModelPainter>();
            for (auto * model_object : m_ordered_model_object_list)
            {
                painter->AddDataObject(model_object);
            }
            for (auto & [class_key, obj_list] : m_ordered_ref_model_object_list_map)
            {
                for (auto * obj : obj_list)
                {
                    painter->AddReferenceDataObject(obj, class_key);
                }
            }
            break;
        case 2:
            painter = std::make_unique<ComparisonPainter>();
            for (auto * model_object : m_ordered_model_object_list)
            {
                painter->AddDataObject(model_object);
            }
            for (auto & [class_key, obj_list] : m_ordered_ref_model_object_list_map)
            {
                for (auto * obj : obj_list)
                {
                    painter->AddReferenceDataObject(obj, class_key);
                }
            }
            break;
        case 3:
            painter = std::make_unique<DemoPainter>();
            for (auto * model_object : m_ordered_model_object_list)
            {
                painter->AddDataObject(model_object);
            }
            for (auto & [class_key, obj_list] : m_ordered_ref_model_object_list_map)
            {
                for (auto * obj : obj_list)
                {
                    painter->AddReferenceDataObject(obj, class_key);
                }
            }
            break;
        default:
            Logger::Log(LogLevel::Warning,
                        "Invalid painter choice input: ["
                        + std::to_string(m_options.painter_choice) + "]");
            Logger::Log(LogLevel::Warning,
                        "Available Painter Choices:\n"
                        "  [0] AtomPainter\n"
                        "  [1] ModelPainter\n"
                        "  [2] ComparisonPainter\n"
                        "  [3] DemoPainter");
            break;
    }
    if (painter)
    {
        painter->SetFolder(m_options.folder_path);
        painter->Painting();
    }
}
