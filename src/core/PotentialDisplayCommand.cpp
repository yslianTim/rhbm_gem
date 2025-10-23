#include "PotentialDisplayCommand.hpp"
#include "DataObjectManager.hpp"
#include "AtomObject.hpp"
#include "ModelObject.hpp"
#include "MapObject.hpp"
#include "PainterBase.hpp"
#include "GausPainter.hpp"
#include "AtomPainter.hpp"
#include "ModelPainter.hpp"
#include "ComparisonPainter.hpp"
#include "DemoPainter.hpp"
#include "StringHelper.hpp"
#include "AtomSelector.hpp"
#include "FilePathHelper.hpp"
#include "Logger.hpp"
#include "ScopeTimer.hpp"
#include "CommandRegistry.hpp"

namespace {
CommandRegistrar<PotentialDisplayCommand> registrar_potential_display{
    "potential_display",
    "Run potential display"};
}

PotentialDisplayCommand::PotentialDisplayCommand(void) :
    CommandBase(), m_options{},
    m_atom_selector{ std::make_unique<AtomSelector>() }
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::PotentialDisplayCommand() called");
}

PotentialDisplayCommand::~PotentialDisplayCommand()
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::~PotentialDisplayCommand() called");
}

void PotentialDisplayCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::RegisterCLIOptionsExtend() called");
    std::map<std::string, PainterType> painter_map
    {
        {"0", PainterType::GAUS},       {"gaus",       PainterType::GAUS},
        {"1", PainterType::MODEL},      {"model",      PainterType::MODEL},
        {"2", PainterType::COMPARISON}, {"comparison", PainterType::COMPARISON},
        {"3", PainterType::DEMO},       {"demo",       PainterType::DEMO},
        {"4", PainterType::ATOM},       {"atom",       PainterType::ATOM},
    };
    cmd->add_option_function<PainterType>("-p,--painter",
        [&](PainterType value) { SetPainterChoice(value); },
        "Painter choice")->required()
        ->transform(CLI::CheckedTransformer(painter_map, CLI::ignore_case));
    cmd->add_option_function<std::string>("-k,--model-keylist",
        [&](const std::string & value) { SetModelKeyTagList(value); },
        "List of model key tag to be display")->required();
    cmd->add_option_function<std::string>("-r,--ref-model-keylist",
        [&](const std::string & value) { SetRefModelKeyTagListMap(value); },
        "List of reference model key tag to be display")
        ->default_val(m_options.ref_model_key_tag_list);
    cmd->add_option_function<std::string>("--pick-chain",
        [&](const std::string & value) { SetPickChainID(value); },
        "Pick chain ID")->default_val(m_options.pick_chain_id);
    cmd->add_option_function<std::string>("--pick-residue",
        [&](const std::string & value) { SetPickResidueType(value); },
        "Pick residue type")->default_val(m_options.pick_residue);
    cmd->add_option_function<std::string>("--pick-element",
        [&](const std::string & value) { SetPickElementType(value); },
        "Pick element type")->default_val(m_options.pick_element);
    cmd->add_option_function<std::string>("--veto-chain",
        [&](const std::string & value) { SetVetoChainID(value); },
        "Veto chain ID")->default_val(m_options.veto_chain_id);
    cmd->add_option_function<std::string>("--veto-residue",
        [&](const std::string & value) { SetVetoResidueType(value); },
        "Veto residue type")->default_val(m_options.veto_residue);
    cmd->add_option_function<std::string>("--veto-element",
        [&](const std::string & value) { SetVetoElementType(value); },
        "Veto element type")->default_val(m_options.veto_element);
}

bool PotentialDisplayCommand::Execute(void)
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::Execute() called");
    if (BuildDataObject() == false) return false;
    RunDataObjectSelection();
    RunDisplay();
    return true;
}

void PotentialDisplayCommand::SetPainterChoice(PainterType value)
{
    m_options.painter_choice = value;
}

void PotentialDisplayCommand::SetModelKeyTagList(const std::string & value)
{
    m_options.model_key_tag_list = StringHelper::ParseListOption<std::string>(value, ',');
    if (m_options.model_key_tag_list.empty())
    {
        Logger::Log(LogLevel::Error, "Model key list cannot be empty");
        m_valiate_options = false;
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

void PotentialDisplayCommand::SetPickChainID(const std::string & value)
{
    m_options.pick_chain_id = value;
}

void PotentialDisplayCommand::SetVetoChainID(const std::string & value)
{
    m_options.veto_chain_id = value;
}

void PotentialDisplayCommand::SetPickResidueType(const std::string & value)
{
    m_options.pick_residue = value;
}

void PotentialDisplayCommand::SetVetoResidueType(const std::string & value)
{
    m_options.veto_residue = value;
}

void PotentialDisplayCommand::SetPickElementType(const std::string & value)
{
    m_options.pick_element = value;
}

void PotentialDisplayCommand::SetVetoElementType(const std::string & value)
{
    m_options.veto_element = value;
}

bool PotentialDisplayCommand::BuildDataObject(void)
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::BuildDataObject() called");
    ScopeTimer timer{ "PotentialDisplayCommand::BuildDataObject" };
    try
    {
        auto data_manager{ GetDataManagerPtr() };
        data_manager->SetDatabaseManager(m_options.database_path);
        auto model_size{ m_options.model_key_tag_list.size() };
        size_t model_count{ 1 };
        Logger::Log(LogLevel::Info, "Load model object list:");
        for (auto & key : m_options.model_key_tag_list)
        {
            Logger::ProgressBar(model_count, model_size);
            data_manager->LoadDataObject(key);
            m_model_object_list.emplace_back(data_manager->GetTypedDataObject<ModelObject>(key));
            model_count++;
        }
        for (auto & [map_key, key_tag_list] : m_ref_model_key_tag_list_map)
        {
            auto ref_model_size{ key_tag_list.size() };
            size_t ref_model_count{ 1 };
            Logger::Log(LogLevel::Info, "Load ["+ map_key +"] reference model object list:");
            for (auto & key_tag : key_tag_list)
            {
                Logger::ProgressBar(ref_model_count, ref_model_size);
                data_manager->LoadDataObject(key_tag);
                m_ref_model_object_list_map[map_key].emplace_back(
                    data_manager->GetTypedDataObject<ModelObject>(key_tag)
                );
                ref_model_count++;
            }
        }
        Logger::Log(LogLevel::Info,
            "Total number of model object sets to be display: "
            + std::to_string(m_options.model_key_tag_list.size()));
    }
    catch(const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "PotentialDisplayCommand::BuildDataObject : " + std::string(e.what()));
        return false;
    }
    return true;
}

void PotentialDisplayCommand::RunDataObjectSelection(void)
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::RunDataObjectSelection() called");
    ScopeTimer timer{ "PotentialDisplayCommand::RunDataObjectSelection" };
    if (m_atom_selector == nullptr) return;
    m_atom_selector->PickChainID(m_options.pick_chain_id);
    m_atom_selector->PickResidueType(m_options.pick_residue);
    m_atom_selector->PickElementType(m_options.pick_element);
    m_atom_selector->VetoChainID(m_options.veto_chain_id);
    m_atom_selector->VetoResidueType(m_options.veto_residue);
    m_atom_selector->VetoElementType(m_options.veto_element);

    for (const auto & model_object : m_model_object_list)
    {
        for (auto & atom : model_object->GetAtomList())
        {
            atom->SetSelectedFlag(
                m_atom_selector->GetSelectionFlag(
                    atom->GetChainID(),
                    atom->GetResidue(),
                    atom->GetElement())
            );
        }
    }
}

void PotentialDisplayCommand::RunDisplay(void)
{
    Logger::Log(LogLevel::Debug, "PotentialDisplayCommand::RunDisplay() called");
    ScopeTimer timer{ "PotentialDisplayCommand::RunDisplay" };
    std::unique_ptr<PainterBase> painter{ nullptr };
    switch (m_options.painter_choice)
    {
        case PainterType::GAUS:
            painter = std::make_unique<GausPainter>();
            for (const auto & model_object : m_model_object_list)
            {
                painter->AddDataObject(model_object.get());
            }
            for (auto & [class_key, model_object_list] : m_ref_model_object_list_map)
            {
                for (const auto & model_object : model_object_list)
                {
                    painter->AddReferenceDataObject(model_object.get(), class_key);
                }
            }
            break;
        case PainterType::MODEL:
            painter = std::make_unique<ModelPainter>();
            for (const auto & model_object : m_model_object_list)
            {
                painter->AddDataObject(model_object.get());
            }
            for (auto & [class_key, model_object_list] : m_ref_model_object_list_map)
            {
                for (const auto & model_object : model_object_list)
                {
                    painter->AddReferenceDataObject(model_object.get(), class_key);
                }
            }
            break;
        case PainterType::COMPARISON:
            painter = std::make_unique<ComparisonPainter>();
            for (const auto & model_object : m_model_object_list)
            {
                painter->AddDataObject(model_object.get());
            }
            for (auto & [class_key, model_object_list] : m_ref_model_object_list_map)
            {
                for (const auto & model_object : model_object_list)
                {
                    painter->AddReferenceDataObject(model_object.get(), class_key);
                }
            }
            break;
        case PainterType::DEMO:
            painter = std::make_unique<DemoPainter>();
            for (const auto & model_object : m_model_object_list)
            {
                painter->AddDataObject(model_object.get());
            }
            for (auto & [class_key, model_object_list] : m_ref_model_object_list_map)
            {
                for (const auto & model_object : model_object_list)
                {
                    painter->AddReferenceDataObject(model_object.get(), class_key);
                }
            }
            break;
        case PainterType::ATOM:
            painter = std::make_unique<AtomPainter>();
            for (const auto & model_object : m_model_object_list)
            {
                for (auto & atom : model_object->GetAtomList())
                {
                    if (atom->GetSelectedFlag() == false) continue;
                    painter->AddDataObject(atom.get());
                }
            }
            m_atom_selector->Print();
            break;
        default:
            Logger::Log(LogLevel::Warning,
                        "Invalid painter choice input: ["
                        + std::to_string(static_cast<int>(m_options.painter_choice)) + "]");
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
        painter->SetFolder(m_options.folder_path.string());
        painter->Painting();
    }
}
