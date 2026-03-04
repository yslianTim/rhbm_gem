#include "PotentialDisplayCommand.hpp"
#include "CommandOptionBinding.hpp"
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
#include "OptionEnumTraits.hpp"
#include "ScopeTimer.hpp"

namespace {
bool ParseReferenceModelKeyTagListMap(
    const std::string & value,
    std::unordered_map<std::string, std::vector<std::string>> & output_map,
    std::string & error_message)
{
    output_map.clear();
    size_t pos{ 0 };
    const size_t len{ value.size() };

    while (pos < len)
    {
        if (value[pos] != '[')
        {
            error_message = "Expected '[' at position " + std::to_string(pos) + ".";
            output_map.clear();
            return false;
        }

        const size_t end_name{ value.find(']', pos + 1) };
        if (end_name == std::string::npos)
        {
            error_message = "Expected ']' after reference group name.";
            output_map.clear();
            return false;
        }

        std::string group_name{ value.substr(pos + 1, end_name - (pos + 1)) };
        const size_t start_members{ end_name + 1 };
        size_t end_block{ value.find(';', start_members) };
        if (end_block == std::string::npos)
        {
            end_block = len;
        }

        std::string members_string{ value.substr(start_members, end_block - start_members) };
        output_map.emplace(
            std::move(group_name),
            StringHelper::SplitStringLineFromDelimiter(members_string, ','));
        pos = end_block + 1;
    }

    return true;
}
}

namespace rhbm_gem {

PotentialDisplayCommand::PotentialDisplayCommand() :
    m_atom_selector{ std::make_unique<AtomSelector>() }
{
}

PotentialDisplayCommand::~PotentialDisplayCommand()
{
}

void PotentialDisplayCommand::RegisterCLIOptionsExtend(CLI::App * cmd)
{
    command_cli::AddEnumOption<PainterType>(
        cmd, "-p,--painter",
        [&](PainterType value) { SetPainterChoice(value); },
        "Painter choice",
        std::nullopt,
        true);
    command_cli::AddStringOption(
        cmd, "-k,--model-keylist",
        [&](const std::string & value) { SetModelKeyTagList(value); },
        "List of model key tag to be display",
        std::nullopt,
        true);
    command_cli::AddStringOption(
        cmd, "-r,--ref-model-keylist",
        [&](const std::string & value) { SetRefModelKeyTagListMap(value); },
        "List of reference model key tag to be display",
        m_options.ref_model_key_tag_list);
    command_cli::AddStringOption(
        cmd, "--pick-chain",
        [&](const std::string & value) { SetPickChainID(value); },
        "Pick chain ID",
        m_options.pick_chain_id);
    command_cli::AddStringOption(
        cmd, "--pick-residue",
        [&](const std::string & value) { SetPickResidueType(value); },
        "Pick residue type",
        m_options.pick_residue);
    command_cli::AddStringOption(
        cmd, "--pick-element",
        [&](const std::string & value) { SetPickElementType(value); },
        "Pick element type",
        m_options.pick_element);
    command_cli::AddStringOption(
        cmd, "--veto-chain",
        [&](const std::string & value) { SetVetoChainID(value); },
        "Veto chain ID",
        m_options.veto_chain_id);
    command_cli::AddStringOption(
        cmd, "--veto-residue",
        [&](const std::string & value) { SetVetoResidueType(value); },
        "Veto residue type",
        m_options.veto_residue);
    command_cli::AddStringOption(
        cmd, "--veto-element",
        [&](const std::string & value) { SetVetoElementType(value); },
        "Veto element type",
        m_options.veto_element);
}

bool PotentialDisplayCommand::ExecuteImpl()
{
    if (BuildDataObject() == false) return false;
    RunDataObjectSelection();
    RunDisplay();
    return true;
}

void PotentialDisplayCommand::ResetRuntimeState()
{
    m_model_object_list.clear();
    m_ref_model_object_list_map.clear();
}

void PotentialDisplayCommand::SetPainterChoice(PainterType value)
{
    SetValidatedEnumOption(
        m_options.painter_choice,
        value,
        "--painter",
        PainterType::MODEL,
        "Painter choice");
}

void PotentialDisplayCommand::SetModelKeyTagList(const std::string & value)
{
    MutateOptions([&]()
    {
        m_options.model_key_tag_list = StringHelper::ParseListOption<std::string>(value, ',');
        ResetParseIssues("--model-keylist");
        if (m_options.model_key_tag_list.empty())
        {
            AddValidationError(
                "--model-keylist",
                "Model key list cannot be empty.",
                ValidationPhase::Parse);
        }
    });
}

void PotentialDisplayCommand::SetRefModelKeyTagListMap(const std::string & value)
{
    MutateOptions([&]()
    {
        m_options.ref_model_key_tag_list = value;
        m_ref_model_key_tag_list_map.clear();
        ResetParseIssues("--ref-model-keylist");
        if (value.empty()) return;

        std::string error_message;
        if (!ParseReferenceModelKeyTagListMap(value, m_ref_model_key_tag_list_map, error_message))
        {
            AddValidationError("--ref-model-keylist", error_message, ValidationPhase::Parse);
            return;
        }

        Logger::Log(
            LogLevel::Debug,
            "Parsed " + std::to_string(m_ref_model_key_tag_list_map.size())
            + " reference model groups.");
    });
}

void PotentialDisplayCommand::SetPickChainID(const std::string & value)
{
    MutateOptions([&]() { m_options.pick_chain_id = value; });
}

void PotentialDisplayCommand::SetVetoChainID(const std::string & value)
{
    MutateOptions([&]() { m_options.veto_chain_id = value; });
}

void PotentialDisplayCommand::SetPickResidueType(const std::string & value)
{
    MutateOptions([&]() { m_options.pick_residue = value; });
}

void PotentialDisplayCommand::SetVetoResidueType(const std::string & value)
{
    MutateOptions([&]() { m_options.veto_residue = value; });
}

void PotentialDisplayCommand::SetPickElementType(const std::string & value)
{
    MutateOptions([&]() { m_options.pick_element = value; });
}

void PotentialDisplayCommand::SetVetoElementType(const std::string & value)
{
    MutateOptions([&]() { m_options.veto_element = value; });
}

bool PotentialDisplayCommand::BuildDataObject()
{
    ScopeTimer timer{ "PotentialDisplayCommand::BuildDataObject" };
    try
    {
        RequireDatabaseManager();
        auto model_size{ m_options.model_key_tag_list.size() };
        size_t model_count{ 1 };
        Logger::Log(LogLevel::Info, "Load model object list:");
        for (auto & key : m_options.model_key_tag_list)
        {
            Logger::ProgressBar(model_count, model_size);
            m_model_object_list.emplace_back(LoadTypedObject<ModelObject>(key, "model object"));
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
                m_ref_model_object_list_map[map_key].emplace_back(
                    LoadTypedObject<ModelObject>(key_tag, "reference model object")
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

void PotentialDisplayCommand::RunDataObjectSelection()
{
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

void PotentialDisplayCommand::RunDisplay()
{
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

} // namespace rhbm_gem
