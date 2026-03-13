#include "PotentialDisplayCommand.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include "CommandDataSupport.hpp"
#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/core/painter/PainterBase.hpp>
#include <rhbm_gem/core/painter/GausPainter.hpp>
#include <rhbm_gem/core/painter/AtomPainter.hpp>
#include <rhbm_gem/core/painter/ModelPainter.hpp>
#include <rhbm_gem/core/painter/ComparisonPainter.hpp>
#include <rhbm_gem/core/painter/DemoPainter.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/AtomSelector.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/core/command/OptionEnumTraits.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>

namespace {
constexpr std::string_view kPainterOption{ "--painter" };
constexpr std::string_view kModelKeyListOption{ "--model-keylist" };
constexpr std::string_view kRefModelKeyListOption{ "--ref-model-keylist" };

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

void IngestModelSetsToPainter(
    rhbm_gem::PainterBase & painter,
    const std::vector<std::shared_ptr<rhbm_gem::ModelObject>> & model_object_list,
    const std::unordered_map<std::string, std::vector<std::shared_ptr<rhbm_gem::ModelObject>>> &
        ref_model_object_list_map)
{
    for (const auto & model_object : model_object_list)
    {
        painter.AddDataObject(model_object.get());
    }
    for (const auto & [class_key, ref_model_list] : ref_model_object_list_map)
    {
        for (const auto & model_object : ref_model_list)
        {
            painter.AddReferenceDataObject(model_object.get(), class_key);
        }
    }
}
}

namespace rhbm_gem {

PotentialDisplayCommand::PotentialDisplayCommand(CommonOptionProfile profile) :
    CommandWithOptions<PotentialDisplayCommandOptions>{
        CommonOptionMaskForProfile(profile) },
    m_atom_selector{ std::make_unique<AtomSelector>() }
{
}

PotentialDisplayCommand::~PotentialDisplayCommand() = default;

void PotentialDisplayCommand::ApplyRequest(const PotentialDisplayRequest & request)
{
    ApplyCommonRequest(request.common);
    SetPainterChoice(request.painter_choice);
    SetModelKeyTagList(request.model_key_tag_list);
    SetRefModelKeyTagListMap(request.ref_model_key_tag_list);
    MutateOptions([&]()
    {
        m_options.pick_chain_id = request.pick_chain_id;
        m_options.veto_chain_id = request.veto_chain_id;
        m_options.pick_residue = request.pick_residue;
        m_options.veto_residue = request.veto_residue;
        m_options.pick_element = request.pick_element;
        m_options.veto_element = request.veto_element;
    });
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
        kPainterOption,
        PainterType::MODEL,
        "Painter choice");
}

void PotentialDisplayCommand::SetModelKeyTagList(const std::string & value)
{
    MutateOptions([&]()
    {
        m_options.model_key_tag_list = StringHelper::ParseListOption<std::string>(value, ',');
        ResetParseIssues(kModelKeyListOption);
        if (m_options.model_key_tag_list.empty())
        {
            AddValidationError(
                kModelKeyListOption,
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
        ResetParseIssues(kRefModelKeyListOption);
        if (value.empty()) return;

        std::string error_message;
        if (!ParseReferenceModelKeyTagListMap(value, m_ref_model_key_tag_list_map, error_message))
        {
            AddValidationError(kRefModelKeyListOption, error_message, ValidationPhase::Parse);
            return;
        }

        Logger::Log(
            LogLevel::Debug,
            "Parsed " + std::to_string(m_ref_model_key_tag_list_map.size())
            + " reference model groups.");
    });
}

bool PotentialDisplayCommand::BuildDataObject()
{
    ScopeTimer timer{ "PotentialDisplayCommand::BuildDataObject" };
    try
    {
        m_data_manager.SetDatabaseManager(m_options.database_path);
        auto model_size{ m_options.model_key_tag_list.size() };
        size_t model_count{ 1 };
        Logger::Log(LogLevel::Info, "Load model object list:");
        for (auto & key : m_options.model_key_tag_list)
        {
            Logger::ProgressBar(model_count, model_size);
            m_model_object_list.emplace_back(
                command_data_loader::LoadModelObject(
                    m_data_manager, key, "model object"));
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
                    command_data_loader::LoadModelObject(
                        m_data_manager, key_tag, "reference model object"));
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
        ApplyModelSelection(*model_object, *m_atom_selector);
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
            IngestModelSetsToPainter(*painter, m_model_object_list, m_ref_model_object_list_map);
            break;
        case PainterType::MODEL:
            painter = std::make_unique<ModelPainter>();
            IngestModelSetsToPainter(*painter, m_model_object_list, m_ref_model_object_list_map);
            break;
        case PainterType::COMPARISON:
            painter = std::make_unique<ComparisonPainter>();
            IngestModelSetsToPainter(*painter, m_model_object_list, m_ref_model_object_list_map);
            break;
        case PainterType::DEMO:
            painter = std::make_unique<DemoPainter>();
            IngestModelSetsToPainter(*painter, m_model_object_list, m_ref_model_object_list_map);
            break;
        case PainterType::ATOM:
            painter = std::make_unique<AtomPainter>();
            for (const auto & model_object : m_model_object_list)
            {
                ModelAtomCollectorOptions collector_options;
                collector_options.selected_only = true;
                auto atom_list{ CollectModelAtoms(*model_object, collector_options) };
                for (auto * atom : atom_list)
                {
                    painter->AddDataObject(atom);
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
