#include "PotentialDisplayCommand.hpp"
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
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>

#include <stdexcept>

namespace {
constexpr std::string_view kPainterOption{ "--painter" };
constexpr std::string_view kModelKeyListOption{ "--model-keylist" };
constexpr std::string_view kRefGroupOption{ "--ref-group" };

void ApplyModelSelection(
    rhbm_gem::ModelObject & model_object,
    ::AtomSelector & selector)
{
    for (auto & atom : model_object.GetAtomList())
    {
        atom->SetSelectedFlag(
            selector.GetSelectionFlag(
                atom->GetChainID(),
                atom->GetResidue(),
                atom->GetElement()));
    }
    model_object.Update();
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

PotentialDisplayCommand::PotentialDisplayCommand() :
    CommandWithRequest<PotentialDisplayRequest>{},
    m_atom_selector{ std::make_unique<AtomSelector>() }
{
}

void PotentialDisplayCommand::NormalizeRequest()
{
    auto & request{ MutableRequest() };
    CoerceEnum(
        request.painter_choice,
        kPainterOption,
        PainterType::MODEL,
        "Painter choice");
    RequireNonEmptyList(request.model_key_tag_list, kModelKeyListOption, "Model key list");
    InvalidatePreparedState();
    ClearParseIssues(kRefGroupOption);
    for (const auto & [group_name, members] : request.reference_model_groups)
    {
        if (group_name.empty())
        {
            AddValidationError(
                kRefGroupOption,
                "Reference group name cannot be empty.",
                ValidationPhase::Parse);
            continue;
        }
        if (members.empty())
        {
            AddValidationError(
                kRefGroupOption,
                "Reference group '" + group_name + "' cannot be empty.",
                ValidationPhase::Parse);
        }
    }
}

PotentialDisplayCommand::~PotentialDisplayCommand() = default;

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

bool PotentialDisplayCommand::BuildDataObject()
{
    const auto & request{ RequestOptions() };
    ScopeTimer timer{ "PotentialDisplayCommand::BuildDataObject" };
    try
    {
        m_data_manager.OpenDatabase(request.database_path);
        auto model_size{ request.model_key_tag_list.size() };
        size_t model_count{ 1 };
        Logger::Log(LogLevel::Info, "Load model object list:");
        for (const auto & key : request.model_key_tag_list)
        {
            Logger::ProgressBar(model_count, model_size);
            m_data_manager.LoadDataObject(key);
            m_model_object_list.emplace_back(
                m_data_manager.GetTypedDataObject<ModelObject>(key));
            model_count++;
        }
        for (const auto & [map_key, key_tag_list] : request.reference_model_groups)
        {
            auto ref_model_size{ key_tag_list.size() };
            size_t ref_model_count{ 1 };
            Logger::Log(LogLevel::Info, "Load ["+ map_key +"] reference model object list:");
            for (auto & key_tag : key_tag_list)
            {
                Logger::ProgressBar(ref_model_count, ref_model_size);
                m_data_manager.LoadDataObject(key_tag);
                m_ref_model_object_list_map[map_key].emplace_back(
                    m_data_manager.GetTypedDataObject<ModelObject>(key_tag));
                ref_model_count++;
            }
        }
        Logger::Log(LogLevel::Info,
            "Total number of model object sets to be display: "
            + std::to_string(request.model_key_tag_list.size()));
    }
    catch(const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "PotentialDisplayCommand::BuildDataObject : Failed to load model data from database: "
                + std::string(e.what()));
        return false;
    }
    return true;
}

void PotentialDisplayCommand::RunDataObjectSelection()
{
    const auto & request{ RequestOptions() };
    ScopeTimer timer{ "PotentialDisplayCommand::RunDataObjectSelection" };
    if (m_atom_selector == nullptr) return;
    m_atom_selector->PickChainID(request.pick_chain_id);
    m_atom_selector->PickResidueType(request.pick_residue);
    m_atom_selector->PickElementType(request.pick_element);
    m_atom_selector->VetoChainID(request.veto_chain_id);
    m_atom_selector->VetoResidueType(request.veto_residue);
    m_atom_selector->VetoElementType(request.veto_element);

    for (const auto & model_object : m_model_object_list)
    {
        ApplyModelSelection(*model_object, *m_atom_selector);
    }
}

void PotentialDisplayCommand::RunDisplay()
{
    const auto & request{ RequestOptions() };
    ScopeTimer timer{ "PotentialDisplayCommand::RunDisplay" };
    std::unique_ptr<PainterBase> painter{ nullptr };
    switch (request.painter_choice)
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
                for (auto * atom : model_object->GetSelectedAtomList())
                {
                    painter->AddDataObject(atom);
                }
            }
            m_atom_selector->Print();
            break;
        default:
            Logger::Log(LogLevel::Warning,
                        "Invalid painter choice input: ["
                        + std::to_string(static_cast<int>(request.painter_choice)) + "]");
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
        painter->SetFolder(OutputFolder().string());
        painter->Painting();
    }
}

} // namespace rhbm_gem
