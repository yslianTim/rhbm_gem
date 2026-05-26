#include "detail/CommandBase.hpp"
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "core/painter/GausPainter.hpp"
#include "core/painter/AtomPainter.hpp"
#include "core/painter/ModelPainter.hpp"
#include "core/painter/ComparisonPainter.hpp"
#include "core/painter/DemoPainter.hpp"
#include <rhbm_gem/utils/domain/AtomSelector.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>

#include <cstddef>
#include <exception>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace rhbm_gem {

class PotentialDisplayCommand final : public CommandBase<PotentialDisplayRequest>
{
public:
    PotentialDisplayCommand();

private:
    void NormalizeAndValidateRequest(PotentialDisplayRequest & request) override;
    bool ExecuteImpl(const PotentialDisplayRequest & request) override;
};

namespace {

struct PotentialDisplayInputs
{
    std::vector<std::unique_ptr<rhbm_gem::ModelObject>> model_objects;
    std::unordered_map<
        std::string,
        std::vector<std::unique_ptr<rhbm_gem::ModelObject>>> reference_model_groups;
};

std::optional<PotentialDisplayInputs> LoadPotentialDisplayInputs(
    const rhbm_gem::PotentialDisplayRequest & request)
{
    ScopeTimer timer{ "PotentialDisplayCommand::BuildDataObject" };
    try
    {
        DataRepository repository{ request.database_path };
        PotentialDisplayInputs inputs;

        auto model_size{ request.model_key_tag_list.size() };
        size_t model_count{ 1 };
        Logger::Log(LogLevel::Info, "Load model object list:");
        inputs.model_objects.reserve(model_size);
        for (const auto & key : request.model_key_tag_list)
        {
            Logger::ProgressBar(model_count, model_size);
            inputs.model_objects.emplace_back(repository.LoadModel(key));
            model_count++;
        }

        for (const auto & [map_key, key_tag_list] : request.reference_model_groups)
        {
            auto ref_model_size{ key_tag_list.size() };
            size_t ref_model_count{ 1 };
            auto & ref_model_objects{ inputs.reference_model_groups[map_key] };
            ref_model_objects.reserve(ref_model_size);
            Logger::Log(LogLevel::Info, "Load ["+ map_key +"] reference model object list:");
            for (const auto & key_tag : key_tag_list)
            {
                Logger::ProgressBar(ref_model_count, ref_model_size);
                ref_model_objects.emplace_back(repository.LoadModel(key_tag));
                ref_model_count++;
            }
        }

        Logger::Log(LogLevel::Info,
            "Total number of model object sets to be display: "
            + std::to_string(request.model_key_tag_list.size()));
        return inputs;
    }
    catch(const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "PotentialDisplayCommand::BuildDataObject : Failed to load model data from database: "
                + std::string(e.what()));
        return std::nullopt;
    }
}

AtomSelector BuildAtomSelector(const rhbm_gem::PotentialDisplayRequest & request)
{
    AtomSelector selector;
    selector.PickChainID(request.pick_chain_id);
    selector.PickResidueType(request.pick_residue);
    selector.PickElementType(request.pick_element);
    selector.VetoChainID(request.veto_chain_id);
    selector.VetoResidueType(request.veto_residue);
    selector.VetoElementType(request.veto_element);
    return selector;
}

void ApplyModelSelection(rhbm_gem::ModelObject & model_object, AtomSelector & selector)
{
    model_object.SelectAtoms(
        [&selector](const rhbm_gem::AtomObject & atom)
        {
            return selector.GetSelectionFlag(
                       atom.GetChainID(),
                       atom.GetResidue(),
                       atom.GetElement())
                && rhbm_gem::AtomLocalPotentialView::For(atom).IsAvailable();
        });
}

void ApplyDataObjectSelection(
    const std::vector<std::unique_ptr<rhbm_gem::ModelObject>> & model_objects,
    AtomSelector & selector)
{
    ScopeTimer timer{ "PotentialDisplayCommand::RunDataObjectSelection" };
    for (const auto & model_object : model_objects)
    {
        ApplyModelSelection(*model_object, selector);
    }
}

} // namespace

PotentialDisplayCommand::PotentialDisplayCommand() : CommandBase<PotentialDisplayRequest>{}
{
}

void PotentialDisplayCommand::NormalizeAndValidateRequest(PotentialDisplayRequest & request)
{
    RequireEnum(request, &PotentialDisplayRequest::painter_choice);
    RequireNonEmptyList(request, &PotentialDisplayRequest::model_key_tag_list);
    for (const auto & [group_name, members] : request.reference_model_groups)
    {
        if (group_name.empty())
        {
            AddFieldValidationError(&PotentialDisplayRequest::reference_model_groups,
                "Reference group name cannot be empty.");
            continue;
        }
        if (members.empty())
        {
            AddFieldValidationError(&PotentialDisplayRequest::reference_model_groups,
                "Reference group '" + group_name + "' cannot be empty.");
        }
    }
}

bool PotentialDisplayCommand::ExecuteImpl(const PotentialDisplayRequest & request)
{
    auto inputs{ LoadPotentialDisplayInputs(request) };
    if (!inputs.has_value()) return false;

    auto selector{ BuildAtomSelector(request) };
    ApplyDataObjectSelection(inputs->model_objects, selector);

    ScopeTimer timer{ "PotentialDisplayCommand::RunDisplay" };
    const auto output_folder{ request.output_dir.string() };
    switch (request.painter_choice)
    {
        case PainterType::GAUS:
        {
            GausPainter painter;
            for (const auto & model_object : inputs->model_objects)
            {
                painter.AddModel(*model_object);
            }
            painter.SetFolder(output_folder);
            painter.Painting();
            break;
        }
        case PainterType::MODEL:
        {
            ModelPainter painter;
            for (const auto & model_object : inputs->model_objects)
            {
                painter.AddModel(*model_object);
            }
            painter.SetFolder(output_folder);
            painter.Painting();
            break;
        }
        case PainterType::COMPARISON:
        {
            ComparisonPainter painter;
            for (const auto & model_object : inputs->model_objects)
            {
                painter.AddModel(*model_object);
            }
            for (const auto & [class_key, ref_model_list] : inputs->reference_model_groups)
            {
                for (const auto & model_object : ref_model_list)
                {
                    painter.AddReferenceModel(*model_object, class_key);
                }
            }
            painter.SetFolder(output_folder);
            painter.Painting();
            break;
        }
        case PainterType::DEMO:
        {
            DemoPainter painter;
            for (const auto & model_object : inputs->model_objects)
            {
                painter.AddModel(*model_object);
            }
            for (const auto & [class_key, ref_model_list] : inputs->reference_model_groups)
            {
                for (const auto & model_object : ref_model_list)
                {
                    painter.AddReferenceModel(*model_object, class_key);
                }
            }
            painter.SetFolder(output_folder);
            painter.Painting();
            break;
        }
        case PainterType::ATOM:
        {
            AtomPainter painter;
            for (const auto & model_object : inputs->model_objects)
            {
                painter.AddModel(*model_object);
            }
            selector.Print();
            painter.SetFolder(output_folder);
            painter.Painting();
            break;
        }
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
    return true;
}

namespace command_internal {

CommandResult ExecutePotentialDisplayCommand(const PotentialDisplayRequest & request)
{
    PotentialDisplayCommand command;
    return command.ExecuteRequest(request);
}

} // namespace command_internal

} // namespace rhbm_gem
