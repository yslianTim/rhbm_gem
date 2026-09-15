#include "detail/CommandRunner.hpp"
#include "detail/MapSimulation.hpp"
#include "detail/SimulationManifest.hpp"

#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace rhbm_gem::core {

namespace {

void NormalizeAndValidateRequest(
    CommandRunner<MapSimulationRequest> & runner,
    MapSimulationRequest & request)
{
    runner.RequireExistingPath(request, &MapSimulationRequest::model_file_path);
    runner.RequireEnum(request, &MapSimulationRequest::potential_model_choice);
    runner.RequireEnum(request, &MapSimulationRequest::partial_charge_choice);
    runner.NormalizeFinitePositiveScalar(request, &MapSimulationRequest::cutoff_distance, 5.0);
    runner.NormalizeFinitePositiveScalar(request, &MapSimulationRequest::grid_spacing, 0.5);

    std::vector<double> filtered_widths;
    filtered_widths.reserve(request.blurring_width_list.size());
    for (const auto width : request.blurring_width_list)
    {
        if (!numeric_validation::IsFinitePositive(width))
        {
            runner.AddFieldNormalizationWarning(&MapSimulationRequest::blurring_width_list,
                "Blurring width must be a finite positive value, dropping current setting: "
                    + std::to_string(width));
            continue;
        }
        filtered_widths.push_back(width);
    }
    request.blurring_width_list = std::move(filtered_widths);
}

bool ExecutePreparedRequest(const MapSimulationRequest & request)
{
    const auto model_sha256{ simulation::FileSha256(request.model_file_path) };
    std::unique_ptr<ModelObject> model_object;
    try
    {
        model_object = ReadModel(request.model_file_path);
        model_object->SetKeyTag("model");
        if (simulation::FileSha256(request.model_file_path) != model_sha256)
            throw std::runtime_error("Model file changed while preparing simulation input.");
    }
    catch(const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "MapSimulationCommand::BuildDataObject : Failed to load model file '"
                + request.model_file_path.string() + "' as ModelObject: " + std::string(e.what()));
        return false;
    }

    model_object->SelectAllAtoms();
    model_object->ApplyElementSelection(Element::HYDROGEN, request.exclude_hydrogen);
    model_object->ApplyBackboneSelection(request.only_backbone);
    //model_object->ApplyComponentIDSelection("HOH", true);
    const auto atom_list{ simulation::PrepareSimulationAtomList(*model_object, request) };
    const simulation::SimulationSource source{ model_object->GetPdbID(), model_sha256 };
    Logger::Log(LogLevel::Info,
        "Total number of blurring width sets to be simulated: "
        + std::to_string(request.blurring_width_list.size()));

    std::vector<std::filesystem::path> outputs;
    std::set<std::filesystem::path> unique_outputs;
    for (const auto blurring_width : request.blurring_width_list)
    {
        auto map_key_tag{
            model_object->GetPdbID() + "_bw" +
            string_helper::ToStringWithPrecision<double>(blurring_width, 2)
        };
        const auto output{ request.output_dir / (request.map_file_name + "_" + map_key_tag + ".map") };
        if (!unique_outputs.insert(output).second)
            throw std::runtime_error("Blurring widths produce the same simulation output filename: " + output.string());
        outputs.push_back(output);
    }
    auto map_object{ simulation::CreateMapObject(request, atom_list) };
    for (size_t i = 0; i < request.blurring_width_list.size(); ++i)
    {
        const auto blurring_width{ request.blurring_width_list[i] };
        const auto actual_job_count{ simulation::PopulateMapValueArray(*map_object, atom_list, request, blurring_width) };
        simulation::WriteSimulationArtifacts(outputs[i], *map_object, atom_list, request,
            blurring_width, actual_job_count, source);
    }
    return true;
}

void ValidatePreparedRequest(
    CommandRunner<MapSimulationRequest> & runner,
    const MapSimulationRequest & request)
{
    runner.RequirePrepareCondition(
        !request.blurring_width_list.empty(),
        "At least one positive blurring width is required.");
}

} // namespace

namespace command_internal {

CommandResult ExecuteMapSimulationCommand(const MapSimulationRequest & request)
{
    try
    {
        return CommandRunner<MapSimulationRequest>{}.Run(
            request,
            NormalizeAndValidateRequest,
            ValidatePreparedRequest,
            ExecutePreparedRequest);
    }
    catch (const std::exception & error)
    {
        Logger::Log(LogLevel::Error, "Map simulation failed: " + std::string(error.what()));
        return { false, { { "request", error.what() } } };
    }
}

} // namespace command_internal

} // namespace rhbm_gem::core
