#include "detail/CommandBase.hpp"

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <memory>
#include <string>

namespace rhbm_gem::core {

class PotentialAnalysisCommand final : public CommandBase<PotentialAnalysisRequest>
{
public:
    PotentialAnalysisCommand();

private:
    void NormalizeAndValidateRequest(PotentialAnalysisRequest & request) override;
    void ValidatePreparedRequest(const PotentialAnalysisRequest & request) override;
    bool ExecuteImpl(const PotentialAnalysisRequest & request) override;
};

PotentialAnalysisCommand::PotentialAnalysisCommand() : CommandBase<PotentialAnalysisRequest>{}
{
}

void PotentialAnalysisCommand::NormalizeAndValidateRequest(PotentialAnalysisRequest & request)
{
    RequireExistingPath(request, &PotentialAnalysisRequest::model_file_path);
    RequireExistingPath(request, &PotentialAnalysisRequest::map_file_path);
    RequireFiniteNonNegativeScalar(request, &PotentialAnalysisRequest::simulated_map_resolution);
    RequireNonEmptyList(request, &PotentialAnalysisRequest::saved_key_tag);
    RequireEnum(request, &PotentialAnalysisRequest::sampling_method);
    RequireFiniteNonNegativeScalar(request, &PotentialAnalysisRequest::fit_range_min);
    RequireFiniteNonNegativeScalar(request, &PotentialAnalysisRequest::fit_range_max);
}

bool PotentialAnalysisCommand::ExecuteImpl(const PotentialAnalysisRequest & request)
{
    std::unique_ptr<ModelObject> model_object;
    std::unique_ptr<MapObject> map_object;
    try
    {
        model_object = ReadModel(request.model_file_path);
        map_object = ReadMap(request.map_file_path);
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error, "PotentialAnalysisCommand : " + std::string(e.what()));
        return false;
    }
    if (model_object == nullptr || map_object == nullptr)
    {
        Logger::Log(LogLevel::Error,
            "PotentialAnalysisCommand : model/map object missing after load.");
        return false;
    }
    if (request.simulation_flag)
    {
        model_object->ApplySimulationMetadata(request.simulated_map_resolution);
    }
    if (!request.simulation_flag && request.map_normalization_flag)
    {
        map_object->MapValueArrayNormalization();
    }

    model_object->SelectAllAtoms();
    model_object->ApplySymmetrySelection(request.asymmetry_flag);
    if (request.exclude_hydrogen) model_object->ApplyElementExclusion(Element::HYDROGEN);
    model_object->LocalPotentialInitialization();
    Logger::Log(LogLevel::Info, model_object->GetAnalysisView().GetAtomCountingSummary());
    Logger::Log(LogLevel::Info, model_object->GetAnalysisView().GetAtomGroupingSummary());
    RunPotentialSamplingWorkflow(*map_object, *model_object, request.sampling_method, request.job_count);

    FitOptions options;
    options.distance_min = request.fit_range_min;
    options.distance_max = request.fit_range_max;
    options.thread_size = request.job_count;
    RunLocalAlphaTraining(*model_object, options);
    RunLocalPotentialFitting(*model_object, options);
    RunGroupAlphaTraining(*model_object, options);
    RunGroupPotentialFitting(*model_object, options);

    DataRepository repository{ request.database_path };
    repository.SaveModel(*model_object, request.saved_key_tag);
    model_object->ClearTransientFitStates();
    return true;
}

void PotentialAnalysisCommand::ValidatePreparedRequest(const PotentialAnalysisRequest & request)
{
    RequirePrepareCondition(
        !request.simulation_flag || request.simulated_map_resolution > 0.0,
        "Expected a positive simulated-map resolution when '--simulation true' is selected.");
    RequirePrepareCondition(
        request.fit_range_min <= request.fit_range_max,
        "Expected --fit-min <= --fit-max.");
}

namespace command_internal {

CommandResult ExecutePotentialAnalysisCommand(const PotentialAnalysisRequest & request)
{
    PotentialAnalysisCommand command;
    return command.ExecuteRequest(request);
}

} // namespace command_internal

} // namespace rhbm_gem::core
