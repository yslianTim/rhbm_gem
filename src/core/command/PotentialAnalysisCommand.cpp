#include "detail/CommandRunner.hpp"

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <map>
#include <memory>
#include <string>
#include <unordered_map>

namespace rhbm_gem::core {

namespace {

std::string BuildAtomCountingSummary(const ModelObject & model_object)
{
    std::map<Element, std::size_t> element_counts;
    for (const auto * atom : model_object.GetSelectedAtoms())
    {
        element_counts[atom->GetElement()]++;
    }

    std::string description{
        "Number of selected atom = " + std::to_string(model_object.GetSelectedAtomCount())
    };
    for (const auto & [element, count] : element_counts)
    {
        description +=
            "\n - Element type: " + ChemicalDataHelper::GetLabel(element) + " include "
            + std::to_string(count) + " atoms.";
    }
    return description;
}

std::string BuildAtomGroupingSummary(const ModelObject & model_object)
{
    return "Atomic model includes "
        + std::to_string(model_object.GetAnalysisView()
            .CollectAtomGroupKeys().size())
        + " atom groups.";
}

void NormalizeAndValidateRequest(
    CommandRunner<PotentialAnalysisRequest> & runner,
    PotentialAnalysisRequest & request)
{
    runner.RequireExistingPath(request, &PotentialAnalysisRequest::model_file_path);
    runner.RequireExistingPath(request, &PotentialAnalysisRequest::map_file_path);
    runner.RequireFiniteNonNegativeScalar(
        request, &PotentialAnalysisRequest::simulated_map_resolution);
    runner.RequireNonEmptyList(request, &PotentialAnalysisRequest::saved_key_tag);
    runner.RequireEnum(request, &PotentialAnalysisRequest::sampling_method);
    runner.RequireFiniteNonNegativeScalar(request, &PotentialAnalysisRequest::fit_range_min);
    runner.RequireFiniteNonNegativeScalar(request, &PotentialAnalysisRequest::fit_range_max);
}

bool ExecutePreparedRequest(const PotentialAnalysisRequest & request)
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

    try
    {
        const auto [reference_height, reference_offset]{
            GetReferenceGaussianParameters(*map_object)
        };
        std::unordered_map<int, double> q_scores_by_serial_id;
        const auto standard_average_qscore{
            CalculateAverageQScores(
                *map_object, *model_object,
                reference_height, reference_offset, q_scores_by_serial_id)
        };
        for (const auto & atom : model_object->GetAtomList())
        {
            atom->SetStandardQScore(atom->GetElement() == Element::HYDROGEN ?
                0.0 : q_scores_by_serial_id.at(atom->GetSerialID()));
        }
        model_object->SetStandardAverageQScore(standard_average_qscore);
        model_object->SetReferenceHeight(reference_height);
        model_object->SetReferenceOffset(reference_offset);
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "PotentialAnalysisCommand : reference Gaussian/Q-score calculation failed: "
                + std::string(e.what()));
        return false;
    }

    model_object->SelectAllAtoms();
    model_object->ApplySymmetrySelection(request.asymmetry_flag);
    model_object->ApplyElementSelection(Element::HYDROGEN, request.exclude_hydrogen);
    model_object->ApplyBackboneSelection(request.only_backbone);
    model_object->EditAnalysis().InitializeFromSelection();
    Logger::Log(LogLevel::Info, BuildAtomCountingSummary(*model_object));
    Logger::Log(LogLevel::Info, BuildAtomGroupingSummary(*model_object));
    RunPotentialSamplingWorkflow(*map_object, *model_object, request.sampling_method, request.job_count);

    FitOptions options;
    options.distance_min = request.fit_range_min;
    options.distance_max = request.fit_range_max;
    options.thread_size = request.job_count;
    options.exclude_hydrogen = request.exclude_hydrogen;
    options.result_csv_path = request.output_dir / ("local_fitting_result_" + path_helper::EnsureSanitizedTag(request.saved_key_tag) + ".csv");
    try
    {
        RunPotentialFittingWorkflow(*model_object, options);
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "PotentialAnalysisCommand : potential fitting failed: " + std::string(e.what()));
        return false;
    }

    DataRepository repository{ request.database_path };
    repository.SaveModel(*model_object, request.saved_key_tag);
    model_object->EditAnalysis().ClearTransientFitStates();
    return true;
}

void ValidatePreparedRequest(
    CommandRunner<PotentialAnalysisRequest> & runner,
    const PotentialAnalysisRequest & request)
{
    runner.RequirePrepareCondition(
        !request.simulation_flag || request.simulated_map_resolution > 0.0,
        "Expected a positive simulated-map resolution when '--simulation true' is selected.");
    runner.RequirePrepareCondition(
        request.fit_range_min <= request.fit_range_max,
        "Expected --fit-min <= --fit-max.");
}

} // namespace

namespace command_internal {

CommandResult ExecutePotentialAnalysisCommand(const PotentialAnalysisRequest & request)
{
    return CommandRunner<PotentialAnalysisRequest>{}.Run(
        request,
        NormalizeAndValidateRequest,
        ValidatePreparedRequest,
        ExecutePreparedRequest);
}

} // namespace command_internal

} // namespace rhbm_gem::core
