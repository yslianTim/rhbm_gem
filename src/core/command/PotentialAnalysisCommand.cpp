#include "detail/CommandRunner.hpp"
#include "utils/domain/FileFingerprint.hpp"

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include "data/io/detail/JointResultJson.hpp"
#include <rhbm_gem/core/MapSampler.hpp>
#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
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
    runner.RequireEnum(request, &PotentialAnalysisRequest::estimator);
    if (request.estimator==PotentialEstimator::JOINT_COMPONENTS)
    {
        if (request.sampling_method!=SphereSamplingMethod::FibonacciDeterministic)
            runner.AddFieldValidationError(&PotentialAnalysisRequest::sampling_method,"Joint initialization requires Fibonacci sampling.");
        if (request.job_count>1) runner.AddFieldNormalizationWarning(&PotentialAnalysisRequest::job_count,"Joint initialization and fitting use one worker.");
        request.exclude_hydrogen=true;
    }
}

bool ExecutePreparedRequest(const PotentialAnalysisRequest & request)
{
    JointAnalysisMetadata metadata;
    const bool joint=request.estimator==PotentialEstimator::JOINT_COMPONENTS;
    std::unique_ptr<ModelObject> model_object;
    std::unique_ptr<MapObject> map_object;
    try
    {
        if (joint)
        {
            metadata.model_sha256=FileSha256(request.model_file_path);
            metadata.map_sha256=FileSha256(request.map_file_path);
        }
        model_object = ReadModel(request.model_file_path);
        map_object = ReadMap(request.map_file_path);
        if (joint && (FileSha256(request.model_file_path)!=*metadata.model_sha256 ||
                      FileSha256(request.map_file_path)!=*metadata.map_sha256))
            throw std::runtime_error("Model/map file changed while preparing joint input.");
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
    JointMapNormalization normalization{request.map_normalization_flag,false,1};
    if (!request.simulation_flag && request.map_normalization_flag)
    {
        const auto divisor=map_object->MapValueArrayNormalization();
        normalization.applied=divisor.has_value();
        normalization.divisor=divisor.value_or(1);
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

    Logger::Log(LogLevel::Info, BuildAtomCountingSummary(*model_object));

    FitOptions options;
    options.thread_size = joint ? 1 : request.job_count;
    options.estimator = request.estimator;
    options.sampling_method = request.sampling_method;
    options.exclude_hydrogen = request.exclude_hydrogen;
    options.enable_second_stage_failed_only_refinement = request.enable_second_stage_failed_only_refinement;
    try
    {
        RunPotentialFittingWorkflow(*map_object, *model_object, options);
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "PotentialAnalysisCommand : potential fitting failed: " + std::string(e.what()));
        return false;
    }

    Logger::Log(LogLevel::Info, BuildAtomGroupingSummary(*model_object));
    if (joint)
    {
        const auto & fit=*model_object->GetAnalysisView().GetJointResult();
        metadata.model_path=request.model_file_path.string(); metadata.map_path=request.map_file_path.string();
        metadata.grid_size=map_object->GetGridSize(); metadata.grid_spacing=map_object->GetGridSpacing();
        metadata.origin=map_object->GetOrigin(); metadata.simulation=request.simulation_flag;
        metadata.map_normalization=normalization;
        metadata.software=fit.metadata.software;
        model_object->EditAnalysis().UpdateJointMetadata(std::move(metadata));

    }

    DataRepository repository{ request.database_path };
    repository.SaveModel(*model_object, request.saved_key_tag);
    if (joint)
    {
        const auto & fit=*model_object->GetAnalysisView().GetJointResult();
        std::size_t available=0;
        for (const auto & component:fit.components)
        {
            available+=component.state.has_value();
            Logger::Log(LogLevel::Info,"Joint component "+component.id+": stop="+component.stop_reason+
                ", runtime_convergence="+std::string(joint_result_io::StatusText(component.runtime_convergence))+
                ", target_runtime_convergence="+std::string(joint_result_io::StatusText(component.target_runtime_convergence)));
        }
        const auto targets=fit.selection_domain->target_indices.size();
        Logger::Log(LogLevel::Info,"Joint targets="+std::to_string(targets)+", halo="+
            std::to_string(fit.atom_ids.size()-targets));
        Logger::Log(LogLevel::Info,"Joint result saved: search_completed="+std::to_string(fit.search_completed)+
            ", available_components="+std::to_string(available)+"/"+std::to_string(fit.components.size())+
            ", runtime_convergence="+std::string(joint_result_io::StatusText(fit.runtime_convergence))+
            ", target_runtime_convergence="+std::string(joint_result_io::StatusText(fit.target_runtime_convergence))+
            ", regular_certificate="+std::string(joint_result_io::StatusText(fit.regular_certificate)));
        if (!fit.initialization.valid) Logger::Log(LogLevel::Warning,"Joint initialization unavailable: "+fit.initialization.reason);
    }
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
}

} // namespace

namespace command_internal {

CommandResult ExecutePotentialAnalysisCommand(const PotentialAnalysisRequest & request)
{
    return CommandRunner<PotentialAnalysisRequest>{}.Run(
        request,
        NormalizeAndValidateRequest,
        ValidatePreparedRequest,
        [](const PotentialAnalysisRequest & prepared)
        {
            try { return ExecutePreparedRequest(prepared); }
            catch (const std::exception & error)
            {
                Logger::Log(LogLevel::Error,"PotentialAnalysisCommand : "+std::string(error.what()));
                return false;
            }
        });
}

} // namespace command_internal

} // namespace rhbm_gem::core
