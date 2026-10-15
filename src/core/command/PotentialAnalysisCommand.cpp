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

#include <algorithm>
#include <map>
#include <cctype>
#include <memory>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace rhbm_gem::core {

namespace {

struct PreparedPotentialAnalysisInputs
{
    std::unique_ptr<ModelObject> model_object;
    JointAnalysisMetadata metadata;
};

std::string TrimChainId(const std::string & chain_id)
{
    const auto first{ std::find_if_not(chain_id.begin(), chain_id.end(), [](unsigned char character)
    {
        return std::isspace(character) != 0;
    }) };
    if (first == chain_id.end())
    {
        return {};
    }
    const auto last{ std::find_if_not(chain_id.rbegin(), chain_id.rend(), [](unsigned char character)
    {
        return std::isspace(character) != 0;
    }).base() };
    return std::string(first, last);
}

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
    runner.RequireEnum(request, &PotentialAnalysisRequest::estimator);
    std::vector<std::string> selected_chain_ids;
    std::unordered_set<std::string> unique_chain_ids;
    bool has_empty_chain_id{};
    for (const auto & chain_id : request.selected_chain_ids)
    {
        auto normalized_chain_id{ TrimChainId(chain_id) };
        if (normalized_chain_id.empty())
        {
            has_empty_chain_id = true;
        }
        else if (unique_chain_ids.insert(normalized_chain_id).second)
        {
            selected_chain_ids.emplace_back(std::move(normalized_chain_id));
        }
    }
    request.selected_chain_ids = std::move(selected_chain_ids);
    if (has_empty_chain_id)
    {
        runner.AddFieldValidationError(&PotentialAnalysisRequest::selected_chain_ids,
            "--select-chain cannot contain an empty chain ID.");
    }
    if (request.estimator==PotentialEstimator::JOINT_COMPONENTS)
    {
        if (request.job_count>1) runner.AddFieldNormalizationWarning(&PotentialAnalysisRequest::job_count,"Joint initialization and fitting use one worker.");
    }
}

bool ExecutePreparedRequest(
    const PotentialAnalysisRequest & request,
    PreparedPotentialAnalysisInputs inputs,
    CommandOutcome & outcome)
{
    JointAnalysisMetadata metadata{ std::move(inputs.metadata) };
    const bool joint=request.estimator==PotentialEstimator::JOINT_COMPONENTS;
    auto model_object{ std::move(inputs.model_object) };
    std::unique_ptr<MapObject> map_object;
    try
    {
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
    JointMapNormalization normalization{!request.simulation_flag,false,1};
    if (!request.simulation_flag)
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

    Logger::Log(LogLevel::Info, BuildAtomCountingSummary(*model_object));

    FitOptions options;
    options.thread_size = joint ? 1 : request.job_count;
    options.estimator = request.estimator;
    options.sampling_method = SphereSamplingMethod::FibonacciDeterministic;
    options.exclude_hydrogen = true;
    options.enable_second_stage_failed_only_refinement = true;
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
        metadata.solver=fit.metadata.solver;
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
        if (!fit.initialization.valid)
        {
            outcome = CommandOutcome::JointInitializationFailed;
            Logger::Log(LogLevel::Warning,"Joint initialization failed; diagnostic result was saved (CLI exit 2).");
        }
        else if (!fit.search_completed || fit.runtime_convergence != JointCheckStatus::Passed ||
            fit.target_runtime_convergence != JointCheckStatus::Passed)
        {
            outcome = CommandOutcome::JointNotConverged;
            Logger::Log(LogLevel::Warning,"Joint fitting did not converge; diagnostic result was saved (CLI exit 3).");
        }
        else
        {
            outcome = CommandOutcome::JointQualified;
            Logger::Log(LogLevel::Info,"Joint fitting qualified.");
        }
    }
    model_object->EditAnalysis().ClearTransientFitStates();
    return true;
}

void ValidatePreparedRequest(
    CommandRunner<PotentialAnalysisRequest> & runner,
    const PotentialAnalysisRequest & request,
    PreparedPotentialAnalysisInputs & inputs)
{
    runner.RequirePrepareCondition(
        !request.simulation_flag || request.simulated_map_resolution > 0.0,
        "Expected a positive simulated-map resolution when '--simulation true' is selected.");
    if (request.simulation_flag && request.simulated_map_resolution <= 0.0)
    {
        return;
    }

    const bool joint{ request.estimator == PotentialEstimator::JOINT_COMPONENTS };
    try
    {
        if (joint)
        {
            inputs.metadata.model_sha256 = FileSha256(request.model_file_path);
            inputs.metadata.map_sha256 = FileSha256(request.map_file_path);
        }
        inputs.model_object = ReadModel(request.model_file_path);
    }
    catch (const std::exception & error)
    {
        runner.AddFieldValidationError(&PotentialAnalysisRequest::model_file_path,
            "Unable to load model for target validation: " + std::string(error.what()));
        return;
    }
    if (inputs.model_object == nullptr)
    {
        runner.AddFieldValidationError(&PotentialAnalysisRequest::model_file_path,
            "Unable to load model for target validation.");
        return;
    }

    std::unordered_set<std::string> available_chain_ids;
    for (const auto & atom : inputs.model_object->GetAtomList())
    {
        available_chain_ids.insert(atom->GetChainID());
    }
    bool has_unknown_chain{};
    for (const auto & chain_id : request.selected_chain_ids)
    {
        if (available_chain_ids.contains(chain_id))
        {
            continue;
        }
        has_unknown_chain = true;
        runner.AddFieldValidationError(&PotentialAnalysisRequest::selected_chain_ids,
            "Chain ID '" + chain_id + "' specified by --select-chain does not exist in the model.");
    }
    if (has_unknown_chain)
    {
        return;
    }

    const std::unordered_set<std::string> selected_chain_ids{
        request.selected_chain_ids.begin(), request.selected_chain_ids.end()
    };
    inputs.model_object->SelectAtoms([&selected_chain_ids](const AtomObject & atom)
    {
        return atom.GetElement() != Element::HYDROGEN
            && (selected_chain_ids.empty() || selected_chain_ids.contains(atom.GetChainID()));
    });
    if (inputs.model_object->GetSelectedAtomCount() == 0)
    {
        runner.AddFieldValidationError(&PotentialAnalysisRequest::selected_chain_ids,
            request.selected_chain_ids.empty()
                ? "The model contains no non-hydrogen target atoms."
                : "--select-chain selected no non-hydrogen target atoms.");
    }
}

} // namespace

namespace command_internal {

CommandResult ExecutePotentialAnalysisCommand(const PotentialAnalysisRequest & request)
{
    PreparedPotentialAnalysisInputs inputs;
    CommandOutcome outcome{CommandOutcome::Completed};
    auto result = CommandRunner<PotentialAnalysisRequest>{}.Run(
        request,
        NormalizeAndValidateRequest,
        [&inputs](auto & runner, const auto & prepared)
        {
            ValidatePreparedRequest(runner, prepared, inputs);
        },
        [&inputs,&outcome](const PotentialAnalysisRequest & prepared)
        {
            try { return ExecutePreparedRequest(prepared, std::move(inputs), outcome); }
            catch (const std::exception & error)
            {
                Logger::Log(LogLevel::Error,"PotentialAnalysisCommand : "+std::string(error.what()));
                return false;
            }
        });
    result.outcome = outcome;
    return result;
}

} // namespace command_internal

} // namespace rhbm_gem::core
