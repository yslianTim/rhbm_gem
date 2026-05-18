#include "detail/CommandBase.hpp"

#include "detail/MapSampling.hpp"

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimator.hpp>

#include <atomic>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem {

class PotentialAnalysisCommand final : public CommandBase<PotentialAnalysisRequest>
{
public:
    PotentialAnalysisCommand();

private:
    void NormalizeAndValidateRequest(PotentialAnalysisRequest & request) override;
    void ValidatePreparedRequest(const PotentialAnalysisRequest & request) override;
    bool ExecuteImpl(const PotentialAnalysisRequest & request) override;
};

namespace {

struct PotentialAnalysisInputs
{
    std::unique_ptr<ModelObject> model_object;
    std::unique_ptr<MapObject> map_object;
};

std::optional<PotentialAnalysisInputs> LoadPotentialAnalysisInputs(
    const PotentialAnalysisRequest & request)
{
    try
    {
        auto model_object{ ReadModel(request.model_file_path) };
        model_object->SetKeyTag("model");
        auto map_object{ ReadMap(request.map_file_path) };
        map_object->SetKeyTag("map");
        if (model_object == nullptr || map_object == nullptr)
        {
            Logger::Log(LogLevel::Error,
                "LoadPotentialAnalysisInputs : model/map object missing after load.");
            return std::nullopt;
        }
        if (request.simulation_flag)
        {
            if (request.simulated_map_resolution == 0.0)
            {
                Logger::Log(LogLevel::Warning,
                    "[Warning] The resolution of input simulated map hasn't been set.\n"
                    "          Please give the corresponding resolution value for this map.\n"
                    "          (-r, --sim-resolution)");
            }
            model_object->SetEmdID("Simulation");
            model_object->SetResolution(request.simulated_map_resolution);
            model_object->SetResolutionMethod("Blurring Width");
        }
        return PotentialAnalysisInputs{ std::move(model_object), std::move(map_object) };
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error, "LoadPotentialAnalysisInputs : " + std::string(e.what()));
        return std::nullopt;
    }
}

std::vector<LocalPotentialEditor> BuildSelectedAtomLocalEditors(ModelObject & model_object)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    std::vector<LocalPotentialEditor> local_editor_list;
    local_editor_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        local_editor_list.emplace_back(analysis.EnsureAtomLocalPotential(*atom));
    }
    return local_editor_list;
}

gaussian_estimator::TrainingOptions MakeGaussianEstimatorOptions(
    const PotentialAnalysisRequest & request)
{
    gaussian_estimator::TrainingOptions options;
    options.alpha_min = request.training_alpha_min;
    options.alpha_max = request.training_alpha_max;
    options.alpha_step = request.training_alpha_step;
    options.fit_range_min = request.fit_range_min;
    options.fit_range_max = request.fit_range_max;
    options.thread_size = request.job_count;
    options.output_progress = true;
    options.output_summary_log = true;
    options.study_plot_dir = request.training_report_dir;
    return options;
}

bool HasEnoughSamplesInFitRange(
    const LocalPotentialSampleList & sample_entries,
    double fit_range_min,
    double fit_range_max,
    std::size_t minimum_sample_count)
{
    std::size_t count{ 0 };
    for (const auto & sample : sample_entries)
    {
        if (sample.distance < fit_range_min || sample.distance > fit_range_max) continue;
        count++;
        if (count >= minimum_sample_count) return true;
    }
    return false;
}

void RunModelObjectPreprocessing(
    ModelObject & model_object, bool asymmetry_flag, double alpha_r, double alpha_g)
{
    auto analysis{ model_object.EditAnalysis() };
    analysis.Clear();
    model_object.SelectAllAtoms();
    model_object.ApplySymmetrySelection(asymmetry_flag);

    for (auto * atom : model_object.GetSelectedAtoms())
    {
        analysis.EnsureAtomLocalPotential(*atom);
    }

    // Establish the model-analysis preprocessing invariant for downstream steps:
    // selection is finalized, local entries exist, atom groups are materialized,
    // and selected atoms carry the initial alpha-r and alpha_g.
    analysis.RebuildAtomGroupsFromSelection();
    for (auto * atom : model_object.GetSelectedAtoms())
    {
        analysis.EnsureAtomLocalPotential(*atom).SetAlphaR(alpha_r);
    }
    const auto analysis_view{ model_object.GetAnalysisView() };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        for (const auto group_key : analysis_view.CollectAtomGroupKeys(class_key))
        {
            analysis.SetAtomGroupAlphaG(group_key, class_key, alpha_g);
        }
    }

    Logger::Log(LogLevel::Info,
        "Number of selected atom = " + std::to_string(model_object.GetSelectedAtomCount()));
    Logger::Log(LogLevel::Info, model_object.GetAnalysisView().DescribeAtomGrouping());
    if (model_object.GetNumberOfAtom() > 0 && model_object.GetSelectedAtomCount() == 0)
    {
        Logger::Log(LogLevel::Warning,
            "No atoms are selected after symmetry filtering. "
            "The input CIF may miss usable _entity/_struct_asym metadata. "
            "Try '--asymmetry true' to bypass symmetry filtering.");
    }
}

void RunLocalPotentialFitting(
    ModelObject & model_object,
    const gaussian_estimator::TrainingOptions & options)
{
    const auto selected_atom_size{ model_object.GetSelectedAtomCount() };
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildSelectedAtomLocalEditors(model_object) };
    std::atomic<size_t> atom_count{ 0 };
    Logger::Log(LogLevel::Info,
        "Run Local atom fitting for " + std::to_string(selected_atom_size) + " atoms.");

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        auto local_editor{ local_editor_list[i] };
        const auto local_view{ LocalPotentialView::RequireFor(*atom_list[i]) };
        const auto result{
            gaussian_estimator::EstimateLocalGaussian(
                local_view.GetSamplingEntries(),
                local_view.GetAlphaR(),
                options)
        };

        local_editor.SetGaussianResult(result);

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            Logger::ProgressPercent(atom_count, selected_atom_size);
        }
    }
}

void RunAtomAlphaTraining(ModelObject & model_object, const PotentialAnalysisRequest & request)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    auto alpha_options{ MakeGaussianEstimatorOptions(request) };
    const auto component_class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    const auto component_group_keys{ analysis_view.CollectAtomGroupKeys(component_class_key) };

    // Alpha_R Training
    size_t count{ 0 };
    for (const auto group_key : component_group_keys)
    {
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(group_key, component_class_key) };
        std::vector<LocalPotentialSampleList> selected_sample_entries_list;
        selected_sample_entries_list.reserve(group_atom_list.size());
        for (auto & atom : group_atom_list)
        {
            analysis.EnsureAtomLocalPotential(*atom);
            const auto local_view{ LocalPotentialView::RequireFor(*atom) };
            const auto & sample_entries{ local_view.GetSamplingEntries() };
            if (!HasEnoughSamplesInFitRange(
                sample_entries, request.fit_range_min, request.fit_range_max, 10)) continue;
            selected_sample_entries_list.emplace_back(sample_entries);
        }
        selected_sample_entries_list.shrink_to_fit();
        if (!selected_sample_entries_list.empty())
        {
            auto alpha_r_options{ alpha_options };
            alpha_r_options.output_progress = false;
            const auto alpha_r{
                gaussian_estimator::TrainAlphaR(selected_sample_entries_list, alpha_r_options)
            };
            for (auto * atom : group_atom_list)
            {
                analysis.EnsureAtomLocalPotential(*atom).SetAlphaR(alpha_r);
            }
        }
        Logger::ProgressPercent(++count, component_group_keys.size());
    }
    
    // Alpha_G Training
    RunLocalPotentialFitting(model_object, alpha_options);

    std::vector<std::vector<LocalPotentialSampleList>> sample_group_list;
    std::vector<std::vector<LocalGaussianResult>> member_result_list;
    sample_group_list.reserve(component_group_keys.size());
    member_result_list.reserve(component_group_keys.size());
    for (const auto group_key : component_group_keys)
    {
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(group_key, component_class_key) };
        if (group_atom_list.size() < 10) continue;
        if (group_atom_list.front()->IsMainChainAtom() == false) continue;

        std::vector<LocalPotentialSampleList> group_samples;
        std::vector<LocalGaussianResult> group_member_results;
        group_samples.reserve(group_atom_list.size());
        group_member_results.reserve(group_atom_list.size());
        for (auto * atom : group_atom_list)
        {
            analysis.EnsureAtomLocalPotential(*atom);
            const auto local_view{ LocalPotentialView::RequireFor(*atom) };
            group_samples.emplace_back(local_view.GetSamplingEntries());
            group_member_results.emplace_back(local_view.GetGaussianResult());
        }
        sample_group_list.emplace_back(std::move(group_samples));
        member_result_list.emplace_back(std::move(group_member_results));
    }

    const auto alpha_g{
        gaussian_estimator::TrainAlphaG(
            sample_group_list, member_result_list, alpha_options, !request.training_report_dir.empty())
    };

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        for (const auto group_key : analysis_view.CollectAtomGroupKeys(class_key))
        {
            analysis.SetAtomGroupAlphaG(group_key, class_key, alpha_g);
        }
    }
}

void RunAtomPotentialFittingWorkflow(ModelObject & model_object, const PotentialAnalysisRequest & request)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomPotentialFittingWorkflow");
    if (request.training_alpha_flag)
    {
        RunAtomAlphaTraining(model_object, request);
    }
    else
    {
        const auto options{ MakeGaussianEstimatorOptions(request) };
        RunLocalPotentialFitting(model_object, options);
    }

    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    for (auto * atom : model_object.GetSelectedAtoms())
    {
        analysis.EnsureAtomLocalPotential(*atom);
    }
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        Logger::Log(LogLevel::Info, "Class type: " + class_key);

        // Group Atom Potential Fitting
        auto group_key_list{ analysis_view.CollectAtomGroupKeys(class_key) };
        auto group_key_size{ group_key_list.size() };
        std::atomic<size_t> key_count{ 0 };

#ifdef USE_OPENMP
        #pragma omp parallel for num_threads(request.job_count)
#endif
        for (size_t k = 0; k < group_key_size; k++)
        {
            auto group_key{ group_key_list[k] };
            const auto & atom_list{ analysis_view.GetAtomObjectList(group_key, class_key) };
            std::vector<LocalPotentialSampleList> member_sample_list;
            std::vector<LocalGaussianResult> member_result_list;
            member_sample_list.reserve(atom_list.size());
            member_result_list.reserve(atom_list.size());
            for (const auto & atom : atom_list)
            {
                const auto local_view{ LocalPotentialView::RequireFor(*atom) };
                member_sample_list.emplace_back(local_view.GetSamplingEntries());
                member_result_list.emplace_back(local_view.GetGaussianResult());
            }
            const auto result{
                gaussian_estimator::EstimateGroupGaussian(
                    member_sample_list,
                    member_result_list,
                    analysis_view.GetAtomAlphaG(group_key, class_key),
                    MakeGaussianEstimatorOptions(request))
            };

#ifdef USE_OPENMP
            #pragma omp critical
#endif
            {
                analysis.ApplyAtomGroupGaussianResult(group_key, class_key, result);
                key_count++;
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

void RunPotentialSamplingWorkflow(
    MapObject & map_object,
    ModelObject & model_object,
    const PotentialAnalysisRequest & request)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunPotentialSamplingWorkflow");
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    const auto atom_size{ atom_list.size() };
    size_t atom_count{ 0 };
    auto local_editor_list{ BuildSelectedAtomLocalEditors(model_object) };
#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(request.job_count)
#endif
    for (size_t i = 0; i < atom_size; i++)
    {
        auto atom{ atom_list[i] };
        auto editor{ local_editor_list[i] };
        auto sampling_entries{
            SampleAtomMapValues(map_object, *atom, request.sampling_profile_choice)
        };
        editor.SetSamplingEntries(std::move(sampling_entries));

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            Logger::ProgressPercent(atom_count, atom_size);
        }
    }
}

} // namespace

PotentialAnalysisCommand::PotentialAnalysisCommand() : CommandBase<PotentialAnalysisRequest>{}
{
}

void PotentialAnalysisCommand::NormalizeAndValidateRequest(PotentialAnalysisRequest & request)
{
    RequireExistingPath(request, &PotentialAnalysisRequest::model_file_path);
    RequireExistingPath(request, &PotentialAnalysisRequest::map_file_path);
    RequireFiniteNonNegativeScalar(request, &PotentialAnalysisRequest::simulated_map_resolution);
    RequireNonEmptyList(request, &PotentialAnalysisRequest::saved_key_tag);
    RequireEnum(request, &PotentialAnalysisRequest::sampling_profile_choice);
    RequireFiniteNonNegativeScalar(request, &PotentialAnalysisRequest::fit_range_min);
    RequireFiniteNonNegativeScalar(request, &PotentialAnalysisRequest::fit_range_max);
    RequireFinitePositiveScalar(request, &PotentialAnalysisRequest::alpha_r);
    RequireFinitePositiveScalar(request, &PotentialAnalysisRequest::alpha_g);
    RequireFiniteNonNegativeScalar(request, &PotentialAnalysisRequest::training_alpha_min);
    RequireFiniteNonNegativeScalar(request, &PotentialAnalysisRequest::training_alpha_max);
    RequireFinitePositiveScalar(request, &PotentialAnalysisRequest::training_alpha_step);
}

bool PotentialAnalysisCommand::ExecuteImpl(const PotentialAnalysisRequest & request)
{
    auto inputs{ LoadPotentialAnalysisInputs(request) };
    if (!inputs.has_value()) return false;

    auto & model_object{ *inputs->model_object };
    auto & map_object{ *inputs->map_object };
    if (!request.simulation_flag && request.map_normalization_flag)
    {
        map_object.MapValueArrayNormalization();
    }
    RunModelObjectPreprocessing(model_object, request.asymmetry_flag, request.alpha_r, request.alpha_g);
    RunPotentialSamplingWorkflow(map_object, model_object, request);
    RunAtomPotentialFittingWorkflow(model_object, request);

    DataRepository repository{ request.database_path };
    repository.SaveModel(model_object, request.saved_key_tag);
    model_object.EditAnalysis().ClearTransientFitStates();
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
    RequirePrepareCondition(
        request.training_alpha_min <= request.training_alpha_max,
        "Expected --training-alpha-min <= --training-alpha-max.");
}

namespace command_internal {

CommandResult ExecutePotentialAnalysisCommand(const PotentialAnalysisRequest & request)
{
    PotentialAnalysisCommand command;
    return command.ExecuteRequest(request);
}

} // namespace command_internal

} // namespace rhbm_gem
