#include "detail/CommandBase.hpp"

#include "detail/MapSampling.hpp"

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/LocalPainter.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

#include <atomic>
#include <cmath>
#include <filesystem>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <Eigen/Dense>

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

RHBMExecutionOptions MakePotentialAnalysisExecutionOptions(int thread_size, bool quiet_mode)
{
    RHBMExecutionOptions execution_options;
    execution_options.quiet_mode = quiet_mode;
    execution_options.thread_size = thread_size;
    return execution_options;
}

bool EmitTrainingReportIfRequested(
    const Eigen::MatrixXd & gaus_bias_matrix,
    const std::vector<double> & alpha_list,
    std::string_view x_label,
    std::string_view y_label,
    const std::filesystem::path & training_report_dir,
    std::string_view report_file_name)
{
    if (training_report_dir.empty()) return false;

    const auto report_path{ training_report_dir / std::filesystem::path{ report_file_name } };
    if (gaus_bias_matrix.rows() < 2 ||
        gaus_bias_matrix.cols() != static_cast<Eigen::Index>(alpha_list.size()))
    {
        Logger::Log(LogLevel::Warning,
            "Skip training report '" + report_path.string() +
                "' because the alpha bias matrix shape does not match the alpha grid.");
        return false;
    }

    std::vector<double> amplitude_values;
    std::vector<double> width_values;
    amplitude_values.reserve(alpha_list.size());
    width_values.reserve(alpha_list.size());
    for (Eigen::Index column = 0; column < gaus_bias_matrix.cols(); ++column)
    {
        amplitude_values.emplace_back(gaus_bias_matrix(0, column));
        width_values.emplace_back(gaus_bias_matrix(1, column));
    }

    LinePlotRequest report_request;
    report_request.output_path = report_path;
    report_request.x_axis.title = std::string(x_label);
    report_request.shared_y_axis_title = std::string(y_label);
    report_request.panels = {
        LinePlotPanel{
            "Amplitude #font[2]{A}",
            AxisSpec{},
            { LineSeries{ "Amplitude", alpha_list, amplitude_values, std::nullopt } }
        },
        LinePlotPanel{
            "Width #tau",
            AxisSpec{},
            { LineSeries{ "Width", alpha_list, width_values, std::nullopt } }
        }
    };

    const auto plot_result{ local_painter::SaveLinePlot(report_request) };
    if (!plot_result.Succeeded())
    {
        Logger::Log(LogLevel::Warning,
            "Failed to emit training report '" + report_path.string() + "': " + plot_result.message);
        return false;
    }

    if (!std::filesystem::exists(report_path))
    {
        Logger::Log(LogLevel::Warning,
            "Training report output was requested but no file was produced: " + report_path.string());
        return false;
    }

    return true;
}

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

std::vector<MutableLocalPotentialView> BuildSelectedAtomLocalEntryViews(ModelObject & model_object)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    std::vector<MutableLocalPotentialView> local_entry_list;
    local_entry_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        local_entry_list.emplace_back(analysis.EnsureAtomLocalPotential(*atom));
    }
    return local_entry_list;
}

void RunModelObjectPreprocessing(ModelObject & model_object, bool asymmetry_flag, double alpha_r)
{
    auto analysis{ model_object.EditAnalysis() };
    analysis.Clear();
    model_object.SelectAllAtoms();
    model_object.SelectAllBonds();
    model_object.ApplySymmetrySelection(asymmetry_flag);

    for (auto * atom : model_object.GetSelectedAtoms())
    {
        analysis.EnsureAtomLocalPotential(*atom);
    }
    for (auto * bond : model_object.GetSelectedBonds())
    {
        analysis.EnsureBondLocalPotential(*bond);
    }

    // Establish the model-analysis preprocessing invariant for downstream steps:
    // selection is finalized, local entries exist, atom groups are materialized,
    // and selected atoms carry the initial alpha-r.
    analysis.RebuildAtomGroupsFromSelection();
    for (auto * atom : model_object.GetSelectedAtoms())
    {
        analysis.EnsureAtomLocalPotential(*atom).SetAlphaR(alpha_r);
    }

    Logger::Log(LogLevel::Info,
        "Number of selected atom = " + std::to_string(model_object.GetSelectedAtomCount()));
    Logger::Log(LogLevel::Info,
        "Number of selected bond = " + std::to_string(model_object.GetSelectedBondCount()));
    Logger::Log(LogLevel::Info, model_object.GetAnalysisView().DescribeAtomGrouping());
    if (model_object.GetNumberOfAtom() > 0 && model_object.GetSelectedAtomCount() == 0)
    {
        Logger::Log(LogLevel::Warning,
            "No atoms are selected after symmetry filtering. "
            "The input CIF may miss usable _entity/_struct_asym metadata. "
            "Try '--asymmetry true' to bypass symmetry filtering.");
    }
}

void RunLocalPotentialFitting(ModelObject & model_object, int thread_size)
{
    const auto selected_atom_size{ model_object.GetSelectedAtomCount() };
    auto local_entry_list{ BuildSelectedAtomLocalEntryViews(model_object) };
    std::atomic<size_t> atom_count{ 0 };
    Logger::Log(LogLevel::Info,
        "Run Local atom fitting for " + std::to_string(selected_atom_size) + " atoms.");

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        auto local_entry{ local_entry_list[i] };
        const auto result{
            rhbm_helper::EstimateBetaMDPDE(
                local_entry.GetAlphaR(),
                local_entry.GetDataset(),
                MakePotentialAnalysisExecutionOptions(thread_size, true))
        };

        local_entry.SetFitResult(result);

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            Logger::ProgressPercent(atom_count, selected_atom_size);
        }
    }
}

void RunAtomAlphaTraining(
    ModelObject & model_object,
    const PotentialAnalysisRequest & request,
    int thread_size)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    rhbm_trainer::AlphaTrainer alpha_trainer(
        request.training_alpha_min, request.training_alpha_max, request.training_alpha_step);
    const auto & alpha_training_list{ alpha_trainer.AlphaGrid() };
    const auto component_class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    const auto component_group_keys{ analysis_view.CollectAtomGroupKeys(component_class_key) };
    Logger::Log(LogLevel::Info, alpha_trainer.GetAlphaGridSummary().str());

    // Alpha_R Training
    size_t count{ 0 };
    for (const auto group_key : component_group_keys)
    {
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(group_key, component_class_key) };
        std::vector<RHBMMemberDataset> selected_atom_dataset_list;
        selected_atom_dataset_list.reserve(group_atom_list.size());
        for (auto & atom : group_atom_list)
        {
            const auto local_entry{ analysis.EnsureAtomLocalPotential(*atom) };
            if (!local_entry.HasDataset()) continue;
            if (local_entry.GetDataset().y.rows() < 10) continue;
            selected_atom_dataset_list.emplace_back(local_entry.GetDataset());
        }
        selected_atom_dataset_list.shrink_to_fit();
        if (!selected_atom_dataset_list.empty())
        {
            rhbm_trainer::AlphaTrainer::AlphaTrainingOptions alpha_r_options;
            alpha_r_options.subset_size = 5;
            alpha_r_options.execution_options = MakePotentialAnalysisExecutionOptions(thread_size, true);
            const auto alpha_r_result{
                alpha_trainer.TrainAlphaR(selected_atom_dataset_list, alpha_r_options)
            };
            for (auto * atom : group_atom_list)
            {
                analysis.EnsureAtomLocalPotential(*atom).SetAlphaR(alpha_r_result.best_alpha);
            }
        }
        Logger::ProgressPercent(++count, component_group_keys.size());
    }
    
    // Alpha_G Training
    RunLocalPotentialFitting(model_object, thread_size);

    std::vector<std::vector<RHBMParameterVector>> beta_group_list;
    beta_group_list.reserve(component_group_keys.size());
    for (const auto group_key : component_group_keys)
    {
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(group_key, component_class_key) };
        if (group_atom_list.size() < 10) continue;
        if (group_atom_list.front()->IsMainChainAtom() == false) continue;

        std::vector<RHBMParameterVector> beta_list;
        beta_list.reserve(group_atom_list.size());
        for (auto * atom : group_atom_list)
        {
            beta_list.emplace_back(
                analysis.EnsureAtomLocalPotential(*atom).GetFitResult().beta_mdpde
            );
        }
        beta_group_list.emplace_back(std::move(beta_list));
    }

    auto selected_group_size{ beta_group_list.size() };
    Logger::Log(LogLevel::Info,
        "Run Alpha_G Training with "+ std::to_string(selected_group_size) +" groups.");
    auto alpha_g{ alpha_training_list.front() };
    if (!beta_group_list.empty())
    {
        rhbm_trainer::AlphaTrainer::AlphaTrainingOptions alpha_g_options;
        alpha_g_options.subset_size = 10;
        alpha_g_options.execution_options = MakePotentialAnalysisExecutionOptions(thread_size, true);
        alpha_g_options.progress_callback =
            [](std::size_t completed, std::size_t total)
            {
                Logger::ProgressPercent(completed, total);
            };
        const auto alpha_g_result{ alpha_trainer.TrainAlphaG(beta_group_list, alpha_g_options) };
        alpha_g = alpha_g_result.best_alpha;
    }
    else
    {
        Logger::Log(LogLevel::Warning,
            "Skip Alpha_G Training because no eligible groups were selected.");
    }
    Logger::Log(LogLevel::Info,
        "Alpha_G Training Results Summary: minimum mu error sum alpha_g = "
        + std::to_string(alpha_g));

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        for (const auto group_key : analysis_view.CollectAtomGroupKeys(class_key))
        {
            analysis.SetAtomGroupAlphaG(group_key, class_key, alpha_g);
        }
    }

    rhbm_trainer::AlphaTrainer::AlphaRunOptions alpha_bias_study_options;
    alpha_bias_study_options.execution_options =
        MakePotentialAnalysisExecutionOptions(thread_size, true);
    alpha_bias_study_options.progress_callback =
        [](std::size_t completed, std::size_t total)
        {
            Logger::ProgressPercent(completed, total);
        };

    if (!beta_group_list.empty())
    {
        const auto alpha_g_bias_matrix{
            alpha_trainer.StudyAlphaGBias(beta_group_list, alpha_bias_study_options)
        };
        const bool report_emitted{
            EmitTrainingReportIfRequested(
                alpha_g_bias_matrix, alpha_training_list, "#alpha_{g}", "Deviation with Mean",
                request.training_report_dir, "alpha_g_bias.pdf")
        };
        if (!report_emitted)
        {
            Logger::Log(LogLevel::Debug,
                "Alpha_G bias report was skipped (set --training-report-dir to emit PDF output).");
        }
    }
    else
    {
        Logger::Log(LogLevel::Warning,
            "Skip Alpha_G bias study because no eligible groups were selected.");
    }
}

void RunAtomPotentialFittingWorkflow(
    ModelObject & model_object,
    const PotentialAnalysisRequest & request,
    int thread_size)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomPotentialFittingWorkflow");
    if (request.training_alpha_flag)
    {
        RunAtomAlphaTraining(model_object, request, thread_size);
    }
    else
    {
        RunLocalPotentialFitting(model_object, thread_size);
    }

    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    std::unordered_map<const AtomObject *, MutableLocalPotentialView> local_entry_map;
    local_entry_map.reserve(model_object.GetSelectedAtomCount());
    for (auto * atom : model_object.GetSelectedAtoms())
    {
        local_entry_map.emplace(atom, analysis.EnsureAtomLocalPotential(*atom));
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
        #pragma omp parallel for num_threads(thread_size)
#endif
        for (size_t idx = 0; idx < group_key_size; idx++)
        {
            auto group_key{ group_key_list[idx] };
            const auto & atom_list{ analysis_view.GetAtomObjectList(group_key, class_key) };
            std::vector<RHBMMemberDataset> member_datasets;
            std::vector<RHBMBetaEstimateResult> member_fit_results;
            member_datasets.reserve(atom_list.size());
            member_fit_results.reserve(atom_list.size());
            for (const auto & atom : atom_list)
            {
                const auto local_entry{ local_entry_map.at(atom) };
                member_datasets.emplace_back(local_entry.GetDataset());
                member_fit_results.emplace_back(local_entry.GetFitResult());
            }
            auto alpha_g{ request.training_alpha_flag ?
                analysis_view.GetAtomAlphaG(group_key, class_key) : request.alpha_g
            };
            const auto result{
                rhbm_helper::EstimateGroup(
                    alpha_g,
                    rhbm_helper::BuildGroupInput(member_datasets, member_fit_results),
                    MakePotentialAnalysisExecutionOptions(thread_size, true))
            };

#ifdef USE_OPENMP
            #pragma omp critical
#endif
            {
                analysis.ApplyAtomGroupEstimateResult(group_key, class_key, result, alpha_g);
                key_count++;
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

void SavePreparedModel(
    ModelObject & model_object,
    const std::filesystem::path & database_path,
    std::string_view saved_key_tag)
{
    DataRepository repository{ database_path };
    repository.SaveModel(model_object, std::string(saved_key_tag));
    model_object.EditAnalysis().ClearTransientFitStates();
}

void RunPotentialSamplingWorkflow(
    MapObject & map_object,
    ModelObject & model_object,
    const PotentialAnalysisRequest & request,
    int thread_size)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunPotentialSamplingWorkflow");
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    const auto atom_size{ atom_list.size() };
    size_t atom_count{ 0 };
    auto local_entry_list{ BuildSelectedAtomLocalEntryViews(model_object) };
#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(thread_size)
#endif
    for (size_t i = 0; i < atom_size; i++)
    {
        auto atom{ atom_list[i] };
        auto entry{ local_entry_list[i] };
        auto sampling_entries{
            SampleAtomMapValues(map_object, *atom, request.sampling_profile_choice)
        };
        auto dataset{
            rhbm_helper::BuildMemberDataset(
                linearization_service::BuildDatasetSeries(
                    sampling_entries, request.fit_range_min, request.fit_range_max))
        };
        entry.SetSamplingEntries(std::move(sampling_entries));
        entry.SetDataset(std::move(dataset));

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
    const auto thread_size{ request.job_count };
    auto inputs{ LoadPotentialAnalysisInputs(request) };
    if (!inputs.has_value()) return false;

    auto & model_object{ *inputs->model_object };
    auto & map_object{ *inputs->map_object };
    if (!request.simulation_flag && request.map_normalization_flag)
    {
        map_object.MapValueArrayNormalization();
    }
    RunModelObjectPreprocessing(model_object, request.asymmetry_flag, request.alpha_r);
    RunPotentialSamplingWorkflow(map_object, model_object, request, thread_size);
    RunAtomPotentialFittingWorkflow(model_object, request, thread_size);
    SavePreparedModel(model_object, request.database_path, request.saved_key_tag);
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
