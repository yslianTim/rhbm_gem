#include "PotentialAnalysisCommand.hpp"

#include "detail/MapSampling.hpp"
#include "experimental/PotentialAnalysisBondWorkflow.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/LocalPainter.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/hrl/GaussianLinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

#include <atomic>
#include <cmath>
#include <iterator>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace {
constexpr std::string_view kModelKey{ "model" };
constexpr std::string_view kMapKey{ "map" };
constexpr std::string_view kModelOption{ "--model" };
constexpr std::string_view kMapOption{ "--map" };
constexpr std::string_view kSimResolutionOption{ "--sim-resolution" };
constexpr std::string_view kSaveKeyOption{ "--save-key" };
constexpr std::string_view kSamplingOption{ "--sampling" };
constexpr std::string_view kSamplingMinOption{ "--sampling-min" };
constexpr std::string_view kSamplingMaxOption{ "--sampling-max" };
constexpr std::string_view kSamplingHeightOption{ "--sampling-height" };
constexpr std::string_view kFitMinOption{ "--fit-min" };
constexpr std::string_view kFitMaxOption{ "--fit-max" };
constexpr std::string_view kAlphaROption{ "--alpha-r" };
constexpr std::string_view kAlphaGOption{ "--alpha-g" };
constexpr std::string_view kTrainingAlphaMinOption{ "--training-alpha-min" };
constexpr std::string_view kTrainingAlphaMaxOption{ "--training-alpha-max" };
constexpr std::string_view kTrainingAlphaStepOption{ "--training-alpha-step" };
constexpr std::string_view kSamplingRangeIssue{ "--sampling-range" };
constexpr std::string_view kFitRangeIssue{ "--fit-range" };
constexpr std::string_view kTrainingAlphaRangeIssue{ "--training-alpha-range" };

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
    if (training_report_dir.empty())
    {
        return false;
    }

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

    rhbm_gem::LinePlotRequest report_request;
    report_request.output_path = report_path;
    report_request.x_axis.title = std::string(x_label);
    report_request.shared_y_axis_title = std::string(y_label);
    report_request.panels = {
        rhbm_gem::LinePlotPanel{
            "Amplitude #font[2]{A}",
            rhbm_gem::AxisSpec{},
            { rhbm_gem::LineSeries{ "Amplitude", alpha_list, amplitude_values, std::nullopt } }
        },
        rhbm_gem::LinePlotPanel{
            "Width #tau",
            rhbm_gem::AxisSpec{},
            { rhbm_gem::LineSeries{ "Width", alpha_list, width_values, std::nullopt } }
        }
    };

    const auto plot_result{ rhbm_gem::local_painter::SaveLinePlot(report_request) };
    if (!plot_result.Succeeded())
    {
        Logger::Log(LogLevel::Warning,
            "Failed to emit training report '" + report_path.string() +
                "': " + plot_result.message);
        return false;
    }

    if (!std::filesystem::exists(report_path))
    {
        Logger::Log(LogLevel::Warning,
            "Training report output was requested but no file was produced: " +
                report_path.string());
        return false;
    }

    return true;
}

rhbm_gem::GaussianLinearizationContext BuildLocalLinearizationContext(
    const rhbm_gem::LocalPotentialView & view)
{
    Eigen::VectorXd gaus_par_init{ Eigen::VectorXd::Zero(3) };
    gaus_par_init(0) = view.GetMomentZeroEstimate();
    gaus_par_init(1) = view.GetMomentTwoEstimate();
    return rhbm_gem::GaussianLinearizationContext::FromModelParameters(gaus_par_init);
}

} // namespace

namespace rhbm_gem {

PotentialAnalysisCommand::PotentialAnalysisCommand() :
    CommandWithRequest<PotentialAnalysisRequest>{},
    m_model_key_tag{ kModelKey }, m_map_key_tag{ kMapKey },
    m_map_object{ nullptr }, m_model_object{ nullptr }
{
}

void PotentialAnalysisCommand::NormalizeRequest()
{
    auto & request{ MutableRequest() };
    ValidateRequiredPath(request.model_file_path, kModelOption, "Model file");
    ValidateRequiredPath(request.map_file_path, kMapOption, "Map file");
    CoerceFiniteNonNegativeScalar(
        request.simulated_map_resolution,
        kSimResolutionOption,
        0.0,
        LogLevel::Error,
        "Simulated map resolution");
    InvalidatePreparedState();
    ClearParseIssues(kSaveKeyOption);
    if (request.saved_key_tag.empty())
    {
        request.saved_key_tag = "model";
        AddValidationError(
            kSaveKeyOption,
            "Saved key tag cannot be empty. Using 'model' instead.",
            ValidationPhase::Parse);
    }
    CoercePositiveScalar(
        request.sampling_size,
        kSamplingOption,
        1500,
        LogLevel::Warning,
        "Sampling size");
    CoerceFiniteNonNegativeScalar(
        request.sampling_range_min,
        kSamplingMinOption,
        0.0,
        LogLevel::Error,
        "Minimum sampling range");
    CoerceFiniteNonNegativeScalar(
        request.sampling_range_max,
        kSamplingMaxOption,
        1.5,
        LogLevel::Error,
        "Maximum sampling range");
    CoerceFinitePositiveScalar(
        request.sampling_height,
        kSamplingHeightOption,
        0.1,
        LogLevel::Error,
        "Sampling height");
    CoerceFiniteNonNegativeScalar(
        request.fit_range_min,
        kFitMinOption,
        0.0,
        LogLevel::Error,
        "Minimum fitting range");
    CoerceFiniteNonNegativeScalar(
        request.fit_range_max,
        kFitMaxOption,
        1.0,
        LogLevel::Error,
        "Maximum fitting range");
    CoerceFinitePositiveScalar(
        request.alpha_r,
        kAlphaROption,
        0.1,
        LogLevel::Error,
        "Alpha-R");
    CoerceFinitePositiveScalar(
        request.alpha_g,
        kAlphaGOption,
        0.2,
        LogLevel::Error,
        "Alpha-G");
    CoerceFiniteNonNegativeScalar(
        request.training_alpha_min,
        kTrainingAlphaMinOption,
        0.0,
        LogLevel::Error,
        "Minimum training alpha");
    CoerceFiniteNonNegativeScalar(
        request.training_alpha_max,
        kTrainingAlphaMaxOption,
        1.0,
        LogLevel::Error,
        "Maximum training alpha");
    CoerceFinitePositiveScalar(
        request.training_alpha_step,
        kTrainingAlphaStepOption,
        0.1,
        LogLevel::Error,
        "Training alpha step");
}

bool PotentialAnalysisCommand::ExecuteImpl()
{
    const auto & request{ RequestOptions() };
    if (!BuildDataObject(request)) return false;
    auto & model_object{ *m_model_object };
    auto & map_object{ *m_map_object };
    RunMapObjectPreprocessing(map_object);
    RunModelObjectPreprocessing(model_object, request.asymmetry_flag);
    RunSamplingWorkflow(
        map_object, model_object,
        request.sampling_size, request.sampling_range_min, request.sampling_range_max);
    RunDatasetPreparationWorkflow(model_object, request.fit_range_min, request.fit_range_max);
    RunAtomPotentialFittingWorkflow(model_object, request);
#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
    experimental::RunPotentialAnalysisBondWorkflow(
        model_object, map_object, request, ThreadSize());
#endif
    SavePreparedModel(model_object, request.saved_key_tag);
    return true;
}

void PotentialAnalysisCommand::ValidateOptions()
{
    const auto & request{ RequestOptions() };
    RequireCondition(
        !request.simulation_flag || request.simulated_map_resolution > 0.0,
        kSimResolutionOption,
        "Expected a positive simulated-map resolution when '--simulation true' is selected.");
    RequireCondition(
        request.sampling_range_min <= request.sampling_range_max,
        kSamplingRangeIssue,
        "Expected --sampling-min <= --sampling-max.");
    RequireCondition(
        request.fit_range_min <= request.fit_range_max,
        kFitRangeIssue,
        "Expected --fit-min <= --fit-max.");
    RequireCondition(
        request.training_alpha_min <= request.training_alpha_max,
        kTrainingAlphaRangeIssue,
        "Expected --training-alpha-min <= --training-alpha-max.");
}

void PotentialAnalysisCommand::ResetRuntimeState()
{
    m_map_object.reset();
    m_model_object.reset();
}

bool PotentialAnalysisCommand::BuildDataObject(const PotentialAnalysisRequest & request)
{
    ScopeTimer timer("PotentialAnalysisCommand::BuildDataObject");
    try
    {
        AttachDataRepository(request.database_path);
        m_model_object = LoadInputFile<ModelObject>(request.model_file_path, m_model_key_tag);
        m_map_object = LoadInputFile<MapObject>(request.map_file_path, m_map_key_tag);
        if (m_model_object == nullptr || m_map_object == nullptr)
        {
            Logger::Log(
                LogLevel::Error,
                "PotentialAnalysisCommand::BuildDataObject : model/map object missing after load.");
            return false;
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
            m_model_object->SetEmdID("Simulation");
            m_model_object->SetResolution(request.simulated_map_resolution);
            m_model_object->SetResolutionMethod("Blurring Width");
        }
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "PotentialAnalysisCommand::BuildDataObject : " + std::string(e.what()));
        return false;
    }
    return true;
}

std::vector<MutableLocalPotentialView>
PotentialAnalysisCommand::BuildSelectedAtomLocalEntryViews(ModelObject & model_object)
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

void PotentialAnalysisCommand::RunMapObjectPreprocessing(MapObject & map_object)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunMapObjectPreprocessing");
    map_object.MapValueArrayNormalization();
}

void PotentialAnalysisCommand::RunModelObjectPreprocessing(
    ModelObject & model_object,
    bool asymmetry_flag)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunModelObjectPreprocessing");
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
    // selection is finalized, local entries exist, and atom groups are materialized.
    analysis.RebuildAtomGroupsFromSelection();

    Logger::Log(LogLevel::Info,
        "Number of selected atom = "
            + std::to_string(model_object.GetSelectedAtomCount()));
    Logger::Log(LogLevel::Info,
        "Number of selected bond = "
            + std::to_string(model_object.GetSelectedBondCount()));
    Logger::Log(LogLevel::Info, model_object.GetAnalysisView().DescribeAtomGrouping());
    if (model_object.GetNumberOfAtom() > 0 &&
        model_object.GetSelectedAtomCount() == 0)
    {
        Logger::Log(LogLevel::Warning,
            "No atoms are selected after symmetry filtering. "
            "The input CIF may miss usable _entity/_struct_asym metadata. "
            "Try '--asymmetry true' to bypass symmetry filtering.");
    }
}

void PotentialAnalysisCommand::RunAtomPotentialFittingWorkflow(
    ModelObject & model_object,
    const PotentialAnalysisRequest & request)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomPotentialFittingWorkflow");
    if (request.training_alpha_flag)
    {
        RunAtomAlphaTraining(model_object, request);
    }
    else
    {
        RunLocalPotentialFitting(model_object, request.alpha_r);
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
        auto group_keys{ analysis_view.CollectAtomGroupKeys(class_key) };
        auto group_key_size{ group_keys.size() };
        std::atomic<size_t> key_count{ 0 };

#ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(ThreadSize())
#endif
        for (size_t idx = 0; idx < group_key_size; idx++)
        {
            auto group_key{ group_keys[idx] };
            const auto & atom_list{ analysis_view.GetAtomObjectList(group_key, class_key) };
            auto group_size{ atom_list.size() };
            std::vector<RHBMMemberDataset> member_datasets;
            std::vector<RHBMBetaEstimateResult> member_fit_results;
            member_datasets.reserve(group_size);
            member_fit_results.reserve(group_size);
            for (const auto & atom : atom_list)
            {
                const auto local_entry{ local_entry_map.at(atom) };
                member_datasets.emplace_back(local_entry.GetDataset());
                member_fit_results.emplace_back(local_entry.GetFitResult());
            }
            auto alpha_g{ request.training_alpha_flag ?
                analysis_view.GetAtomAlphaG(group_key, class_key) : request.alpha_g
            };
            const auto input{
                rhbm_gem::rhbm_helper::BuildGroupInput(member_datasets, member_fit_results)
            };
            const auto result{
                rhbm_gem::rhbm_helper::EstimateGroup(
                    alpha_g,
                    input,
                    MakePotentialAnalysisExecutionOptions(ThreadSize(), true))
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

void PotentialAnalysisCommand::SavePreparedModel(
    ModelObject & model_object, std::string_view saved_key_tag)
{
    ScopeTimer timer("PotentialAnalysisCommand::SavePreparedModel");
    SaveStoredObject(m_model_key_tag, std::string(saved_key_tag));
    model_object.EditAnalysis().ClearTransientFitStates();
}

void PotentialAnalysisCommand::RunSamplingWorkflow(
    MapObject & map_object,
    ModelObject & model_object,
    int sampling_size,
    double sampling_range_min,
    double sampling_range_max)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunSamplingWorkflow");
    auto thread_size{ ThreadSize() };
    SphereSampler sampler;
    sampler.SetSamplingProfile(
        //SphereSamplingProfile::RadiusUniformRandom(
        //    SphereDistanceRange{ sampling_range_min, sampling_range_max },
        //    static_cast<unsigned int>(sampling_size))
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ sampling_range_min, sampling_range_max },
            0.1,
            static_cast<unsigned int>(sampling_size))
    );
    sampler.Print();

    const auto & atom_list{ model_object.GetSelectedAtoms() };
    const auto atom_size{ atom_list.size() };
    size_t atom_count{ 0 };
    auto local_entry_list{ BuildSelectedAtomLocalEntryViews(model_object) };

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(thread_size)
    {
        #pragma omp for
        for (size_t i = 0; i < atom_size; i++)
        {
            auto atom{ atom_list[i] };
            auto entry{ local_entry_list[i] };
            auto sampling_entries{
                SampleMapValues(map_object, sampler, *atom, sampling_range_max)
            };
            entry.SetSamplingEntries(sampling_entries);
            #pragma omp critical
            {
                atom_count++;
                Logger::ProgressPercent(atom_count, atom_size);
            }
        }
    }
#else
    for (size_t i = 0; i < atom_size; i++)
    {
        auto atom{ atom_list[i] };
        auto entry{ local_entry_list[i] };
        auto sampling_entries{
            SampleMapValues(map_object, sampler, *atom, sampling_range_max)
        };
        entry.SetSamplingEntries(sampling_entries);
        atom_count++;
        Logger::ProgressPercent(atom_count, atom_size);
    }
#endif
}

void PotentialAnalysisCommand::RunDatasetPreparationWorkflow(
    ModelObject & model_object,
    double fit_range_min,
    double fit_range_max)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunDatasetPreparationWorkflow");
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    const auto atom_size{ atom_list.size() };
    size_t atom_count{ 0 };
    auto local_entry_list{ BuildSelectedAtomLocalEntryViews(model_object) };
    const auto dataset_service{
        rhbm_gem::GaussianLinearizationService{
            rhbm_gem::GaussianLinearizationSpec::DefaultDataset()
        }
    };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(ThreadSize())
#endif
    for (size_t i = 0; i < atom_size; i++)
    {
        auto entry{ local_entry_list[i] };
        const auto local_view{ LocalPotentialView::RequireFor(*atom_list[i]) };
        entry.SetDataset(
            dataset_service.BuildDataset(
                local_view.GetSamplingEntries(),
                fit_range_min,
                fit_range_max,
                BuildLocalLinearizationContext(local_view))
        );

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            Logger::ProgressPercent(atom_count, atom_size);
        }
    }
}

void PotentialAnalysisCommand::RunAtomAlphaTraining(
    ModelObject & model_object,
    const PotentialAnalysisRequest & request)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomAlphaTraining");
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    rhbm_gem::rhbm_trainer::AlphaTrainer alpha_trainer(
        request.training_alpha_min,
        request.training_alpha_max,
        request.training_alpha_step);
    const auto & alpha_training_list{ alpha_trainer.AlphaGrid() };
    Logger::Log(LogLevel::Info, alpha_trainer.GetAlphaGridSummary().str());

    // Alpha_R Training
    std::vector<RHBMMemberDataset> selected_atom_dataset_list;
    selected_atom_dataset_list.reserve(model_object.GetSelectedAtomCount());
    for (auto & atom : model_object.GetSelectedAtoms())
    {
        if (atom->IsMainChainAtom() == false) continue;
        const auto local_entry{ analysis.EnsureAtomLocalPotential(*atom) };
        if (!local_entry.HasDataset()) continue;
        selected_atom_dataset_list.emplace_back(local_entry.GetDataset());
    }
    selected_atom_dataset_list.shrink_to_fit();
    
    auto selected_atom_size{ selected_atom_dataset_list.size() };
    Logger::Log(LogLevel::Info,
        "Run Alpha_R Training with "+ std::to_string(selected_atom_size) +" atoms.");
    auto alpha_r{ alpha_training_list.front() };
    if (!selected_atom_dataset_list.empty())
    {
        rhbm_gem::rhbm_trainer::AlphaTrainer::AlphaTrainingOptions alpha_r_options;
        alpha_r_options.subset_size = 5;
        alpha_r_options.execution_options =
            MakePotentialAnalysisExecutionOptions(ThreadSize(), true);
        alpha_r_options.progress_callback =
            [](std::size_t completed, std::size_t total)
            {
                Logger::ProgressPercent(completed, total);
            };
        const auto alpha_r_result{
            alpha_trainer.TrainAlphaR(selected_atom_dataset_list, alpha_r_options)
        };
        alpha_r = alpha_r_result.best_alpha;
    }
    else
    {
        Logger::Log(LogLevel::Warning,
            "Skip Alpha_R Training because no eligible atoms were selected.");
    }
    Logger::Log(LogLevel::Info,
        "Alpha_R Training Results Summary: best alpha_r = "+ std::to_string(alpha_r));
    
    // Alpha_G Training
    RunLocalPotentialFitting(model_object, alpha_r);
    
    std::vector<std::vector<RHBMBetaVector>> beta_group_list;
    const auto component_class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    const auto component_group_keys{ analysis_view.CollectAtomGroupKeys(component_class_key) };
    beta_group_list.reserve(component_group_keys.size());
    for (const auto group_key : component_group_keys)
    {
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(group_key, component_class_key) };
        if (group_atom_list.size() < 10) continue;
        if (group_atom_list.front()->IsMainChainAtom() == false) continue;

        std::vector<RHBMBetaVector> beta_list;
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
        rhbm_gem::rhbm_trainer::AlphaTrainer::AlphaTrainingOptions alpha_g_options;
        alpha_g_options.subset_size = 10;
        alpha_g_options.execution_options =
            MakePotentialAnalysisExecutionOptions(ThreadSize(), true);
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

    rhbm_gem::rhbm_trainer::AlphaTrainer::AlphaRunOptions alpha_bias_study_options;
    alpha_bias_study_options.execution_options =
        MakePotentialAnalysisExecutionOptions(ThreadSize(), true);
    alpha_bias_study_options.progress_callback =
        [](std::size_t completed, std::size_t total)
        {
            Logger::ProgressPercent(completed, total);
        };

    if (!selected_atom_dataset_list.empty())
    {
        const auto alpha_r_bias_matrix{
            alpha_trainer.StudyAlphaRBias(selected_atom_dataset_list, alpha_bias_study_options)
        };
        const bool report_emitted{
            EmitTrainingReportIfRequested(
                alpha_r_bias_matrix,
                alpha_training_list,
                "#alpha_{r}",
                "Deviation with OLS",
                request.training_report_dir,
                "alpha_r_bias.pdf")
        };
        if (!report_emitted)
        {
            Logger::Log(LogLevel::Debug,
                "Alpha_R bias report was skipped (set --training-report-dir to emit PDF output).");
        }
    }
    else
    {
        Logger::Log(LogLevel::Warning,
            "Skip Alpha_R bias study because no eligible atoms were selected.");
    }

    if (!beta_group_list.empty())
    {
        const auto alpha_g_bias_matrix{
            alpha_trainer.StudyAlphaGBias(beta_group_list, alpha_bias_study_options)
        };
        const bool report_emitted{
            EmitTrainingReportIfRequested(
                alpha_g_bias_matrix,
                alpha_training_list,
                "#alpha_{g}",
                "Deviation with Mean",
                request.training_report_dir,
                "alpha_g_bias.pdf")
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

void PotentialAnalysisCommand::RunLocalPotentialFitting(
    ModelObject & model_object,
    double alpha_r)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunLocalPotentialFitting");
    auto thread_size{ ThreadSize() };
    std::atomic<size_t> atom_count{ 0 };
    const auto & selected_atom_list{ model_object.GetSelectedAtoms() };
    const auto selected_atom_size{ selected_atom_list.size() };
    auto local_entry_list{ BuildSelectedAtomLocalEntryViews(model_object) };
    Logger::Log(LogLevel::Info,
        "Run Local atom fitting for " + std::to_string(selected_atom_size) + " atoms.");

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        auto local_entry{ local_entry_list[i] };
        local_entry.SetAlphaR(alpha_r);
        const auto result{
            rhbm_gem::rhbm_helper::EstimateBetaMDPDE(
                alpha_r,
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

} // namespace rhbm_gem
