#include "internal/command/PotentialAnalysisCommand.hpp"
#include "internal/command/CommandCliSupport.hpp"
#include "internal/command/CommandDataSupport.hpp"
#include "internal/command/MapSampling.hpp"
#include "experimental/PotentialAnalysisBondWorkflow.hpp"
#include <rhbm_gem/core/command/CommandApi.hpp>
#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/GroupPotentialEntry.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/LocalPainter.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/hrl/HRLAlphaTrainer.hpp>
#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>
#include <rhbm_gem/utils/hrl/HRLGroupEstimator.hpp>
#include <rhbm_gem/utils/hrl/HRLModelAlgorithms.hpp>
#include <rhbm_gem/utils/hrl/HRLModelTypes.hpp>
#include <rhbm_gem/utils/math/GausLinearTransformHelper.hpp>
#include <rhbm_gem/utils/math/SphereSampler.hpp>

#include <algorithm>
#include <atomic>
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
constexpr std::string_view kSamplingRangeIssue{ "--sampling-range" };
constexpr std::string_view kFitRangeIssue{ "--fit-range" };

HRLExecutionOptions MakePotentialAnalysisExecutionOptions(
    int thread_size,
    bool quiet_mode)
{
    HRLExecutionOptions execution_options;
    execution_options.quiet_mode = quiet_mode;
    execution_options.thread_size = thread_size;
    return execution_options;
}
}

namespace rhbm_gem {

PotentialAnalysisCommand::PotentialAnalysisCommand(CommonOptionProfile profile) :
    CommandBase{ profile },
    m_model_key_tag{ kModelKey }, m_map_key_tag{ kMapKey },
    m_map_object{ nullptr }, m_model_object{ nullptr }
{
}

void BindPotentialAnalysisRequestOptions(
    CLI::App * command,
    PotentialAnalysisRequest & request)
{
    command_cli::AddPathOption(
        command,
        "-a,--model",
        [&](const std::filesystem::path & value) { request.model_file_path = value; },
        "Model file path",
        std::nullopt,
        true);
    command_cli::AddPathOption(
        command,
        "-m,--map",
        [&](const std::filesystem::path & value) { request.map_file_path = value; },
        "Map file path",
        std::nullopt,
        true);
    command_cli::AddPathOption(
        command,
        "--training-report-dir",
        [&](const std::filesystem::path & value) { request.training_report_dir = value; },
        "Optional output directory for alpha training reports");
    command_cli::AddScalarOption<bool>(
        command,
        "--simulation",
        [&](bool value) { request.simulation_flag = value; },
        "Simulation flag",
        request.simulation_flag);
    command_cli::AddScalarOption<double>(
        command,
        "-r,--sim-resolution",
        [&](double value) { request.simulated_map_resolution = value; },
        "Set simulated map's resolution (blurring width)",
        request.simulated_map_resolution);
    command_cli::AddStringOption(
        command,
        "-k,--save-key",
        [&](const std::string & value) { request.saved_key_tag = value; },
        "New key tag for saving ModelObject results into database",
        request.saved_key_tag);
    command_cli::AddScalarOption<bool>(
        command,
        "--training-alpha",
        [&](bool value) { request.training_alpha_flag = value; },
        "Turn On/Off alpha training flag",
        request.training_alpha_flag);
    command_cli::AddScalarOption<bool>(
        command,
        "--asymmetry",
        [&](bool value) { request.asymmetry_flag = value; },
        "Turn On/Off asymmetry flag",
        request.asymmetry_flag);
    command_cli::AddScalarOption<int>(
        command,
        "-s,--sampling",
        [&](int value) { request.sampling_size = value; },
        "Number of sampling points per atom",
        request.sampling_size);
    command_cli::AddScalarOption<double>(
        command,
        "--sampling-min",
        [&](double value) { request.sampling_range_min = value; },
        "Minimum sampling range",
        request.sampling_range_min);
    command_cli::AddScalarOption<double>(
        command,
        "--sampling-max",
        [&](double value) { request.sampling_range_max = value; },
        "Maximum sampling range",
        request.sampling_range_max);
    command_cli::AddScalarOption<double>(
        command,
        "--sampling-height",
        [&](double value) { request.sampling_height = value; },
        "Maximum sampling height",
        request.sampling_height);
    command_cli::AddScalarOption<double>(
        command,
        "--fit-min",
        [&](double value) { request.fit_range_min = value; },
        "Minimum fitting range",
        request.fit_range_min);
    command_cli::AddScalarOption<double>(
        command,
        "--fit-max",
        [&](double value) { request.fit_range_max = value; },
        "Maximum fitting range",
        request.fit_range_max);
    command_cli::AddScalarOption<double>(
        command,
        "--alpha-r",
        [&](double value) { request.alpha_r = value; },
        "Alpha value for R",
        request.alpha_r);
    command_cli::AddScalarOption<double>(
        command,
        "--alpha-g",
        [&](double value) { request.alpha_g = value; },
        "Alpha value for G",
        request.alpha_g);
}

void PotentialAnalysisCommand::ApplyRequest(const PotentialAnalysisRequest & request)
{
    ApplyCommonRequest(request.common);
    SetModelFilePath(request.model_file_path);
    SetMapFilePath(request.map_file_path);
    AssignOption(m_options.is_simulation, request.simulation_flag);
    SetSimulatedMapResolution(request.simulated_map_resolution);
    SetSavedKeyTag(request.saved_key_tag);
    AssignOption(m_options.training_report_dir, request.training_report_dir);
    AssignOption(m_options.use_training_alpha, request.training_alpha_flag);
    AssignOption(m_options.is_asymmetry, request.asymmetry_flag);
    SetSamplingSize(request.sampling_size);
    SetSamplingRangeMinimum(request.sampling_range_min);
    SetSamplingRangeMaximum(request.sampling_range_max);
    SetSamplingHeight(request.sampling_height);
    SetFitRangeMinimum(request.fit_range_min);
    SetFitRangeMaximum(request.fit_range_max);
    SetAlphaR(request.alpha_r);
    SetAlphaG(request.alpha_g);
}

bool PotentialAnalysisCommand::ExecuteImpl()
{
    if (BuildDataObject() == false) return false;
    RunMapObjectPreprocessing();
    RunModelObjectPreprocessing();

    RunAtomMapValueSampling();
    RunAtomGroupClassification();
    if (m_options.use_training_alpha) RunAtomAlphaTraining();
    else RunLocalAtomFitting(m_options.alpha_r);
    RunAtomPotentialFitting();
    RunExperimentalBondWorkflowIfEnabled();
    SaveDataObject();
    return true;
}

void PotentialAnalysisCommand::ValidateOptions()
{
    ResetPrepareIssues(kSimResolutionOption);
    if (m_options.is_simulation && m_options.resolution_simulation <= 0.0)
    {
        AddValidationError(
            kSimResolutionOption,
            "Expected a positive simulated-map resolution when '--simulation true' is selected.");
    }

    ResetPrepareIssues(kSamplingRangeIssue);
    if (m_options.sampling_range_min > m_options.sampling_range_max)
    {
        AddValidationError(kSamplingRangeIssue,
            "Expected --sampling-min <= --sampling-max.");
    }

    ResetPrepareIssues(kFitRangeIssue);
    if (m_options.fit_range_min > m_options.fit_range_max)
    {
        AddValidationError(kFitRangeIssue,
            "Expected --fit-min <= --fit-max.");
    }
}

void PotentialAnalysisCommand::ResetRuntimeState()
{
    m_map_object.reset();
    m_model_object.reset();
}

void PotentialAnalysisCommand::SetSimulatedMapResolution(double value)
{
    SetFiniteNonNegativeScalarOption(
        m_options.resolution_simulation,
        value,
        kSimResolutionOption,
        0.0,
        "Simulated-map resolution must be a finite non-negative value.");
}

void PotentialAnalysisCommand::SetSavedKeyTag(const std::string & tag)
{
    MutateOptions([&]()
    {
        ResetParseIssues(kSaveKeyOption);
        if (!tag.empty())
        {
            m_options.saved_key_tag = tag;
            return;
        }

        m_options.saved_key_tag = "model";
        AddValidationError(
            kSaveKeyOption,
            "Saved key tag cannot be empty.",
            ValidationPhase::Parse);
    });
}

void PotentialAnalysisCommand::SetSamplingSize(int value)
{
    SetNormalizedScalarOption(
        m_options.sampling_size,
        value,
        kSamplingOption,
        [](int candidate) { return candidate > 0; },
        1500,
        "Sampling size must be positive, reset to default value = 1500");
}

void PotentialAnalysisCommand::SetSamplingRangeMinimum(double value)
{
    SetFiniteNonNegativeScalarOption(
        m_options.sampling_range_min,
        value,
        kSamplingMinOption,
        0.0,
        "Minimum sampling range must be a finite non-negative value.");
}

void PotentialAnalysisCommand::SetSamplingRangeMaximum(double value)
{
    SetFiniteNonNegativeScalarOption(
        m_options.sampling_range_max,
        value,
        kSamplingMaxOption,
        1.5,
        "Maximum sampling range must be a finite non-negative value.");
}

void PotentialAnalysisCommand::SetSamplingHeight(double value)
{
    SetFinitePositiveScalarOption(
        m_options.sampling_height,
        value,
        kSamplingHeightOption,
        0.1,
        "Sampling height must be a finite positive value.");
}

void PotentialAnalysisCommand::SetFitRangeMinimum(double value)
{
    SetFiniteNonNegativeScalarOption(
        m_options.fit_range_min,
        value,
        kFitMinOption,
        0.0,
        "Minimum fitting range must be a finite non-negative value.");
}

void PotentialAnalysisCommand::SetFitRangeMaximum(double value)
{
    SetFiniteNonNegativeScalarOption(
        m_options.fit_range_max,
        value,
        kFitMaxOption,
        1.0,
        "Maximum fitting range must be a finite non-negative value.");
}

void PotentialAnalysisCommand::SetAlphaR(double value)
{
    SetFinitePositiveScalarOption(
        m_options.alpha_r,
        value,
        kAlphaROption,
        0.1,
        "Alpha-R must be a finite positive value.");
}

void PotentialAnalysisCommand::SetAlphaG(double value)
{
    SetFinitePositiveScalarOption(
        m_options.alpha_g,
        value,
        kAlphaGOption,
        0.2,
        "Alpha-G must be a finite positive value.");
}

void PotentialAnalysisCommand::SetModelFilePath(const std::filesystem::path & path)
{
    SetRequiredExistingPathOption(m_options.model_file_path, path, kModelOption, "Model file");
}

void PotentialAnalysisCommand::SetMapFilePath(const std::filesystem::path & path)
{
    SetRequiredExistingPathOption(m_options.map_file_path, path, kMapOption, "Map file");
}

} // namespace rhbm_gem

namespace rhbm_gem::detail {

std::vector<double> BuildOrderedAlphaRTrainingList()
{
    std::vector<double> alpha_list{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    std::sort(alpha_list.begin(), alpha_list.end());
    return alpha_list;
}

std::vector<double> BuildOrderedAlphaGTrainingList()
{
    std::vector<double> alpha_list{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    std::sort(alpha_list.begin(), alpha_list.end());
    return alpha_list;
}

std::filesystem::path ResolveTrainingReportPath(
    const PotentialAnalysisCommandOptions & options,
    std::string_view report_file_name)
{
    if (options.training_report_dir.empty())
    {
        return {};
    }
    return options.training_report_dir / std::filesystem::path{ report_file_name };
}

bool EmitTrainingReportIfRequested(
    const Eigen::MatrixXd & gaus_bias_matrix,
    const std::vector<double> & alpha_list,
    std::string_view x_label,
    std::string_view y_label,
    const PotentialAnalysisCommandOptions & options,
    std::string_view report_file_name)
{
    const auto report_path{ ResolveTrainingReportPath(options, report_file_name) };
    if (report_path.empty())
    {
        return false;
    }

    std::error_code ec;
    std::filesystem::create_directories(report_path.parent_path(), ec);
    if (ec)
    {
        Logger::Log(
            LogLevel::Warning,
            "Failed to create training report directory '" +
                report_path.parent_path().string() + "': " + ec.message());
        return false;
    }

    try
    {
        LocalPainter::PaintTemplate1(
            gaus_bias_matrix,
            alpha_list,
            std::string(x_label),
            std::string(y_label),
            report_path.string());
    }
    catch (const std::exception & ex)
    {
        Logger::Log(
            LogLevel::Warning,
            "Failed to emit training report '" + report_path.string() + "': " + ex.what());
        return false;
    }

    if (!std::filesystem::exists(report_path))
    {
        Logger::Log(
            LogLevel::Warning,
            "Training report output was requested but no file was produced: " +
                report_path.string());
        return false;
    }

    return true;
}

} // namespace rhbm_gem::detail

namespace rhbm_gem {
bool PotentialAnalysisCommand::BuildDataObject()
{
    ScopeTimer timer("PotentialAnalysisCommand::BuildDataObject");
    try
    {
        m_data_manager.SetDatabaseManager(DatabasePath());
        m_model_object = command_data_loader::ProcessModelFile(
            m_data_manager, m_options.model_file_path, m_model_key_tag, "model file");
        m_map_object = command_data_loader::ProcessMapFile(
            m_data_manager, m_options.map_file_path, m_map_key_tag, "map file");
        if (m_options.is_simulation == true)
        {
            UpdateModelObjectForSimulation(m_model_object.get());
        }
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error,
            "PotentialAnalysisCommand::Execute() : " + std::string(e.what()));
        return false;
    }
    return true;
}

void PotentialAnalysisCommand::UpdateModelObjectForSimulation(ModelObject * model_object)
{
    if (model_object == nullptr) return;
    if (m_options.resolution_simulation == 0.0)
    {
        Logger::Log(LogLevel::Warning,
            "[Warning] The resolution of input simulated map hasn't been set.\n"
            "          Please give the corresponding resolution value for this map.\n"
            "          (-r, --sim-resolution)");
    }
    model_object->SetEmdID("Simulation");
    model_object->SetResolution(m_options.resolution_simulation);
    model_object->SetResolutionMethod("Blurring Width");
}

void PotentialAnalysisCommand::RunMapObjectPreprocessing()
{
    ScopeTimer timer("PotentialAnalysisCommand::RunMapObjectPreprocessing");
    NormalizeMapObject(*m_map_object);
}

void PotentialAnalysisCommand::RunModelObjectPreprocessing()
{
    ScopeTimer timer("PotentialAnalysisCommand::RunModelObjectPreprocessing");
    ModelPreparationOptions options;
    options.select_all_atoms = true;
    options.select_all_bonds = true;
    options.apply_atom_symmetry_filter = true;
    options.apply_bond_symmetry_filter = true;
    options.asymmetry_flag = m_options.is_asymmetry;
    options.update_model = true;
    options.initialize_atom_local_entries = true;
    options.initialize_bond_local_entries = true;
    PrepareModelObject(*m_model_object, options);
    Logger::Log(LogLevel::Info,
        "Number of selected atom = " + std::to_string(m_model_object->GetNumberOfSelectedAtom()));
    Logger::Log(LogLevel::Info,
        "Number of selected bond = " + std::to_string(m_model_object->GetNumberOfSelectedBond()));
    if (m_model_object->GetNumberOfAtom() > 0 &&
        m_model_object->GetNumberOfSelectedAtom() == 0)
    {
        Logger::Log(LogLevel::Warning,
            "No atoms are selected after symmetry filtering. "
            "The input CIF may miss usable _entity/_struct_asym metadata. "
            "Try '--asymmetry true' to bypass symmetry filtering.");
    }
}

void PotentialAnalysisCommand::RunAtomPotentialFitting()
{
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomPotentialFitting");
    if (m_model_object == nullptr) return;
    const int basis_size{ 2 };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        Logger::Log(LogLevel::Info, "Class type: " + class_key);

        // Group Atom Potential Fitting
        auto group_potential_entry{ m_model_object->GetAtomGroupPotentialEntry(class_key) };
        const auto & key_set{ group_potential_entry->GetGroupKeySet() };
        std::vector<GroupKey> group_keys(key_set.begin(), key_set.end());
        auto group_key_size{ group_keys.size() };
        std::atomic<size_t> key_count{ 0 };

#ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(ThreadSize())
#endif
        for (size_t idx = 0; idx < group_key_size; idx++)
        {
            auto group_key{ group_keys[idx] };
            const auto & atom_list{ group_potential_entry->GetAtomObjectPtrList(group_key) };
            auto group_size{ atom_list.size() };
            std::vector<std::vector<Eigen::VectorXd>> data_entry_list;
            std::vector<Eigen::VectorXd> beta_mdpde_list;
            std::vector<double> sigma_square_list;
            std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> data_weight_list;
            std::vector<Eigen::DiagonalMatrix<double, Eigen::Dynamic>> data_covariance_list;
            data_entry_list.reserve(group_size);
            beta_mdpde_list.reserve(group_size);
            sigma_square_list.reserve(group_size);
            data_weight_list.reserve(group_size);
            data_covariance_list.reserve(group_size);
            for (const auto & atom : atom_list)
            {
                auto entry{ atom->GetLocalPotentialEntry() };
                data_entry_list.emplace_back(entry->GetBasisAndResponseEntryList());
                beta_mdpde_list.emplace_back(entry->GetBetaEstimateMDPDE());
                sigma_square_list.emplace_back(entry->GetSigmaSquare());
                data_weight_list.emplace_back(entry->GetDataWeight());
                data_covariance_list.emplace_back(entry->GetDataCovariance());
            }
            auto alpha_g{ (m_options.use_training_alpha) ?
                group_potential_entry->GetAlphaG(group_key) : m_options.alpha_g
            };
            const auto input = HRLDataTransform::BuildGroupInput(
                basis_size,
                data_entry_list,
                beta_mdpde_list,
                sigma_square_list,
                data_weight_list,
                data_covariance_list
            );
            HRLGroupEstimator estimator(
                MakePotentialAnalysisExecutionOptions(ThreadSize(), true));
            const auto result = estimator.Estimate(input, alpha_g);

            auto gaus_group_mean{
                GausLinearTransformHelper::BuildGaus3DModel(result.mu_mean)
            };

            auto gaus_group_mdpde{
                GausLinearTransformHelper::BuildGaus3DModel(result.mu_mdpde)
            };

            auto gaus_prior{
                GausLinearTransformHelper::BuildGaus3DModelWithVariance(
                    result.mu_prior, result.capital_lambda)
            };
            auto prior_estimate{ std::get<0>(gaus_prior) };
            auto prior_variance{ std::get<1>(gaus_prior) };

            auto count{ 0 };
            for (const auto & atom : atom_list)
            {
                auto atom_entry{ atom->GetLocalPotentialEntry() };
                const auto beta_vector_posterior{
                    result.beta_posterior_array.col(static_cast<Eigen::Index>(count))
                };
                const auto & sigma_matrix_posterior{
                    result.capital_sigma_posterior_list.at(static_cast<std::size_t>(count))
                };
                auto gaus_posterior{
                    GausLinearTransformHelper::BuildGaus3DModelWithVariance(
                        beta_vector_posterior, sigma_matrix_posterior)
                };
                auto posterior_estimate{ std::get<0>(gaus_posterior) };
                auto posterior_variance{ std::get<1>(gaus_posterior) };
                atom_entry->AddGausEstimatePosterior(class_key, posterior_estimate(0), posterior_estimate(1));
                atom_entry->AddGausVariancePosterior(class_key, posterior_variance(0), posterior_variance(1));
                atom_entry->AddOutlierTag(class_key, result.outlier_flag_array(count));
                atom_entry->AddStatisticalDistance(class_key, result.statistical_distance_array(count));
                count++;
            }

#ifdef USE_OPENMP
            #pragma omp critical
#endif
            {
                group_potential_entry->AddGausEstimateMean(
                    group_key, gaus_group_mean(0), gaus_group_mean(1)
                );
                group_potential_entry->AddGausEstimateMDPDE(
                    group_key, gaus_group_mdpde(0), gaus_group_mdpde(1)
                );
                group_potential_entry->AddGausEstimatePrior(
                    group_key, prior_estimate(0), prior_estimate(1)
                );
                group_potential_entry->AddGausVariancePrior(
                    group_key, prior_variance(0), prior_variance(1)
                );
                group_potential_entry->AddAlphaG(group_key, alpha_g);
                key_count++;
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

void PotentialAnalysisCommand::RunExperimentalBondWorkflowIfEnabled()
{
#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
    if (m_model_object == nullptr || m_map_object == nullptr) return;
    experimental::RunPotentialAnalysisBondWorkflow(
        *m_model_object, *m_map_object, m_options, ThreadSize());
#endif
}


} // namespace rhbm_gem

namespace rhbm_gem {
void PotentialAnalysisCommand::StudyAtomLocalFittingViaAlphaR(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<double> & alpha_list)
{
    ScopeTimer timer("PotentialAnalysisCommand::StudyAtomLocalFittingViaAlphaR");
    if (m_model_object == nullptr) return;

    auto atom_size{ atom_list.size() };
    if (atom_size == 0)
    {
        Logger::Log(
            LogLevel::Warning,
            "Skip Alpha_R bias study because no eligible atoms were selected.");
        return;
    }

    auto alpha_size{ static_cast<int>(alpha_list.size()) };
    std::atomic<size_t> atom_count{ 0 };
    Eigen::MatrixXd gaus_bias_matrix{ Eigen::MatrixXd::Zero(3, alpha_size) };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(ThreadSize())
#endif
    for (size_t i = 0; i < atom_size; i++)
    {
        auto local_entry{ atom_list[i]->GetLocalPotentialEntry() };
        const auto & data_entry_list{ local_entry->GetBasisAndResponseEntryList() };
        const auto dataset{ HRLDataTransform::BuildMemberDataset(data_entry_list) };
        const auto algorithm_options{
            MakePotentialAnalysisExecutionOptions(ThreadSize(), true)
        };

        Eigen::MatrixXd local_bias_array{ Eigen::MatrixXd::Zero(3, alpha_size) };
        for (int j = 0; j < alpha_size; j++)
        {
            auto alpha_r{ alpha_list[static_cast<size_t>(j)] };
            const auto result = HRLModelAlgorithms::EstimateBetaMDPDE(
                alpha_r,
                dataset.X,
                dataset.y,
                algorithm_options
            );
            Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
            auto gaus_ols{
                GausLinearTransformHelper::BuildGaus3DModel(result.beta_ols, model_par_init)
            };
            auto gaus_mdpde{
                GausLinearTransformHelper::BuildGaus3DModel(result.beta_mdpde, model_par_init)
            };
            local_bias_array.col(j) = (gaus_mdpde - gaus_ols).array().abs();
        }
        
#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            gaus_bias_matrix += local_bias_array;
            atom_count++;
            Logger::ProgressPercent(atom_count, atom_size);
        }
    }
    gaus_bias_matrix /= static_cast<double>(atom_size);

    const bool report_emitted{
        detail::EmitTrainingReportIfRequested(
            gaus_bias_matrix,
            alpha_list,
            "#alpha_{r}",
            "Deviation with OLS",
            m_options,
            "alpha_r_bias.pdf")
    };
    if (!report_emitted)
    {
        Logger::Log(
            LogLevel::Debug,
            "Alpha_R bias report was skipped (set --training-report-dir to emit PDF output).");
    }

    std::ostringstream alpha_r_bias_stream;
    alpha_r_bias_stream << "Alpha_R bias matrix:\n" << gaus_bias_matrix;
    Logger::Log(LogLevel::Debug, alpha_r_bias_stream.str());
}

void PotentialAnalysisCommand::StudyAtomGroupFittingViaAlphaG(
    const std::vector<std::vector<AtomObject *>> & atom_list_set,
    const std::vector<double> & alpha_list)
{
    ScopeTimer timer("PotentialAnalysisCommand::StudyAtomGroupFittingViaAlphaG");
    if (m_model_object == nullptr) return;

    auto alpha_size{ static_cast<int>(alpha_list.size()) };
    auto group_size{ atom_list_set.size() };
    if (group_size == 0)
    {
        Logger::Log(
            LogLevel::Warning,
            "Skip Alpha_G bias study because no eligible groups were selected.");
        return;
    }

    std::atomic<size_t> group_count{ 0 };
    Eigen::MatrixXd gaus_bias_matrix{ Eigen::MatrixXd::Zero(3, alpha_size) };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(ThreadSize())
#endif
    for (size_t i = 0; i < group_size; i++)
    {
        auto & group_atom_list{ atom_list_set[i] };
        std::vector<Eigen::VectorXd> data_entry_list;
        data_entry_list.reserve(group_atom_list.size());
        for (auto atom : group_atom_list)
        {
            data_entry_list.emplace_back(
                atom->GetLocalPotentialEntry()->GetBetaEstimateMDPDE()
            );
        }
        const auto beta_matrix{ HRLDataTransform::BuildBetaMatrix(data_entry_list, true) };
        const auto algorithm_options{
            MakePotentialAnalysisExecutionOptions(ThreadSize(), true)
        };

        Eigen::MatrixXd local_bias_array{ Eigen::MatrixXd::Zero(3, alpha_size) };
        for (int j = 0; j < alpha_size; j++)
        {
            auto alpha_g{ alpha_list[static_cast<size_t>(j)] };
            const auto result =
                HRLModelAlgorithms::EstimateMuMDPDE(alpha_g, beta_matrix, algorithm_options);
            Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
            auto gaus_mean{
                GausLinearTransformHelper::BuildGaus3DModel(result.mu_mean, model_par_init)
            };
            auto gaus_mdpde{
                GausLinearTransformHelper::BuildGaus3DModel(result.mu_mdpde, model_par_init)
            };
            local_bias_array.col(j) = (gaus_mdpde - gaus_mean).array().abs();
        }
        
#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            gaus_bias_matrix += local_bias_array;
            group_count++;
            Logger::ProgressPercent(group_count, group_size);
        }
    }
    gaus_bias_matrix /= static_cast<double>(group_size);

    const bool report_emitted{
        detail::EmitTrainingReportIfRequested(
            gaus_bias_matrix,
            alpha_list,
            "#alpha_{g}",
            "Deviation with Mean",
            m_options,
            "alpha_g_bias.pdf")
    };
    if (!report_emitted)
    {
        Logger::Log(
            LogLevel::Debug,
            "Alpha_G bias report was skipped (set --training-report-dir to emit PDF output).");
    }
    std::ostringstream alpha_g_bias_stream;
    alpha_g_bias_stream << "Alpha_G bias matrix:\n" << gaus_bias_matrix;
    Logger::Log(LogLevel::Debug, alpha_g_bias_stream.str());
}

void PotentialAnalysisCommand::SaveDataObject()
{
    ScopeTimer timer("PotentialAnalysisCommand::SaveDataObject");
    if (m_model_object == nullptr) return;

    m_data_manager.SaveDataObject(m_model_key_tag, m_options.saved_key_tag);

    for (auto atom : m_model_object->GetSelectedAtomList())
    {
        auto entry{ atom->GetLocalPotentialEntry() };
        if (entry != nullptr)
        {
            entry->ClearDistanceAndMapValueList();
        }
    }
}


} // namespace rhbm_gem

namespace rhbm_gem {
namespace detail {

void RunAtomSamplingWorkflow(
    ModelObject & model_object,
    const MapObject & map_object,
    const PotentialAnalysisCommandOptions & options,
    int thread_size)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomMapValueSampling");
    auto sampler{ std::make_unique<SphereSampler>() };
    sampler->SetSamplingSize(static_cast<unsigned int>(options.sampling_size));
    sampler->SetDistanceRangeMinimum(options.sampling_range_min);
    sampler->SetDistanceRangeMaximum(options.sampling_range_max);
    sampler->Print();

    const auto & atom_list{ model_object.GetSelectedAtomList() };
    const auto atom_size{ atom_list.size() };
    size_t atom_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(thread_size)
    {
        #pragma omp for
        for (size_t i = 0; i < atom_size; i++)
        {
            auto atom{ atom_list[i] };
            auto entry{ atom->GetLocalPotentialEntry() };
            entry->AddDistanceAndMapValueList(
                SampleMapValues(map_object, *sampler, atom->GetPosition()));
            entry->AddBasisAndResponseEntryList(
                GausLinearTransformHelper::MapValueTransform(
                    entry->GetDistanceAndMapValueList(),
                    options.fit_range_min,
                    options.fit_range_max));
            if (!options.use_training_alpha)
            {
                entry->SetAlphaR(options.alpha_r);
            }
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
        auto entry{ atom->GetLocalPotentialEntry() };
        entry->AddDistanceAndMapValueList(
            SampleMapValues(map_object, *sampler, atom->GetPosition()));
        entry->AddBasisAndResponseEntryList(
            GausLinearTransformHelper::MapValueTransform(
                entry->GetDistanceAndMapValueList(),
                options.fit_range_min,
                options.fit_range_max));
        if (!options.use_training_alpha)
        {
            entry->SetAlphaR(options.alpha_r);
        }
        atom_count++;
        Logger::ProgressPercent(atom_count, atom_size);
    }
#endif
}

void RunAtomGroupingWorkflow(ModelObject & model_object)
{
    ScopeTimer timer("RunAtomGroupClassification");
    Logger::Log(LogLevel::Info, "Atom Classification Summary:");
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto group_potential_entry{ std::make_unique<GroupPotentialEntry>() };
        for (auto atom : model_object.GetSelectedAtomList())
        {
            auto group_key{ AtomClassifier::GetGroupKeyInClass(atom, class_key) };
            group_potential_entry->AddAtomObjectPtr(group_key, atom);
            group_potential_entry->InsertGroupKey(group_key);
        }
        const auto group_size{ group_potential_entry->GetGroupKeySet().size() };
        model_object.AddAtomGroupPotentialEntry(class_key, group_potential_entry);
        Logger::Log(
            LogLevel::Info,
            " - Class type: " + class_key + " include " + std::to_string(group_size)
                + " groups.");
    }
}

void RunLocalAtomFittingWorkflow(
    ModelObject & model_object,
    int thread_size,
    double universal_alpha_r)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunLocalAtomFitting");
    std::atomic<size_t> atom_count{ 0 };
    auto & selected_atom_list{ model_object.GetSelectedAtomList() };
    const auto selected_atom_size{ selected_atom_list.size() };
    Logger::Log(
        LogLevel::Info,
        "Run Local atom fitting for " + std::to_string(selected_atom_size) + " atoms.");
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        auto local_entry{ selected_atom_list[i]->GetLocalPotentialEntry() };
        auto & data_entry_list{ local_entry->GetBasisAndResponseEntryList() };
        const auto dataset{ HRLDataTransform::BuildMemberDataset(data_entry_list) };
        const auto result{
            HRLModelAlgorithms::EstimateBetaMDPDE(
                universal_alpha_r,
                dataset.X,
                dataset.y,
                MakePotentialAnalysisExecutionOptions(thread_size, true))
        };

        local_entry->SetBetaEstimateOLS(result.beta_ols);
        local_entry->SetBetaEstimateMDPDE(result.beta_mdpde);
        local_entry->SetSigmaSquare(result.sigma_square);
        local_entry->SetDataWeight(result.data_weight);
        local_entry->SetDataCovariance(result.data_covariance);

        Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
        model_par_init(0) = local_entry->GetMomentZeroEstimate();
        model_par_init(1) = local_entry->GetMomentTwoEstimate();
        auto gaus_ols{
            GausLinearTransformHelper::BuildGaus3DModel(result.beta_ols, model_par_init)
        };
        auto gaus_mdpde{
            GausLinearTransformHelper::BuildGaus3DModel(result.beta_mdpde, model_par_init)
        };
        local_entry->AddGausEstimateOLS(gaus_ols(0), gaus_ols(1));
        local_entry->AddGausEstimateMDPDE(gaus_mdpde(0), gaus_mdpde(1));

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            Logger::ProgressPercent(atom_count, selected_atom_size);
        }
    }
}

} // namespace detail
void PotentialAnalysisCommand::RunAtomMapValueSampling()
{
    if (m_map_object == nullptr || m_model_object == nullptr) return;
    detail::RunAtomSamplingWorkflow(*m_model_object, *m_map_object, m_options, ThreadSize());
}

void PotentialAnalysisCommand::RunAtomGroupClassification()
{
    if (m_map_object == nullptr || m_model_object == nullptr) return;
    detail::RunAtomGroupingWorkflow(*m_model_object);
}

void PotentialAnalysisCommand::RunAtomAlphaTraining()
{
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomAlphaTraining");
    if (m_map_object == nullptr) return;

    const size_t subset_size_alpha_r{ 5 };
    const auto ordered_alpha_r_list{ detail::BuildOrderedAlphaRTrainingList() };
    
    // Alpha_R Training
    auto & atom_list{ m_model_object->GetSelectedAtomList() };
    std::vector<AtomObject *> selected_atom_list;
    selected_atom_list.reserve(atom_list.size());
    for (auto & atom : atom_list)
    {
        if (atom->IsMainChainAtom() == false) continue;
        if (atom->GetLocalPotentialEntry()->GetBasisAndResponseEntryListSize() < 500) continue;
        selected_atom_list.emplace_back(atom);
    }
    selected_atom_list.shrink_to_fit();
    
    auto selected_atom_size{ selected_atom_list.size() };
    Logger::Log(LogLevel::Info,
        "Run Alpha_R Training with "+ std::to_string(selected_atom_size) +" atoms.");
    auto alpha_r{ TrainUniversalAlphaR(
        selected_atom_list, subset_size_alpha_r, ordered_alpha_r_list)
    };
    
    for (auto & atom : atom_list)
    {
        atom->GetLocalPotentialEntry()->SetAlphaR(alpha_r);
    }
    
    
    // Alpha_G Training
    RunLocalAtomFitting(alpha_r);
    const size_t subset_size_alpha_g{ 10 };
    const auto ordered_alpha_g_list{ detail::BuildOrderedAlphaGTrainingList() };
    
    auto component_group_potential_entry{
        m_model_object->GetAtomGroupPotentialEntry(ChemicalDataHelper::GetComponentAtomClassKey())
    };
    std::vector<std::vector<AtomObject *>> atom_list_set;
    atom_list_set.reserve(component_group_potential_entry->GetGroupKeySet().size());
    for (auto group_key : component_group_potential_entry->GetGroupKeySet())
    {
        auto & group_atom_list{ component_group_potential_entry->GetAtomObjectPtrList(group_key) };
        if (group_atom_list.size() < 10) continue;
        if (group_atom_list.front()->IsMainChainAtom() == false) continue;
        atom_list_set.emplace_back(group_atom_list);
    }

    auto selected_group_size{ atom_list_set.size() };
    Logger::Log(LogLevel::Info,
        "Run Alpha_G Training with "+ std::to_string(selected_group_size) +" groups.");
    auto alpha_g{ TrainUniversalAlphaG(
        atom_list_set, subset_size_alpha_g, ordered_alpha_g_list)
    };

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto group_potential_entry{ m_model_object->GetAtomGroupPotentialEntry(class_key) };
        const auto & group_key_set{ group_potential_entry->GetGroupKeySet() };
        for (auto group_key : group_key_set)
        {
            group_potential_entry->AddAlphaG(group_key, alpha_g);
        }
    }

    StudyAtomLocalFittingViaAlphaR(selected_atom_list, ordered_alpha_r_list);
    StudyAtomGroupFittingViaAlphaG(atom_list_set, ordered_alpha_g_list);
}

double PotentialAnalysisCommand::TrainUniversalAlphaR(
    const std::vector<AtomObject *> & atom_list,
    const size_t subset_size,
    const std::vector<double> & alpha_list)
{
    auto atom_size{ atom_list.size() };
    auto alpha_size{ static_cast<int>(alpha_list.size()) };
    std::atomic<size_t> atom_count{ 0 };
    Eigen::ArrayXd beta_error_sum_array{ Eigen::ArrayXd::Zero(alpha_size) };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(ThreadSize())
#endif
    for (size_t i = 0; i < atom_size; i++)
    {
        auto local_entry{ atom_list[i]->GetLocalPotentialEntry() };
        const auto & data_entry_list{ local_entry->GetBasisAndResponseEntryList() };
        auto error_array{
            HRLAlphaTrainer::EvaluateAlphaR(
                data_entry_list,
                subset_size,
                alpha_list,
                MakePotentialAnalysisExecutionOptions(ThreadSize(), true)
            )
        };
        
#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            beta_error_sum_array += error_array.array();
            atom_count++;
            Logger::ProgressPercent(atom_count, atom_size);
        }
    }

    int error_min_id;
    beta_error_sum_array.minCoeff(&error_min_id);
    auto alpha_r_error{ alpha_list.at(static_cast<size_t>(error_min_id)) };

    Logger::Log(
        LogLevel::Info,
        "Alpha_R Training Results Summary: minimum beta error sum alpha_r = "
        + std::to_string(alpha_r_error));

    return alpha_r_error;
}

double PotentialAnalysisCommand::TrainUniversalAlphaG(
    const std::vector<std::vector<AtomObject *>> & atom_list_set,
    const size_t subset_size,
    const std::vector<double> & alpha_list)
{
    auto alpha_size{ static_cast<int>(alpha_list.size()) };
    auto group_size{ atom_list_set.size() };
    std::atomic<size_t> group_count{ 0 };
    Eigen::ArrayXd mu_error_sum_array{ Eigen::ArrayXd::Zero(alpha_size) };

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(ThreadSize())
#endif
    for (size_t i = 0; i < group_size; i++)
    {
        auto & group_atom_list{ atom_list_set[i] };
        std::vector<Eigen::VectorXd> data_entry_list;
        data_entry_list.reserve(group_atom_list.size());
        for (auto atom : group_atom_list)
        {
            data_entry_list.emplace_back(
                atom->GetLocalPotentialEntry()->GetBetaEstimateMDPDE()
            );
        }

        auto error_array{
            HRLAlphaTrainer::EvaluateAlphaG(
                data_entry_list,
                subset_size,
                alpha_list,
                MakePotentialAnalysisExecutionOptions(ThreadSize(), true)
            )
        };
        
#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            mu_error_sum_array += error_array.array();
            group_count++;
            Logger::ProgressPercent(group_count, group_size);
        }
    }

    int error_min_id;
    mu_error_sum_array.minCoeff(&error_min_id);
    auto alpha_g_error{ alpha_list.at(static_cast<size_t>(error_min_id)) };

    Logger::Log(
        LogLevel::Info,
        "Alpha_G Training Results Summary: minimum mu error sum alpha_g = "
        + std::to_string(alpha_g_error));

    return alpha_g_error;
}

void PotentialAnalysisCommand::RunLocalAtomFitting(double universal_alpha_r)
{
    if (m_model_object == nullptr) return;
    detail::RunLocalAtomFittingWorkflow(
        *m_model_object,
        ThreadSize(),
        universal_alpha_r);
}


} // namespace rhbm_gem
