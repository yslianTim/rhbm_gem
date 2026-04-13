#include "PotentialAnalysisCommand.hpp"

#include "MapSampling.hpp"
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

HRLExecutionOptions MakePotentialAnalysisExecutionOptions(int thread_size, bool quiet_mode)
{
    HRLExecutionOptions execution_options;
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
    const auto report_path{ training_report_dir / std::filesystem::path{ report_file_name } };
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

} // namespace

namespace rhbm_gem {

struct PotentialAnalysisCommand::AtomSamplingOptions
{
    int sampling_size;
    double sampling_range_min;
    double sampling_range_max;
    double fit_range_min;
    double fit_range_max;
};

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
    CoerceScalar(
        request.simulated_map_resolution,
        kSimResolutionOption,
        [](double candidate)
        {
            return std::isfinite(candidate) && candidate >= 0.0;
        },
        0.0,
        LogLevel::Error,
        "Simulated map resolution must be a finite non-negative value.");
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
    CoerceScalar(
        request.sampling_size,
        kSamplingOption,
        [](int candidate) { return candidate > 0; },
        1500,
        LogLevel::Warning,
        "Sampling size must be positive, reset to default value = 1500");
    CoerceScalar(
        request.sampling_range_min,
        kSamplingMinOption,
        [](double candidate)
        {
            return std::isfinite(candidate) && candidate >= 0.0;
        },
        0.0,
        LogLevel::Error,
        "Minimum sampling range must be a finite non-negative value.");
    CoerceScalar(
        request.sampling_range_max,
        kSamplingMaxOption,
        [](double candidate)
        {
            return std::isfinite(candidate) && candidate >= 0.0;
        },
        1.5,
        LogLevel::Error,
        "Maximum sampling range must be a finite non-negative value.");
    CoerceScalar(
        request.sampling_height,
        kSamplingHeightOption,
        [](double candidate)
        {
            return std::isfinite(candidate) && candidate > 0.0;
        },
        0.1,
        LogLevel::Error,
        "Sampling height must be a finite positive value.");
    CoerceScalar(
        request.fit_range_min,
        kFitMinOption,
        [](double candidate)
        {
            return std::isfinite(candidate) && candidate >= 0.0;
        },
        0.0,
        LogLevel::Error,
        "Minimum fitting range must be a finite non-negative value.");
    CoerceScalar(
        request.fit_range_max,
        kFitMaxOption,
        [](double candidate)
        {
            return std::isfinite(candidate) && candidate >= 0.0;
        },
        1.0,
        LogLevel::Error,
        "Maximum fitting range must be a finite non-negative value.");
    CoerceScalar(
        request.alpha_r,
        kAlphaROption,
        [](double candidate)
        {
            return std::isfinite(candidate) && candidate > 0.0;
        },
        0.1,
        LogLevel::Error,
        "Alpha-R must be a finite positive value.");
    CoerceScalar(
        request.alpha_g,
        kAlphaGOption,
        [](double candidate)
        {
            return std::isfinite(candidate) && candidate > 0.0;
        },
        0.2,
        LogLevel::Error,
        "Alpha-G must be a finite positive value.");
}

bool PotentialAnalysisCommand::ExecuteImpl()
{
    const auto & request{ RequestOptions() };
    const AtomSamplingOptions atom_sampling_options{
        request.sampling_size,
        request.sampling_range_min,
        request.sampling_range_max,
        request.fit_range_min,
        request.fit_range_max
    };

    if (!BuildDataObject(request)) return false;
    if (m_model_object == nullptr || m_map_object == nullptr)
    {
        Logger::Log(
            LogLevel::Error,
            "PotentialAnalysisCommand invariant violated: model/map object missing after load.");
        return false;
    }

    auto & model_object{ *m_model_object };
    auto & map_object{ *m_map_object };
    RunMapObjectPreprocessing(map_object);
    RunModelObjectPreprocessing(model_object, request.asymmetry_flag);
    RunSamplingWorkflow(map_object, model_object, atom_sampling_options);
    RunGroupingWorkflow(model_object);
    RunFittingWorkflow(model_object, map_object, request);
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
        if (request.simulation_flag)
        {
            UpdateModelObjectForSimulation(
                *m_model_object,
                request.simulated_map_resolution);
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

void PotentialAnalysisCommand::RunFittingWorkflow(
    ModelObject & model_object,
    MapObject & map_object,
    const PotentialAnalysisRequest & request)
{
    if (request.training_alpha_flag)
    {
        RunAtomAlphaTraining(model_object, request.training_report_dir);
        return;
    }

    RunLocalFitting(model_object, request.alpha_r);
    RunAtomPotentialFitting(model_object, request.training_alpha_flag, request.alpha_g);
    RunExperimentalBondWorkflowIfEnabled(model_object, map_object, request);
}

void PotentialAnalysisCommand::UpdateModelObjectForSimulation(
    ModelObject & model_object,
    double simulated_map_resolution)
{
    if (simulated_map_resolution == 0.0)
    {
        Logger::Log(LogLevel::Warning,
            "[Warning] The resolution of input simulated map hasn't been set.\n"
            "          Please give the corresponding resolution value for this map.\n"
            "          (-r, --sim-resolution)");
    }
    model_object.SetEmdID("Simulation");
    model_object.SetResolution(simulated_map_resolution);
    model_object.SetResolutionMethod("Blurring Width");
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

    Logger::Log(LogLevel::Info,
        "Number of selected atom = "
            + std::to_string(model_object.GetSelectedAtomCount()));
    Logger::Log(LogLevel::Info,
        "Number of selected bond = "
            + std::to_string(model_object.GetSelectedBondCount()));
    if (model_object.GetNumberOfAtom() > 0 &&
        model_object.GetSelectedAtomCount() == 0)
    {
        Logger::Log(LogLevel::Warning,
            "No atoms are selected after symmetry filtering. "
            "The input CIF may miss usable _entity/_struct_asym metadata. "
            "Try '--asymmetry true' to bypass symmetry filtering.");
    }
}

void PotentialAnalysisCommand::RunAtomPotentialFitting(
    ModelObject & model_object,
    bool training_alpha_flag,
    double fallback_alpha_g)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomPotentialFitting");
    const int basis_size{ 2 };
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
                const auto local_entry{ local_entry_map.at(atom) };
                const auto & dataset{ local_entry.GetDataset() };
                const auto & fit_result{ local_entry.GetFitResult() };
                data_entry_list.emplace_back(dataset.basis_and_response_entry_list);
                beta_mdpde_list.emplace_back(fit_result.beta_mdpde);
                sigma_square_list.emplace_back(fit_result.sigma_square);
                data_weight_list.emplace_back(fit_result.data_weight);
                data_covariance_list.emplace_back(fit_result.data_covariance);
            }
            auto alpha_g{ training_alpha_flag ?
                model_object.GetAnalysisView().GetAtomAlphaG(group_key, class_key) :
                fallback_alpha_g
            };
            const auto input = HRLDataTransform::BuildGroupInput(
                basis_size,
                data_entry_list,
                beta_mdpde_list,
                sigma_square_list,
                data_weight_list,
                data_covariance_list
            );
            HRLGroupEstimator estimator(MakePotentialAnalysisExecutionOptions(ThreadSize(), true));
            const auto result{ estimator.Estimate(input, alpha_g) };

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
                auto atom_entry{ local_entry_map.at(atom) };
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
                GaussianPosterior posterior;
                posterior.estimate = GaussianEstimate{ posterior_estimate(0), posterior_estimate(1) };
                posterior.variance = GaussianEstimate{ posterior_variance(0), posterior_variance(1) };
                atom_entry.SetAnnotation(
                    class_key,
                    LocalPotentialAnnotationData{
                        posterior,
                        static_cast<bool>(result.outlier_flag_array(count)),
                        result.statistical_distance_array(count)
                    });
                count++;
            }

#ifdef USE_OPENMP
            #pragma omp critical
#endif
            {
                analysis.SetAtomGroupStatistics(
                    group_key,
                    class_key,
                    GroupPotentialStatistics{
                        GaussianEstimate{ gaus_group_mean(0), gaus_group_mean(1) },
                        GaussianEstimate{ gaus_group_mdpde(0), gaus_group_mdpde(1) },
                        GaussianEstimate{ prior_estimate(0), prior_estimate(1) },
                        GaussianEstimate{ prior_variance(0), prior_variance(1) },
                        alpha_g
                    });
                key_count++;
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

void PotentialAnalysisCommand::RunExperimentalBondWorkflowIfEnabled(
    ModelObject & model_object,
    MapObject & map_object,
    const PotentialAnalysisRequest & request)
{
#ifdef RHBM_GEM_ENABLE_EXPERIMENTAL_FEATURE
    experimental::RunPotentialAnalysisBondWorkflow(
        model_object, map_object, request, ThreadSize());
#else
    (void)model_object;
    (void)map_object;
    (void)request;
#endif
}

void PotentialAnalysisCommand::StudyAtomLocalFittingViaAlphaR(
    ModelObject & model_object,
    const std::vector<AtomObject *> & atom_list,
    const std::vector<double> & alpha_list,
    const std::filesystem::path & training_report_dir)
{
    ScopeTimer timer("PotentialAnalysisCommand::StudyAtomLocalFittingViaAlphaR");
    auto analysis{ model_object.EditAnalysis() };
    std::unordered_map<const AtomObject *, MutableLocalPotentialView> local_entry_map;
    local_entry_map.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        local_entry_map.emplace(atom, analysis.EnsureAtomLocalPotential(*atom));
    }

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
        const auto local_entry{ local_entry_map.at(atom_list[i]) };
        const auto & data_entry_list{ local_entry.GetDataset().basis_and_response_entry_list };
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
        EmitTrainingReportIfRequested(
            gaus_bias_matrix,
            alpha_list,
            "#alpha_{r}",
            "Deviation with OLS",
            training_report_dir,
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
    ModelObject & model_object,
    const std::vector<std::vector<AtomObject *>> & atom_list_set,
    const std::vector<double> & alpha_list,
    const std::filesystem::path & training_report_dir)
{
    ScopeTimer timer("PotentialAnalysisCommand::StudyAtomGroupFittingViaAlphaG");
    auto analysis{ model_object.EditAnalysis() };
    std::unordered_map<const AtomObject *, MutableLocalPotentialView> local_entry_map;
    for (const auto & group_atom_list : atom_list_set)
    {
        for (auto * atom : group_atom_list)
        {
            local_entry_map.try_emplace(atom, analysis.EnsureAtomLocalPotential(*atom));
        }
    }

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
                local_entry_map.at(atom).GetFitResult().beta_mdpde
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
        EmitTrainingReportIfRequested(
            gaus_bias_matrix,
            alpha_list,
            "#alpha_{g}",
            "Deviation with Mean",
            training_report_dir,
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
    const AtomSamplingOptions & options)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunSamplingWorkflow");
    auto thread_size{ ThreadSize() };
    SphereSampler sampler;
    sampler.SetSamplingProfile(
        //SphereSamplingProfile::RadiusUniformRandom(
        //    SphereDistanceRange{ options.sampling_range_min, options.sampling_range_max },
        //    static_cast<unsigned int>(options.sampling_size))
        SphereSamplingProfile::FibonacciDeterministic(
            SphereDistanceRange{ options.sampling_range_min, options.sampling_range_max },
            0.05,
            static_cast<unsigned int>(options.sampling_size))
    );
    sampler.Print();

    const auto & atom_list{ model_object.GetSelectedAtoms() };
    const auto atom_size{ atom_list.size() };
    size_t atom_count{ 0 };
    auto analysis{ model_object.EditAnalysis() };
    std::vector<MutableLocalPotentialView> local_entry_list;
    local_entry_list.reserve(atom_size);
    for (auto * atom : atom_list)
    {
        local_entry_list.emplace_back(analysis.EnsureAtomLocalPotential(*atom));
    }

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(thread_size)
    {
        #pragma omp for
        for (size_t i = 0; i < atom_size; i++)
        {
            auto atom{ atom_list[i] };
            auto entry{ local_entry_list[i] };
            auto sampling_entries{
                SampleMapValues(map_object, sampler, atom->GetPosition())
            };
            entry.SetSamplingEntries(sampling_entries);
            entry.SetDataset(LocalPotentialDataset{
                GausLinearTransformHelper::MapValueTransform(
                    sampling_entries,
                    options.fit_range_min,
                    options.fit_range_max)
            });
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
            SampleMapValues(map_object, sampler, atom->GetPosition())
        };
        entry.SetSamplingEntries(sampling_entries);
        entry.SetDataset(LocalPotentialDataset{
            GausLinearTransformHelper::MapValueTransform(
                sampling_entries,
                options.fit_range_min,
                options.fit_range_max)
        });
        atom_count++;
        Logger::ProgressPercent(atom_count, atom_size);
    }
#endif
}

void PotentialAnalysisCommand::RunGroupingWorkflow(ModelObject & model_object)
{
    ScopeTimer timer("RunGroupingWorkflow");
    Logger::Log(LogLevel::Info, "Atom Grouping Summary:");
    auto analysis{ model_object.EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    const auto analysis_view{ model_object.GetAnalysisView() };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        const auto group_size{ analysis_view.CollectAtomGroupKeys(class_key).size() };
        Logger::Log(
            LogLevel::Info,
            " - Class type: " + class_key + " include " + std::to_string(group_size) + " groups.");
    }
}

void PotentialAnalysisCommand::RunAtomAlphaTraining(
    ModelObject & model_object,
    const std::filesystem::path & training_report_dir)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunAtomAlphaTraining");
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };

    const size_t subset_size_alpha_r{ 5 };
    const auto ordered_alpha_r_list{ BuildOrderedAlphaRTrainingList() };
    
    // Alpha_R Training
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    std::vector<AtomObject *> selected_atom_list;
    selected_atom_list.reserve(atom_list.size());
    for (auto & atom : atom_list)
    {
        if (atom->IsMainChainAtom() == false) continue;
        const auto local_entry{ analysis.EnsureAtomLocalPotential(*atom) };
        if (!local_entry.HasDataset() ||
            local_entry.GetDataset().basis_and_response_entry_list.size() < 500) continue;
        selected_atom_list.emplace_back(atom);
    }
    selected_atom_list.shrink_to_fit();
    
    auto selected_atom_size{ selected_atom_list.size() };
    Logger::Log(LogLevel::Info,
        "Run Alpha_R Training with "+ std::to_string(selected_atom_size) +" atoms.");
    auto alpha_r{ TrainUniversalAlphaR(
        model_object, selected_atom_list, subset_size_alpha_r, ordered_alpha_r_list)
    };
    
    // Alpha_G Training
    RunLocalFitting(model_object, alpha_r);
    const size_t subset_size_alpha_g{ 10 };
    const auto ordered_alpha_g_list{ BuildOrderedAlphaGTrainingList() };
    
    std::vector<std::vector<AtomObject *>> atom_list_set;
    const auto component_class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    const auto component_group_keys{ analysis_view.CollectAtomGroupKeys(component_class_key) };
    atom_list_set.reserve(component_group_keys.size());
    for (const auto group_key : component_group_keys)
    {
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(group_key, component_class_key) };
        if (group_atom_list.size() < 10) continue;
        if (group_atom_list.front()->IsMainChainAtom() == false) continue;
        atom_list_set.emplace_back(group_atom_list);
    }

    auto selected_group_size{ atom_list_set.size() };
    Logger::Log(LogLevel::Info,
        "Run Alpha_G Training with "+ std::to_string(selected_group_size) +" groups.");
    auto alpha_g{ TrainUniversalAlphaG(
        model_object,
        atom_list_set, subset_size_alpha_g, ordered_alpha_g_list)
    };

    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        for (const auto group_key : analysis_view.CollectAtomGroupKeys(class_key))
        {
            analysis.SetAtomGroupAlphaG(group_key, class_key, alpha_g);
        }
    }

    StudyAtomLocalFittingViaAlphaR(
        model_object,
        selected_atom_list,
        ordered_alpha_r_list,
        training_report_dir);
    StudyAtomGroupFittingViaAlphaG(
        model_object,
        atom_list_set,
        ordered_alpha_g_list,
        training_report_dir);
}

double PotentialAnalysisCommand::TrainUniversalAlphaR(
    ModelObject & model_object,
    const std::vector<AtomObject *> & atom_list,
    const size_t subset_size,
    const std::vector<double> & alpha_list)
{
    auto atom_size{ atom_list.size() };
    auto alpha_size{ static_cast<int>(alpha_list.size()) };
    std::atomic<size_t> atom_count{ 0 };
    Eigen::ArrayXd beta_error_sum_array{ Eigen::ArrayXd::Zero(alpha_size) };
    auto analysis{ model_object.EditAnalysis() };
    std::unordered_map<const AtomObject *, MutableLocalPotentialView> local_entry_map;
    local_entry_map.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        local_entry_map.emplace(atom, analysis.EnsureAtomLocalPotential(*atom));
    }

#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(ThreadSize())
#endif
    for (size_t i = 0; i < atom_size; i++)
    {
        const auto local_entry{ local_entry_map.at(atom_list[i]) };
        const auto & data_entry_list{ local_entry.GetDataset().basis_and_response_entry_list };
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
    ModelObject & model_object,
    const std::vector<std::vector<AtomObject *>> & atom_list_set,
    const size_t subset_size,
    const std::vector<double> & alpha_list)
{
    auto alpha_size{ static_cast<int>(alpha_list.size()) };
    auto group_size{ atom_list_set.size() };
    std::atomic<size_t> group_count{ 0 };
    Eigen::ArrayXd mu_error_sum_array{ Eigen::ArrayXd::Zero(alpha_size) };
    auto analysis{ model_object.EditAnalysis() };
    std::unordered_map<const AtomObject *, MutableLocalPotentialView> local_entry_map;
    for (const auto & group_atom_list : atom_list_set)
    {
        for (auto * atom : group_atom_list)
        {
            local_entry_map.try_emplace(atom, analysis.EnsureAtomLocalPotential(*atom));
        }
    }

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
                local_entry_map.at(atom).GetFitResult().beta_mdpde
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

void PotentialAnalysisCommand::RunLocalFitting(
    ModelObject & model_object,
    double universal_alpha_r)
{
    ScopeTimer timer("PotentialAnalysisCommand::RunLocalAtomFitting");
    auto thread_size{ ThreadSize() };
    std::atomic<size_t> atom_count{ 0 };
    const auto & selected_atom_list{ model_object.GetSelectedAtoms() };
    const auto selected_atom_size{ selected_atom_list.size() };
    auto analysis{ model_object.EditAnalysis() };
    std::vector<MutableLocalPotentialView> local_entry_list;
    local_entry_list.reserve(selected_atom_size);
    for (auto * atom : selected_atom_list)
    {
        local_entry_list.emplace_back(analysis.EnsureAtomLocalPotential(*atom));
        analysis.EnsureAtomLocalPotential(*atom).SetAlphaR(universal_alpha_r);
    }
    Logger::Log(
        LogLevel::Info,
        "Run Local atom fitting for " + std::to_string(selected_atom_size) + " atoms.");
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        auto local_entry{ local_entry_list[i] };
        const auto & data_entry_list{ local_entry.GetDataset().basis_and_response_entry_list };
        const auto dataset{ HRLDataTransform::BuildMemberDataset(data_entry_list) };
        const auto result{
            HRLModelAlgorithms::EstimateBetaMDPDE(
                universal_alpha_r,
                dataset.X,
                dataset.y,
                MakePotentialAnalysisExecutionOptions(thread_size, true))
        };

        local_entry.SetFitResult(LocalPotentialFitResult{
            result.beta_ols,
            result.beta_mdpde,
            result.sigma_square,
            result.data_weight,
            result.data_covariance
        });

        Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
        model_par_init(0) = LocalPotentialView::RequireFor(*selected_atom_list[i]).GetMomentZeroEstimate();
        model_par_init(1) = LocalPotentialView::RequireFor(*selected_atom_list[i]).GetMomentTwoEstimate();
        auto gaus_ols{
            GausLinearTransformHelper::BuildGaus3DModel(result.beta_ols, model_par_init)
        };
        auto gaus_mdpde{
            GausLinearTransformHelper::BuildGaus3DModel(result.beta_mdpde, model_par_init)
        };
        local_entry.SetEstimates(LocalPotentialEstimates{
            GaussianEstimate{ gaus_ols(0), gaus_ols(1) },
            GaussianEstimate{ gaus_mdpde(0), gaus_mdpde(1) }
        });

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
