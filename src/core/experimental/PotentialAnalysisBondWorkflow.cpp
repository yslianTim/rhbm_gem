#include "experimental/PotentialAnalysisBondWorkflow.hpp"

#include "command/MapSampling.hpp"
#include "command/PotentialAnalysisCommand.hpp"
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/hrl/HRLDataTransform.hpp>
#include <rhbm_gem/utils/hrl/HRLGroupEstimator.hpp>
#include <rhbm_gem/utils/hrl/HRLModelAlgorithms.hpp>
#include <rhbm_gem/utils/math/CylinderSampler.hpp>
#include <rhbm_gem/utils/math/GausLinearTransformHelper.hpp>

#include <array>
#include <atomic>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::experimental {

namespace {

struct PotentialAnalysisBondWorkflowContext
{
    ModelObject & model_object;
    MapObject & map_object;
    const PotentialAnalysisRequest & options;
    int thread_size;
};

HRLExecutionOptions MakePotentialAnalysisExecutionOptions(
    int thread_size,
    bool quiet_mode)
{
    HRLExecutionOptions execution_options;
    execution_options.quiet_mode = quiet_mode;
    execution_options.thread_size = thread_size;
    return execution_options;
}

void RunBondSampling(
    ModelObject & model_object,
    const MapObject & map_object,
    const PotentialAnalysisRequest & options,
    int thread_size)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunBondMapValueSampling");
    CylinderSampler sampler;
    sampler.SetSampleCount(static_cast<unsigned int>(options.sampling_size));
    sampler.SetDistanceRangeMinimum(options.sampling_range_min);
    sampler.SetDistanceRangeMaximum(options.sampling_range_max);
    sampler.SetHeight(options.sampling_height);
    sampler.Print();

    const auto & bond_list{ model_object.GetSelectedBonds() };
    const auto bond_size{ bond_list.size() };
    size_t bond_count{ 0 };
    auto analysis{ model_object.EditAnalysis() };
    std::vector<MutableLocalPotentialView> local_entry_list;
    local_entry_list.reserve(bond_size);
    for (auto * bond : bond_list)
    {
        local_entry_list.emplace_back(analysis.EnsureBondLocalPotential(*bond));
    }

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(thread_size)
    {
        #pragma omp for
        for (size_t i = 0; i < bond_size; i++)
        {
            auto bond{ bond_list[i] };
            auto entry{ local_entry_list[i] };
            auto bond_vector{ bond->GetBondVector() };
            auto bond_position{ bond->GetPosition() };
            constexpr float adjusted_rate{ 0.0f };
            std::array<float, 3> adjusted_position{
                bond_position[0] + 0.5f * bond_vector[0] * adjusted_rate,
                bond_position[1] + 0.5f * bond_vector[1] * adjusted_rate,
                bond_position[2] + 0.5f * bond_vector[2] * adjusted_rate
            };
            auto distance_and_map_value_list{
                SampleMapValues(
                    map_object,
                    sampler,
                    adjusted_position,
                    bond_vector)
            };
            entry.SetDistanceAndMapValueList(distance_and_map_value_list);
            entry.SetDataset(LocalPotentialDataset{
                GausLinearTransformHelper::MapValueTransform(
                    distance_and_map_value_list,
                    options.fit_range_min,
                    options.fit_range_max)
            });
            entry.SetAlphaR(options.alpha_r);
            #pragma omp critical
            {
                bond_count++;
                Logger::ProgressPercent(bond_count, bond_size);
            }
        }
    }
#else
    for (size_t i = 0; i < bond_size; i++)
    {
        auto bond{ bond_list[i] };
        auto entry{ local_entry_list[i] };
        auto distance_and_map_value_list{
            SampleMapValues(
                map_object,
                sampler,
                bond->GetPosition(),
                bond->GetBondVector())
        };
        entry.SetDistanceAndMapValueList(distance_and_map_value_list);
        entry.SetDataset(LocalPotentialDataset{
            GausLinearTransformHelper::MapValueTransform(
                distance_and_map_value_list,
                options.fit_range_min,
                options.fit_range_max)
        });
        entry.SetAlphaR(options.alpha_r);
        bond_count++;
        Logger::ProgressPercent(bond_count, bond_size);
    }
#endif
}

void RunBondGrouping(ModelObject & model_object)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunBondGroupClassification");
    Logger::Log(LogLevel::Info, "Bond Classification Summary:");
    auto analysis{ model_object.EditAnalysis() };
    analysis.RebuildBondGroupsFromSelection();
    const auto analysis_view{ model_object.GetAnalysisView() };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        const auto group_size{ analysis_view.CollectBondGroupKeys(class_key).size() };
        Logger::Log(
            LogLevel::Info,
            " - Class type: " + class_key + " include "
                + std::to_string(group_size) + " groups.");
    }
}

void RunBondMapValueSampling(const PotentialAnalysisBondWorkflowContext & context)
{
    RunBondSampling(context.model_object, context.map_object, context.options, context.thread_size);
}

void RunBondGroupClassification(const PotentialAnalysisBondWorkflowContext & context)
{
    RunBondGrouping(context.model_object);
}

void RunLocalBondFitting(
    const PotentialAnalysisBondWorkflowContext & context,
    double universal_alpha_r)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunLocalBondFitting");

    std::atomic<size_t> bond_count{ 0 };
    const auto & selected_bond_list{ context.model_object.GetSelectedBonds() };
    const auto selected_bond_size{ selected_bond_list.size() };
    auto analysis{ context.model_object.EditAnalysis() };
    std::vector<MutableLocalPotentialView> local_entry_list;
    local_entry_list.reserve(selected_bond_size);
    for (auto * bond : selected_bond_list)
    {
        local_entry_list.emplace_back(analysis.EnsureBondLocalPotential(*bond));
    }
    Logger::Log(
        LogLevel::Info,
        "Run Local bond fitting for " + std::to_string(selected_bond_size) + " bonds.");
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(context.thread_size)
#endif
    for (size_t i = 0; i < selected_bond_size; i++)
    {
        auto local_entry{ local_entry_list[i] };
        const auto & data_entry_list{ local_entry.GetDataset().basis_and_response_entry_list };
        const auto dataset{ HRLDataTransform::BuildMemberDataset(data_entry_list) };
        const auto result{
            HRLModelAlgorithms::EstimateBetaMDPDE(
                universal_alpha_r,
                dataset.X,
                dataset.y,
                MakePotentialAnalysisExecutionOptions(context.thread_size, true))
        };

        local_entry.SetFitResult(LocalPotentialFitResult{
            result.beta_ols,
            result.beta_mdpde,
            result.sigma_square,
            result.data_weight,
            result.data_covariance
        });

        Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
        model_par_init(0) =
            LocalPotentialView::RequireFor(*selected_bond_list[i]).GetMomentZeroEstimate();
        model_par_init(1) =
            LocalPotentialView::RequireFor(*selected_bond_list[i]).GetMomentTwoEstimate();
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
            bond_count++;
            Logger::ProgressPercent(bond_count, selected_bond_size);
        }
    }
}

void RunBondPotentialFitting(const PotentialAnalysisBondWorkflowContext & context)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunBondPotentialFitting");
    constexpr int basis_size{ 2 };
    auto analysis{ context.model_object.EditAnalysis() };
    const auto analysis_view{ context.model_object.GetAnalysisView() };
    std::unordered_map<const BondObject *, MutableLocalPotentialView> local_entry_map;
    local_entry_map.reserve(context.model_object.GetSelectedBondCount());
    for (auto * bond : context.model_object.GetSelectedBonds())
    {
        local_entry_map.emplace(bond, analysis.EnsureBondLocalPotential(*bond));
    }
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        Logger::Log(LogLevel::Info, "Class type: " + class_key);

        auto group_keys{ analysis_view.CollectBondGroupKeys(class_key) };
        const auto group_key_size{ group_keys.size() };
        std::atomic<size_t> key_count{ 0 };

#ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(context.thread_size)
#endif
        for (size_t idx = 0; idx < group_key_size; idx++)
        {
            auto group_key{ group_keys[idx] };
            const auto & bond_list{ analysis_view.GetBondObjectList(group_key, class_key) };
            const auto group_size{ bond_list.size() };
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
            for (const auto & bond : bond_list)
            {
                const auto local_entry{ local_entry_map.at(bond) };
                const auto & dataset{ local_entry.GetDataset() };
                const auto & fit_result{ local_entry.GetFitResult() };
                data_entry_list.emplace_back(dataset.basis_and_response_entry_list);
                beta_mdpde_list.emplace_back(fit_result.beta_mdpde);
                sigma_square_list.emplace_back(fit_result.sigma_square);
                data_weight_list.emplace_back(fit_result.data_weight);
                data_covariance_list.emplace_back(fit_result.data_covariance);
            }
            const auto input{
                HRLDataTransform::BuildGroupInput(
                    basis_size,
                    data_entry_list,
                    beta_mdpde_list,
                    sigma_square_list,
                    data_weight_list,
                    data_covariance_list)
            };
            HRLGroupEstimator estimator(
                MakePotentialAnalysisExecutionOptions(context.thread_size, true));
            const auto result{ estimator.Estimate(input, context.options.alpha_g) };

            auto gaus_group_mean{
                GausLinearTransformHelper::BuildGaus2DModel(result.mu_mean)
            };
            auto gaus_group_mdpde{
                GausLinearTransformHelper::BuildGaus2DModel(result.mu_mdpde)
            };
            auto gaus_prior{
                GausLinearTransformHelper::BuildGaus2DModelWithVariance(
                    result.mu_prior,
                    result.capital_lambda)
            };
            auto prior_estimate{ std::get<0>(gaus_prior) };
            auto prior_variance{ std::get<1>(gaus_prior) };

            auto count{ 0 };
            for (const auto & bond : bond_list)
            {
                auto bond_entry{ local_entry_map.at(bond) };
                const auto beta_vector_posterior{
                    result.beta_posterior_array.col(static_cast<Eigen::Index>(count))
                };
                const auto & sigma_matrix_posterior{
                    result.capital_sigma_posterior_list.at(static_cast<std::size_t>(count))
                };
                auto gaus_posterior{
                    GausLinearTransformHelper::BuildGaus2DModelWithVariance(
                        beta_vector_posterior,
                        sigma_matrix_posterior)
                };
                auto posterior_estimate{ std::get<0>(gaus_posterior) };
                auto posterior_variance{ std::get<1>(gaus_posterior) };
                GaussianPosterior posterior;
                posterior.estimate = GaussianEstimate{ posterior_estimate(0), posterior_estimate(1) };
                posterior.variance = GaussianEstimate{ posterior_variance(0), posterior_variance(1) };
                bond_entry.SetAnnotation(
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
                analysis.SetBondGroupStatistics(
                    group_key,
                    class_key,
                    GroupPotentialStatistics{
                        GaussianEstimate{ gaus_group_mean(0), gaus_group_mean(1) },
                        GaussianEstimate{ gaus_group_mdpde(0), gaus_group_mdpde(1) },
                        GaussianEstimate{ prior_estimate(0), prior_estimate(1) },
                        GaussianEstimate{ prior_variance(0), prior_variance(1) },
                        context.options.alpha_g
                    });
                key_count++;
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

void RunPotentialAnalysisBondWorkflowImpl(
    const PotentialAnalysisBondWorkflowContext & context)
{
    RunBondMapValueSampling(context);
    RunBondGroupClassification(context);
    RunLocalBondFitting(context, context.options.alpha_r);
    RunBondPotentialFitting(context);
}

} // namespace

void RunPotentialAnalysisBondWorkflow(
    ModelObject & model_object,
    MapObject & map_object,
    const PotentialAnalysisRequest & options,
    int thread_size)
{
    PotentialAnalysisBondWorkflowContext context{
        model_object,
        map_object,
        options,
        thread_size,
    };
    RunPotentialAnalysisBondWorkflowImpl(context);
}

} // namespace rhbm_gem::experimental
