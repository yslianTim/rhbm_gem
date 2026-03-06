#include "PotentialAnalysisBondWorkflow.hpp"
#include "PotentialAnalysisExecutionOptions.hpp"

#include "BondClassifier.hpp"
#include "BondObject.hpp"
#include "ChemicalDataHelper.hpp"
#include "CylinderSampler.hpp"
#include "GausLinearTransformHelper.hpp"
#include "GroupPotentialEntry.hpp"
#include "HRLDataTransform.hpp"
#include "HRLGroupEstimator.hpp"
#include "HRLModelAlgorithms.hpp"
#include "LocalPotentialEntry.hpp"
#include "Logger.hpp"
#include "MapObject.hpp"
#include "MapSampling.hpp"
#include "ModelObject.hpp"
#include "ScopeTimer.hpp"

#include <array>
#include <atomic>
#include <memory>
#include <tuple>
#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::detail {

namespace {

void RunBondSampling(
    ModelObject & model_object,
    const MapObject & map_object,
    const PotentialAnalysisCommandOptions & options)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunBondMapValueSampling");
    auto sampler{ std::make_unique<CylinderSampler>() };
    sampler->SetSamplingSize(static_cast<unsigned int>(options.sampling_size));
    sampler->SetDistanceRangeMinimum(options.sampling_range_min);
    sampler->SetDistanceRangeMaximum(options.sampling_range_max);
    sampler->SetHeight(options.sampling_height);
    sampler->Print();

    const auto & bond_list{ model_object.GetSelectedBondList() };
    const auto bond_size{ bond_list.size() };
    size_t bond_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel num_threads(options.thread_size)
    {
        #pragma omp for
        for (size_t i = 0; i < bond_size; i++)
        {
            auto bond{ bond_list[i] };
            auto entry{ bond->GetLocalPotentialEntry() };
            auto bond_vector{ bond->GetBondVector() };
            auto bond_position{ bond->GetPosition() };
            constexpr float adjusted_rate{ 0.0f };
            std::array<float, 3> adjusted_position{
                bond_position[0] + 0.5f * bond_vector[0] * adjusted_rate,
                bond_position[1] + 0.5f * bond_vector[1] * adjusted_rate,
                bond_position[2] + 0.5f * bond_vector[2] * adjusted_rate
            };
            entry->AddDistanceAndMapValueList(
                SampleMapValues(
                    map_object,
                    *sampler,
                    adjusted_position,
                    bond_vector));
            entry->AddBasisAndResponseEntryList(
                GausLinearTransformHelper::MapValueTransform(
                    entry->GetDistanceAndMapValueList(),
                    options.fit_range_min,
                    options.fit_range_max));
            entry->SetAlphaR(options.alpha_r);
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
        auto entry{ bond->GetLocalPotentialEntry() };
        entry->AddDistanceAndMapValueList(
            SampleMapValues(
                map_object,
                *sampler,
                bond->GetPosition(),
                bond->GetBondVector()));
        entry->AddBasisAndResponseEntryList(
            GausLinearTransformHelper::MapValueTransform(
                entry->GetDistanceAndMapValueList(),
                options.fit_range_min,
                options.fit_range_max));
        entry->SetAlphaR(options.alpha_r);
        bond_count++;
        Logger::ProgressPercent(bond_count, bond_size);
    }
#endif
}

void RunBondGrouping(ModelObject & model_object)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunBondGroupClassification");
    Logger::Log(LogLevel::Info, "Bond Classification Summary:");
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        auto group_potential_entry{ std::make_unique<GroupPotentialEntry>() };
        for (auto bond : model_object.GetSelectedBondList())
        {
            auto group_key{ BondClassifier::GetGroupKeyInClass(bond, class_key) };
            group_potential_entry->AddBondObjectPtr(group_key, bond);
            group_potential_entry->InsertGroupKey(group_key);
        }
        const auto group_size{ group_potential_entry->GetGroupKeySet().size() };
        model_object.AddBondGroupPotentialEntry(class_key, group_potential_entry);
        Logger::Log(
            LogLevel::Info,
            " - Class type: " + class_key + " include "
                + std::to_string(group_size) + " groups.");
    }
}

void RunBondMapValueSampling(const PotentialAnalysisBondWorkflowContext & context)
{
    RunBondSampling(context.model_object, context.map_object, context.options);
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
    auto & selected_bond_list{ context.model_object.GetSelectedBondList() };
    const auto selected_bond_size{ selected_bond_list.size() };
    Logger::Log(
        LogLevel::Info,
        "Run Local bond fitting for " + std::to_string(selected_bond_size) + " bonds.");
#ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(context.options.thread_size)
#endif
    for (size_t i = 0; i < selected_bond_size; i++)
    {
        auto local_entry{ selected_bond_list[i]->GetLocalPotentialEntry() };
        auto & data_entry_list{ local_entry->GetBasisAndResponseEntryList() };
        const auto dataset{ HRLDataTransform::BuildMemberDataset(data_entry_list) };
        const auto result{
            HRLModelAlgorithms::EstimateBetaMDPDE(
                universal_alpha_r,
                dataset.X,
                dataset.y,
                MakePotentialAnalysisExecutionOptions(context.options, true))
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
            bond_count++;
            Logger::ProgressPercent(bond_count, selected_bond_size);
        }
    }
}

void RunBondPotentialFitting(const PotentialAnalysisBondWorkflowContext & context)
{
    ScopeTimer timer("PotentialAnalysisBondWorkflow::RunBondPotentialFitting");
    constexpr int basis_size{ 2 };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupBondClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupBondClassKey(i) };
        Logger::Log(LogLevel::Info, "Class type: " + class_key);

        auto group_potential_entry{ context.model_object.GetBondGroupPotentialEntry(class_key) };
        const auto & key_set{ group_potential_entry->GetGroupKeySet() };
        std::vector<GroupKey> group_keys(key_set.begin(), key_set.end());
        const auto group_key_size{ group_keys.size() };
        std::atomic<size_t> key_count{ 0 };

#ifdef USE_OPENMP
        #pragma omp parallel for schedule(dynamic) num_threads(context.options.thread_size)
#endif
        for (size_t idx = 0; idx < group_key_size; idx++)
        {
            auto group_key{ group_keys[idx] };
            const auto & bond_list{ group_potential_entry->GetBondObjectPtrList(group_key) };
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
                auto entry{ bond->GetLocalPotentialEntry() };
                data_entry_list.emplace_back(entry->GetBasisAndResponseEntryList());
                beta_mdpde_list.emplace_back(entry->GetBetaEstimateMDPDE());
                sigma_square_list.emplace_back(entry->GetSigmaSquare());
                data_weight_list.emplace_back(entry->GetDataWeight());
                data_covariance_list.emplace_back(entry->GetDataCovariance());
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
                MakePotentialAnalysisExecutionOptions(context.options, true));
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
                auto bond_entry{ bond->GetLocalPotentialEntry() };
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
                bond_entry->AddGausEstimatePosterior(
                    class_key, posterior_estimate(0), posterior_estimate(1));
                bond_entry->AddGausVariancePosterior(
                    class_key, posterior_variance(0), posterior_variance(1));
                bond_entry->AddOutlierTag(class_key, result.outlier_flag_array(count));
                bond_entry->AddStatisticalDistance(
                    class_key,
                    result.statistical_distance_array(count));
                count++;
            }

#ifdef USE_OPENMP
            #pragma omp critical
#endif
            {
                group_potential_entry->AddGausEstimateMean(
                    group_key, gaus_group_mean(0), gaus_group_mean(1));
                group_potential_entry->AddGausEstimateMDPDE(
                    group_key, gaus_group_mdpde(0), gaus_group_mdpde(1));
                group_potential_entry->AddGausEstimatePrior(
                    group_key, prior_estimate(0), prior_estimate(1));
                group_potential_entry->AddGausVariancePrior(
                    group_key, prior_variance(0), prior_variance(1));
                group_potential_entry->AddAlphaG(group_key, context.options.alpha_g);
                key_count++;
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

} // namespace

void RunPotentialAnalysisBondWorkflow(const PotentialAnalysisBondWorkflowContext & context)
{
    RunBondMapValueSampling(context);
    RunBondGroupClassification(context);
    RunLocalBondFitting(context, context.options.alpha_r);
    RunBondPotentialFitting(context);
}

} // namespace rhbm_gem::detail
