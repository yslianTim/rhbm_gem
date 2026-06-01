#include <rhbm_gem/core/GaussianEstimator.hpp>

#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem::core {
namespace {
RHBMExecutionOptions MakeExecutionOptions(const FitOptions & options)
{
    RHBMExecutionOptions execution_options;
    execution_options.quiet_mode = true;
    execution_options.thread_size = options.thread_size;
    return execution_options;
}

double EstimateInitialIntercept(const LocalPotentialSampleList & sample_entries)
{
    float maximum_distance{ 0.0f };
    for (const auto & sample : sample_entries)
    {
        if (sample.point.distance > maximum_distance) maximum_distance = sample.point.distance;
    }

    std::vector<double> response_list;
    response_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        if (sample.point.distance != maximum_distance) continue;
        const auto response{ static_cast<double>(sample.response) };
        numeric_validation::RequireFinite(response, "maximum distance shell response");
        response_list.emplace_back(response);
    }
    if (response_list.empty()) return 0.0;

    const auto lowest_response_list{
        array_helper::ComputeSmallestProportionValues(response_list, 0.10)
    };
    return array_helper::ComputeMedian(lowest_response_list);
}

double EstimateResidualIntercept(const LocalPotentialSampleList & residual_sample_entries)
{
    std::vector<double> response_list;
    response_list.reserve(residual_sample_entries.size());
    for (const auto & sample : residual_sample_entries)
    {
        response_list.emplace_back(static_cast<double>(sample.response));
    }
    if (response_list.empty()) return 0.0;
    return array_helper::ComputeMean(response_list.data(), response_list.size());
}

LocalPotentialSampleList BuildShiftedSampleEntries(
    const LocalPotentialSampleList & sample_entries,
    double intercept)
{
    LocalPotentialSampleList shifted_sample_entries;
    shifted_sample_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        shifted_sample_entries.emplace_back(
            LocalPotentialSample{
                static_cast<float>(static_cast<double>(sample.response) - intercept),
                sample.point
            }
        );
    }
    return shifted_sample_entries;
}

LocalPotentialSampleList BuildResidualSampleEntries(
    const LocalPotentialSampleList & sample_entries,
    const RHBMParameterVector & beta,
    LocalGaussianFitModel fit_model,
    double range_min,
    double range_max)
{
    const auto signal_model{ linearization_service::DecodeParameterVector(beta, fit_model) };
    LocalPotentialSampleList residual_sample_entries;
    residual_sample_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ sample.point.distance };
        if (distance < static_cast<float>(range_min)) continue;
        if (distance > static_cast<float>(range_max)) continue;

        const auto residual{
            static_cast<double>(sample.response) -
                signal_model.ResponseAtDistance(static_cast<double>(distance))
        };
        residual_sample_entries.emplace_back(
            LocalPotentialSample{ static_cast<float>(residual), sample.point }
        );
    }
    return residual_sample_entries;
}

rhbm_trainer::RHBMTrainingOptions MakeTrainingOptions(const FitOptions & options)
{
    rhbm_trainer::RHBMTrainingOptions training_options;
    training_options.execution_options = MakeExecutionOptions(options);
    return training_options;
}

std::vector<RHBMMemberDataset> BuildMemberDatasetList(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options)
{
    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (const auto & sample_entries : sample_entries_list)
    {
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(
                sample_entries,
                options.distance_min,
                options.distance_max,
                options.local_fit_model)
        );
    }
    return dataset_list;
}

std::vector<RHBMMemberDataset> BuildMemberDatasetList(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    const FitOptions & options)
{
    if (sample_entries_list.size() != member_result_list.size())
    {
        throw std::invalid_argument("sample_entries_list and member_result_list sizes are inconsistent.");
    }

    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (std::size_t i = 0; i < sample_entries_list.size(); i++)
    {
        const auto intercept{ member_result_list.at(i).mdpde.GetModel().GetIntercept() };
        const auto shifted_sample_entries{
            BuildShiftedSampleEntries(sample_entries_list.at(i), intercept)
        };
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(
                shifted_sample_entries,
                options.distance_min,
                options.distance_max,
                options.local_fit_model)
        );
    }
    return dataset_list;
}

LocalGaussianResult DecodeLocalGaussianResult(
    double alpha_r,
    const RHBMBetaEstimateResult & fit_result,
    LocalGaussianFitModel fit_model,
    double intercept = 0.0)
{
    const auto ols_model{
        linearization_service::DecodeParameterVector(fit_result.beta_ols, fit_model)
            .WithIntercept(intercept)
    };
    const auto mdpde_model{
        linearization_service::DecodeParameterVector(fit_result.beta_mdpde, fit_model)
            .WithIntercept(intercept)
    };
    return LocalGaussianResult{
        alpha_r,
        GaussianModel3DWithUncertainty{ ols_model, GaussianModel3DUncertainty{} },
        GaussianModel3DWithUncertainty{ mdpde_model, GaussianModel3DUncertainty{} },
        false,
        0.0,
        fit_model,
        fit_result
    };
}

GroupGaussianResult DecodeGroupGaussianResult(
    double alpha_g,
    const RHBMGroupEstimationResult & result,
    LocalGaussianFitModel fit_model)
{
    return GroupGaussianResult{
        alpha_g,
        linearization_service::DecodeParameterVector(result.mu_mean, fit_model),
        linearization_service::DecodeParameterVector(result.mu_mdpde, fit_model),
        linearization_service::DecodeParameterVector(
            result.mu_prior, result.capital_lambda, fit_model)
    };
}

std::vector<LocalGaussianResult> DecodeMemberGaussianResults(
    const RHBMGroupEstimationResult & result,
    LocalGaussianFitModel fit_model)
{
    const auto member_count{ static_cast<std::size_t>(result.beta_posterior_matrix.cols()) };
    if (result.capital_sigma_posterior_list.size() != member_count)
    {
        throw std::invalid_argument("Group Gaussian member result count is inconsistent.");
    }
    eigen_validation::RequireVectorSize(
        result.outlier_flag_array, result.beta_posterior_matrix.cols(),
        "outlier_flag_array", "Group Gaussian member result count is inconsistent.");
    eigen_validation::RequireVectorSize(
        result.statistical_distance_array, result.beta_posterior_matrix.cols(),
        "statistical_distance_array", "Group Gaussian member result count is inconsistent.");

    std::vector<LocalGaussianResult> member_results;
    member_results.reserve(member_count);
    for (Eigen::Index i = 0; i < result.beta_posterior_matrix.cols(); i++)
    {
        const auto gaussian{
            linearization_service::DecodeParameterVector(
                result.beta_posterior_matrix.col(i),
                result.capital_sigma_posterior_list.at(static_cast<std::size_t>(i)),
                fit_model)
        };
        member_results.emplace_back(LocalGaussianResult{
            0.0,
            gaussian,
            gaussian,
            static_cast<bool>(result.outlier_flag_array(i)),
            result.statistical_distance_array(i),
            fit_model
        });
    }
    return member_results;
}

std::vector<RHBMBetaEstimateResult> BuildMemberFitResultList(
    const std::vector<LocalGaussianResult> & member_result_list,
    LocalGaussianFitModel fit_model)
{
    std::vector<RHBMBetaEstimateResult> fit_result_list;
    fit_result_list.reserve(member_result_list.size());
    for (const auto & member_result : member_result_list)
    {
        if (!member_result.fit_result.has_value())
        {
            throw std::invalid_argument(
                "member_result_list contains a result without transient fit state.");
        }
        if (member_result.fit_model != fit_model)
        {
            throw std::invalid_argument(
                "group Gaussian estimation requires member results to match the requested local fit model.");
        }
        fit_result_list.emplace_back(*member_result.fit_result);
    }
    return fit_result_list;
}

} // namespace

double TrainAlphaR(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");

    const auto dataset_list{ BuildMemberDatasetList(sample_entries_list, options) };
    const auto training_options{ MakeTrainingOptions(options) };
    return rhbm_trainer::CrossValidationAlphaR(dataset_list, training_options).best_alpha;
}

double TrainAlphaG(
    const std::vector<std::vector<LocalGaussianResult>> & member_result_list,
    const FitOptions & options)
{
    std::vector<std::vector<RHBMParameterVector>> beta_group_list;
    beta_group_list.reserve(member_result_list.size());
    for (const auto & member_results : member_result_list)
    {
        std::vector<RHBMParameterVector> beta_list;
        beta_list.reserve(member_results.size());
        for (const auto & member_result : member_results)
        {
            if (member_result.fit_model != options.local_fit_model)
            {
                throw std::invalid_argument(
                    "Alpha_G training requires member results to match the requested local fit model.");
            }
            beta_list.emplace_back(
                linearization_service::EncodeGaussianToParameterVector(
                    member_result.mdpde.GetModel(),
                    member_result.fit_model));
        }
        beta_group_list.emplace_back(std::move(beta_list));
    }

    const auto training_options{ MakeTrainingOptions(options) };
    if (beta_group_list.empty())
    {
        return training_options.alpha_min;
    }

    return rhbm_trainer::CrossValidationAlphaG(beta_group_list, training_options).best_alpha;
}

LocalGaussianResult EstimateLocalGaussian(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");

    auto execution_options{ MakeExecutionOptions(options) };
    auto dataset{
        rhbm_helper::BuildMemberDataset(
            sample_entries,
            options.distance_min,
            options.distance_max,
            options.local_fit_model)
    };
    const auto result{ rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, execution_options) };
    return DecodeLocalGaussianResult(alpha_r, result, options.local_fit_model);
}

LocalGaussianResult EstimateLocalGaussianWithIntercept(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");

    auto execution_options{ MakeExecutionOptions(options) };
    double intercept{ EstimateInitialIntercept(sample_entries) };
    RHBMBetaEstimateResult result;
    bool has_previous_intercept{ false };

    for (int t = 0; t < execution_options.max_iterations; t++)
    {
        auto shifted_sample_entries{ BuildShiftedSampleEntries(sample_entries, intercept) };
        auto dataset{
            rhbm_helper::BuildMemberDataset(
                shifted_sample_entries,
                options.distance_min,
                options.distance_max,
                options.local_fit_model)
        };
        result = rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, execution_options);
        const auto residual_sample_entries{
            BuildResidualSampleEntries(
                sample_entries,
                result.beta_mdpde,
                options.local_fit_model,
                options.distance_min,
                options.distance_max)
        };
        const auto next_intercept{ EstimateResidualIntercept(residual_sample_entries) };
        if (has_previous_intercept && std::fabs(next_intercept - intercept) < execution_options.tolerance)
        {
            break;
        }
        has_previous_intercept = true;
        if (t + 1 < execution_options.max_iterations)
        {
            intercept = next_intercept;
        }
    }
    return DecodeLocalGaussianResult(alpha_r, result, options.local_fit_model, intercept);
}

GroupGaussianResult EstimateGroupGaussian(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    double alpha_g,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_g, "alpha_g");

    if (sample_entries_list.size() != member_result_list.size())
    {
        throw std::invalid_argument("sample_entries_list and member_result_list sizes are inconsistent.");
    }

    auto execution_options{ MakeExecutionOptions(options) };
    const auto dataset_list{ BuildMemberDatasetList(sample_entries_list, member_result_list, options) };
    const auto fit_result_list{ BuildMemberFitResultList(
        member_result_list,
        options.local_fit_model) };
    const auto group_input{ rhbm_helper::BuildGroupInput(dataset_list, fit_result_list) };
    const auto raw_result{ rhbm_helper::EstimateGroup(alpha_g, group_input, execution_options) };
    auto result{ DecodeGroupGaussianResult(alpha_g, raw_result, options.local_fit_model) };
    result.member_results = DecodeMemberGaussianResults(raw_result, options.local_fit_model);
    return result;
}

} // namespace rhbm_gem::core
