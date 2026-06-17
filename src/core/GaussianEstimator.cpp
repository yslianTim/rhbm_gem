#include <cstddef>
#include <rhbm_gem/core/GaussianEstimator.hpp>

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Constants.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/SampleFilter.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
#include <rhbm_gem/utils/math/EigenValidation.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::core {
namespace {
constexpr std::size_t kMinimumAlphaRTrainingSampleCount{ 10 };
constexpr std::size_t kMinimumAlphaGTrainingMemberCount{ 10 };
constexpr double kResidualInterceptRangeMin{ 1.0 };
constexpr double kResidualInterceptRangeMax{ 2.0 };
constexpr double kEstimatedInterceptMin{ -1.0 };
constexpr double kEstimatedInterceptMax{ 1.0 };
constexpr double kInterceptDampingFactor{ 0.5 };
constexpr double kNeighborContributionCutoffStart{ 2.0 };
constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr std::size_t kLocalFittingMaximumIterations{ 100 };
constexpr double kLocalFittingParameterChangeTolerance{ 1.0e-6 };
constexpr double kLocalFittingChangePercentile{ 0.95 };
constexpr int kHuberSlopeMaximumIterations{ 50 };
constexpr double kHuberSlopeTolerance{ 1.0e-8 };
constexpr double kHuberScaleMultiplier{ 1.4826 };
constexpr double kHuberScaleMin{ 1.0e-12 };
constexpr double kHuberCutoffMultiplier{ 1.345 };

struct ResidualInterceptSample
{
    double basis{ 0.0 };
    double residual{ 0.0 };
};

struct LocalFittingParameterChangeStats
{
    double amplitude_change_percentile{ 0.0 };
    double width_change_percentile{ 0.0 };
    double intercept_change_percentile{ 0.0 };
};

struct LocalFittingIterationResult
{
    std::vector<LocalPotentialSampleList> sample_entries_list;
    std::vector<LocalGaussianResult> result_list;
    std::vector<Eigen::VectorXd> estimation_list;
};

double ClampEstimatedIntercept(double intercept)
{
    if (intercept < kEstimatedInterceptMin) return kEstimatedInterceptMin;
    if (intercept > kEstimatedInterceptMax) return kEstimatedInterceptMax;
    return intercept;
}

double CalculateNeighborContributionCutoffWeight(double distance)
{
    if (distance <= kNeighborContributionCutoffStart) return 1.0;
    if (distance >= kNeighborContributionDistanceMax) return 0.0;

    const auto transition_fraction{
        (distance - kNeighborContributionCutoffStart) /
        (kNeighborContributionDistanceMax - kNeighborContributionCutoffStart)
    };
    return 0.5 * (1.0 + std::cos(Constants::pi * transition_fraction));
}

std::vector<AtomLocalPotentialEditor> BuildSelectedAtomLocalEditors(ModelObject & model_object)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    std::vector<AtomLocalPotentialEditor> local_editor_list;
    local_editor_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        local_editor_list.emplace_back(analysis.EnsureAtomLocalPotential(*atom));
    }
    return local_editor_list;
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
        if (sample.point.distance < fit_range_min || sample.point.distance > fit_range_max) continue;
        count++;
        if (count >= minimum_sample_count) return true;
    }
    return false;
}

RHBMExecutionOptions MakeExecutionOptions(const FitOptions & options)
{
    RHBMExecutionOptions execution_options;
    execution_options.quiet_mode = true;
    execution_options.thread_size = options.thread_size;
    return execution_options;
}

bool CanBuildFiniteZeroInterceptSamples(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model)
{
    for (const auto & sample : sample_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto model_offset{ model.ResponseAtDistance(distance) - model.SignalAtDistance(distance) };
        const auto response{ static_cast<double>(sample.response) - model_offset };
        if (!std::isfinite(response)) return false;
        if (std::abs(response) > static_cast<double>(std::numeric_limits<float>::max()))
        {
            return false;
        }
    }
    return true;
}

double CalculateRelativeScalarChange(double previous_value, double current_value)
{
    if (!std::isfinite(previous_value) || !std::isfinite(current_value))
    {
        return std::numeric_limits<double>::infinity();
    }
    const auto scale{
        std::max({ 1.0, std::abs(previous_value), std::abs(current_value) })
    };
    return std::abs(current_value - previous_value) / scale;
}

double CalculateMaxRelativeParameterChange(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & current_model)
{
    double max_change{ 0.0 };
    for (int i = 0; i < GaussianModel3D::ParameterSize(); i++)
    {
        const auto parameter_change{
            CalculateRelativeScalarChange(
                previous_model.GetModelParameter(i),
                current_model.GetModelParameter(i))
        };
        if (!std::isfinite(parameter_change))
        {
            return std::numeric_limits<double>::infinity();
        }
        max_change = std::max(max_change, parameter_change);
    }
    return max_change;
}

bool EstimateOrdinarySlopeThroughOrigin(
    const std::vector<ResidualInterceptSample> & sample_list,
    double & slope)
{
    double numerator{ 0.0 };
    double denominator{ 0.0 };
    for (const auto & sample : sample_list)
    {
        numerator += sample.basis * sample.residual;
        denominator += sample.basis * sample.basis;
    }
    if (!std::isfinite(numerator) ||
        !std::isfinite(denominator) ||
        denominator <= std::numeric_limits<double>::epsilon())
    {
        return false;
    }
    slope = numerator / denominator;
    return std::isfinite(slope);
}

double ComputeHuberResidualScale(
    const std::vector<ResidualInterceptSample> & sample_list,
    double slope)
{
    std::vector<double> residual_list;
    residual_list.reserve(sample_list.size());
    for (const auto & sample : sample_list)
    {
        residual_list.emplace_back(sample.residual - slope * sample.basis);
    }

    const auto median_residual{ array_helper::ComputeMedian(residual_list) };
    std::vector<double> deviation_list;
    deviation_list.reserve(residual_list.size());
    for (const auto residual : residual_list)
    {
        deviation_list.emplace_back(std::abs(residual - median_residual));
    }

    return std::max(
        kHuberScaleMultiplier * array_helper::ComputeMedian(deviation_list),
        kHuberScaleMin);
}

bool EstimateHuberSlopeThroughOrigin(
    const std::vector<ResidualInterceptSample> & sample_list,
    double & slope)
{
    if (sample_list.empty())
    {
        return false;
    }
    if (!EstimateOrdinarySlopeThroughOrigin(sample_list, slope))
    {
        return false;
    }

    for (int t = 0; t < kHuberSlopeMaximumIterations; t++)
    {
        const auto scale{ ComputeHuberResidualScale(sample_list, slope) };
        const auto cutoff{ kHuberCutoffMultiplier * scale };
        double numerator{ 0.0 };
        double denominator{ 0.0 };
        for (const auto & sample : sample_list)
        {
            const auto error{ sample.residual - slope * sample.basis };
            const auto abs_error{ std::abs(error) };
            const auto weight{ abs_error <= cutoff ? 1.0 : cutoff / abs_error };
            numerator += weight * sample.basis * sample.residual;
            denominator += weight * sample.basis * sample.basis;
        }
        if (!std::isfinite(numerator) ||
            !std::isfinite(denominator) ||
            denominator <= std::numeric_limits<double>::epsilon())
        {
            return false;
        }

        const auto updated_slope{ numerator / denominator };
        if (!std::isfinite(updated_slope))
        {
            return false;
        }
        if (std::abs(updated_slope - slope) < kHuberSlopeTolerance)
        {
            slope = updated_slope;
            return true;
        }
        slope = updated_slope;
    }
    return true;
}

double EstimateResidualInterceptParameter(
    const LocalPotentialSampleList & sample_entries,
    const RHBMBetaEstimateResult & fit_result,
    double current_intercept)
{
    const auto signal_model{ linearization_service::DecodeParameterVector(fit_result.beta_mdpde) };
    const auto width{ signal_model.GetWidth() };
    if (!std::isfinite(width) || width <= 0.0) return current_intercept;

    std::vector<ResidualInterceptSample> residual_sample_list;
    residual_sample_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        if (distance < kResidualInterceptRangeMin) continue;
        if (distance > kResidualInterceptRangeMax) continue;

        const auto basis{ signal_model.InterceptBasisAtDistance(distance) };
        if (!std::isfinite(basis) || std::abs(basis) <= std::numeric_limits<double>::epsilon())
        {
            continue;
        }
        const auto residual{
            static_cast<double>(sample.response) - signal_model.SignalAtDistance(distance)
        };
        if (!std::isfinite(residual))
        {
            continue;
        }
        residual_sample_list.emplace_back(ResidualInterceptSample{ basis, residual });
    }
    double candidate_intercept{ current_intercept };
    if (!EstimateHuberSlopeThroughOrigin(residual_sample_list, candidate_intercept))
    {
        return current_intercept;
    }
    const auto candidate_model{ signal_model.WithIntercept(candidate_intercept) };
    if (!CanBuildFiniteZeroInterceptSamples(sample_entries, candidate_model))
    {
        return current_intercept;
    }
    return candidate_intercept;
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
                options.distance_max)
        );
    }
    return dataset_list;
}

std::size_t GetMinimumDatasetResponseCount(const std::vector<RHBMMemberDataset> & dataset_list)
{
    std::size_t minimum_response_count{ std::numeric_limits<std::size_t>::max() };
    for (const auto & dataset : dataset_list)
    {
        const auto response_count{ static_cast<std::size_t>(dataset.y.size()) };
        if (response_count < minimum_response_count)
        {
            minimum_response_count = response_count;
        }
    }
    return minimum_response_count;
}

LocalPotentialSampleList BuildSamplesForZeroInterceptGaussianFit(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model)
{
    LocalPotentialSampleList updated_sample_entries;
    updated_sample_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto model_offset{ model.ResponseAtDistance(distance) - model.SignalAtDistance(distance) };
        const auto response{ static_cast<float>(static_cast<double>(sample.response) - model_offset)};
        updated_sample_entries.emplace_back(LocalPotentialSample{ response, sample.point });
    }
    return updated_sample_entries;
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
    auto range_min{ options.distance_min };
    auto range_max{ options.distance_max };
    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (std::size_t i = 0; i < sample_entries_list.size(); i++)
    {
        const auto sampling_entries{
            BuildSamplesForZeroInterceptGaussianFit(
                sample_entries_list.at(i),
                member_result_list.at(i).mdpde.GetModel())
        };
        dataset_list.emplace_back(rhbm_helper::BuildMemberDataset(sampling_entries, range_min, range_max));
    }
    return dataset_list;
}

LocalGaussianResult DecodeLocalGaussianResult(
    double alpha_r,
    const RHBMBetaEstimateResult & fit_result,
    double intercept = 0.0)
{
    const auto ols_model{
        linearization_service::DecodeParameterVector(fit_result.beta_ols)
            .WithIntercept(intercept)
    };
    const auto mdpde_model{
        linearization_service::DecodeParameterVector(fit_result.beta_mdpde)
            .WithIntercept(intercept)
    };
    return LocalGaussianResult{
        alpha_r,
        GaussianModel3DWithUncertainty{ ols_model, GaussianModel3DUncertainty{} },
        GaussianModel3DWithUncertainty{ mdpde_model, GaussianModel3DUncertainty{} },
        std::nullopt,
        false,
        0.0,
        fit_result
    };
}

GaussianModel3DWithUncertainty WithModelIntercept(
    const GaussianModel3DWithUncertainty & gaussian,
    double intercept)
{
    return GaussianModel3DWithUncertainty{
        gaussian.GetModel().WithIntercept(intercept),
        gaussian.GetStandardDeviationModel()
    };
}

GroupGaussianResult DecodeGroupGaussianResult(
    double alpha_g,
    const RHBMGroupEstimationResult & result,
    double intercept)
{
    return GroupGaussianResult{
        alpha_g,
        linearization_service::DecodeParameterVector(result.mu_mean).WithIntercept(intercept),
        linearization_service::DecodeParameterVector(result.mu_mdpde).WithIntercept(intercept),
        WithModelIntercept(
            linearization_service::DecodeParameterVector(result.mu_prior, result.capital_lambda),
            intercept)
    };
}

std::vector<LocalGaussianResult> DecodeMemberGaussianResults(
    const RHBMGroupEstimationResult & result,
    const std::vector<LocalGaussianResult> & member_result_list)
{
    const auto member_count{ static_cast<std::size_t>(result.beta_posterior_matrix.cols()) };
    if (member_result_list.size() != member_count)
    {
        throw std::invalid_argument("Group Gaussian member result count is inconsistent.");
    }
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
        const auto member_index{ static_cast<std::size_t>(i) };
        const auto intercept{
            member_result_list.at(member_index).mdpde.GetModel().GetIntercept()
        };
        const auto gaussian{
            WithModelIntercept(
                linearization_service::DecodeParameterVector(
                    result.beta_posterior_matrix.col(i),
                    result.capital_sigma_posterior_list.at(member_index)),
                intercept)
        };
        member_results.emplace_back(LocalGaussianResult{
            0.0,
            gaussian,
            gaussian,
            gaussian,
            static_cast<bool>(result.outlier_flag_array(i)),
            result.statistical_distance_array(i)
        });
    }
    return member_results;
}

std::vector<RHBMBetaEstimateResult> BuildMemberFitResultList(
    const std::vector<RHBMMemberDataset> & dataset_list,
    const std::vector<LocalGaussianResult> & member_result_list,
    const FitOptions & options)
{
    if (dataset_list.size() != member_result_list.size())
    {
        throw std::invalid_argument("dataset_list and member_result_list sizes are inconsistent.");
    }
    const auto execution_options{ MakeExecutionOptions(options) };
    std::vector<RHBMBetaEstimateResult> fit_result_list;
    fit_result_list.reserve(member_result_list.size());
    for (std::size_t i = 0; i < member_result_list.size(); i++)
    {
        fit_result_list.emplace_back(
            rhbm_helper::EstimateBetaMDPDE(
                member_result_list.at(i).alpha_r,
                dataset_list.at(i),
                execution_options));
    }
    return fit_result_list;
}

using FittedGaussianSnapshot = std::unordered_map<const AtomObject *, GaussianModel3D>;

FittedGaussianSnapshot BuildFittedGaussianSnapshot(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<Eigen::VectorXd> & estimation_list)
{
    if (atom_list.size() != estimation_list.size())
    {
        throw std::invalid_argument("atom_list and estimation_list sizes are inconsistent.");
    }

    FittedGaussianSnapshot snapshot;
    snapshot.reserve(atom_list.size());
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        snapshot.emplace(atom_list.at(i), GaussianModel3D::FromVector(estimation_list.at(i)));
    }
    return snapshot;
}

LocalPotentialSampleList UpdateSampleListWithFittedGaussian(
    const AtomObject & atom,
    const FittedGaussianSnapshot & snapshot)
{
    const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
    const auto sample_entries{ local_view.GetSamplingEntries(false) };
    const auto & neighbor_atom_list{ atom.FindNeighborAtoms() };
    LocalPotentialSampleList updated_list;
    updated_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        auto sample_position{ sample.point.position };
        auto response_value{ sample.response };
        for (const auto * neighbor_atom : neighbor_atom_list)
        {
            const auto gaussian_iter{ snapshot.find(neighbor_atom) };
            if (gaussian_iter == snapshot.end()) continue;

            auto neighbor_position{ neighbor_atom->GetPosition() };
            auto distance{
                static_cast<double>(
                    array_helper::ComputeNorm<float>(sample_position, neighbor_position))
            };
            const auto cutoff_weight{ CalculateNeighborContributionCutoffWeight(distance) };
            if (cutoff_weight == 0.0) continue;
            response_value -= static_cast<float>(
                cutoff_weight * gaussian_iter->second.ResponseAtDistance(distance));
        }
        updated_list.emplace_back(LocalPotentialSample{response_value, sample.point });
    }
    return updated_list;
}

bool HasFiniteSampleResponses(const LocalPotentialSampleList & sample_entries)
{
    for (const auto & sample : sample_entries)
    {
        if (!std::isfinite(static_cast<double>(sample.response))) return false;
    }
    return true;
}

bool CanRefreshLocalFittingSampleEntries(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<Eigen::VectorXd> & estimation_list)
{
    try
    {
        const auto snapshot{ BuildFittedGaussianSnapshot(atom_list, estimation_list) };
        for (const auto * atom : atom_list)
        {
            const auto sample_entries{ UpdateSampleListWithFittedGaussian(*atom, snapshot) };
            if (!HasFiniteSampleResponses(sample_entries))
            {
                return false;
            }
            const auto model_iter{ snapshot.find(atom) };
            if (model_iter == snapshot.end()) return false;
            if (!CanBuildFiniteZeroInterceptSamples(sample_entries, model_iter->second))
            {
                return false;
            }
        }
    }
    catch (const std::exception &)
    {
        return false;
    }
    return true;
}

LocalFittingParameterChangeStats CalculateLocalFittingParameterChangeStats(
    const std::vector<Eigen::VectorXd> & current_estimation_list,
    const std::vector<Eigen::VectorXd> & previous_estimation_list)
{
    LocalFittingParameterChangeStats stats;
    std::vector<double> amplitude_change_list(current_estimation_list.size());
    std::vector<double> width_change_list(current_estimation_list.size());
    std::vector<double> intercept_change_list(current_estimation_list.size());
    for (size_t i = 0; i < current_estimation_list.size(); i++)
    {
        const auto parameter_delta{ current_estimation_list[i] - previous_estimation_list[i] };
        const auto amplitude_change{
            std::abs(parameter_delta(GaussianModel3D::AmplitudeIndex()))
        };
        const auto width_change{
            std::abs(parameter_delta(GaussianModel3D::WidthIndex()))
        };
        const auto intercept_change{
            std::abs(parameter_delta(GaussianModel3D::InterceptIndex()))
        };
        amplitude_change_list[i] = amplitude_change;
        width_change_list[i] = width_change;
        intercept_change_list[i] = intercept_change;
    }

    stats.amplitude_change_percentile = array_helper::ComputePercentile(
        amplitude_change_list,
        kLocalFittingChangePercentile);
    stats.width_change_percentile = array_helper::ComputePercentile(
        width_change_list,
        kLocalFittingChangePercentile);
    stats.intercept_change_percentile = array_helper::ComputePercentile(
        intercept_change_list,
        kLocalFittingChangePercentile);
    return stats;
}

bool IsLocalFittingParameterChangeConverged(const LocalFittingParameterChangeStats & stats)
{
    return
        (std::pow(stats.amplitude_change_percentile, 2) < kLocalFittingParameterChangeTolerance) &&
        (std::pow(stats.width_change_percentile, 2) < kLocalFittingParameterChangeTolerance) &&
        (std::pow(stats.intercept_change_percentile, 2) < kLocalFittingParameterChangeTolerance);
}

double GetLocalFittingParameterChange(const LocalFittingParameterChangeStats & stats)
{
    const auto shape_change{
        stats.amplitude_change_percentile > stats.width_change_percentile ?
            stats.amplitude_change_percentile :
            stats.width_change_percentile
    };
    return shape_change > stats.intercept_change_percentile ?
        shape_change :
        stats.intercept_change_percentile;
}

bool IsBetterLocalFittingCandidate(
    const LocalFittingParameterChangeStats & stats,
    const LocalFittingParameterChangeStats & best_stats)
{
    return GetLocalFittingParameterChange(stats) < GetLocalFittingParameterChange(best_stats);
}

void ApplyLocalFittingUnderRelaxation(
    LocalFittingIterationResult & iteration_result,
    const std::vector<Eigen::VectorXd> & previous_estimation_list,
    double beta)
{
    if (iteration_result.estimation_list.size() != previous_estimation_list.size() ||
        iteration_result.result_list.size() != previous_estimation_list.size())
    {
        throw std::invalid_argument("Local fitting relaxation input sizes are inconsistent.");
    }
    for (std::size_t i = 0; i < iteration_result.estimation_list.size(); i++)
    {
        const auto relaxed_estimation{
            beta * iteration_result.estimation_list.at(i) +
            (1.0 - beta) * previous_estimation_list.at(i)
        };
        const auto relaxed_model{ GaussianModel3D::FromVector(relaxed_estimation) };
        auto & result{ iteration_result.result_list.at(i) };
        result.mdpde = GaussianModel3DWithUncertainty{
            relaxed_model,
            result.mdpde.GetStandardDeviationModel()
        };
        iteration_result.estimation_list.at(i) = relaxed_estimation;
    }
}

LocalGaussianResult EstimateLocalGaussianWithOffsetModel(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options,
    const GaussianModel3D & offset_model)
{
    auto range_min{ options.distance_min };
    auto range_max{ options.distance_max };
    numeric_validation::RequireFiniteNonNegativeRange(range_min, range_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");
    numeric_validation::RequireFinite(offset_model.GetIntercept(), "intercept");

    auto execution_options{ MakeExecutionOptions(options) };
    const auto updated_sample_entries{
        BuildSamplesForZeroInterceptGaussianFit(sample_entries, offset_model)
    };
    auto dataset{
        rhbm_helper::BuildMemberDataset(updated_sample_entries, range_min, range_max)
    };
    const auto result{ rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, execution_options) };
    return DecodeLocalGaussianResult(alpha_r, result, offset_model.GetIntercept());
}

LocalFittingIterationResult RunLocalFittingIteration(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<Eigen::VectorXd> & input_estimation_list,
    const std::vector<LocalGaussianResult> & input_result_list,
    const FitOptions & options)
{
    const auto selected_atom_size{ atom_list.size() };
    if (input_result_list.size() != selected_atom_size)
    {
        throw std::invalid_argument("Local fitting iteration input sizes are inconsistent.");
    }
    const auto snapshot{ BuildFittedGaussianSnapshot(atom_list, input_estimation_list) };
    LocalFittingIterationResult iteration_result{
        std::vector<LocalPotentialSampleList>(selected_atom_size),
        std::vector<LocalGaussianResult>(selected_atom_size),
        std::vector<Eigen::VectorXd>(selected_atom_size)
    };

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        const auto & atom{ *atom_list[i] };
        const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
        auto sample_entries{
            UpdateSampleListWithFittedGaussian(atom, snapshot)
        };
        auto result{ input_result_list.at(i) };
        const auto intercept{ snapshot.at(&atom).GetIntercept() };
        try
        {
            auto candidate_result{
                EstimateLocalGaussianWithIntercept(
                    sample_entries, local_view.GetAlphaR(), options, intercept)
            };
            if (CanBuildFiniteZeroInterceptSamples(sample_entries, candidate_result.mdpde.GetModel()))
            {
                result = std::move(candidate_result);
            }
        }
        catch (const std::exception &)
        {
            result = input_result_list.at(i);
        }
        const auto fitted_model{ result.mdpde.GetModel() };
        iteration_result.sample_entries_list[i] = std::move(sample_entries);
        iteration_result.estimation_list[i] = fitted_model.ToVector();
        iteration_result.result_list[i] = std::move(result);
    }
    if (!CanRefreshLocalFittingSampleEntries(atom_list, iteration_result.estimation_list))
    {
        const auto input_snapshot{ BuildFittedGaussianSnapshot(atom_list, input_estimation_list) };
        for (size_t i = 0; i < selected_atom_size; i++)
        {
            iteration_result.sample_entries_list[i] =
                UpdateSampleListWithFittedGaussian(*atom_list.at(i), input_snapshot);
            iteration_result.result_list[i] = input_result_list.at(i);
            iteration_result.estimation_list[i] = input_estimation_list.at(i);
        }
    }
    return iteration_result;
}

void RefreshLocalFittingIterationSampleEntries(
    LocalFittingIterationResult & iteration_result,
    const std::vector<AtomObject *> & atom_list)
{
    const auto snapshot{ BuildFittedGaussianSnapshot(atom_list, iteration_result.estimation_list) };
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        iteration_result.sample_entries_list.at(i) =
            UpdateSampleListWithFittedGaussian(*atom_list.at(i), snapshot);
    }
}

void ApplyLocalFittingIterationResult(
    LocalFittingIterationResult & iteration_result,
    const std::vector<AtomObject *> & atom_list,
    std::vector<AtomLocalPotentialEditor> & local_editor_list)
{
    if (local_editor_list.size() != iteration_result.result_list.size())
    {
        throw std::invalid_argument(
            "local_editor_list and iteration_result sizes are inconsistent.");
    }
    if (atom_list.size() != iteration_result.result_list.size())
    {
        throw std::invalid_argument("atom_list and iteration_result sizes are inconsistent.");
    }

    RefreshLocalFittingIterationSampleEntries(iteration_result, atom_list);
    for (std::size_t i = 0; i < local_editor_list.size(); i++)
    {
        local_editor_list.at(i).SetGaussianResult(iteration_result.result_list.at(i));
        local_editor_list.at(i).SetSamplingEntries(
            std::move(iteration_result.sample_entries_list.at(i)));
    }
}

} // namespace

double TrainAlphaR(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");

    const auto dataset_list{ BuildMemberDatasetList(sample_entries_list, options) };
    auto training_options{ MakeTrainingOptions(options) };
    if (!dataset_list.empty())
    {
        const auto minimum_response_count{ GetMinimumDatasetResponseCount(dataset_list) };
        if (minimum_response_count < 2)
        {
            return training_options.alpha_min;
        }
        if (training_options.subset_size > minimum_response_count)
        {
            training_options.subset_size = minimum_response_count;
        }
    }
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
            beta_list.emplace_back(
                linearization_service::EncodeGaussianToParameterVector(member_result.mdpde.GetModel()));
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
    const FitOptions & options,
    double intercept)
{
    numeric_validation::RequireFinite(intercept, "intercept");
    const auto zero_intercept_result{
        EstimateLocalGaussianWithOffsetModel(
            sample_entries,
            alpha_r,
            options,
            GaussianModel3D{ 0.0, 1.0, 0.0 })
    };
    if (intercept == 0.0)
    {
        return zero_intercept_result;
    }
    const auto offset_model{ zero_intercept_result.mdpde.GetModel().WithIntercept(intercept) };
    return EstimateLocalGaussianWithOffsetModel(sample_entries, alpha_r, options, offset_model);
}

LocalGaussianResult EstimateLocalGaussianWithIntercept(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options,
    double intercept_initial)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");
    numeric_validation::RequireFinite(intercept_initial, "intercept_initial");

    auto execution_options{ MakeExecutionOptions(options) };
    auto result{
        EstimateLocalGaussianWithOffsetModel(
            sample_entries,
            alpha_r,
            options,
            GaussianModel3D{ 0.0, 1.0, 0.0 })
    };
    auto current_model{
        result.mdpde.GetModel().WithIntercept(ClampEstimatedIntercept(intercept_initial))
    };
    auto previous_convergence_model{ current_model };
    double previous_robust_objective{ std::numeric_limits<double>::infinity() };
    bool has_previous_convergence_state{ false };
    double best_robust_objective{ std::numeric_limits<double>::infinity() };
    auto best_result{ result };
    bool has_best_result{ false };
    auto max_iterations{ execution_options.max_iterations };
    auto tolerance{ execution_options.tolerance };
    for (int t = 0; t < max_iterations; t++)
    {
        const auto intercept{ current_model.GetIntercept() };
        result = EstimateLocalGaussianWithOffsetModel(sample_entries, alpha_r, options, current_model);
        const auto fitted_model{ result.mdpde.GetModel() };
        const auto robust_objective{
            result.fit_result.has_value() ?
                result.fit_result->mdpde_objective :
                std::numeric_limits<double>::infinity()
        };
        if (std::isfinite(robust_objective) && robust_objective < best_robust_objective)
        {
            best_robust_objective = robust_objective;
            best_result = result;
            has_best_result = true;
        }
        double parameter_change{ std::numeric_limits<double>::infinity() };
        double robust_objective_change{ std::numeric_limits<double>::infinity() };
        if (has_previous_convergence_state)
        {
            parameter_change =
                CalculateMaxRelativeParameterChange(previous_convergence_model, fitted_model);
            robust_objective_change =
                CalculateRelativeScalarChange(previous_robust_objective, robust_objective);
            if (parameter_change < tolerance ||
                robust_objective_change < tolerance)
            {
                break;
            }
        }

        if (t + 1 == max_iterations)
        {
            result = has_best_result ?
                best_result :
                result;
            if (!options.quiet_mode)
            {
                Logger::Log(LogLevel::Debug,
                    "Maximum iterations reached in local Gaussian estimation with intercept; "
                    "returning best robust-objective candidate with objective = " +
                    std::to_string(best_robust_objective) +
                    " and latest parameter relative change = " +
                    std::to_string(parameter_change) +
                    ", robust objective relative change = " +
                    std::to_string(robust_objective_change) + ".");
            }
            break;
        }

        const auto raw_intercept{
            ClampEstimatedIntercept(
                EstimateResidualInterceptParameter(sample_entries, *result.fit_result, intercept))
        };
        const auto damped_intercept{
            ClampEstimatedIntercept(intercept + kInterceptDampingFactor * (raw_intercept - intercept))
        };
        current_model = result.mdpde.GetModel().WithIntercept(damped_intercept);
        previous_convergence_model = fitted_model;
        previous_robust_objective = robust_objective;
        has_previous_convergence_state = true;
    }
    return result;
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
    const auto fit_result_list{ BuildMemberFitResultList(dataset_list, member_result_list, options) };
    const auto group_input{ rhbm_helper::BuildGroupInput(dataset_list, fit_result_list) };
    const auto raw_result{ rhbm_helper::EstimateGroup(alpha_g, group_input, execution_options) };
    std::vector<double> member_intercept_list;
    member_intercept_list.reserve(member_result_list.size());
    for (const auto & member_result : member_result_list)
    {
        member_intercept_list.emplace_back(member_result.mdpde.GetModel().GetIntercept());
    }
    const auto group_intercept{ array_helper::ComputeMedian(member_intercept_list) };
    auto result{ DecodeGroupGaussianResult(alpha_g, raw_result, group_intercept) };
    result.member_results = DecodeMemberGaussianResults(raw_result, member_result_list);
    return result;
}

void RunLocalAlphaTraining(ModelObject & model_object, const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys() };

    size_t count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run local alpha training for " + std::to_string(group_key_list.size()) + " groups.");
    }
    for (const auto group_key : group_key_list)
    {
        const auto & group_atom_list{
            analysis_view.GetAtomObjectList(group_key)
        };
        std::vector<LocalPotentialSampleList> sample_entries_list;
        sample_entries_list.reserve(group_atom_list.size());
        for (auto * atom : group_atom_list)
        {
            analysis.EnsureAtomLocalPotential(*atom);
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            const auto & sample_entries{ local_view.GetSamplingEntries() };
            if (!HasEnoughSamplesInFitRange(
                    sample_entries,
                    options.distance_min,
                    options.distance_max,
                    kMinimumAlphaRTrainingSampleCount)) continue;
            sample_entries_list.emplace_back(sample_entries);
        }
        sample_entries_list.shrink_to_fit();
        if (!sample_entries_list.empty())
        {
            const auto alpha_r{ TrainAlphaR(sample_entries_list, options) };
            for (auto * atom : group_atom_list)
            {
                analysis.EnsureAtomLocalPotential(*atom).SetAlphaR(alpha_r);
            }
        }
        count++;
        if (!options.quiet_mode)
        {
            Logger::ProgressPercent(count, group_key_list.size());
        }
    }
}

void RunGroupAlphaTraining(ModelObject & model_object, const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys() };

    std::vector<std::vector<LocalGaussianResult>> member_result_list;
    member_result_list.reserve(group_key_list.size());
    for (const auto group_key : group_key_list)
    {
        const auto & group_atom_list{
            analysis_view.GetAtomObjectList(group_key)
        };
        if (group_atom_list.size() < kMinimumAlphaGTrainingMemberCount) continue;
        if (group_atom_list.front()->IsMainChainAtom() == false) continue;

        std::vector<LocalGaussianResult> group_member_results;
        group_member_results.reserve(group_atom_list.size());
        for (auto * atom : group_atom_list)
        {
            analysis.EnsureAtomLocalPotential(*atom);
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            group_member_results.emplace_back(local_view.GetGaussianResult());
        }
        member_result_list.emplace_back(std::move(group_member_results));
    }

    const auto alpha_g{ TrainAlphaG(member_result_list, options) };
    for (const auto group_key : analysis_view.CollectAtomGroupKeys())
    {
        analysis.SetAtomGroupAlphaG(group_key, alpha_g);
    }
}

void RunFirstStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    const auto selected_atom_size{ model_object.GetSelectedAtomCount() };
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildSelectedAtomLocalEditors(model_object) };
    std::atomic<size_t> atom_count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info,
            "Run first-stage local atom fitting for " +
            std::to_string(selected_atom_size) + " atoms.");
    }

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list[i]) };
        auto sample_entries{ local_view.GetSamplingEntries() };
        const auto result{
            EstimateLocalGaussianWithIntercept(
                sample_entries, local_view.GetAlphaR(), options, 0.0)
        };
        local_editor_list[i].SetGaussianResult(result);

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            atom_count++;
            if (!options.quiet_mode)
            {
                Logger::ProgressPercent(atom_count, selected_atom_size);
            }
        }
    }
}

void RunSecondStageLocalFitting(
    ModelObject & model_object,
    const std::vector<AtomObject *> & atom_list,
    const FitOptions & options)
{
    const auto under_relaxation_factor{
        numeric_validation::RequireFiniteExclusiveInclusiveRange(
            options.local_fitting_under_relaxation_factor,
            0.0,
            1.0,
            "local_fitting_under_relaxation_factor")
    };
    auto analysis{ model_object.EditAnalysis() };
    const auto atom_size{ atom_list.size() };
    std::vector<AtomLocalPotentialEditor> local_editor_list;
    local_editor_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        local_editor_list.emplace_back(analysis.EnsureAtomLocalPotential(*atom));
    }

    std::vector<Eigen::VectorXd> previous_estimation_list(atom_size);
    std::vector<LocalGaussianResult> previous_result_list(atom_size);
    for (size_t i = 0; i < atom_size; i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list[i]) };
        previous_result_list[i] = local_view.GetGaussianResult();
        previous_estimation_list[i] = previous_result_list[i].mdpde.GetModel().ToVector();
    }

    LocalFittingIterationResult best_iteration_result;
    LocalFittingParameterChangeStats best_change_stats;
    bool has_best_iteration_result{ false };
    for (size_t iter = 0; iter < kLocalFittingMaximumIterations; iter++)
    {
        auto iteration_result{
            RunLocalFittingIteration(
                atom_list,
                previous_estimation_list,
                previous_result_list,
                options)
        };
        ApplyLocalFittingUnderRelaxation(
            iteration_result,
            previous_estimation_list,
            under_relaxation_factor);
        const auto change_stats{
            CalculateLocalFittingParameterChangeStats(
                iteration_result.estimation_list,
                previous_estimation_list)
        };
        if (!has_best_iteration_result ||
            IsBetterLocalFittingCandidate(change_stats, best_change_stats))
        {
            best_iteration_result = iteration_result;
            best_change_stats = change_stats;
            has_best_iteration_result = true;
        }

        if (!options.quiet_mode)
        {
            std::ostringstream progress_message;
            progress_message << "Local fitting iteration " << iter + 1 << '/'
                << kLocalFittingMaximumIterations
                << std::fixed << std::setprecision(5)
                << ", percentile amplitude change = "
                << change_stats.amplitude_change_percentile
                << ", percentile width change = "
                << change_stats.width_change_percentile
                << ", percentile intercept change = "
                << change_stats.intercept_change_percentile;
            Logger::ProgressLine(progress_message.str());
        }

        const auto converged{ IsLocalFittingParameterChangeConverged(change_stats) };
        if (converged)
        {
            ApplyLocalFittingIterationResult(
                iteration_result,
                atom_list,
                local_editor_list);
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                Logger::Log(LogLevel::Info,
                    "Converged after " + std::to_string(iter + 1) +
                    " iterations with percentile amplitude change = " +
                    std::to_string(change_stats.amplitude_change_percentile) +
                    ", percentile width change = " +
                    std::to_string(change_stats.width_change_percentile) +
                    ", and percentile intercept change = " +
                    std::to_string(change_stats.intercept_change_percentile) + ".");
            }
            break;
        }

        if (iter + 1 == kLocalFittingMaximumIterations)
        {
            ApplyLocalFittingIterationResult(
                best_iteration_result,
                atom_list,
                local_editor_list);
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                Logger::Log(LogLevel::Warning,
                    "Reached maximum iteration size; refitting at best fixed-point candidate "
                    "with percentile amplitude change = " +
                    std::to_string(best_change_stats.amplitude_change_percentile) +
                    ", percentile width change = " +
                    std::to_string(best_change_stats.width_change_percentile) +
                    ", and percentile intercept change = " +
                    std::to_string(best_change_stats.intercept_change_percentile));
            }
        }
        previous_estimation_list = std::move(iteration_result.estimation_list);
        previous_result_list = std::move(iteration_result.result_list);
    }
}

void RunLocalPotentialFitting(ModelObject & model_object, const FitOptions & options)
{
    RunFirstStageLocalFitting(model_object, options);

    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run updated local atom fitting with iterations...");
    }

    const auto & atom_list{ model_object.GetSelectedAtoms() };
    RunSecondStageLocalFitting(model_object, atom_list, options);
}

void RunGroupPotentialFitting(ModelObject & model_object, const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    for (auto * atom : model_object.GetSelectedAtoms())
    {
        analysis.EnsureAtomLocalPotential(*atom);
    }
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run component atom group fitting.");
    }

    auto group_key_list{ analysis_view.CollectAtomGroupKeys() };
    auto group_key_size{ group_key_list.size() };
    std::atomic<size_t> key_count{ 0 };

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t k = 0; k < group_key_size; k++)
    {
        auto group_key{ group_key_list[k] };
        const auto & atom_list{ analysis_view.GetAtomObjectList(group_key) };
        std::vector<LocalPotentialSampleList> member_sample_list;
        std::vector<LocalGaussianResult> member_result_list;
        member_sample_list.reserve(atom_list.size());
        member_result_list.reserve(atom_list.size());
        for (const auto & atom : atom_list)
        {
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            member_sample_list.emplace_back(local_view.GetSamplingEntries(false));
            member_result_list.emplace_back(local_view.GetGaussianResult());
        }
        const auto result{
            EstimateGroupGaussian(
                member_sample_list,
                member_result_list,
                analysis_view.GetAtomAlphaG(group_key),
                options)
        };

#ifdef USE_OPENMP
        #pragma omp critical
#endif
        {
            analysis.ApplyAtomGroupGaussianResult(group_key, result);
            key_count++;
            if (!options.quiet_mode)
            {
                Logger::ProgressBar(key_count, group_key_size);
            }
        }
    }
}

void RunPotentialFittingWorkflow(ModelObject & model_object, const FitOptions & options)
{
    RunLocalAlphaTraining(model_object, options);
    RunLocalPotentialFitting(model_object, options);
    RunGroupAlphaTraining(model_object, options);
    RunGroupPotentialFitting(model_object, options);
}

} // namespace rhbm_gem::core
