#include <cstddef>
#include <rhbm_gem/core/GaussianEstimator.hpp>

#include "detail/LocalFittingThawHysteresisTracker.hpp"

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/AdaptiveRelaxationController.hpp>
#include <rhbm_gem/utils/algorithm/ConvergenceFreezeTracker.hpp>
#include <rhbm_gem/utils/algorithm/IterationState.hpp>
#include <rhbm_gem/utils/algorithm/LinearRegressionSample.hpp>
#include <rhbm_gem/utils/algorithm/NormalizedChange.hpp>
#include <rhbm_gem/utils/algorithm/ParameterChangeStats.hpp>
#include <rhbm_gem/utils/algorithm/RobustSlopeEstimator.hpp>
#include <rhbm_gem/utils/algorithm/RobustSlopeOptions.hpp>
#include <rhbm_gem/utils/algorithm/SparseRegressionRow.hpp>
#include <rhbm_gem/utils/algorithm/WeightedRidgeSolver.hpp>
#include <rhbm_gem/utils/algorithm/WeightedRidgeSystem.hpp>
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
#include <map>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace rhbm_gem::core {
namespace {
constexpr double kLocalFittingNormalizedChangeTolerance{ 1.0e-4 };
constexpr double kLocalFittingNormalizedChangeScaleFloor{ 1.0 };
constexpr std::size_t kMinimumAlphaRTrainingSampleCount{ 10 };
constexpr std::size_t kMinimumAlphaGTrainingMemberCount{ 10 };
constexpr double kResidualOffsetRangeMin{ 1.0 };
constexpr double kResidualOffsetRangeMax{ 2.0 };
constexpr double kOffsetDampingFactor{ 0.5 };
constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr std::size_t kLocalFittingMaximumIterations{ 200 };
constexpr double kLocalFittingParameterChangeTolerance{ 1.0e-6 };
constexpr double kLocalFittingChangePercentile{ 0.95 };
constexpr int kHuberSlopeMaximumIterations{ 50 };
constexpr double kHuberSlopeTolerance{ 1.0e-8 };
constexpr double kHuberScaleMultiplier{ 1.4826 };
constexpr double kHuberScaleMin{ 1.0e-12 };
constexpr double kHuberCutoffMultiplier{ 1.345 };
constexpr double kOffsetRegularizationAmplitudeRatio{ 0.1 };
constexpr double kOffsetRegularizationPriorScaleMin{ 1.0e-12 };
constexpr double kJointOffsetRidgeRatio{ 1.0e-3 };
constexpr double kJointOffsetRidgeRatioMin{ 1.0e-4 };
constexpr double kJointOffsetRidgeRatioMax{ 1.0 };
constexpr double kJointOffsetRidgeGrowth{ 2.0 };
constexpr double kJointOffsetRidgeShrink{ 0.8 };
constexpr double kSuspiciousJointOffsetRidgeMultiplier{ 10.0 };
constexpr double kJointOffsetCollinearityOverlapThreshold{ 0.98 };
constexpr double kCollinearJointOffsetRidgeMultiplier{ 10.0 };
constexpr double kJointOffsetIrlsScaleFloor{ 1.0e-2 };
constexpr double kJointOffsetIrlsNormalizedChangeTolerance{ 1.0e-6 };
constexpr double kJointOffsetIrlsObjectiveRelativeTolerance{ 1.0e-10 };
constexpr double kAdaptiveRelaxationInitialBeta{ 0.5 };
constexpr double kAdaptiveRelaxationMin{ 0.05 };
constexpr double kAdaptiveRelaxationMax{ 1.0 };
constexpr double kAdaptiveRelaxationGrowth{ 1.2 };
constexpr double kAdaptiveRelaxationShrink{ 0.5 };
constexpr double kAdaptiveRelaxationImprovementRatio{ 0.01 };
constexpr int kAdaptiveRelaxationIncreaseStreak{ 2 };
constexpr int kAdaptiveRelaxationShrinkStreak{ 3 };
constexpr double kLocalFittingFreezeChangeRatio{ 0.1 };
constexpr int kLocalFittingFreezeStableIterations{ 3 };
constexpr double kLocalFittingDependencyThawHysteresisGrowth{ 2.0 };
constexpr double kLocalFittingDependencyThawHysteresisMax{ 8.0 };
constexpr double kLocalFittingDependencyThawHysteresisFrozenDecay{ 0.9 };
constexpr double kLocalFittingObjectiveTieRelativeTolerance{ 1.0e-8 };
constexpr double kLocalFittingConvergenceObjectiveRelativeTolerance{ 1.0e-3 };
constexpr int kLocalFittingObjectiveBacktrackingMaximumAttempts{ 3 };
constexpr std::size_t kLocalFittingObjectiveScaleWarmupCount{ 5 };
constexpr double kLocalFittingObjectiveResidualScaleFloorRatio{ 1.0e-6 };

using GaussianFittingState = algorithm::IterationState<LocalGaussianResult, Eigen::VectorXd>;

struct JointOffsetSolveResult
{
    Eigen::VectorXd offset{};
    bool used_fallback{ false };
    std::vector<std::vector<std::size_t>> active_coupling_adjacency{};
};

struct JointOffsetBuildResult
{
    algorithm::WeightedRidgeSystem system{};
    std::vector<std::vector<std::size_t>> active_coupling_adjacency{};
};

double IncreaseLocalFittingRidgeRatio(double ridge_ratio)
{
    return std::min(kJointOffsetRidgeRatioMax, ridge_ratio * kJointOffsetRidgeGrowth);
}

double DecreaseLocalFittingRidgeRatio(double ridge_ratio)
{
    return std::max(kJointOffsetRidgeRatioMin, ridge_ratio * kJointOffsetRidgeShrink);
}

struct LocalRefitResult
{
    LocalGaussianResult result{};
    bool used_fallback{ false };
    bool suspicious_offset_fallback{ false };
};

struct LocalFittingIterationResult
{
    GaussianFittingState state{};
    bool joint_offset_used_fallback{ false };
    std::vector<std::size_t> refit_fallback_state_index_list{};
    std::vector<std::size_t> suspicious_offset_state_index_list{};
};

struct LocalFittingFallbackStats
{
    std::size_t joint_offset_fallback_iterations{ 0 };
    std::size_t refit_fallback_atom_events{ 0 };
    std::size_t suspicious_offset_atom_events{ 0 };
    std::vector<bool> refit_fallback_atom_seen{};
    std::vector<bool> suspicious_offset_atom_seen{};

    explicit LocalFittingFallbackStats(std::size_t atom_size)
        : refit_fallback_atom_seen(atom_size, false),
          suspicious_offset_atom_seen(atom_size, false)
    {
    }

    void Accumulate(const LocalFittingIterationResult & iteration_result)
    {
        if (iteration_result.joint_offset_used_fallback)
        {
            joint_offset_fallback_iterations++;
        }
        refit_fallback_atom_events += iteration_result.refit_fallback_state_index_list.size();
        for (const auto state_index : iteration_result.refit_fallback_state_index_list)
        {
            if (state_index >= refit_fallback_atom_seen.size())
            {
                throw std::invalid_argument("Local fitting fallback atom index is out of range.");
            }
            refit_fallback_atom_seen.at(state_index) = true;
        }
        suspicious_offset_atom_events += iteration_result.suspicious_offset_state_index_list.size();
        for (const auto state_index : iteration_result.suspicious_offset_state_index_list)
        {
            if (state_index >= suspicious_offset_atom_seen.size())
            {
                throw std::invalid_argument(
                    "Local fitting suspicious offset atom index is out of range.");
            }
            suspicious_offset_atom_seen.at(state_index) = true;
        }
    }

    std::size_t GetDistinctRefitFallbackAtomCount() const
    {
        return static_cast<std::size_t>(
            std::count(
                refit_fallback_atom_seen.begin(),
                refit_fallback_atom_seen.end(),
                true));
    }

    std::size_t GetDistinctSuspiciousOffsetAtomCount() const
    {
        return static_cast<std::size_t>(
            std::count(
                suspicious_offset_atom_seen.begin(),
                suspicious_offset_atom_seen.end(),
                true));
    }

    bool HasFallback() const
    {
        return joint_offset_fallback_iterations > 0 ||
            refit_fallback_atom_events > 0 ||
            suspicious_offset_atom_events > 0;
    }
};

struct LocalFittingObjectiveStats
{
    bool has_quality_objective{ false };
    double quality_objective{ std::numeric_limits<double>::infinity() };
    double residual_scale_sample{ std::numeric_limits<double>::infinity() };
};

struct LocalFittingObjectiveReference
{
    bool has_reference{ false };
    double residual_scale{ std::numeric_limits<double>::infinity() };
};

struct LocalFittingObjectiveSamples
{
    std::vector<double> residual_list{};
    std::vector<double> response_list{};
};

struct ParameterSummaryStats
{
    double mean{ 0.0 };
    double standard_deviation{ 0.0 };
};

struct GroupPriorSpotSampleList
{
    std::vector<double> amplitude_list{};
    std::vector<double> width_list{};
    std::vector<double> offset_list{};
};

std::vector<AtomLocalPotentialEditor> BuildAtomLocalEditors(
    ModelObject & model_object,
    const std::vector<AtomObject *> & atom_list)
{
    auto analysis{ model_object.EditAnalysis() };
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
    execution_options.quiet_mode = false;
    execution_options.thread_size = options.thread_size;
    return execution_options;
}

bool CanBuildFiniteZeroOffsetSamples(
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

algorithm::RobustSlopeOptions MakeResidualOffsetSlopeOptions(double amplitude)
{
    algorithm::RobustSlopeOptions options;
    options.maximum_iterations = kHuberSlopeMaximumIterations;
    options.tolerance = kHuberSlopeTolerance;
    options.scale_multiplier = kHuberScaleMultiplier;
    options.scale_min = kHuberScaleMin;
    options.cutoff_multiplier = kHuberCutoffMultiplier;
    options.regularization_prior_scale =
        std::max(
            std::abs(amplitude) * kOffsetRegularizationAmplitudeRatio,
            kOffsetRegularizationPriorScaleMin);
    return options;
}

double EstimateResidualOffsetParameter(
    const LocalPotentialSampleList & sample_entries,
    const RHBMBetaEstimateResult & fit_result,
    double current_offset)
{
    const auto signal_model{
        linearization_service::DecodeParameterVector(fit_result.beta_mdpde)
    };
    const auto width{ signal_model.GetWidth() };
    if (!std::isfinite(width) || width <= 0.0) return current_offset;

    std::vector<algorithm::LinearRegressionSample> residual_sample_list;
    residual_sample_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        if (distance < kResidualOffsetRangeMin || distance > kResidualOffsetRangeMax)
        {
            continue;
        }

        const auto basis{ signal_model.OffsetBasisAtDistance(distance) };
        if (!std::isfinite(basis) ||
            std::abs(basis) <= std::numeric_limits<double>::epsilon())
        {
            continue;
        }
        const auto residual{
            static_cast<double>(sample.response) -
            signal_model.SignalAtDistance(distance)
        };
        if (!std::isfinite(residual)) continue;
        residual_sample_list.emplace_back(algorithm::LinearRegressionSample{ basis, residual });
    }

    double candidate_offset{ current_offset };
    if (!algorithm::RobustSlopeEstimator::EstimateHuberSlopeThroughOrigin(
            residual_sample_list,
            MakeResidualOffsetSlopeOptions(signal_model.GetAmplitude()),
            candidate_offset))
    {
        return current_offset;
    }
    const auto candidate_model{ signal_model.WithOffset(candidate_offset) };
    if (!CanBuildFiniteZeroOffsetSamples(sample_entries, candidate_model))
    {
        return current_offset;
    }
    return candidate_offset;
}

bool IsSuspiciousJointOffset(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & previous_model,
    const GaussianModel3D & offset_model)
{
    return CanBuildFiniteZeroOffsetSamples(sample_entries, previous_model) &&
        !CanBuildFiniteZeroOffsetSamples(sample_entries, offset_model);
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

LocalPotentialSampleList BuildSamplesForZeroOffsetGaussianFit(
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
            BuildSamplesForZeroOffsetGaussianFit(
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
    double offset = 0.0)
{
    const auto ols_model{
        linearization_service::DecodeParameterVector(fit_result.beta_ols)
            .WithOffset(offset)
    };
    const auto mdpde_model{
        linearization_service::DecodeParameterVector(fit_result.beta_mdpde)
            .WithOffset(offset)
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

GaussianModel3DWithUncertainty WithModelOffset(
    const GaussianModel3DWithUncertainty & gaussian,
    double offset)
{
    return GaussianModel3DWithUncertainty{
        gaussian.GetModel().WithOffset(offset),
        gaussian.GetStandardDeviationModel()
    };
}

GroupGaussianResult DecodeGroupGaussianResult(
    double alpha_g,
    const RHBMGroupEstimationResult & result,
    double offset)
{
    return GroupGaussianResult{
        alpha_g,
        linearization_service::DecodeParameterVector(result.mu_mean).WithOffset(offset),
        linearization_service::DecodeParameterVector(result.mu_mdpde).WithOffset(offset),
        WithModelOffset(
            linearization_service::DecodeParameterVector(result.mu_prior, result.capital_lambda),
            offset)
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
        const auto offset{
            member_result_list.at(member_index).mdpde.GetModel().GetOffset()
        };
        const auto gaussian{
            WithModelOffset(
                linearization_service::DecodeParameterVector(
                    result.beta_posterior_matrix.col(i),
                    result.capital_sigma_posterior_list.at(member_index)),
                offset)
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
using SpotMedianModelMap = std::unordered_map<Spot, GaussianModel3D>;

struct GaussianModelParameterSamples
{
    std::vector<double> amplitude_list{};
    std::vector<double> width_list{};
    std::vector<double> offset_list{};
};

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

SpotMedianModelMap BuildSpotMedianMDPDEModelMap(const std::vector<AtomObject *> & atom_list)
{
    std::unordered_map<Spot, GaussianModelParameterSamples> parameter_samples_by_spot;
    parameter_samples_by_spot.reserve(atom_list.size());
    for (const auto * atom : atom_list)
    {
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        if (!local_view.IsAvailable()) continue;

        const auto & model{ local_view.GetEstimateMDPDE() };
        auto & parameter_samples{ parameter_samples_by_spot[atom->GetSpot()] };
        parameter_samples.amplitude_list.emplace_back(model.GetAmplitude());
        parameter_samples.width_list.emplace_back(model.GetWidth());
        parameter_samples.offset_list.emplace_back(model.GetOffset());
    }

    SpotMedianModelMap median_model_by_spot;
    median_model_by_spot.reserve(parameter_samples_by_spot.size());
    for (const auto & [spot, parameter_samples] : parameter_samples_by_spot)
    {
        if (parameter_samples.amplitude_list.empty()) continue;

        median_model_by_spot.emplace(
            spot,
            GaussianModel3D{
                array_helper::ComputeMedian(parameter_samples.amplitude_list),
                array_helper::ComputeMedian(parameter_samples.width_list),
                array_helper::ComputeMedian(parameter_samples.offset_list)
            });
    }
    return median_model_by_spot;
}

JointOffsetBuildResult BuildJointOffsetSystem(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<std::size_t> & active_index_list,
    const FittedGaussianSnapshot & snapshot,
    double ridge_ratio,
    const std::vector<double> & ridge_multiplier_list)
{
    if (atom_list.size() != ridge_multiplier_list.size())
    {
        throw std::invalid_argument("Joint offset ridge multiplier size is inconsistent.");
    }

    std::unordered_map<const AtomObject *, Eigen::Index> atom_column_map;
    atom_column_map.reserve(active_index_list.size());
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        if (atom_index >= atom_list.size())
        {
            throw std::invalid_argument("Joint offset active atom index is out of range.");
        }
        atom_column_map.emplace(atom_list.at(atom_index), static_cast<Eigen::Index>(i));
    }

    std::vector<algorithm::SparseRegressionRow> row_list;
    for (const auto active_index : active_index_list)
    {
        if (active_index >= atom_list.size())
        {
            throw std::invalid_argument("Joint offset active atom index is out of range.");
        }
        const auto * atom{ atom_list.at(active_index) };
        const auto model_iter{ snapshot.find(atom) };
        if (model_iter == snapshot.end())
        {
            throw std::invalid_argument("Joint offset snapshot is missing an atom.");
        }
        const auto target_column{ atom_column_map.at(atom) };
        const auto & target_model{ model_iter->second };
        const auto sample_entries{
            AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false)
        };
        const auto neighbor_atom_list{ atom->FindNeighborAtoms() };
        for (const auto & sample : sample_entries)
        {
            if (!std::isfinite(static_cast<double>(sample.response)))
            {
                throw std::runtime_error("Joint offset sample response is not finite.");
            }
            const auto target_distance{ static_cast<double>(sample.point.distance) };
            const auto target_signal{ target_model.SignalAtDistance(target_distance) };
            const auto target_basis{ target_model.OffsetBasisAtDistance(target_distance) };
            if (!std::isfinite(target_signal) || !std::isfinite(target_basis))
            {
                throw std::runtime_error("Joint offset target model evaluation is not finite.");
            }
            auto residual{ static_cast<double>(sample.response) - target_signal };
            algorithm::SparseRegressionRow row;
            if (std::abs(target_basis) > std::numeric_limits<double>::epsilon())
            {
                row.basis_entries.emplace_back(target_column, target_basis);
            }

            for (const auto * neighbor_atom : neighbor_atom_list)
            {
                const auto neighbor_model_iter{ snapshot.find(neighbor_atom) };
                if (neighbor_model_iter == snapshot.end()) continue;

                const auto distance{
                    static_cast<double>(
                        array_helper::ComputeNorm<float>(sample.point.position, neighbor_atom->GetPosition()))
                };
                if (distance > kNeighborContributionDistanceMax) continue;
                const auto & neighbor_model{ neighbor_model_iter->second };
                const auto column_iter{ atom_column_map.find(neighbor_atom) };
                if (column_iter == atom_column_map.end())
                {
                    const auto response{ neighbor_model.ResponseAtDistance(distance) };
                    if (!std::isfinite(response))
                    {
                        throw std::runtime_error(
                            "Joint offset fixed neighbor model evaluation is not finite.");
                    }
                    residual -= response;
                    continue;
                }

                const auto signal{ neighbor_model.SignalAtDistance(distance) };
                const auto basis{ neighbor_model.OffsetBasisAtDistance(distance) };
                if (!std::isfinite(signal) || !std::isfinite(basis))
                {
                    throw std::runtime_error(
                        "Joint offset active neighbor model evaluation is not finite.");
                }
                residual -= signal;
                if (std::abs(basis) > std::numeric_limits<double>::epsilon())
                {
                    row.basis_entries.emplace_back(column_iter->second, basis);
                }
            }
            if (!std::isfinite(residual))
            {
                throw std::runtime_error("Joint offset residual is not finite.");
            }
            if (row.basis_entries.empty()) continue;
            row.response = residual;
            row_list.emplace_back(std::move(row));
        }
    }

    const auto row_count{ static_cast<Eigen::Index>(row_list.size()) };
    const auto column_count{ static_cast<Eigen::Index>(active_index_list.size()) };
    std::vector<Eigen::Triplet<double>> triplet_list;
    Eigen::VectorXd response{ Eigen::VectorXd::Zero(row_count) };
    Eigen::VectorXd column_square_sum{ Eigen::VectorXd::Zero(column_count) };
    std::map<std::pair<Eigen::Index, Eigen::Index>, double> column_cross_sum_map;
    std::vector<std::vector<std::size_t>> active_coupling_adjacency(active_index_list.size());
    for (Eigen::Index row_index = 0; row_index < row_count; row_index++)
    {
        const auto & row{ row_list.at(static_cast<std::size_t>(row_index)) };
        response(row_index) = row.response;
        for (const auto & [column_index, basis] : row.basis_entries)
        {
            triplet_list.emplace_back(row_index, column_index, basis);
            column_square_sum(column_index) += basis * basis;
        }
        for (std::size_t i = 0; i < row.basis_entries.size(); i++)
        {
            const auto [left_column, left_basis]{ row.basis_entries.at(i) };
            for (std::size_t j = i + 1; j < row.basis_entries.size(); j++)
            {
                const auto [right_column, right_basis]{ row.basis_entries.at(j) };
                if (left_column == right_column) continue;
                const auto column_pair{
                    std::minmax(left_column, right_column)
                };
                column_cross_sum_map[column_pair] += left_basis * right_basis;
                const auto left_index{ static_cast<std::size_t>(column_pair.first) };
                const auto right_index{ static_cast<std::size_t>(column_pair.second) };
                active_coupling_adjacency.at(left_index).emplace_back(right_index);
                active_coupling_adjacency.at(right_index).emplace_back(left_index);
            }
        }
    }
    for (auto & neighbor_list : active_coupling_adjacency)
    {
        std::sort(neighbor_list.begin(), neighbor_list.end());
        neighbor_list.erase(
            std::unique(neighbor_list.begin(), neighbor_list.end()),
            neighbor_list.end());
    }

    Eigen::VectorXd proactive_ridge_multiplier{
        Eigen::VectorXd::Ones(column_count)
    };
    for (const auto & [column_pair, cross_sum] : column_cross_sum_map)
    {
        const auto left_column{ column_pair.first };
        const auto right_column{ column_pair.second };
        const auto left_square_sum{ column_square_sum(left_column) };
        const auto right_square_sum{ column_square_sum(right_column) };
        if (left_square_sum <= std::numeric_limits<double>::epsilon() ||
            right_square_sum <= std::numeric_limits<double>::epsilon())
        {
            continue;
        }
        const auto overlap{
            std::abs(cross_sum) / std::sqrt(left_square_sum * right_square_sum)
        };
        if (!std::isfinite(overlap) ||
            overlap < kJointOffsetCollinearityOverlapThreshold)
        {
            continue;
        }
        proactive_ridge_multiplier(left_column) = std::max(
            proactive_ridge_multiplier(left_column),
            kCollinearJointOffsetRidgeMultiplier);
        proactive_ridge_multiplier(right_column) = std::max(
            proactive_ridge_multiplier(right_column),
            kCollinearJointOffsetRidgeMultiplier);
    }

    algorithm::WeightedRidgeSystem system;
    system.design_matrix.resize(row_count, column_count);
    system.design_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
    system.response = std::move(response);
    system.previous_parameter = Eigen::VectorXd::Zero(column_count);
    system.ridge_diagonal = Eigen::VectorXd::Zero(column_count);
    for (Eigen::Index column_index = 0; column_index < column_count; column_index++)
    {
        const auto atom_index{ active_index_list.at(static_cast<std::size_t>(column_index)) };
        const auto & model{ snapshot.at(atom_list.at(atom_index)) };
        system.previous_parameter(column_index) = model.GetOffset();
        const auto square_sum{ column_square_sum(column_index) };
        const auto multiplier{
            ridge_multiplier_list.at(atom_index)
        };
        if (!std::isfinite(multiplier) || multiplier <= 0.0)
        {
            throw std::invalid_argument("Joint offset ridge multiplier must be positive and finite.");
        }
        const auto combined_multiplier{
            std::max(multiplier, proactive_ridge_multiplier(column_index))
        };
        const auto base_ridge{
            square_sum > std::numeric_limits<double>::epsilon()
                ? ridge_ratio * square_sum
                : ridge_ratio / kJointOffsetRidgeRatio
        };
        system.ridge_diagonal(column_index) = combined_multiplier * base_ridge;
    }
    return JointOffsetBuildResult{
        std::move(system),
        std::move(active_coupling_adjacency)
    };
}

double CalculateWeightedRidgeSurrogateObjective(
    const algorithm::WeightedRidgeSystem & system,
    const Eigen::VectorXd & weight,
    const Eigen::VectorXd & offset)
{
    if (system.response.size() != weight.size())
    {
        throw std::invalid_argument("Weighted ridge objective input sizes are inconsistent.");
    }
    if (system.previous_parameter.size() != offset.size() ||
        system.ridge_diagonal.size() != offset.size())
    {
        throw std::invalid_argument("Weighted ridge objective parameter sizes are inconsistent.");
    }
    if (system.response.size() == 0)
    {
        return std::numeric_limits<double>::infinity();
    }

    const Eigen::VectorXd residual{ system.response - system.design_matrix * offset };
    const auto weighted_residual_loss{
        weight.cwiseProduct(residual.cwiseAbs2()).sum()
    };
    const Eigen::VectorXd offset_delta{ offset - system.previous_parameter };
    const auto ridge_loss{
        system.ridge_diagonal.cwiseProduct(offset_delta.cwiseAbs2()).sum()
    };
    const auto objective{
        (weighted_residual_loss + ridge_loss) /
        static_cast<double>(system.response.size())
    };
    return std::isfinite(objective) ? objective : std::numeric_limits<double>::infinity();
}

bool IsJointOffsetObjectiveDeteriorated(
    double updated_objective,
    double current_objective)
{
    if (!std::isfinite(updated_objective))
    {
        return true;
    }
    if (!std::isfinite(current_objective))
    {
        return false;
    }
    const auto scale{
        std::max({ std::abs(updated_objective), std::abs(current_objective), 1.0 })
    };
    return updated_objective >
        current_objective + kJointOffsetIrlsObjectiveRelativeTolerance * scale;
}

JointOffsetSolveResult EstimateJointOffsets(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<std::size_t> & active_index_list,
    const FittedGaussianSnapshot & snapshot,
    double ridge_ratio,
    const std::vector<double> & ridge_multiplier_list)
{
    Eigen::VectorXd previous_offset{
        Eigen::VectorXd::Zero(static_cast<Eigen::Index>(active_index_list.size()))
    };
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        if (atom_index >= atom_list.size())
        {
            throw std::invalid_argument("Joint offset active atom index is out of range.");
        }
        previous_offset(static_cast<Eigen::Index>(i)) =
            snapshot.at(atom_list.at(atom_index)).GetOffset();
    }
    JointOffsetBuildResult build_result;
    try
    {
        build_result = BuildJointOffsetSystem(
            atom_list,
            active_index_list,
            snapshot,
            ridge_ratio,
            ridge_multiplier_list);
    }
    catch (const std::exception &)
    {
        return JointOffsetSolveResult{
            previous_offset,
            true,
            {}
        };
    }
    auto system{ std::move(build_result.system) };
    auto active_coupling_adjacency{ std::move(build_result.active_coupling_adjacency) };
    if (system.response.size() == 0 || system.previous_parameter.size() == 0)
    {
        return JointOffsetSolveResult{
            previous_offset,
            true,
            std::move(active_coupling_adjacency)
        };
    }

    Eigen::VectorXd weight{ Eigen::VectorXd::Ones(system.response.size()) };
    algorithm::WeightedRidgeSolver solver{ system };
    Eigen::VectorXd offset;
    if (!solver.Solve(system, weight, offset))
    {
        return JointOffsetSolveResult{
            previous_offset,
            true,
            std::move(active_coupling_adjacency)
        };
    }

    for (int iteration = 0; iteration < kHuberSlopeMaximumIterations; iteration++)
    {
        const Eigen::VectorXd residual{ system.response - system.design_matrix * offset };
        std::vector<double> residual_list(residual.data(), residual.data() + residual.size());
        const auto median_residual{ array_helper::ComputeMedian(residual_list) };
        std::vector<double> deviation_list;
        deviation_list.reserve(residual_list.size());
        for (const auto value : residual_list)
        {
            deviation_list.emplace_back(std::abs(value - median_residual));
        }
        const auto residual_scale{ std::max(
            kHuberScaleMultiplier * array_helper::ComputeMedian(deviation_list),
            kHuberScaleMin)
        };
        const auto cutoff{ kHuberCutoffMultiplier * residual_scale };
        for (Eigen::Index i = 0; i < residual.size(); i++)
        {
            const auto absolute_residual{ std::abs(residual(i)) };
            weight(i) = absolute_residual <= cutoff ? 1.0 : cutoff / absolute_residual;
        }

        Eigen::VectorXd updated_offset;
        if (!solver.Solve(system, weight, updated_offset))
        {
            return JointOffsetSolveResult{
                system.previous_parameter,
                true,
                std::move(active_coupling_adjacency)
            };
        }
        const auto current_objective{
            CalculateWeightedRidgeSurrogateObjective(system, weight, offset)
        };
        const auto updated_objective{
            CalculateWeightedRidgeSurrogateObjective(system, weight, updated_offset)
        };
        if (IsJointOffsetObjectiveDeteriorated(updated_objective, current_objective))
        {
            break;
        }
        const auto maximum_change{
            algorithm::CalculateMaximumNormalizedVectorChange(
                updated_offset,
                offset,
                kJointOffsetIrlsScaleFloor)
        };
        offset = std::move(updated_offset);
        if (maximum_change < kJointOffsetIrlsNormalizedChangeTolerance) break;
    }

    return JointOffsetSolveResult{
        offset,
        false,
        std::move(active_coupling_adjacency)
    };
}

void ApplyJointOffsetsToSnapshot(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<std::size_t> & active_index_list,
    const Eigen::VectorXd & offset,
    FittedGaussianSnapshot & snapshot)
{
    if (active_index_list.size() != static_cast<std::size_t>(offset.size()))
    {
        throw std::invalid_argument("Joint offset result size is inconsistent.");
    }
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        if (atom_index >= atom_list.size())
        {
            throw std::invalid_argument("Joint offset active atom index is out of range.");
        }
        snapshot.at(atom_list.at(atom_index)) =
            snapshot.at(atom_list.at(atom_index)).WithOffset(offset(static_cast<Eigen::Index>(i)));
    }
}

template <typename GaussianLookup>
LocalPotentialSampleList UpdateSampleListWithGaussianLookup(
    const AtomObject & atom,
    GaussianLookup lookup_gaussian)
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
            const auto * gaussian{ lookup_gaussian(*neighbor_atom) };
            if (gaussian == nullptr) continue;

            auto neighbor_position{ neighbor_atom->GetPosition() };
            auto distance{
                static_cast<double>(
                    array_helper::ComputeNorm<float>(sample_position, neighbor_position))
            };
            if (distance > kNeighborContributionDistanceMax) continue;
            response_value -= static_cast<float>(gaussian->ResponseAtDistance(distance));
        }
        updated_list.emplace_back(LocalPotentialSample{response_value, sample.point });
    }
    return updated_list;
}

LocalPotentialSampleList UpdateSampleListWithFittedGaussian(
    const AtomObject & atom,
    const FittedGaussianSnapshot & snapshot)
{
    return UpdateSampleListWithGaussianLookup(
        atom,
        [&snapshot](const AtomObject & neighbor_atom) -> const GaussianModel3D *
        {
            const auto gaussian_iter{ snapshot.find(&neighbor_atom) };
            return gaussian_iter == snapshot.end() ? nullptr : &gaussian_iter->second;
        });
}

LocalPotentialSampleList UpdateSampleListWithFittedGaussian(const AtomObject & atom)
{
    return UpdateSampleListWithGaussianLookup(
        atom,
        [](const AtomObject & neighbor_atom) -> const GaussianModel3D *
        {
            const auto local_view{ AtomLocalPotentialView::For(neighbor_atom) };
            return local_view.IsAvailable() ? &local_view.GetEstimateMDPDE() : nullptr;
        });
}

LocalPotentialSampleList UpdateSampleListWithSpotMedianGaussian(
    const AtomObject & atom,
    const SpotMedianModelMap & median_model_by_spot)
{
    return UpdateSampleListWithGaussianLookup(
        atom,
        [&median_model_by_spot](const AtomObject & neighbor_atom) -> const GaussianModel3D *
        {
            const auto median_model_iter{ median_model_by_spot.find(neighbor_atom.GetSpot()) };
            if (median_model_iter != median_model_by_spot.end())
            {
                return &median_model_iter->second;
            }

            const auto local_view{ AtomLocalPotentialView::For(neighbor_atom) };
            return local_view.IsAvailable() ? &local_view.GetEstimateMDPDE() : nullptr;
        });
}

double CalculateHuberLoss(double residual, double cutoff)
{
    const auto absolute_residual{ std::abs(residual) };
    if (absolute_residual <= cutoff)
    {
        return 0.5 * residual * residual;
    }
    return cutoff * (absolute_residual - 0.5 * cutoff);
}

std::optional<LocalFittingObjectiveSamples> CollectLocalFittingObjectiveSamples(
    const std::vector<AtomObject *> & atom_list,
    const GaussianFittingState & fitting_state)
{
    const auto snapshot{
        BuildFittedGaussianSnapshot(atom_list, fitting_state.estimation_list)
    };
    LocalFittingObjectiveSamples objective_samples;
    for (const auto * atom : atom_list)
    {
        const auto model_iter{ snapshot.find(atom) };
        if (model_iter == snapshot.end())
        {
            throw std::invalid_argument("Local fitting objective snapshot is missing an atom.");
        }
        const auto sample_entries{ UpdateSampleListWithFittedGaussian(*atom, snapshot) };
        objective_samples.residual_list.reserve(
            objective_samples.residual_list.size() + sample_entries.size());
        objective_samples.response_list.reserve(
            objective_samples.response_list.size() + sample_entries.size());
        const auto & target_model{ model_iter->second };
        for (const auto & sample : sample_entries)
        {
            const auto distance{ static_cast<double>(sample.point.distance) };
            const auto expected_response{ target_model.ResponseAtDistance(distance) };
            const auto response{ static_cast<double>(sample.response) };
            const auto residual{ response - expected_response };
            if (!std::isfinite(response) || !std::isfinite(residual))
            {
                return std::nullopt;
            }
            objective_samples.residual_list.emplace_back(residual);
            objective_samples.response_list.emplace_back(response);
        }
    }
    return objective_samples;
}

double CalculateMedianAbsoluteDeviationScale(const std::vector<double> & value_list)
{
    if (value_list.empty())
    {
        return std::numeric_limits<double>::infinity();
    }

    const auto median_value{ array_helper::ComputeMedian(value_list) };
    std::vector<double> deviation_list;
    deviation_list.reserve(value_list.size());
    for (const auto value : value_list)
    {
        deviation_list.emplace_back(std::abs(value - median_value));
    }
    return kHuberScaleMultiplier * array_helper::ComputeMedian(deviation_list);
}

std::optional<double> CalculateLocalFittingResidualScaleSample(
    const LocalFittingObjectiveSamples & objective_samples)
{
    if (objective_samples.residual_list.empty() ||
        objective_samples.residual_list.size() != objective_samples.response_list.size())
    {
        return std::nullopt;
    }

    const auto residual_scale{
        CalculateMedianAbsoluteDeviationScale(objective_samples.residual_list)
    };
    const auto response_scale_floor{
        kLocalFittingObjectiveResidualScaleFloorRatio *
        CalculateMedianAbsoluteDeviationScale(objective_samples.response_list)
    };
    const auto scale_sample{
        std::max({ residual_scale, response_scale_floor, kHuberScaleMin })
    };
    if (!std::isfinite(scale_sample))
    {
        return std::nullopt;
    }
    return scale_sample;
}

LocalFittingObjectiveStats CalculateLocalFittingObjectiveStats(
    const LocalFittingObjectiveSamples & objective_samples,
    const LocalFittingObjectiveReference & objective_reference)
{
    LocalFittingObjectiveStats stats;
    if (!objective_reference.has_reference)
    {
        return stats;
    }

    const auto residual_scale_sample{
        CalculateLocalFittingResidualScaleSample(objective_samples)
    };
    if (!residual_scale_sample.has_value())
    {
        return stats;
    }

    double loss_sum{ 0.0 };
    for (const auto residual : objective_samples.residual_list)
    {
        const auto normalized_residual{ residual / objective_reference.residual_scale };
        loss_sum += CalculateHuberLoss(normalized_residual, kHuberCutoffMultiplier);
    }
    const auto quality_objective{
        loss_sum / static_cast<double>(objective_samples.residual_list.size())
    };
    if (!std::isfinite(quality_objective))
    {
        return stats;
    }

    stats.has_quality_objective = true;
    stats.quality_objective = quality_objective;
    stats.residual_scale_sample = *residual_scale_sample;
    return stats;
}

LocalFittingObjectiveStats CalculateLocalFittingObjectiveStats(
    const std::vector<AtomObject *> & atom_list,
    const GaussianFittingState & fitting_state,
    const LocalFittingObjectiveReference & objective_reference)
{
    const auto objective_samples{
        CollectLocalFittingObjectiveSamples(atom_list, fitting_state)
    };
    if (!objective_samples.has_value())
    {
        return LocalFittingObjectiveStats{};
    }
    return CalculateLocalFittingObjectiveStats(*objective_samples, objective_reference);
}

class LocalFittingObjectiveScaleTracker
{
    std::vector<double> m_scale_sample_list{};
    std::size_t m_accepted_scale_sample_count{ 0 };
    bool m_locked{ false };

    static void AddScaleSample(
        std::vector<double> & scale_sample_list,
        double scale_sample)
    {
        if (!std::isfinite(scale_sample)) return;

        scale_sample_list.emplace_back(scale_sample);
        while (scale_sample_list.size() > kLocalFittingObjectiveScaleWarmupCount)
        {
            scale_sample_list.erase(scale_sample_list.begin());
        }
    }

    static LocalFittingObjectiveReference BuildReference(
        const std::vector<double> & scale_sample_list)
    {
        LocalFittingObjectiveReference reference;
        if (scale_sample_list.empty())
        {
            return reference;
        }

        const auto scale_sum{
            std::accumulate(scale_sample_list.begin(), scale_sample_list.end(), 0.0)
        };
        const auto residual_scale{
            scale_sum / static_cast<double>(scale_sample_list.size())
        };
        if (!std::isfinite(residual_scale))
        {
            return reference;
        }

        reference.has_reference = true;
        reference.residual_scale = residual_scale;
        return reference;
    }

public:
    explicit LocalFittingObjectiveScaleTracker(
        std::optional<double> initial_scale_sample)
    {
        if (initial_scale_sample.has_value())
        {
            AddScaleSample(m_scale_sample_list, *initial_scale_sample);
        }
    }

    bool HasReference() const
    {
        return !m_scale_sample_list.empty();
    }

    bool IsLocked() const
    {
        return m_locked;
    }

    LocalFittingObjectiveReference GetCommittedReference() const
    {
        return BuildReference(m_scale_sample_list);
    }

    LocalFittingObjectiveReference GetProvisionalReference(double scale_sample) const
    {
        if (m_locked)
        {
            return GetCommittedReference();
        }

        auto scale_sample_list{ m_scale_sample_list };
        AddScaleSample(scale_sample_list, scale_sample);
        return BuildReference(scale_sample_list);
    }

    void CommitScaleSample(double scale_sample)
    {
        if (m_locked) return;

        AddScaleSample(m_scale_sample_list, scale_sample);
        m_accepted_scale_sample_count++;
        if (m_accepted_scale_sample_count >= kLocalFittingObjectiveScaleWarmupCount)
        {
            m_locked = true;
        }
    }
};

ParameterSummaryStats SummarizeParameterValues(const std::vector<double> & value_list)
{
    if (value_list.empty())
    {
        return {};
    }

    double sum{ 0.0 };
    for (const auto value : value_list)
    {
        sum += value;
    }
    const auto mean{ sum / static_cast<double>(value_list.size()) };
    if (value_list.size() < 2)
    {
        return ParameterSummaryStats{ mean, 0.0 };
    }

    double squared_error_sum{ 0.0 };
    for (const auto value : value_list)
    {
        const auto error{ value - mean };
        squared_error_sum += error * error;
    }
    return ParameterSummaryStats{
        mean,
        std::sqrt(squared_error_sum / static_cast<double>(value_list.size() - 1))
    };
}

std::vector<std::string> BuildGroupPriorSpotSummaryLines(const ModelObject & model_object)
{
    const auto analysis_view{ model_object.GetAnalysisView() };
    std::map<std::string, GroupPriorSpotSampleList> spot_sample_map;
    for (const auto group_key : analysis_view.CollectAtomGroupKeys())
    {
        const auto & atom_list{ analysis_view.GetAtomObjectList(group_key) };
        if (atom_list.empty()) continue;

        const auto & prior{ analysis_view.GetAtomGroupPrior(group_key) };
        auto & sample_list{
            spot_sample_map["Spot::" + ChemicalDataHelper::GetLabel(atom_list.front()->GetSpot())]
        };
        sample_list.amplitude_list.emplace_back(prior.GetAmplitude());
        sample_list.width_list.emplace_back(prior.GetWidth());
        sample_list.offset_list.emplace_back(prior.GetOffset());
    }

    std::vector<std::string> summary_lines;
    summary_lines.reserve(spot_sample_map.size());
    for (const auto & [spot_label, sample_list] : spot_sample_map)
    {
        if (spot_label != "Spot::CA" && spot_label != "Spot::C" && spot_label != "Spot::N" && spot_label != "Spot::O")
        {
            continue;
        }
        const auto amplitude_stats{ SummarizeParameterValues(sample_list.amplitude_list) };
        const auto width_stats{ SummarizeParameterValues(sample_list.width_list) };
        const auto offset_stats{ SummarizeParameterValues(sample_list.offset_list) };

        std::ostringstream stream;
        stream << spot_label << std::fixed << std::setprecision(2)
            << " , amplitude mean = " << amplitude_stats.mean
            << ", amplitude s.d. = " << amplitude_stats.standard_deviation
            << ", width mean = " << width_stats.mean
            << ", width s.d. = " << width_stats.standard_deviation
            << ", offset mean = " << offset_stats.mean
            << ", offset s.d. = " << offset_stats.standard_deviation;
        summary_lines.emplace_back(stream.str());
    }
    return summary_lines;
}

void LogGroupPriorSpotSummary(const ModelObject & model_object)
{
    const auto summary_lines{ BuildGroupPriorSpotSummaryLines(model_object) };
    if (summary_lines.empty())
    {
        Logger::Log(LogLevel::Info,
            "Group fitting prior summary by Spot: no atom groups available.");
        return;
    }

    Logger::Log(LogLevel::Info, "Group fitting prior summary by Spot:");
    for (const auto & line : summary_lines)
    {
        Logger::Log(LogLevel::Info, line);
    }
}

std::vector<algorithm::ParameterChange> CalculateLocalFittingParameterChanges(
    const std::vector<Eigen::VectorXd> & current_estimation_list,
    const std::vector<Eigen::VectorXd> & previous_estimation_list)
{
    if (current_estimation_list.size() != previous_estimation_list.size())
    {
        throw std::invalid_argument("Local fitting parameter change input sizes are inconsistent.");
    }

    std::vector<algorithm::ParameterChange> change_list(current_estimation_list.size());
    for (size_t i = 0; i < current_estimation_list.size(); i++)
    {
        const auto parameter_delta{ current_estimation_list[i] - previous_estimation_list[i] };
        change_list.at(i).value_list = {
            std::abs(parameter_delta(GaussianModel3D::AmplitudeIndex())),
            std::abs(parameter_delta(GaussianModel3D::WidthIndex())),
            std::abs(parameter_delta(GaussianModel3D::OffsetIndex()))
        };
    }
    return change_list;
}

std::vector<algorithm::ParameterChange> CalculateLocalFittingNormalizedParameterChanges(
    const std::vector<Eigen::VectorXd> & current_estimation_list,
    const std::vector<Eigen::VectorXd> & previous_estimation_list)
{
    if (current_estimation_list.size() != previous_estimation_list.size())
    {
        throw std::invalid_argument(
            "Local fitting normalized parameter change input sizes are inconsistent.");
    }

    std::vector<algorithm::ParameterChange> change_list(current_estimation_list.size());
    for (size_t i = 0; i < current_estimation_list.size(); i++)
    {
        change_list.at(i).value_list = {
            algorithm::CalculateNormalizedChange(
                current_estimation_list[i](GaussianModel3D::AmplitudeIndex()),
                previous_estimation_list[i](GaussianModel3D::AmplitudeIndex()),
                kLocalFittingNormalizedChangeScaleFloor),
            algorithm::CalculateNormalizedChange(
                current_estimation_list[i](GaussianModel3D::WidthIndex()),
                previous_estimation_list[i](GaussianModel3D::WidthIndex()),
                kLocalFittingNormalizedChangeScaleFloor),
            algorithm::CalculateNormalizedChange(
                current_estimation_list[i](GaussianModel3D::OffsetIndex()),
                previous_estimation_list[i](GaussianModel3D::OffsetIndex()),
                kLocalFittingNormalizedChangeScaleFloor)
        };
    }
    return change_list;
}

bool IsLocalFittingNormalizedParameterChangeConverged(
    const algorithm::ParameterChangeStats & stats)
{
    if (stats.percentile_list.size() !=
        static_cast<std::size_t>(GaussianModel3D::ParameterSize()))
    {
        throw std::invalid_argument(
            "Local fitting normalized parameter change stats size is inconsistent.");
    }
    for (std::size_t i = 0; i < stats.percentile_list.size(); i++)
    {
        if (stats.percentile_list.at(i) >= kLocalFittingNormalizedChangeTolerance)
        {
            return false;
        }
    }
    return true;
}

void ApplyLocalFittingUnderRelaxation(
    GaussianFittingState & current_state,
    const GaussianFittingState & previous_state,
    double beta)
{
    if (current_state.estimation_list.size() != previous_state.estimation_list.size() ||
        current_state.result_list.size() != previous_state.result_list.size() ||
        current_state.result_list.size() != previous_state.estimation_list.size())
    {
        throw std::invalid_argument("Local fitting relaxation input sizes are inconsistent.");
    }
    for (std::size_t i = 0; i < current_state.estimation_list.size(); i++)
    {
        auto relaxed_estimation{
            (beta * current_state.estimation_list.at(i) +
            (1.0 - beta) * previous_state.estimation_list.at(i)).eval()
        };
        const auto relaxed_model{ GaussianModel3D::FromVector(relaxed_estimation) };
        auto & result{ current_state.result_list.at(i) };
        result.mdpde = GaussianModel3DWithUncertainty{
            relaxed_model,
            result.mdpde.GetStandardDeviationModel()
        };
        current_state.estimation_list.at(i) = relaxed_estimation;
    }
}

struct LocalFittingAttemptResult
{
    GaussianFittingState state{};
    std::vector<algorithm::ParameterChange> change_list{};
    algorithm::ParameterChangeStats change_stats{};
    algorithm::ParameterChangeStats normalized_change_stats{};
    algorithm::FittingQualityCandidateStats candidate_stats{};
    algorithm::FittingQualityCandidateStats previous_candidate_stats{};
    algorithm::FittingQualityCandidateStats best_candidate_stats{};
    double objective_scale_sample{ std::numeric_limits<double>::infinity() };
    algorithm::FittingQualityBacktrackingDecision backtracking_decision{ true, false, false };
};

LocalFittingAttemptResult EvaluateLocalFittingAttempt(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<std::size_t> & active_index_list,
    const GaussianFittingState & raw_state,
    const GaussianFittingState & previous_state,
    const GaussianFittingState & best_state,
    const algorithm::FittingQualityCandidateStats & previous_candidate_stats,
    const algorithm::FittingQualityCandidateStats & best_candidate_stats,
    bool has_best_candidate,
    const LocalFittingObjectiveScaleTracker & objective_scale_tracker,
    double beta,
    int attempt)
{
    LocalFittingAttemptResult result;
    result.state = raw_state;
    ApplyLocalFittingUnderRelaxation(result.state, previous_state, beta);
    result.change_list = CalculateLocalFittingParameterChanges(
        result.state.estimation_list,
        previous_state.estimation_list);
    result.change_stats = algorithm::SummarizeParameterChangeStats(
        result.change_list,
        active_index_list,
        kLocalFittingChangePercentile);

    const auto normalized_change_list{
        CalculateLocalFittingNormalizedParameterChanges(
            result.state.estimation_list,
            previous_state.estimation_list)
    };
    result.normalized_change_stats = algorithm::SummarizeParameterChangeStats(
        normalized_change_list,
        active_index_list,
        kLocalFittingChangePercentile);

    auto objective_reference{ objective_scale_tracker.GetCommittedReference() };
    auto objective_stats{ LocalFittingObjectiveStats{} };
    result.previous_candidate_stats = previous_candidate_stats;
    result.best_candidate_stats = best_candidate_stats;

    const auto objective_samples{
        CollectLocalFittingObjectiveSamples(atom_list, result.state)
    };
    if (objective_scale_tracker.HasReference() && objective_samples.has_value())
    {
        const auto scale_sample{
            CalculateLocalFittingResidualScaleSample(*objective_samples)
        };
        if (scale_sample.has_value())
        {
            objective_reference = objective_scale_tracker.GetProvisionalReference(*scale_sample);
            objective_stats = CalculateLocalFittingObjectiveStats(
                *objective_samples,
                objective_reference);
            if (objective_stats.has_quality_objective)
            {
                const auto recalculated_previous_objective_stats{
                    CalculateLocalFittingObjectiveStats(
                        atom_list,
                        previous_state,
                        objective_reference)
                };
                result.previous_candidate_stats.has_quality_objective =
                    recalculated_previous_objective_stats.has_quality_objective;
                result.previous_candidate_stats.quality_objective =
                    recalculated_previous_objective_stats.quality_objective;
                if (has_best_candidate)
                {
                    const auto recalculated_best_objective_stats{
                        CalculateLocalFittingObjectiveStats(
                            atom_list,
                            best_state,
                            objective_reference)
                    };
                    result.best_candidate_stats.has_quality_objective =
                        recalculated_best_objective_stats.has_quality_objective;
                    result.best_candidate_stats.quality_objective =
                        recalculated_best_objective_stats.quality_objective;
                }
            }
        }
    }

    result.candidate_stats = algorithm::FittingQualityCandidateStats{
        objective_stats.has_quality_objective,
        objective_stats.quality_objective,
        result.normalized_change_stats
    };
    result.objective_scale_sample = objective_stats.residual_scale_sample;
    if (objective_reference.has_reference)
    {
        result.backtracking_decision = algorithm::EvaluateFittingQualityBacktracking(
            result.candidate_stats,
            result.previous_candidate_stats,
            has_best_candidate,
            result.best_candidate_stats,
            kLocalFittingConvergenceObjectiveRelativeTolerance,
            attempt,
            kLocalFittingObjectiveBacktrackingMaximumAttempts);
    }
    return result;
}

LocalRefitResult FitAtomWithJointOffsetFallback(
    const AtomObject & atom,
    const LocalGaussianResult & previous_result,
    const FittedGaussianSnapshot & offset_snapshot,
    const FitOptions & options)
{
    const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
    auto sample_entries{ UpdateSampleListWithFittedGaussian(atom, offset_snapshot) };
    const auto & offset_model{ offset_snapshot.at(&atom) };
    try
    {
        auto candidate_result{
            EstimateLocalGaussian(sample_entries, local_view.GetAlphaR(), options, offset_model)
        };
        if (CanBuildFiniteZeroOffsetSamples(sample_entries, candidate_result.mdpde.GetModel()))
        {
            return LocalRefitResult{ candidate_result, false, false };
        }
    }
    catch (const std::exception &)
    {
    }

    auto result{ previous_result };
    result.ols = WithModelOffset(result.ols, offset_model.GetOffset());
    result.mdpde = WithModelOffset(result.mdpde, offset_model.GetOffset());
    if (IsSuspiciousJointOffset(
            sample_entries,
            previous_result.mdpde.GetModel(),
            result.mdpde.GetModel()))
    {
        return LocalRefitResult{ previous_result, true, true };
    }
    return LocalRefitResult{ result, true, false };
}

std::unordered_map<const AtomObject *, std::size_t> BuildSelectedAtomIndexMap(
    const std::vector<AtomObject *> & atom_list)
{
    std::unordered_map<const AtomObject *, std::size_t> atom_index_map;
    atom_index_map.reserve(atom_list.size());
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        atom_index_map.emplace(atom_list.at(i), i);
    }
    return atom_index_map;
}

std::size_t ThawChangedActiveAtomNeighbors(
    const std::vector<AtomObject *> & atom_list,
    const std::unordered_map<const AtomObject *, std::size_t> & atom_index_map,
    const std::vector<algorithm::ParameterChange> & change_list,
    const std::vector<std::size_t> & active_index_list,
    algorithm::ConvergenceFreezeTracker & freeze_tracker,
    detail::LocalFittingThawHysteresisTracker & thaw_hysteresis_tracker)
{
    if (change_list.size() != atom_list.size())
    {
        throw std::invalid_argument("Local fitting dependency thaw input size is inconsistent.");
    }

    const double thaw_threshold{ std::sqrt(kLocalFittingParameterChangeTolerance) };
    std::size_t thaw_count{ 0 };
    for (const auto active_index : active_index_list)
    {
        if (active_index >= atom_list.size())
        {
            throw std::invalid_argument("Local fitting dependency thaw active index is out of range.");
        }
        const auto active_change{
            algorithm::GetMaximumParameterChange(change_list.at(active_index))
        };

        for (const auto * neighbor_atom : atom_list.at(active_index)->FindNeighborAtoms())
        {
            const auto neighbor_iter{ atom_index_map.find(neighbor_atom) };
            if (neighbor_iter == atom_index_map.end()) continue;

            const auto neighbor_index{ neighbor_iter->second };
            if (!freeze_tracker.IsFrozen(neighbor_index)) continue;
            if (!thaw_hysteresis_tracker.ShouldThaw(
                    neighbor_index,
                    active_change,
                    thaw_threshold))
            {
                continue;
            }
            if (freeze_tracker.Thaw(neighbor_index))
            {
                thaw_hysteresis_tracker.RecordDependencyThaw(neighbor_index);
                thaw_count++;
            }
        }
    }
    return thaw_count;
}

void ExpandSuspiciousOffsetClusters(
    const std::vector<std::vector<std::size_t>> & active_coupling_adjacency,
    std::vector<int> & suspicious_offset_flag_list,
    std::vector<int> & refit_fallback_flag_list)
{
    if (active_coupling_adjacency.empty()) return;
    if (active_coupling_adjacency.size() != suspicious_offset_flag_list.size() ||
        refit_fallback_flag_list.size() != suspicious_offset_flag_list.size())
    {
        throw std::invalid_argument("Suspicious offset cluster input sizes are inconsistent.");
    }

    std::vector<char> visited(suspicious_offset_flag_list.size(), 0);
    for (std::size_t seed_index = 0; seed_index < suspicious_offset_flag_list.size(); seed_index++)
    {
        if (suspicious_offset_flag_list.at(seed_index) == 0 || visited.at(seed_index) != 0)
        {
            continue;
        }

        std::vector<std::size_t> stack{ seed_index };
        visited.at(seed_index) = 1;
        while (!stack.empty())
        {
            const auto active_index{ stack.back() };
            stack.pop_back();
            suspicious_offset_flag_list.at(active_index) = 1;
            refit_fallback_flag_list.at(active_index) = 1;

            for (const auto neighbor_index : active_coupling_adjacency.at(active_index))
            {
                if (neighbor_index >= active_coupling_adjacency.size())
                {
                    throw std::invalid_argument(
                        "Suspicious offset cluster adjacency index is out of range.");
                }
                if (visited.at(neighbor_index) != 0) continue;
                visited.at(neighbor_index) = 1;
                stack.emplace_back(neighbor_index);
            }
        }
    }
}

void RollBackSuspiciousOffsetClusters(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<std::size_t> & active_index_list,
    const GaussianFittingState & previous_state,
    const std::vector<int> & suspicious_offset_flag_list,
    FittedGaussianSnapshot & current_snapshot,
    GaussianFittingState & iteration_state)
{
    if (active_index_list.size() != suspicious_offset_flag_list.size())
    {
        throw std::invalid_argument("Suspicious offset rollback input sizes are inconsistent.");
    }
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        if (suspicious_offset_flag_list.at(i) == 0) continue;

        const auto state_index{ active_index_list.at(i) };
        if (state_index >= atom_list.size() ||
            state_index >= previous_state.estimation_list.size() ||
            state_index >= previous_state.result_list.size() ||
            state_index >= iteration_state.estimation_list.size() ||
            state_index >= iteration_state.result_list.size())
        {
            throw std::invalid_argument("Suspicious offset rollback state index is out of range.");
        }
        const auto previous_model{
            GaussianModel3D::FromVector(previous_state.estimation_list.at(state_index))
        };
        current_snapshot.at(atom_list.at(state_index)) = previous_model;
        iteration_state.estimation_list.at(state_index) =
            previous_state.estimation_list.at(state_index);
        iteration_state.result_list.at(state_index) =
            previous_state.result_list.at(state_index);
    }
}

LocalFittingIterationResult RunLocalFittingIteration(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<std::size_t> & active_index_list,
    const GaussianFittingState & previous_state,
    const FitOptions & options,
    double ridge_ratio,
    const std::vector<double> & ridge_multiplier_list)
{
    const auto selected_atom_size{ atom_list.size() };
    if (previous_state.result_list.size() != selected_atom_size ||
        previous_state.estimation_list.size() != selected_atom_size ||
        ridge_multiplier_list.size() != selected_atom_size)
    {
        throw std::invalid_argument("Local fitting iteration input sizes are inconsistent.");
    }
    auto current_snapshot{
        BuildFittedGaussianSnapshot(atom_list, previous_state.estimation_list)
    };
    const auto joint_offset_result{
        EstimateJointOffsets(
            atom_list,
            active_index_list,
            current_snapshot,
            ridge_ratio,
            ridge_multiplier_list)
    };
    ApplyJointOffsetsToSnapshot(
        atom_list,
        active_index_list,
        joint_offset_result.offset,
        current_snapshot);
    auto iteration_state{ previous_state };
    std::vector<int> refit_fallback_flag_list(active_index_list.size(), 0);
    std::vector<int> suspicious_offset_flag_list(active_index_list.size(), 0);

    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto state_index{ active_index_list.at(i) };
        const auto * atom{ atom_list.at(state_index) };
        const auto sample_entries{
            AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false)
        };
        const auto previous_model{
            GaussianModel3D::FromVector(previous_state.estimation_list.at(state_index))
        };
        const auto & offset_model{ current_snapshot.at(atom) };
        if (!IsSuspiciousJointOffset(sample_entries, previous_model, offset_model))
        {
            continue;
        }
        suspicious_offset_flag_list.at(i) = 1;
    }
    ExpandSuspiciousOffsetClusters(
        joint_offset_result.active_coupling_adjacency,
        suspicious_offset_flag_list,
        refit_fallback_flag_list);
    RollBackSuspiciousOffsetClusters(
        atom_list,
        active_index_list,
        previous_state,
        suspicious_offset_flag_list,
        current_snapshot,
        iteration_state);

    for (size_t i = 0; i < active_index_list.size(); i++)
    {
        if (suspicious_offset_flag_list.at(i) != 0) continue;

        const auto state_index{ active_index_list.at(i) };
        const auto & atom{ *atom_list.at(state_index) };
        auto refit_result{
            FitAtomWithJointOffsetFallback(
                atom,
                previous_state.result_list.at(state_index),
                current_snapshot,
                options)
        };
        if (refit_result.used_fallback)
        {
            refit_fallback_flag_list.at(i) = 1;
        }
        if (refit_result.suspicious_offset_fallback)
        {
            suspicious_offset_flag_list.at(i) = 1;
        }
        auto result{ std::move(refit_result.result) };
        const auto fitted_model{ result.mdpde.GetModel() };
        iteration_state.estimation_list.at(state_index) = fitted_model.ToVector();
        iteration_state.result_list.at(state_index) = std::move(result);
    }
    ExpandSuspiciousOffsetClusters(
        joint_offset_result.active_coupling_adjacency,
        suspicious_offset_flag_list,
        refit_fallback_flag_list);
    RollBackSuspiciousOffsetClusters(
        atom_list,
        active_index_list,
        previous_state,
        suspicious_offset_flag_list,
        current_snapshot,
        iteration_state);

    LocalFittingIterationResult iteration_result;
    iteration_result.state = std::move(iteration_state);
    iteration_result.joint_offset_used_fallback = joint_offset_result.used_fallback;
    for (std::size_t i = 0; i < refit_fallback_flag_list.size(); i++)
    {
        if (refit_fallback_flag_list.at(i) == 0) continue;
        iteration_result.refit_fallback_state_index_list.emplace_back(active_index_list.at(i));
    }
    for (std::size_t i = 0; i < suspicious_offset_flag_list.size(); i++)
    {
        if (suspicious_offset_flag_list.at(i) == 0) continue;
        iteration_result.suspicious_offset_state_index_list.emplace_back(active_index_list.at(i));
    }
    return iteration_result;
}

void ApplyLocalFittingState(
    const GaussianFittingState & iteration_state,
    std::vector<AtomLocalPotentialEditor> & local_editor_list)
{
    if (local_editor_list.size() != iteration_state.result_list.size())
    {
        throw std::invalid_argument(
            "local_editor_list and local fitting state sizes are inconsistent.");
    }

    for (std::size_t i = 0; i < local_editor_list.size(); i++)
    {
        local_editor_list.at(i).SetGaussianResult(iteration_state.result_list.at(i));
    }
}

void LogLocalFittingFallbackSummary(const LocalFittingFallbackStats & fallback_stats)
{
    if (!fallback_stats.HasFallback()) return;

    std::ostringstream message;
    message << "Second-stage local fitting fallback summary: "
        << "joint offset fallback iterations = "
        << fallback_stats.joint_offset_fallback_iterations
        << ", refit fallback atom-events = "
        << fallback_stats.refit_fallback_atom_events
        << ", refit fallback distinct atoms = "
        << fallback_stats.GetDistinctRefitFallbackAtomCount()
        << ", suspicious offset atom-events = "
        << fallback_stats.suspicious_offset_atom_events
        << ", suspicious offset distinct atoms = "
        << fallback_stats.GetDistinctSuspiciousOffsetAtomCount()
        << ".";
    Logger::Log(LogLevel::Warning, message.str());
}

void InitializeLocalFittingSeedModels(ModelObject & model_object)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    const auto seed_model{ GaussianModel3D{ 0.0, 1.0, 0.0 } };
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list[i]) };
        auto result{ local_view.GetGaussianResult() };
        result.ols = GaussianModel3DWithUncertainty{
            seed_model,
            GaussianModel3DUncertainty{}
        };
        result.mdpde = GaussianModel3DWithUncertainty{
            seed_model,
            GaussianModel3DUncertainty{}
        };
        result.posterior.reset();
        result.is_outlier = false;
        result.statistical_distance = 0.0;
        result.fit_result.reset();
        local_editor_list[i].SetGaussianResult(result);
    }
}

[[maybe_unused]] LocalGaussianResult EstimateLocalGaussianWithOffset(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options,
    double offset_initial)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min,
        options.distance_max,
        "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");
    numeric_validation::RequireFinite(offset_initial, "offset_initial");

    const auto execution_options{ MakeExecutionOptions(options) };
    auto result{ EstimateLocalGaussian(sample_entries, alpha_r, options) };
    auto current_model{ result.mdpde.GetModel().WithOffset(offset_initial) };
    double best_error{ std::numeric_limits<double>::infinity() };
    auto best_result{ result };
    for (int iteration = 0; iteration < execution_options.max_iterations; iteration++)
    {
        const auto offset{ current_model.GetOffset() };
        result = EstimateLocalGaussian(sample_entries, alpha_r, options, current_model);
        const auto raw_offset{
            EstimateResidualOffsetParameter(
                sample_entries,
                *result.fit_result,
                offset)
        };
        const auto candidate_model{ result.mdpde.GetModel().WithOffset(raw_offset) };
        const auto error{ (candidate_model.ToVector() - current_model.ToVector()).norm() };
        if (error < best_error)
        {
            best_error = error;
            best_result = result;
        }
        if (error < execution_options.tolerance)
        {
            break;
        }

        if (iteration + 1 == execution_options.max_iterations)
        {
            result = best_result;
            if (!options.quiet_mode)
            {
                Logger::Log(LogLevel::Debug,
                    "Maximum iterations reached in local Gaussian estimation with offset; "
                    "refitting at best fixed-point candidate with error = " +
                    std::to_string(best_error) + ".");
            }
            break;
        }

        const auto damped_offset{
            offset + kOffsetDampingFactor * (raw_offset - offset)
        };
        current_model = result.mdpde.GetModel().WithOffset(damped_offset);
    }
    return result;
}

} // namespace

void ApplySpotMedianMDPDEOffsets(
    ModelObject & model_object,
    const std::vector<AtomObject *> & atom_list)
{
    std::unordered_map<Spot, std::vector<double>> offset_list_by_spot;
    offset_list_by_spot.reserve(atom_list.size());
    for (const auto * atom : atom_list)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        offset_list_by_spot[atom->GetSpot()].emplace_back(
            local_view.GetGaussianResult().mdpde.GetModel().GetOffset());
    }

    std::unordered_map<Spot, double> median_offset_by_spot;
    median_offset_by_spot.reserve(offset_list_by_spot.size());
    for (const auto & [spot, offset_list] : offset_list_by_spot)
    {
        median_offset_by_spot.emplace(spot, array_helper::ComputeMedian(offset_list));
    }

    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list.at(i)) };
        auto result{ local_view.GetGaussianResult() };
        result.mdpde = WithModelOffset(
            result.mdpde,
            median_offset_by_spot.at(atom_list.at(i)->GetSpot()));
        local_editor_list.at(i).SetGaussianResult(std::move(result));
    }
}

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
    const GaussianModel3D & offset_model)
{
    auto range_min{ options.distance_min };
    auto range_max{ options.distance_max };
    numeric_validation::RequireFiniteNonNegativeRange(range_min, range_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");
    numeric_validation::RequireFinite(offset_model.GetOffset(), "offset");

    auto execution_options{ MakeExecutionOptions(options) };
    const auto updated_sample_entries{
        BuildSamplesForZeroOffsetGaussianFit(sample_entries, offset_model)
    };
    auto dataset{
        rhbm_helper::BuildMemberDataset(updated_sample_entries, range_min, range_max)
    };
    const auto result{ rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset, execution_options) };
    return DecodeLocalGaussianResult(alpha_r, result, offset_model.GetOffset());
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
    std::vector<double> member_offset_list;
    member_offset_list.reserve(member_result_list.size());
    for (const auto & member_result : member_result_list)
    {
        member_offset_list.emplace_back(member_result.mdpde.GetModel().GetOffset());
    }
    const auto group_offset{ array_helper::ComputeMedian(member_offset_list) };
    auto result{ DecodeGroupGaussianResult(alpha_g, raw_result, group_offset) };
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
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    const auto selected_atom_size{ atom_list.size() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    std::atomic<size_t> atom_count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info,
            "Run 1st-stage local atom fitting for " +
            std::to_string(selected_atom_size) + " atoms.");
    }

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list[i]) };
        auto sample_entries{ local_view.GetSamplingEntries() };
        const auto offset_model{ local_view.GetGaussianResult().mdpde.GetModel() };
        auto result{
            EstimateLocalGaussian(sample_entries, local_view.GetAlphaR(), options, offset_model)
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

void RunSecondStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    const auto atom_size{ atom_list.size() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run 2nd-stage local atom fitting with iterations...");
    }

    GaussianFittingState previous_state{
        std::vector<LocalGaussianResult>(atom_size),
        std::vector<Eigen::VectorXd>(atom_size)
    };
    for (size_t i = 0; i < atom_size; i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list[i]) };
        previous_state.result_list[i] = local_view.GetGaussianResult();
        previous_state.estimation_list[i] = previous_state.result_list[i].mdpde.GetModel().ToVector();
    }

    std::optional<double> initial_objective_scale_sample;
    const auto initial_objective_samples{
        CollectLocalFittingObjectiveSamples(atom_list, previous_state)
    };
    if (initial_objective_samples.has_value())
    {
        initial_objective_scale_sample =
            CalculateLocalFittingResidualScaleSample(*initial_objective_samples);
    }
    LocalFittingObjectiveScaleTracker objective_scale_tracker{
        initial_objective_scale_sample
    };
    const auto objective_reference{ objective_scale_tracker.GetCommittedReference() };
    LocalFittingObjectiveStats previous_objective_stats;
    if (initial_objective_samples.has_value())
    {
        previous_objective_stats = CalculateLocalFittingObjectiveStats(
            *initial_objective_samples,
            objective_reference);
    }
    algorithm::FittingQualityCandidateStats previous_candidate_stats{
        previous_objective_stats.has_quality_objective,
        previous_objective_stats.quality_objective,
        algorithm::ParameterChangeStats{
            std::vector<double>(
                static_cast<std::size_t>(GaussianModel3D::ParameterSize()),
                0.0)
        }
    };
    GaussianFittingState best_state;
    algorithm::FittingQualityCandidateStats best_candidate_stats;
    bool has_best_candidate{ false };
    if (previous_objective_stats.has_quality_objective)
    {
        best_state = previous_state;
        best_candidate_stats = previous_candidate_stats;
        has_best_candidate = true;
    }
    algorithm::AdaptiveRelaxationController relaxation_controller{
        kAdaptiveRelaxationInitialBeta,
        kAdaptiveRelaxationMin,
        kAdaptiveRelaxationMax,
        kAdaptiveRelaxationGrowth,
        kAdaptiveRelaxationShrink,
        kAdaptiveRelaxationImprovementRatio,
        kAdaptiveRelaxationIncreaseStreak,
        kAdaptiveRelaxationShrinkStreak
    };
    algorithm::ConvergenceFreezeTracker freeze_tracker{
        atom_size,
        kLocalFittingParameterChangeTolerance,
        kLocalFittingFreezeChangeRatio,
        kLocalFittingFreezeStableIterations
    };
    LocalFittingFallbackStats fallback_stats{ atom_size };
    const auto atom_index_map{ BuildSelectedAtomIndexMap(atom_list) };
    detail::LocalFittingThawHysteresisTracker thaw_hysteresis_tracker{
        atom_size,
        kLocalFittingDependencyThawHysteresisGrowth,
        kLocalFittingDependencyThawHysteresisMax,
        kLocalFittingDependencyThawHysteresisFrozenDecay
    };
    double ridge_ratio{ kJointOffsetRidgeRatio };
    std::vector<double> joint_offset_ridge_multiplier_list(atom_size, 1.0);
    for (size_t iter = 0; iter < kLocalFittingMaximumIterations; iter++)
    {
        const auto active_index_list{ freeze_tracker.BuildActiveIndexList() };
        if (active_index_list.empty())
        {
            ApplyLocalFittingState(previous_state, local_editor_list);
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                Logger::Log(LogLevel::Info,
                    "Converged after " + std::to_string(iter) +
                    " iterations because all local atoms are frozen.");
            }
            break;
        }

        const auto iteration_ridge_ratio{ ridge_ratio };
        auto iteration_result{
            RunLocalFittingIteration(
                atom_list,
                active_index_list,
                previous_state,
                options,
                iteration_ridge_ratio,
                joint_offset_ridge_multiplier_list)
        };
        fallback_stats.Accumulate(iteration_result);
        const auto suspicious_offset_state_index_list{
            iteration_result.suspicious_offset_state_index_list
        };
        const auto has_suspicious_offset_fallback{
            !suspicious_offset_state_index_list.empty()
        };
        std::vector<bool> suspicious_offset_atom_seen(atom_size, false);
        for (const auto state_index : suspicious_offset_state_index_list)
        {
            if (state_index >= suspicious_offset_atom_seen.size())
            {
                throw std::invalid_argument(
                    "Local fitting suspicious offset atom index is out of range.");
            }
            suspicious_offset_atom_seen.at(state_index) = true;
            joint_offset_ridge_multiplier_list.at(state_index) =
                kSuspiciousJointOffsetRidgeMultiplier;
        }
        const auto raw_state{ std::move(iteration_result.state) };
        GaussianFittingState current_state;
        std::vector<algorithm::ParameterChange> change_list;
        algorithm::ParameterChangeStats change_stats;
        algorithm::ParameterChangeStats normalized_change_stats;
        algorithm::FittingQualityCandidateStats current_candidate_stats;
        double beta{ relaxation_controller.GetBeta() };
        bool has_current_candidate{ false };
        bool has_backtracking_rejection{ false };
        double current_objective_scale_sample{ std::numeric_limits<double>::infinity() };
        for (int attempt = 0; attempt < kLocalFittingObjectiveBacktrackingMaximumAttempts; attempt++)
        {
            beta = relaxation_controller.GetBeta();
            auto attempt_result{
                EvaluateLocalFittingAttempt(
                    atom_list,
                    active_index_list,
                    raw_state,
                    previous_state,
                    best_state,
                    previous_candidate_stats,
                    best_candidate_stats,
                    has_best_candidate,
                    objective_scale_tracker,
                    beta,
                    attempt)
            };

            if (attempt_result.backtracking_decision.accepted)
            {
                current_state = std::move(attempt_result.state);
                change_list = std::move(attempt_result.change_list);
                change_stats = std::move(attempt_result.change_stats);
                normalized_change_stats = std::move(attempt_result.normalized_change_stats);
                current_candidate_stats = std::move(attempt_result.candidate_stats);
                previous_candidate_stats = std::move(attempt_result.previous_candidate_stats);
                if (has_best_candidate)
                {
                    best_candidate_stats = std::move(attempt_result.best_candidate_stats);
                }
                current_objective_scale_sample = attempt_result.objective_scale_sample;
                has_current_candidate = true;
                break;
            }
            has_backtracking_rejection = true;
            ridge_ratio = IncreaseLocalFittingRidgeRatio(ridge_ratio);
            if (attempt_result.backtracking_decision.should_shrink_beta &&
                !relaxation_controller.IsAtMinimum())
            {
                relaxation_controller.Shrink();
                continue;
            }
            break;
        }

        if (!has_current_candidate)
        {
            if (!relaxation_controller.IsAtMinimum() && iter + 1 < kLocalFittingMaximumIterations)
            {
                const auto next_beta{ relaxation_controller.Shrink() };
                if (!options.quiet_mode)
                {
                    std::ostringstream progress_message;
                    progress_message << "Local fitting iteration " << iter + 1 << '/'
                        << kLocalFittingMaximumIterations
                        << " rejected by objective backtracking; retry with beta = "
                        << std::fixed << std::setprecision(5)
                        << next_beta
                        << ", next ridge ratio = "
                        << ridge_ratio;
                    Logger::ProgressLine(progress_message.str());
                }
                continue;
            }

            ApplyLocalFittingState(
                has_best_candidate ? best_state : previous_state,
                local_editor_list);
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                std::ostringstream warning_message;
                warning_message
                    << "Stopped local fitting because objective backtracking rejected the candidate at "
                    << (relaxation_controller.IsAtMinimum() ? "minimum beta" : "the maximum iteration limit")
                    << "; applying "
                    << (has_best_candidate ? "best fixed-point candidate." : "previous state.");
                Logger::Log(LogLevel::Warning, warning_message.str());
            }
            break;
        }

        if (current_candidate_stats.has_quality_objective)
        {
            objective_scale_tracker.CommitScaleSample(current_objective_scale_sample);
        }
        for (const auto active_index : active_index_list)
        {
            if (!suspicious_offset_atom_seen.at(active_index))
            {
                joint_offset_ridge_multiplier_list.at(active_index) = 1.0;
            }
        }
        if (!has_backtracking_rejection)
        {
            ridge_ratio = DecreaseLocalFittingRidgeRatio(ridge_ratio);
        }
        relaxation_controller.Update(normalized_change_stats);
        if (!has_best_candidate ||
            algorithm::IsBetterFittingQualityCandidate(
                current_candidate_stats,
                best_candidate_stats,
                kLocalFittingObjectiveTieRelativeTolerance))
        {
            best_state = current_state;
            best_candidate_stats = current_candidate_stats;
            has_best_candidate = true;
        }
        freeze_tracker.Update(change_list, active_index_list);
        auto thaw_count{
            ThawChangedActiveAtomNeighbors(
                atom_list,
                atom_index_map,
                change_list,
                active_index_list,
                freeze_tracker,
                thaw_hysteresis_tracker)
        };
        for (const auto state_index : suspicious_offset_state_index_list)
        {
            if (freeze_tracker.Thaw(state_index))
            {
                thaw_count++;
            }
        }
        for (std::size_t state_index = 0; state_index < atom_size; state_index++)
        {
            if (freeze_tracker.IsFrozen(state_index))
            {
                thaw_hysteresis_tracker.DecayFrozen(state_index);
            }
        }

        if (!options.quiet_mode)
        {
            std::ostringstream progress_message;
            progress_message << "Iter. " << iter + 1 << '/' << kLocalFittingMaximumIterations
                << std::fixed << std::setprecision(4)
                << ", d_amplitude = "<< change_stats.percentile_list.at(GaussianModel3D::AmplitudeIndex())
                << ", d_width = "<< change_stats.percentile_list.at(GaussianModel3D::WidthIndex())
                << ", d_offset = "<< change_stats.percentile_list.at(GaussianModel3D::OffsetIndex())
                << ", objective = "<< current_candidate_stats.quality_objective
                << ", beta = "<< beta
                //<< ", ridge ratio = "<< iteration_ridge_ratio
                //<< ", next ridge ratio = "<< ridge_ratio
                << ", active/frozen/thawed atoms = "<< freeze_tracker.GetActiveCount()
                << "/" << freeze_tracker.GetFrozenCount() << "/" << thaw_count;
            Logger::ProgressLine(progress_message.str());
        }

        if (freeze_tracker.GetActiveCount() == 0)
        {
            ApplyLocalFittingState(current_state, local_editor_list);
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                Logger::Log(LogLevel::Info,
                    "Converged after " + std::to_string(iter + 1) +
                    " iterations because all local atoms are frozen.");
            }
            break;
        }

        const auto converged{
            !has_suspicious_offset_fallback &&
            (!objective_scale_tracker.HasReference() || objective_scale_tracker.IsLocked()) &&
            IsLocalFittingNormalizedParameterChangeConverged(normalized_change_stats)
        };
        if (converged)
        {
            ApplyLocalFittingState(current_state, local_editor_list);
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                Logger::Log(LogLevel::Info,
                    "Converged after " + std::to_string(iter + 1) +
                    " iterations with normalized percentile amplitude change = " +
                    std::to_string(
                        normalized_change_stats.percentile_list.at(
                            GaussianModel3D::AmplitudeIndex())) +
                    ", normalized percentile width change = " +
                    std::to_string(
                        normalized_change_stats.percentile_list.at(
                            GaussianModel3D::WidthIndex())) +
                    ", and normalized percentile offset change = " +
                    std::to_string(
                        normalized_change_stats.percentile_list.at(
                            GaussianModel3D::OffsetIndex())) +
                    ", objective = " +
                    std::to_string(current_candidate_stats.quality_objective) + ".");
            }
            break;
        }

        if (iter + 1 == kLocalFittingMaximumIterations)
        {
            ApplyLocalFittingState(best_state, local_editor_list);
            if (!options.quiet_mode)
            {
                Logger::FinishProgressLine();
                Logger::Log(LogLevel::Warning,
                    "Reached maximum iteration size; refitting at best fixed-point candidate "
                    "with normalized percentile amplitude change = " +
                    std::to_string(
                        best_candidate_stats.parameter_change_stats.percentile_list.at(
                            GaussianModel3D::AmplitudeIndex())) +
                    ", normalized percentile width change = " +
                    std::to_string(
                        best_candidate_stats.parameter_change_stats.percentile_list.at(
                            GaussianModel3D::WidthIndex())) +
                    ", and normalized percentile offset change = " +
                    std::to_string(
                        best_candidate_stats.parameter_change_stats.percentile_list.at(
                            GaussianModel3D::OffsetIndex())) +
                    ", objective = " +
                    std::to_string(best_candidate_stats.quality_objective));
            }
        }
        previous_state = std::move(current_state);
        previous_candidate_stats = current_candidate_stats;
    }
    if (!options.quiet_mode)
    {
        LogLocalFittingFallbackSummary(fallback_stats);
    }
}

void RunThirdStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    const auto selected_atom_size{ atom_list.size() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    const auto median_model_by_spot{ BuildSpotMedianMDPDEModelMap(atom_list) };
    std::atomic<size_t> atom_count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info,
            "Run 3rd-stage local atom fitting for " +
            std::to_string(selected_atom_size) + " atoms.");
    }

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom_list[i]) };
        auto sample_entries{
            UpdateSampleListWithSpotMedianGaussian(*atom_list[i], median_model_by_spot)
        };
        auto result{
            EstimateLocalGaussian(sample_entries, local_view.GetAlphaR(), options,
                median_model_by_spot.at(atom_list[i]->GetSpot()))
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

void RunLocalPotentialFitting(ModelObject & model_object, const FitOptions & options)
{
    InitializeLocalFittingSeedModels(model_object);
    RunFirstStageLocalFitting(model_object, options);
    RunSecondStageLocalFitting(model_object, options);
    RunThirdStageLocalFitting(model_object, options);
}

void RunGroupPotentialFitting(ModelObject & model_object, const FitOptions & options)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto & selected_atom_list{ model_object.GetSelectedAtoms() };
    for (auto * atom : selected_atom_list)
    {
        analysis.EnsureAtomLocalPotential(*atom);
    }
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run atom group fitting.");
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
        const auto alpha_g{ analysis_view.GetAtomAlphaG(group_key) };
        std::vector<LocalPotentialSampleList> sample_entries_list;
        std::vector<LocalGaussianResult> member_result_list;
        sample_entries_list.reserve(atom_list.size());
        member_result_list.reserve(atom_list.size());
        for (const auto & atom : atom_list)
        {
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            sample_entries_list.emplace_back(UpdateSampleListWithFittedGaussian(*atom));
            member_result_list.emplace_back(local_view.GetGaussianResult());
        }
        const auto result{
            EstimateGroupGaussian(sample_entries_list, member_result_list, alpha_g, options)
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
    if (!options.quiet_mode)
    {
        LogGroupPriorSpotSummary(model_object);
    }
}

} // namespace rhbm_gem::core
