#include <cstddef>
#include <rhbm_gem/core/GaussianEstimator.hpp>

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
#include <array>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
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
constexpr std::array<double, 3> kLocalFittingNormalizedChangeToleranceList{
    1.0e-5,
    1.0e-5,
    1.0e-5
};
constexpr std::array<double, 3> kLocalFittingNormalizedChangeScaleFloorList{
    1.0,
    1.0,
    1.0
};
constexpr std::size_t kMinimumAlphaRTrainingSampleCount{ 10 };
constexpr std::size_t kMinimumAlphaGTrainingMemberCount{ 10 };
constexpr double kResidualOffsetRangeMin{ 1.0 };
constexpr double kResidualOffsetRangeMax{ 2.0 };
constexpr double kOffsetDampingFactor{ 0.5 };
constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr std::size_t kLocalFittingMaximumIterations{ 1000 };
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
constexpr double kJointOffsetIrlsScaleFloor{ 1.0e-2 };
constexpr double kJointOffsetIrlsNormalizedChangeTolerance{ 1.0e-6 };
constexpr double kJointOffsetIrlsObjectiveRelativeTolerance{ 1.0e-10 };
constexpr double kAdaptiveRelaxationMin{ 0.05 };
constexpr double kAdaptiveRelaxationMax{ 1.0 };
constexpr double kAdaptiveRelaxationGrowth{ 1.2 };
constexpr double kAdaptiveRelaxationShrink{ 0.5 };
constexpr double kAdaptiveRelaxationImprovementRatio{ 0.01 };
constexpr int kAdaptiveRelaxationIncreaseStreak{ 2 };
constexpr double kLocalFittingFreezeChangeRatio{ 0.1 };
constexpr int kLocalFittingFreezeStableIterations{ 3 };
constexpr double kLocalFittingObjectiveTieRelativeTolerance{ 1.0e-8 };
constexpr double kLocalFittingConvergenceObjectiveRelativeTolerance{ 1.0e-3 };
constexpr int kLocalFittingObjectiveBacktrackingMaximumAttempts{ 3 };
constexpr std::size_t kAmplitudeChangeIndex{ 0 };
constexpr std::size_t kWidthChangeIndex{ 1 };
constexpr std::size_t kOffsetChangeIndex{ 2 };

using GaussianFittingState = algorithm::IterationState<LocalGaussianResult, Eigen::VectorXd>;

enum class JointOffsetFallbackReason
{
    None,
    BuildSystemFailed,
    EmptySystem,
    InitialSolveFailed,
    RobustSolveFailed
};

struct JointOffsetSolveResult
{
    Eigen::VectorXd offset{};
    bool used_fallback{ false };
    JointOffsetFallbackReason fallback_reason{ JointOffsetFallbackReason::None };
};

struct LocalRefitResult
{
    LocalGaussianResult result{};
    bool used_fallback{ false };
};

struct LocalFittingIterationDiagnostics
{
    bool joint_offset_used_fallback{ false };
    std::vector<std::size_t> refit_fallback_state_index_list{};
};

struct LocalFittingIterationResult
{
    GaussianFittingState state{};
    LocalFittingIterationDiagnostics diagnostics{};
};

struct LocalFittingFallbackStats
{
    std::size_t joint_offset_fallback_iterations{ 0 };
    std::size_t refit_fallback_atom_events{ 0 };
    std::vector<bool> refit_fallback_atom_seen{};

    explicit LocalFittingFallbackStats(std::size_t atom_size)
        : refit_fallback_atom_seen(atom_size, false)
    {
    }

    void Accumulate(const LocalFittingIterationDiagnostics & diagnostics)
    {
        if (diagnostics.joint_offset_used_fallback)
        {
            joint_offset_fallback_iterations++;
        }
        refit_fallback_atom_events += diagnostics.refit_fallback_state_index_list.size();
        for (const auto state_index : diagnostics.refit_fallback_state_index_list)
        {
            if (state_index >= refit_fallback_atom_seen.size())
            {
                throw std::invalid_argument("Local fitting fallback atom index is out of range.");
            }
            refit_fallback_atom_seen.at(state_index) = true;
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

    bool HasFallback() const
    {
        return joint_offset_fallback_iterations > 0 || refit_fallback_atom_events > 0;
    }
};

struct LocalFittingObjectiveStats
{
    bool has_quality_objective{ false };
    double quality_objective{ std::numeric_limits<double>::infinity() };
    std::size_t sample_count{ 0 };
};

struct LocalFittingObjectiveReference
{
    bool has_reference{ false };
    double residual_scale{ std::numeric_limits<double>::infinity() };
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

double ComputeOffsetRegularizationPriorScale(double amplitude)
{
    return std::max(
        std::abs(amplitude) * kOffsetRegularizationAmplitudeRatio,
        kOffsetRegularizationPriorScaleMin);
}

double EstimateResidualOffsetParameter(
    const LocalPotentialSampleList & sample_entries,
    const RHBMBetaEstimateResult & fit_result,
    double current_offset)
{
    const auto signal_model{ linearization_service::DecodeParameterVector(fit_result.beta_mdpde) };
    const auto width{ signal_model.GetWidth() };
    if (!std::isfinite(width) || width <= 0.0) return current_offset;

    std::vector<algorithm::LinearRegressionSample> residual_sample_list;
    residual_sample_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        if (distance < kResidualOffsetRangeMin) continue;
        if (distance > kResidualOffsetRangeMax) continue;

        const auto basis{ signal_model.OffsetBasisAtDistance(distance) };
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
        residual_sample_list.emplace_back(algorithm::LinearRegressionSample{ basis, residual });
    }
    double candidate_offset{ current_offset };
    algorithm::RobustSlopeOptions slope_options;
    slope_options.maximum_iterations = kHuberSlopeMaximumIterations;
    slope_options.tolerance = kHuberSlopeTolerance;
    slope_options.scale_multiplier = kHuberScaleMultiplier;
    slope_options.scale_min = kHuberScaleMin;
    slope_options.cutoff_multiplier = kHuberCutoffMultiplier;
    slope_options.regularization_prior_scale =
        ComputeOffsetRegularizationPriorScale(signal_model.GetAmplitude());
    if (!algorithm::RobustSlopeEstimator::EstimateHuberSlopeThroughOrigin(
            residual_sample_list,
            slope_options,
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

algorithm::WeightedRidgeSystem BuildJointOffsetSystem(
    const std::vector<AtomObject *> & atom_list,
    const FittedGaussianSnapshot & snapshot)
{
    std::unordered_map<const AtomObject *, Eigen::Index> atom_column_map;
    atom_column_map.reserve(atom_list.size());
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        atom_column_map.emplace(atom_list.at(i), static_cast<Eigen::Index>(i));
    }

    std::vector<algorithm::SparseRegressionRow> row_list;
    for (const auto * atom : atom_list)
    {
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
    const auto column_count{ static_cast<Eigen::Index>(atom_list.size()) };
    std::vector<Eigen::Triplet<double>> triplet_list;
    Eigen::VectorXd response{ Eigen::VectorXd::Zero(row_count) };
    Eigen::VectorXd column_square_sum{ Eigen::VectorXd::Zero(column_count) };
    for (Eigen::Index row_index = 0; row_index < row_count; row_index++)
    {
        const auto & row{ row_list.at(static_cast<std::size_t>(row_index)) };
        response(row_index) = row.response;
        for (const auto & [column_index, basis] : row.basis_entries)
        {
            triplet_list.emplace_back(row_index, column_index, basis);
            column_square_sum(column_index) += basis * basis;
        }
    }

    algorithm::WeightedRidgeSystem system;
    system.design_matrix.resize(row_count, column_count);
    system.design_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
    system.response = std::move(response);
    system.previous_parameter = Eigen::VectorXd::Zero(column_count);
    system.ridge_diagonal = Eigen::VectorXd::Zero(column_count);
    for (Eigen::Index column_index = 0; column_index < column_count; column_index++)
    {
        const auto & model{ snapshot.at(atom_list.at(static_cast<std::size_t>(column_index))) };
        system.previous_parameter(column_index) = model.GetOffset();
        const auto square_sum{ column_square_sum(column_index) };
        system.ridge_diagonal(column_index) =
            square_sum > std::numeric_limits<double>::epsilon() ? kJointOffsetRidgeRatio * square_sum : 1.0;
    }
    return system;
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
    const FittedGaussianSnapshot & snapshot)
{
    Eigen::VectorXd previous_offset{
        Eigen::VectorXd::Zero(static_cast<Eigen::Index>(atom_list.size()))
    };
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        previous_offset(static_cast<Eigen::Index>(i)) = snapshot.at(atom_list.at(i)).GetOffset();
    }
    algorithm::WeightedRidgeSystem system;
    try
    {
        system = BuildJointOffsetSystem(atom_list, snapshot);
    }
    catch (const std::exception &)
    {
        return JointOffsetSolveResult{
            previous_offset,
            true,
            JointOffsetFallbackReason::BuildSystemFailed
        };
    }
    if (system.response.size() == 0 || system.previous_parameter.size() == 0)
    {
        return JointOffsetSolveResult{
            previous_offset,
            true,
            JointOffsetFallbackReason::EmptySystem
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
            JointOffsetFallbackReason::InitialSolveFailed
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
                JointOffsetFallbackReason::RobustSolveFailed
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
        JointOffsetFallbackReason::None
    };
}

void ApplyJointOffsetsToSnapshot(
    const std::vector<AtomObject *> & atom_list,
    const Eigen::VectorXd & offset,
    FittedGaussianSnapshot & snapshot)
{
    if (atom_list.size() != static_cast<std::size_t>(offset.size()))
    {
        throw std::invalid_argument("Joint offset result size is inconsistent.");
    }
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        snapshot.at(atom_list.at(i)) =
            snapshot.at(atom_list.at(i)).WithOffset(offset(static_cast<Eigen::Index>(i)));
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

double CalculateHuberLoss(double residual, double cutoff)
{
    const auto absolute_residual{ std::abs(residual) };
    if (absolute_residual <= cutoff)
    {
        return 0.5 * residual * residual;
    }
    return cutoff * (absolute_residual - 0.5 * cutoff);
}

std::optional<std::vector<double>> CollectLocalFittingResiduals(
    const std::vector<AtomObject *> & atom_list,
    const GaussianFittingState & fitting_state)
{
    const auto snapshot{
        BuildFittedGaussianSnapshot(atom_list, fitting_state.estimation_list)
    };
    std::vector<double> residual_list;
    for (const auto * atom : atom_list)
    {
        const auto model_iter{ snapshot.find(atom) };
        if (model_iter == snapshot.end())
        {
            throw std::invalid_argument("Local fitting objective snapshot is missing an atom.");
        }
        const auto sample_entries{ UpdateSampleListWithFittedGaussian(*atom, snapshot) };
        residual_list.reserve(residual_list.size() + sample_entries.size());
        const auto & target_model{ model_iter->second };
        for (const auto & sample : sample_entries)
        {
            const auto distance{ static_cast<double>(sample.point.distance) };
            const auto expected_response{ target_model.ResponseAtDistance(distance) };
            const auto residual{ static_cast<double>(sample.response) - expected_response };
            if (!std::isfinite(residual))
            {
                return std::nullopt;
            }
            residual_list.emplace_back(residual);
        }
    }
    return residual_list;
}

LocalFittingObjectiveReference BuildLocalFittingObjectiveReference(
    const std::vector<AtomObject *> & atom_list,
    const GaussianFittingState & fitting_state)
{
    LocalFittingObjectiveReference reference;
    const auto residual_list{ CollectLocalFittingResiduals(atom_list, fitting_state) };
    if (!residual_list.has_value() || residual_list->empty())
    {
        return reference;
    }

    const auto median_residual{ array_helper::ComputeMedian(*residual_list) };
    std::vector<double> deviation_list;
    deviation_list.reserve(residual_list->size());
    for (const auto residual : *residual_list)
    {
        deviation_list.emplace_back(std::abs(residual - median_residual));
    }
    const auto residual_scale{ std::max(
        kHuberScaleMultiplier * array_helper::ComputeMedian(deviation_list),
        kHuberScaleMin)
    };
    if (!std::isfinite(residual_scale))
    {
        return reference;
    }

    reference.has_reference = true;
    reference.residual_scale = residual_scale;
    return reference;
}

LocalFittingObjectiveStats CalculateLocalFittingObjectiveStats(
    const std::vector<AtomObject *> & atom_list,
    const GaussianFittingState & fitting_state,
    const LocalFittingObjectiveReference & objective_reference)
{
    LocalFittingObjectiveStats stats;
    if (!objective_reference.has_reference)
    {
        return stats;
    }

    const auto residual_list{ CollectLocalFittingResiduals(atom_list, fitting_state) };
    if (!residual_list.has_value() || residual_list->empty())
    {
        return stats;
    }

    double loss_sum{ 0.0 };
    for (const auto residual : *residual_list)
    {
        const auto normalized_residual{ residual / objective_reference.residual_scale };
        loss_sum += CalculateHuberLoss(normalized_residual, kHuberCutoffMultiplier);
    }
    const auto quality_objective{
        loss_sum / static_cast<double>(residual_list->size())
    };
    if (!std::isfinite(quality_objective))
    {
        return stats;
    }

    stats.has_quality_objective = true;
    stats.quality_objective = quality_objective;
    stats.sample_count = residual_list->size();
    return stats;
}

algorithm::FittingQualityCandidateStats BuildLocalFittingCandidateStats(
    const LocalFittingObjectiveStats & objective_stats,
    const algorithm::ParameterChangeStats & change_stats)
{
    return algorithm::FittingQualityCandidateStats{
        objective_stats.has_quality_objective,
        objective_stats.quality_objective,
        change_stats
    };
}

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
                kLocalFittingNormalizedChangeScaleFloorList.at(kAmplitudeChangeIndex)),
            algorithm::CalculateNormalizedChange(
                current_estimation_list[i](GaussianModel3D::WidthIndex()),
                previous_estimation_list[i](GaussianModel3D::WidthIndex()),
                kLocalFittingNormalizedChangeScaleFloorList.at(kWidthChangeIndex)),
            algorithm::CalculateNormalizedChange(
                current_estimation_list[i](GaussianModel3D::OffsetIndex()),
                previous_estimation_list[i](GaussianModel3D::OffsetIndex()),
                kLocalFittingNormalizedChangeScaleFloorList.at(kOffsetChangeIndex))
        };
    }
    return change_list;
}

algorithm::ParameterChangeStats SummarizeLocalFittingParameterChangeStats(
    const std::vector<algorithm::ParameterChange> & change_list,
    const std::vector<std::size_t> & index_list)
{
    return algorithm::SummarizeParameterChangeStats(
        change_list,
        index_list,
        kLocalFittingChangePercentile);
}

double GetLocalFittingParameterChangePercentile(
    const algorithm::ParameterChangeStats & stats,
    std::size_t index)
{
    return stats.percentile_list.at(index);
}

bool IsLocalFittingNormalizedParameterChangeConverged(
    const algorithm::ParameterChangeStats & stats)
{
    if (stats.percentile_list.size() != kLocalFittingNormalizedChangeToleranceList.size())
    {
        throw std::invalid_argument(
            "Local fitting normalized parameter change stats size is inconsistent.");
    }
    for (std::size_t i = 0; i < stats.percentile_list.size(); i++)
    {
        if (stats.percentile_list.at(i) >= kLocalFittingNormalizedChangeToleranceList.at(i))
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
            EstimateLocalGaussianWithOffsetModel(
                sample_entries,
                local_view.GetAlphaR(),
                options,
                offset_model)
        };
        if (CanBuildFiniteZeroOffsetSamples(sample_entries, candidate_result.mdpde.GetModel()))
        {
            return LocalRefitResult{
                candidate_result,
                false
            };
        }
    }
    catch (const std::exception &)
    {
    }

    auto result{ previous_result };
    result.ols = WithModelOffset(result.ols, offset_model.GetOffset());
    result.mdpde = WithModelOffset(result.mdpde, offset_model.GetOffset());
    return LocalRefitResult{
        result,
        true
    };
}

std::vector<AtomObject *> BuildActiveAtomList(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<std::size_t> & active_index_list)
{
    std::vector<AtomObject *> active_atom_list;
    active_atom_list.reserve(active_index_list.size());
    for (const auto index : active_index_list)
    {
        if (index >= atom_list.size())
        {
            throw std::invalid_argument("Local fitting active atom index is out of range.");
        }
        active_atom_list.emplace_back(atom_list.at(index));
    }
    return active_atom_list;
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
    algorithm::ConvergenceFreezeTracker & freeze_tracker)
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
        if (algorithm::GetMaximumParameterChange(change_list.at(active_index)) < thaw_threshold)
        {
            continue;
        }

        for (const auto * neighbor_atom : atom_list.at(active_index)->FindNeighborAtoms())
        {
            const auto neighbor_iter{ atom_index_map.find(neighbor_atom) };
            if (neighbor_iter == atom_index_map.end()) continue;

            const auto neighbor_index{ neighbor_iter->second };
            if (!freeze_tracker.IsFrozen(neighbor_index)) continue;
            if (freeze_tracker.Thaw(neighbor_index))
            {
                thaw_count++;
            }
        }
    }
    return thaw_count;
}

LocalFittingIterationResult RunLocalFittingIteration(
    const std::vector<AtomObject *> & atom_list,
    const std::vector<std::size_t> & active_index_list,
    const GaussianFittingState & previous_state,
    const FitOptions & options)
{
    const auto selected_atom_size{ atom_list.size() };
    if (previous_state.result_list.size() != selected_atom_size ||
        previous_state.estimation_list.size() != selected_atom_size)
    {
        throw std::invalid_argument("Local fitting iteration input sizes are inconsistent.");
    }
    auto current_snapshot{
        BuildFittedGaussianSnapshot(atom_list, previous_state.estimation_list)
    };
    const auto active_atom_list{ BuildActiveAtomList(atom_list, active_index_list) };
    const auto joint_offset_result{
        EstimateJointOffsets(active_atom_list, current_snapshot)
    };
    ApplyJointOffsetsToSnapshot(active_atom_list, joint_offset_result.offset, current_snapshot);
    auto iteration_state{ previous_state };
    std::vector<int> refit_fallback_flag_list(active_index_list.size(), 0);

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < active_index_list.size(); i++)
    {
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
        auto result{ std::move(refit_result.result) };
        const auto fitted_model{ result.mdpde.GetModel() };
        iteration_state.estimation_list.at(state_index) = fitted_model.ToVector();
        iteration_state.result_list.at(state_index) = std::move(result);
    }

    LocalFittingIterationDiagnostics diagnostics;
    diagnostics.joint_offset_used_fallback = joint_offset_result.used_fallback;
    for (std::size_t i = 0; i < refit_fallback_flag_list.size(); i++)
    {
        if (refit_fallback_flag_list.at(i) == 0) continue;
        diagnostics.refit_fallback_state_index_list.emplace_back(active_index_list.at(i));
    }
    return LocalFittingIterationResult{
        std::move(iteration_state),
        std::move(diagnostics)
    };
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
        << ".";
    Logger::Log(LogLevel::Warning, message.str());
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
    double offset)
{
    numeric_validation::RequireFinite(offset, "offset");
    const auto zero_offset_result{
        EstimateLocalGaussianWithOffsetModel(
            sample_entries,
            alpha_r,
            options,
            GaussianModel3D{ 0.0, 1.0, 0.0 })
    };
    if (offset == 0.0)
    {
        return zero_offset_result;
    }
    const auto offset_model{ zero_offset_result.mdpde.GetModel().WithOffset(offset) };
    return EstimateLocalGaussianWithOffsetModel(sample_entries, alpha_r, options, offset_model);
}

LocalGaussianResult EstimateLocalGaussianWithOffset(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const FitOptions & options,
    double offset_initial)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");
    numeric_validation::RequireFiniteNonNegative(alpha_r, "alpha_r");
    numeric_validation::RequireFinite(offset_initial, "offset_initial");

    auto execution_options{ MakeExecutionOptions(options) };
    auto result{
        EstimateLocalGaussianWithOffsetModel(
            sample_entries,
            alpha_r,
            options,
            GaussianModel3D{ 0.0, 1.0, 0.0 })
    };
    auto current_model{ result.mdpde.GetModel().WithOffset(offset_initial) };
    double best_error{ std::numeric_limits<double>::infinity() };
    auto best_result{ result };
    auto max_iterations{ execution_options.max_iterations };
    auto tolerance{ execution_options.tolerance };
    for (int t = 0; t < max_iterations; t++)
    {
        const auto offset{ current_model.GetOffset() };
        result = EstimateLocalGaussianWithOffsetModel(sample_entries, alpha_r, options, current_model);
        const auto raw_offset{
            EstimateResidualOffsetParameter(sample_entries, *result.fit_result, offset)
        };
        const auto candidate_model{ result.mdpde.GetModel().WithOffset(raw_offset) };
        const auto error{ (candidate_model.ToVector() - current_model.ToVector()).norm() };
        if (error < best_error)
        {
            best_error = error;
            best_result = result;
        }
        if (error < tolerance)
        {
            break;
        }

        if (t + 1 == max_iterations)
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

        const auto damped_offset{ offset + kOffsetDampingFactor * (raw_offset - offset) };
        current_model = result.mdpde.GetModel().WithOffset(damped_offset);
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
        auto result{
            EstimateLocalGaussianWithOffset(
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
    auto analysis{ model_object.EditAnalysis() };
    const auto atom_size{ atom_list.size() };
    std::vector<AtomLocalPotentialEditor> local_editor_list;
    local_editor_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        local_editor_list.emplace_back(analysis.EnsureAtomLocalPotential(*atom));
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

    const auto objective_reference{
        BuildLocalFittingObjectiveReference(atom_list, previous_state)
    };
    const auto previous_objective_stats{
        CalculateLocalFittingObjectiveStats(
            atom_list,
            previous_state,
            objective_reference)
    };
    algorithm::FittingQualityCandidateStats previous_candidate_stats{
        previous_objective_stats.has_quality_objective,
        previous_objective_stats.quality_objective,
        algorithm::ParameterChangeStats{ std::vector<double>{
            0.0,
            0.0,
            0.0
        } }
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
        options.relaxation_factor,
        kAdaptiveRelaxationMin,
        kAdaptiveRelaxationMax,
        kAdaptiveRelaxationGrowth,
        kAdaptiveRelaxationShrink,
        kAdaptiveRelaxationImprovementRatio,
        kAdaptiveRelaxationIncreaseStreak
    };
    algorithm::ConvergenceFreezeTracker freeze_tracker{
        atom_size,
        kLocalFittingParameterChangeTolerance,
        kLocalFittingFreezeChangeRatio,
        kLocalFittingFreezeStableIterations
    };
    LocalFittingFallbackStats fallback_stats{ atom_size };
    const auto atom_index_map{ BuildSelectedAtomIndexMap(atom_list) };
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

        auto iteration_result{
            RunLocalFittingIteration(
                atom_list,
                active_index_list,
                previous_state,
                options)
        };
        fallback_stats.Accumulate(iteration_result.diagnostics);
        const auto raw_state{ std::move(iteration_result.state) };
        GaussianFittingState current_state;
        std::vector<algorithm::ParameterChange> change_list;
        algorithm::ParameterChangeStats change_stats;
        algorithm::ParameterChangeStats normalized_change_stats;
        algorithm::FittingQualityCandidateStats current_candidate_stats;
        double beta{ relaxation_controller.GetBeta() };
        bool has_current_candidate{ false };
        for (int attempt = 0; attempt < kLocalFittingObjectiveBacktrackingMaximumAttempts; attempt++)
        {
            beta = relaxation_controller.GetBeta();
            auto attempt_state{ raw_state };
            ApplyLocalFittingUnderRelaxation(attempt_state, previous_state, beta);
            auto attempt_change_list{
                CalculateLocalFittingParameterChanges(
                    attempt_state.estimation_list,
                    previous_state.estimation_list)
            };
            auto attempt_change_stats{
                SummarizeLocalFittingParameterChangeStats(attempt_change_list, active_index_list)
            };
            auto attempt_normalized_change_list{
                CalculateLocalFittingNormalizedParameterChanges(
                    attempt_state.estimation_list,
                    previous_state.estimation_list)
            };
            auto attempt_normalized_change_stats{
                SummarizeLocalFittingParameterChangeStats(
                    attempt_normalized_change_list,
                    active_index_list)
            };
            const auto attempt_objective_stats{
                CalculateLocalFittingObjectiveStats(
                    atom_list,
                    attempt_state,
                    objective_reference)
            };
            auto attempt_candidate_stats{
                BuildLocalFittingCandidateStats(
                    attempt_objective_stats,
                    attempt_normalized_change_stats)
            };
            algorithm::FittingQualityBacktrackingDecision backtracking_decision{
                true,
                false,
                false
            };
            if (objective_reference.has_reference)
            {
                backtracking_decision = algorithm::EvaluateFittingQualityBacktracking(
                    attempt_candidate_stats,
                    previous_candidate_stats,
                    has_best_candidate,
                    best_candidate_stats,
                    kLocalFittingConvergenceObjectiveRelativeTolerance,
                    attempt,
                    kLocalFittingObjectiveBacktrackingMaximumAttempts);
            }

            if (backtracking_decision.accepted)
            {
                current_state = std::move(attempt_state);
                change_list = std::move(attempt_change_list);
                change_stats = std::move(attempt_change_stats);
                normalized_change_stats = std::move(attempt_normalized_change_stats);
                current_candidate_stats = std::move(attempt_candidate_stats);
                has_current_candidate = true;
                break;
            }
            if (backtracking_decision.should_shrink_beta && !relaxation_controller.IsAtMinimum())
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
                        << next_beta;
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
        const auto thaw_count{
            ThawChangedActiveAtomNeighbors(
                atom_list,
                atom_index_map,
                change_list,
                active_index_list,
                freeze_tracker)
        };

        if (!options.quiet_mode)
        {
            std::ostringstream progress_message;
            progress_message << "Local fitting iteration " << iter + 1 << '/'
                << kLocalFittingMaximumIterations
                << std::fixed << std::setprecision(5)
                << ", percentile amplitude change = "
                << GetLocalFittingParameterChangePercentile(change_stats, kAmplitudeChangeIndex)
                << ", percentile width change = "
                << GetLocalFittingParameterChangePercentile(change_stats, kWidthChangeIndex)
                << ", percentile offset change = "
                << GetLocalFittingParameterChangePercentile(change_stats, kOffsetChangeIndex)
                << ", objective = "
                << current_candidate_stats.quality_objective
                << ", beta = "
                << beta
                << ", active atoms = "
                << freeze_tracker.GetActiveCount()
                << ", frozen atoms = "
                << freeze_tracker.GetFrozenCount()
                << ", thawed atoms = "
                << thaw_count;
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
                        GetLocalFittingParameterChangePercentile(
                            normalized_change_stats,
                            kAmplitudeChangeIndex)) +
                    ", normalized percentile width change = " +
                    std::to_string(
                        GetLocalFittingParameterChangePercentile(
                            normalized_change_stats,
                            kWidthChangeIndex)) +
                    ", and normalized percentile offset change = " +
                    std::to_string(
                        GetLocalFittingParameterChangePercentile(
                            normalized_change_stats,
                            kOffsetChangeIndex)) +
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
                        GetLocalFittingParameterChangePercentile(
                            best_candidate_stats.parameter_change_stats,
                            kAmplitudeChangeIndex)) +
                    ", normalized percentile width change = " +
                    std::to_string(
                        GetLocalFittingParameterChangePercentile(
                            best_candidate_stats.parameter_change_stats,
                            kWidthChangeIndex)) +
                    ", and normalized percentile offset change = " +
                    std::to_string(
                        GetLocalFittingParameterChangePercentile(
                            best_candidate_stats.parameter_change_stats,
                            kOffsetChangeIndex)) +
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
