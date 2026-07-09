#include <cstddef>
#include <rhbm_gem/core/GaussianEstimator.hpp>

#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/algorithm/AndersonAcceleration.hpp>
#include <rhbm_gem/utils/algorithm/ClusteredAndersonAcceleration.hpp>
#include <rhbm_gem/utils/algorithm/ClusteredFittingQualityState.hpp>
#include <rhbm_gem/utils/algorithm/ConvergenceFreezeTracker.hpp>
#include <rhbm_gem/utils/algorithm/DependencyThawHysteresisTracker.hpp>
#include <rhbm_gem/utils/algorithm/IterationState.hpp>
#include <rhbm_gem/utils/algorithm/NormalizedChange.hpp>
#include <rhbm_gem/utils/algorithm/ParameterChangeStats.hpp>
#include <rhbm_gem/utils/algorithm/ScaleReferenceTracker.hpp>
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
constexpr double kLocalFittingNormalizedChangeTolerance{ 1.0e-4 };
constexpr double kLocalFittingNormalizedChangeScaleFloor{ 1.0 };
constexpr std::size_t kMinimumAlphaRTrainingSampleCount{ 10 };
constexpr std::size_t kMinimumAlphaGTrainingMemberCount{ 10 };
constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr double kNeighborAtomSearchRange{ 2.0 * kNeighborContributionDistanceMax };
constexpr std::size_t kLocalFittingMaximumIterations{ 200 };
constexpr double kLocalFittingParameterChangeTolerance{ 1.0e-6 };
constexpr double kLocalFittingChangePercentile{ 0.99 };
constexpr int kHuberSlopeMaximumIterations{ 50 };
constexpr double kHuberScaleMultiplier{ 1.4826 };
constexpr double kHuberScaleMin{ 1.0e-12 };
constexpr double kHuberCutoffMultiplier{ 1.345 };
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
constexpr std::size_t kLocalFittingAndersonHistoryDepth{ 5 };
constexpr double kLocalFittingAndersonCoefficientL1Limit{ 10.0 };
constexpr double kLocalFittingAndersonRegularization{ 1.0e-12 };
constexpr std::array<double, 3> kLocalFittingAccelerationDampingList{ 1.0, 0.5, 0.25 };
constexpr double kLocalFittingFreezeChangeRatio{ 0.1 };
constexpr int kLocalFittingFreezeStableIterations{ 3 };
constexpr double kLocalFittingDependencyThawHysteresisGrowth{ 2.0 };
constexpr double kLocalFittingDependencyThawHysteresisMax{ 8.0 };
constexpr double kLocalFittingDependencyThawHysteresisFrozenDecay{ 0.9 };
constexpr int kLocalFittingDependencyThawMaximumCount{ 5 };
constexpr double kLocalFittingObjectiveTieRelativeTolerance{ 1.0e-8 };
constexpr double kLocalFittingConvergenceObjectiveRelativeTolerance{ 1.0e-3 };
constexpr std::size_t kLocalFittingObjectiveScaleWarmupCount{ 5 };
constexpr double kLocalFittingObjectiveResidualScaleFloorRatio{ 1.0e-6 };
constexpr std::size_t kSuspiciousOffsetClusterMaxDepth{ 2 };
constexpr double kSuspiciousOffsetClusterMinimumOverlap{ 0.05 };

using GaussianFittingState = algorithm::IterationState<LocalGaussianResult, Eigen::VectorXd>;

enum class LocalFittingAccelerationKind
{
    Anderson,
    DampedAnderson,
    DampedFixedPoint
};

struct LocalFittingAccelerationAttempt
{
    LocalFittingAccelerationKind kind{ LocalFittingAccelerationKind::DampedFixedPoint };
    double damping{ 1.0 };
};

struct ActiveCouplingEdge
{
    std::size_t neighbor_index{ 0 };
    double overlap{ 0.0 };
};

using ActiveCouplingGraph = std::vector<std::vector<ActiveCouplingEdge>>;

struct JointOffsetSolveResult
{
    Eigen::VectorXd offset{};
    bool used_fallback{ false };
    ActiveCouplingGraph active_coupling_graph{};
};

struct JointOffsetBuildResult
{
    algorithm::WeightedRidgeSystem system{};
    ActiveCouplingGraph active_coupling_graph{};
};

std::vector<double> CombineLocalFittingRidgeMultiplierLists(
    const std::vector<double> & left_multiplier_list,
    const std::vector<double> & right_multiplier_list)
{
    if (left_multiplier_list.size() != right_multiplier_list.size())
    {
        throw std::invalid_argument("Local fitting ridge multiplier sizes are inconsistent.");
    }
    std::vector<double> multiplier_list(left_multiplier_list.size(), 1.0);
    for (std::size_t i = 0; i < multiplier_list.size(); i++)
    {
        multiplier_list.at(i) = std::max(left_multiplier_list.at(i), right_multiplier_list.at(i));
    }
    return multiplier_list;
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
    std::vector<std::size_t> suspicious_offset_state_index_list{};
};

struct LocalFittingObjectiveStats
{
    bool has_quality_objective{ false };
    double quality_objective{ std::numeric_limits<double>::infinity() };
};

struct LocalFittingObjectiveSamples
{
    std::vector<double> residual_list{};
    std::vector<double> response_list{};
};

struct LocalFittingObjectiveSampleRef
{
    std::size_t atom_index{ 0 };
    std::size_t sample_index{ 0 };
};

using LocalFittingAndersonClusterKey = algorithm::ClusterKey;

struct LocalFittingAndersonCluster
{
    LocalFittingAndersonClusterKey active_index_list{};
    std::vector<LocalFittingObjectiveSampleRef> objective_sample_ref_list{};
};

struct LocalFittingClusterLayout
{
    std::vector<LocalFittingAndersonCluster> cluster_list{};
    std::vector<LocalFittingAndersonClusterKey> key_list{};
    std::map<LocalFittingAndersonClusterKey, std::size_t> index_by_key{};
};

struct LocalFittingClusterAttemptStatus
{
    bool accepted{ false };
    bool stopped{ false };
    bool accepted_by_anderson{ false };
    bool rejected{ false };
    bool rejected_anderson{ false };
};

struct SecondStageNeighborSample
{
    std::size_t atom_index{ 0 };
    double distance{ 0.0 };
};

struct SecondStageAtomContext
{
    LocalPotentialSampleList sample_entries{};
    std::vector<std::size_t> selected_neighbor_index_list{};
    std::vector<std::vector<SecondStageNeighborSample>> sample_neighbor_list{};
    double alpha_r{ 0.0 };
};

struct SecondStageLocalFittingContext
{
    std::vector<AtomObject *> atom_list{};
    std::unordered_map<const AtomObject *, std::size_t> atom_index_map{};
    std::vector<SecondStageAtomContext> atom_context_list{};

    std::size_t AtomSize() const { return atom_list.size(); }
};

struct ParameterSummaryStats
{
    double mean{ 0.0 };
    double standard_deviation{ 0.0 };
};

struct GaussianModelParameterSamples
{
    std::vector<double> amplitude_list{};
    std::vector<double> width_list{};
    std::vector<double> offset_list{};
};

enum class LocalFittingPass
{
    FirstStage,
    ThirdStage
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

double CalculateZeroOffsetResponse(
    const LocalPotentialSample & sample,
    const GaussianModel3D & model)
{
    const auto distance{ static_cast<double>(sample.point.distance) };
    const auto model_offset{ model.ResponseAtDistance(distance) - model.SignalAtDistance(distance) };
    return static_cast<double>(sample.response) - model_offset;
}

bool CanBuildFiniteZeroOffsetSamples(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model)
{
    for (const auto & sample : sample_entries)
    {
        const auto response{ CalculateZeroOffsetResponse(sample, model) };
        if (!std::isfinite(response)) return false;
        if (std::abs(response) > static_cast<double>(std::numeric_limits<float>::max())) return false;
    }
    return true;
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
        const auto response{ static_cast<float>(CalculateZeroOffsetResponse(sample, model)) };
        updated_sample_entries.emplace_back(LocalPotentialSample{ response, sample.point });
    }
    return updated_sample_entries;
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

GroupGaussianResult DecodeGroupGaussianResult(
    double alpha_g,
    const RHBMGroupEstimationResult & result,
    double offset)
{
    const auto prior{
        linearization_service::DecodeParameterVector(result.mu_prior, result.capital_lambda)
    };
    return GroupGaussianResult{
        alpha_g,
        linearization_service::DecodeParameterVector(result.mu_mean).WithOffset(offset),
        linearization_service::DecodeParameterVector(result.mu_mdpde).WithOffset(offset),
        GaussianModel3DWithUncertainty{
            prior.GetModel().WithOffset(offset),
            prior.GetStandardDeviationModel()
        }
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
            linearization_service::DecodeParameterVector(
                result.beta_posterior_matrix.col(i),
                result.capital_sigma_posterior_list.at(member_index))
        };
        const auto gaussian_with_offset{
            GaussianModel3DWithUncertainty{
                gaussian.GetModel().WithOffset(offset),
                gaussian.GetStandardDeviationModel()
            }
        };
        member_results.emplace_back(LocalGaussianResult{
            0.0,
            gaussian_with_offset,
            gaussian_with_offset,
            gaussian_with_offset,
            static_cast<bool>(result.outlier_flag_array(i)),
            result.statistical_distance_array(i)
        });
    }
    return member_results;
}

using FittedGaussianSnapshot = std::vector<GaussianModel3D>;
using GroupMedianModelMap = std::unordered_map<GroupKey, GaussianModel3D>;

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

SecondStageLocalFittingContext BuildSecondStageLocalFittingContext(ModelObject & model_object)
{
    SecondStageLocalFittingContext context;
    context.atom_list = model_object.GetSelectedAtoms();
    context.atom_index_map = BuildSelectedAtomIndexMap(context.atom_list);
    context.atom_context_list.resize(context.atom_list.size());

    for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
    {
        const auto * atom{ context.atom_list.at(atom_index) };
        auto & atom_context{ context.atom_context_list.at(atom_index) };
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        atom_context.sample_entries = local_view.GetSamplingEntries(false);
        atom_context.alpha_r = local_view.GetAlphaR();
    }

    for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
    {
        const auto * atom{ context.atom_list.at(atom_index) };
        auto & atom_context{ context.atom_context_list.at(atom_index) };
        const auto neighbor_atom_list{ atom->FindNeighborAtoms(kNeighborAtomSearchRange) };
        atom_context.selected_neighbor_index_list.reserve(neighbor_atom_list.size());
        for (const auto * neighbor_atom : neighbor_atom_list)
        {
            const auto neighbor_iter{ context.atom_index_map.find(neighbor_atom) };
            if (neighbor_iter == context.atom_index_map.end()) continue;

            atom_context.selected_neighbor_index_list.emplace_back(neighbor_iter->second);
        }

        atom_context.sample_neighbor_list.resize(atom_context.sample_entries.size());
        for (std::size_t sample_index = 0; sample_index < atom_context.sample_entries.size(); sample_index++)
        {
            const auto & sample{ atom_context.sample_entries.at(sample_index) };
            auto & sample_neighbor_list{ atom_context.sample_neighbor_list.at(sample_index) };
            sample_neighbor_list.reserve(atom_context.selected_neighbor_index_list.size());
            for (const auto neighbor_index : atom_context.selected_neighbor_index_list)
            {
                const auto distance{
                    static_cast<double>(
                        array_helper::ComputeNorm<float>(
                            sample.point.position,
                            context.atom_list.at(neighbor_index)->GetPositionRef()))
                };
                if (distance > kNeighborContributionDistanceMax) continue;

                sample_neighbor_list.emplace_back(SecondStageNeighborSample{ neighbor_index, distance });
            }
        }
    }

    return context;
}

std::size_t FindLocalFittingAndersonClusterRoot(std::vector<std::size_t> & parent_list, std::size_t index)
{
    if (index >= parent_list.size())
    {
        throw std::invalid_argument("Local fitting Anderson cluster index is out of range.");
    }
    auto root{ index };
    while (parent_list.at(root) != root)
    {
        root = parent_list.at(root);
    }
    while (parent_list.at(index) != index)
    {
        const auto parent{ parent_list.at(index) };
        parent_list.at(index) = root;
        index = parent;
    }
    return root;
}

void MergeLocalFittingAndersonClusters(
    std::vector<std::size_t> & parent_list,
    std::size_t left_index,
    std::size_t right_index)
{
    const auto left_root{ FindLocalFittingAndersonClusterRoot(parent_list, left_index) };
    const auto right_root{ FindLocalFittingAndersonClusterRoot(parent_list, right_index) };
    if (left_root == right_root) return;

    parent_list.at(right_root) = left_root;
}

std::vector<LocalFittingAndersonCluster> BuildLocalFittingAndersonClusters(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list)
{
    std::vector<int> active_position_by_atom_index(context.AtomSize(), -1);
    for (std::size_t active_position = 0; active_position < active_index_list.size(); active_position++)
    {
        const auto active_index{ active_index_list.at(active_position) };
        if (active_index >= context.AtomSize())
        {
            throw std::invalid_argument("Local fitting Anderson active index is out of range.");
        }
        active_position_by_atom_index.at(active_index) = static_cast<int>(active_position);
    }

    std::vector<std::size_t> parent_list(active_index_list.size());
    for (std::size_t i = 0; i < parent_list.size(); i++)
    {
        parent_list.at(i) = i;
    }

    for (std::size_t atom_index = 0; atom_index < context.AtomSize(); atom_index++)
    {
        const auto & atom_context{ context.atom_context_list.at(atom_index) };
        for (std::size_t sample_index = 0; sample_index < atom_context.sample_entries.size(); sample_index++)
        {
            std::vector<std::size_t> contributor_position_list;
            const auto target_active_position{ active_position_by_atom_index.at(atom_index) };
            if (target_active_position >= 0)
            {
                contributor_position_list.emplace_back(static_cast<std::size_t>(target_active_position));
            }
            for (const auto & neighbor_sample : atom_context.sample_neighbor_list.at(sample_index))
            {
                if (neighbor_sample.atom_index >= active_position_by_atom_index.size())
                {
                    throw std::invalid_argument("Local fitting Anderson neighbor index is out of range.");
                }
                const auto neighbor_active_position{
                    active_position_by_atom_index.at(neighbor_sample.atom_index)
                };
                if (neighbor_active_position >= 0)
                {
                    contributor_position_list.emplace_back(static_cast<std::size_t>(neighbor_active_position));
                }
            }

            if (contributor_position_list.size() < 2) continue;

            std::sort(contributor_position_list.begin(), contributor_position_list.end());
            contributor_position_list.erase(
                std::unique(contributor_position_list.begin(), contributor_position_list.end()),
                contributor_position_list.end());
            for (std::size_t i = 1; i < contributor_position_list.size(); i++)
            {
                MergeLocalFittingAndersonClusters(
                    parent_list,
                    contributor_position_list.front(),
                    contributor_position_list.at(i));
            }
        }
    }

    std::map<std::size_t, LocalFittingAndersonClusterKey> active_index_list_by_root;
    for (std::size_t active_position = 0; active_position < active_index_list.size(); active_position++)
    {
        const auto root{ FindLocalFittingAndersonClusterRoot(parent_list, active_position) };
        active_index_list_by_root[root].emplace_back(active_index_list.at(active_position));
    }

    std::map<std::size_t, std::vector<LocalFittingObjectiveSampleRef>> sample_ref_list_by_root;
    for (std::size_t atom_index = 0; atom_index < context.AtomSize(); atom_index++)
    {
        const auto & atom_context{ context.atom_context_list.at(atom_index) };
        for (std::size_t sample_index = 0; sample_index < atom_context.sample_entries.size(); sample_index++)
        {
            std::vector<std::size_t> contributor_root_list;
            const auto target_active_position{ active_position_by_atom_index.at(atom_index) };
            if (target_active_position >= 0)
            {
                contributor_root_list.emplace_back(
                    FindLocalFittingAndersonClusterRoot(
                        parent_list, static_cast<std::size_t>(target_active_position)));
            }
            for (const auto & neighbor_sample : atom_context.sample_neighbor_list.at(sample_index))
            {
                const auto neighbor_active_position{
                    active_position_by_atom_index.at(neighbor_sample.atom_index)
                };
                if (neighbor_active_position >= 0)
                {
                    contributor_root_list.emplace_back(
                        FindLocalFittingAndersonClusterRoot(
                            parent_list,
                            static_cast<std::size_t>(neighbor_active_position)));
                }
            }
            if (contributor_root_list.empty()) continue;

            std::sort(contributor_root_list.begin(), contributor_root_list.end());
            contributor_root_list.erase(
                std::unique(contributor_root_list.begin(), contributor_root_list.end()),
                contributor_root_list.end());
            if (contributor_root_list.size() != 1)
            {
                throw std::logic_error(
                    "Local fitting objective sample spans multiple Anderson clusters.");
            }
            sample_ref_list_by_root[contributor_root_list.front()].emplace_back(
                LocalFittingObjectiveSampleRef{ atom_index, sample_index });
        }
    }

    std::vector<LocalFittingAndersonCluster> cluster_list;
    cluster_list.reserve(active_index_list_by_root.size());
    for (auto & [root, cluster_active_index_list] : active_index_list_by_root)
    {
        std::sort(cluster_active_index_list.begin(), cluster_active_index_list.end());
        cluster_list.emplace_back(LocalFittingAndersonCluster{
            std::move(cluster_active_index_list),
            std::move(sample_ref_list_by_root[root])
        });
    }
    return cluster_list;
}

GaussianFittingState BuildInitialLocalFittingState(const SecondStageLocalFittingContext & context)
{
    GaussianFittingState state{
        std::vector<LocalGaussianResult>(context.AtomSize()),
        std::vector<Eigen::VectorXd>(context.AtomSize())
    };
    for (std::size_t i = 0; i < context.AtomSize(); i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*context.atom_list.at(i)) };
        state.result_list.at(i) = local_view.GetGaussianResult();
        state.estimation_list.at(i) = state.result_list.at(i).mdpde.GetModel().ToVector();
    }
    return state;
}

FittedGaussianSnapshot BuildFittedGaussianSnapshot(
    const SecondStageLocalFittingContext & context,
    const std::vector<Eigen::VectorXd> & estimation_list)
{
    if (context.AtomSize() != estimation_list.size())
    {
        throw std::invalid_argument("atom context and estimation_list sizes are inconsistent.");
    }

    FittedGaussianSnapshot snapshot;
    snapshot.reserve(context.AtomSize());
    for (const auto & estimation : estimation_list)
    {
        snapshot.emplace_back(GaussianModel3D::FromVector(estimation));
    }
    return snapshot;
}

GroupMedianModelMap BuildGroupMedianMDPDEModelMap(const std::vector<AtomObject *> & atom_list)
{
    std::unordered_map<GroupKey, GaussianModelParameterSamples> parameter_samples_by_group;
    parameter_samples_by_group.reserve(atom_list.size());
    for (const auto * atom : atom_list)
    {
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        if (!local_view.IsAvailable()) continue;

        const auto & model{ local_view.GetEstimateMDPDE() };
        auto & parameter_samples{
            parameter_samples_by_group[data_internal::GetGroupKey(atom)]
        };
        parameter_samples.amplitude_list.emplace_back(model.GetAmplitude());
        parameter_samples.width_list.emplace_back(model.GetWidth());
        parameter_samples.offset_list.emplace_back(model.GetOffset());
    }

    GroupMedianModelMap median_model_by_group;
    median_model_by_group.reserve(parameter_samples_by_group.size());
    for (const auto & [group_key, parameter_samples] : parameter_samples_by_group)
    {
        if (parameter_samples.amplitude_list.empty()) continue;

        median_model_by_group.emplace(
            group_key,
            GaussianModel3D{
                array_helper::ComputeMedian(parameter_samples.amplitude_list),
                array_helper::ComputeMedian(parameter_samples.width_list),
                array_helper::ComputeMedian(parameter_samples.offset_list)
            });
    }
    return median_model_by_group;
}

JointOffsetBuildResult BuildJointOffsetSystem(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const FittedGaussianSnapshot & snapshot,
    double ridge_ratio,
    const std::vector<double> & ridge_multiplier_list)
{
    const auto atom_size{ context.AtomSize() };
    if (snapshot.size() != atom_size || ridge_multiplier_list.size() != atom_size)
    {
        throw std::invalid_argument("Joint offset input sizes are inconsistent.");
    }

    std::vector<int> active_column_by_atom_index(atom_size, -1);
    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto atom_index{ active_index_list.at(i) };
        if (atom_index >= atom_size)
        {
            throw std::invalid_argument("Joint offset active atom index is out of range.");
        }
        active_column_by_atom_index.at(atom_index) = static_cast<int>(i);
    }

    const auto column_count{ static_cast<Eigen::Index>(active_index_list.size()) };
    std::vector<Eigen::Triplet<double>> triplet_list;
    std::vector<double> response_list;
    Eigen::VectorXd column_square_sum{ Eigen::VectorXd::Zero(column_count) };
    std::map<std::pair<Eigen::Index, Eigen::Index>, double> column_cross_sum_map;
    ActiveCouplingGraph active_coupling_graph(active_index_list.size());
    std::vector<std::pair<Eigen::Index, double>> row_basis_entries;
    for (const auto active_index : active_index_list)
    {
        if (active_index >= atom_size)
        {
            throw std::invalid_argument("Joint offset active atom index is out of range.");
        }
        const auto target_column{ active_column_by_atom_index.at(active_index) };
        if (target_column < 0)
        {
            throw std::invalid_argument("Joint offset active column is missing.");
        }

        const auto & atom_context{ context.atom_context_list.at(active_index) };
        const auto & target_model{ snapshot.at(active_index) };
        row_basis_entries.reserve(atom_context.selected_neighbor_index_list.size() + 1);
        for (std::size_t sample_index = 0; sample_index < atom_context.sample_entries.size(); sample_index++)
        {
            const auto & sample{ atom_context.sample_entries.at(sample_index) };
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
            row_basis_entries.clear();
            if (std::abs(target_basis) > std::numeric_limits<double>::epsilon())
            {
                row_basis_entries.emplace_back(static_cast<Eigen::Index>(target_column), target_basis);
            }

            for (const auto & neighbor_sample : atom_context.sample_neighbor_list.at(sample_index))
            {
                const auto & neighbor_model{ snapshot.at(neighbor_sample.atom_index) };
                const auto neighbor_column{ active_column_by_atom_index.at(neighbor_sample.atom_index) };
                if (neighbor_column < 0)
                {
                    const auto response{ neighbor_model.ResponseAtDistance(neighbor_sample.distance) };
                    if (!std::isfinite(response))
                    {
                        throw std::runtime_error("Joint offset fixed neighbor model evaluation is not finite.");
                    }
                    residual -= response;
                    continue;
                }

                const auto signal{ neighbor_model.SignalAtDistance(neighbor_sample.distance) };
                const auto basis{ neighbor_model.OffsetBasisAtDistance(neighbor_sample.distance) };
                if (!std::isfinite(signal) || !std::isfinite(basis))
                {
                    throw std::runtime_error("Joint offset active neighbor model evaluation is not finite.");
                }
                residual -= signal;
                if (std::abs(basis) > std::numeric_limits<double>::epsilon())
                {
                    row_basis_entries.emplace_back(static_cast<Eigen::Index>(neighbor_column), basis);
                }
            }
            if (!std::isfinite(residual))
            {
                throw std::runtime_error("Joint offset residual is not finite.");
            }
            if (row_basis_entries.empty()) continue;

            const auto row_index{ static_cast<Eigen::Index>(response_list.size()) };
            response_list.emplace_back(residual);
            for (const auto & [column_index, basis] : row_basis_entries)
            {
                triplet_list.emplace_back(row_index, column_index, basis);
                column_square_sum(column_index) += basis * basis;
            }
            for (std::size_t i = 0; i < row_basis_entries.size(); i++)
            {
                const auto [left_column, left_basis]{ row_basis_entries.at(i) };
                for (std::size_t j = i + 1; j < row_basis_entries.size(); j++)
                {
                    const auto [right_column, right_basis]{ row_basis_entries.at(j) };
                    if (left_column == right_column) continue;
                    const auto column_pair{ std::minmax(left_column, right_column) };
                    column_cross_sum_map[column_pair] += left_basis * right_basis;
                }
            }
        }
    }

    Eigen::VectorXd proactive_ridge_multiplier{ Eigen::VectorXd::Ones(column_count) };
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
        const auto overlap{ std::abs(cross_sum) / std::sqrt(left_square_sum * right_square_sum) };
        if (!std::isfinite(overlap))
        {
            continue;
        }
        active_coupling_graph.at(static_cast<std::size_t>(left_column)).emplace_back(
            ActiveCouplingEdge{ static_cast<std::size_t>(right_column), overlap });
        active_coupling_graph.at(static_cast<std::size_t>(right_column)).emplace_back(
            ActiveCouplingEdge{ static_cast<std::size_t>(left_column), overlap });
        if (overlap < kJointOffsetCollinearityOverlapThreshold) continue;

        proactive_ridge_multiplier(left_column) = std::max(
            proactive_ridge_multiplier(left_column),
            kCollinearJointOffsetRidgeMultiplier);
        proactive_ridge_multiplier(right_column) = std::max(
            proactive_ridge_multiplier(right_column),
            kCollinearJointOffsetRidgeMultiplier);
    }

    const auto row_count{ static_cast<Eigen::Index>(response_list.size()) };
    Eigen::VectorXd response{ Eigen::VectorXd::Zero(row_count) };
    for (Eigen::Index row_index = 0; row_index < row_count; row_index++)
    {
        response(row_index) = response_list.at(static_cast<std::size_t>(row_index));
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
        const auto & model{ snapshot.at(atom_index) };
        system.previous_parameter(column_index) = model.GetOffset();
        const auto square_sum{ column_square_sum(column_index) };
        const auto multiplier{ ridge_multiplier_list.at(atom_index) };
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
        std::move(active_coupling_graph)
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
    if (system.previous_parameter.size() != offset.size() || system.ridge_diagonal.size() != offset.size())
    {
        throw std::invalid_argument("Weighted ridge objective parameter sizes are inconsistent.");
    }
    if (system.response.size() == 0)
    {
        return std::numeric_limits<double>::infinity();
    }

    const Eigen::VectorXd residual{ system.response - system.design_matrix * offset };
    const auto weighted_residual_loss{ weight.cwiseProduct(residual.cwiseAbs2()).sum() };
    const Eigen::VectorXd offset_delta{ offset - system.previous_parameter };
    const auto ridge_loss{ system.ridge_diagonal.cwiseProduct(offset_delta.cwiseAbs2()).sum() };
    const auto objective{ (weighted_residual_loss + ridge_loss) / static_cast<double>(system.response.size()) };
    return std::isfinite(objective) ? objective : std::numeric_limits<double>::infinity();
}

bool IsJointOffsetObjectiveDeteriorated(double updated_objective, double current_objective)
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
    return updated_objective > current_objective + kJointOffsetIrlsObjectiveRelativeTolerance * scale;
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

JointOffsetSolveResult EstimateJointOffsets(
    const SecondStageLocalFittingContext & context,
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
        if (atom_index >= context.AtomSize())
        {
            throw std::invalid_argument("Joint offset active atom index is out of range.");
        }
        previous_offset(static_cast<Eigen::Index>(i)) = snapshot.at(atom_index).GetOffset();
    }
    JointOffsetBuildResult build_result;
    try
    {
        build_result = BuildJointOffsetSystem(context, active_index_list, snapshot, ridge_ratio, ridge_multiplier_list);
    }
    catch (const std::exception &)
    {
        return JointOffsetSolveResult{previous_offset, true, {} };
    }
    auto system{ std::move(build_result.system) };
    auto active_coupling_graph{ std::move(build_result.active_coupling_graph) };
    if (system.response.size() == 0 || system.previous_parameter.size() == 0)
    {
        return JointOffsetSolveResult{
            previous_offset,
            true,
            std::move(active_coupling_graph)
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
            std::move(active_coupling_graph)
        };
    }

    for (int iteration = 0; iteration < kHuberSlopeMaximumIterations; iteration++)
    {
        const Eigen::VectorXd residual{ system.response - system.design_matrix * offset };
        std::vector<double> residual_list(residual.data(), residual.data() + residual.size());
        const auto residual_scale{
            std::max(CalculateMedianAbsoluteDeviationScale(residual_list), kHuberScaleMin)
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
                std::move(active_coupling_graph)
            };
        }
        const auto current_objective{
            CalculateWeightedRidgeSurrogateObjective(system, weight, offset)
        };
        const auto updated_objective{
            CalculateWeightedRidgeSurrogateObjective(system, weight, updated_offset)
        };
        if (IsJointOffsetObjectiveDeteriorated(updated_objective, current_objective)) break;
        const auto maximum_change{
            algorithm::CalculateMaximumNormalizedVectorChange(updated_offset, offset, kJointOffsetIrlsScaleFloor)
        };
        offset = std::move(updated_offset);
        if (maximum_change < kJointOffsetIrlsNormalizedChangeTolerance) break;
    }

    return JointOffsetSolveResult{
        offset,
        false,
        std::move(active_coupling_graph)
    };
}

void ApplyJointOffsetsToSnapshot(
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
        if (atom_index >= snapshot.size())
        {
            throw std::invalid_argument("Joint offset active atom index is out of range.");
        }
        snapshot.at(atom_index) = snapshot.at(atom_index).WithOffset(offset(static_cast<Eigen::Index>(i)));
    }
}

template <typename GaussianLookup>
LocalPotentialSampleList UpdateSampleListWithGaussianLookup(
    const AtomObject & atom,
    GaussianLookup lookup_gaussian)
{
    const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
    const auto sample_entries{ local_view.GetSamplingEntries(false) };
    const auto & neighbor_atom_list{ atom.FindNeighborAtoms(kNeighborAtomSearchRange) };
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
                static_cast<double>(array_helper::ComputeNorm<float>(sample_position, neighbor_position))
            };
            if (distance > kNeighborContributionDistanceMax) continue;
            response_value -= static_cast<float>(gaussian->ResponseAtDistance(distance));
        }
        updated_list.emplace_back(LocalPotentialSample{response_value, sample.point });
    }
    return updated_list;
}

double CalculateSecondStageAdjustedResponse(
    const SecondStageLocalFittingContext & context,
    std::size_t atom_index,
    std::size_t sample_index,
    const FittedGaussianSnapshot & snapshot)
{
    const auto & atom_context{ context.atom_context_list.at(atom_index) };
    auto response_value{ static_cast<double>(atom_context.sample_entries.at(sample_index).response) };
    for (const auto & neighbor_sample : atom_context.sample_neighbor_list.at(sample_index))
    {
        response_value -= snapshot.at(neighbor_sample.atom_index).ResponseAtDistance(neighbor_sample.distance);
    }
    return response_value;
}

LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const SecondStageLocalFittingContext & context,
    std::size_t atom_index,
    const FittedGaussianSnapshot & snapshot)
{
    const auto & atom_context{ context.atom_context_list.at(atom_index) };
    LocalPotentialSampleList updated_list;
    updated_list.reserve(atom_context.sample_entries.size());
    for (std::size_t sample_index = 0; sample_index < atom_context.sample_entries.size(); sample_index++)
    {
        auto sample{ atom_context.sample_entries.at(sample_index) };
        sample.response = static_cast<float>(CalculateSecondStageAdjustedResponse(context, atom_index, sample_index, snapshot));
        updated_list.emplace_back(sample);
    }
    return updated_list;
}

LocalPotentialSampleList UpdateSampleListWithGroupMedianGaussian(
    const AtomObject & atom,
    const GroupMedianModelMap & median_model_by_group)
{
    return UpdateSampleListWithGaussianLookup(
        atom,
        [&median_model_by_group](const AtomObject & neighbor_atom) -> const GaussianModel3D *
        {
            const auto median_model_iter{
                median_model_by_group.find(data_internal::GetGroupKey(&neighbor_atom))
            };
            if (median_model_iter != median_model_by_group.end())
            {
                return &median_model_iter->second;
            }

            const auto local_view{ AtomLocalPotentialView::For(neighbor_atom) };
            return local_view.IsAvailable() ? &local_view.GetEstimateMDPDE() : nullptr;
        });
}

LocalPotentialSampleList UpdateSampleListWithFittedGroupGaussian(
    const AtomObject & atom,
    const ModelAnalysisView & analysis_view)
{
    return UpdateSampleListWithGaussianLookup(
        atom,
        [&analysis_view](const AtomObject & neighbor_atom) -> const GaussianModel3D *
        {
            const auto group_key{ data_internal::GetGroupKey(&neighbor_atom) };
            if (!analysis_view.HasAtomGroup(group_key))
            {
                return nullptr;
            }
            return &analysis_view.GetAtomGroupPrior(group_key);
        });
}

void SetUpdatedSamplingEntriesFromGroupMedianGaussian(ModelObject & model_object)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    const auto median_model_by_group{ BuildGroupMedianMDPDEModelMap(atom_list) };
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        local_editor_list[i].SetUpdatedSamplingEntries(
            UpdateSampleListWithGroupMedianGaussian(*atom_list[i], median_model_by_group));
    }
}

void SetUpdatedSamplingEntriesFromFittedGroupGaussian(ModelObject & model_object)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    const auto analysis_view{ model_object.GetAnalysisView() };
    for (size_t i = 0; i < atom_list.size(); i++)
    {
        local_editor_list[i].SetUpdatedSamplingEntries(
            UpdateSampleListWithFittedGroupGaussian(*atom_list[i], analysis_view));
    }
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
    const SecondStageLocalFittingContext & context,
    const FittedGaussianSnapshot & snapshot,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list)
{
    if (snapshot.size() != context.AtomSize())
    {
        throw std::invalid_argument("Local fitting objective snapshot size is inconsistent.");
    }

    LocalFittingObjectiveSamples objective_samples;
    objective_samples.residual_list.reserve(sample_ref_list.size());
    objective_samples.response_list.reserve(sample_ref_list.size());
    for (const auto & sample_ref : sample_ref_list)
    {
        if (sample_ref.atom_index >= context.AtomSize())
        {
            throw std::invalid_argument("Local fitting objective sample atom index is out of range.");
        }
        const auto & atom_context{ context.atom_context_list.at(sample_ref.atom_index) };
        if (sample_ref.sample_index >= atom_context.sample_entries.size())
        {
            throw std::invalid_argument("Local fitting objective sample index is out of range.");
        }

        const auto & sample{ atom_context.sample_entries.at(sample_ref.sample_index) };
        const auto & target_model{ snapshot.at(sample_ref.atom_index) };
        const auto distance{ static_cast<double>(sample.point.distance) };
        const auto expected_response{ target_model.ResponseAtDistance(distance) };
        const auto response{
            CalculateSecondStageAdjustedResponse(
                context, sample_ref.atom_index, sample_ref.sample_index, snapshot)
        };
        const auto residual{ response - expected_response };
        if (!std::isfinite(response) || !std::isfinite(residual)) return std::nullopt;
        objective_samples.residual_list.emplace_back(residual);
        objective_samples.response_list.emplace_back(response);
    }
    return objective_samples;
}

std::optional<LocalFittingObjectiveSamples> CollectLocalFittingObjectiveSamples(
    const SecondStageLocalFittingContext & context,
    const GaussianFittingState & fitting_state,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list)
{
    const auto snapshot{
        BuildFittedGaussianSnapshot(context, fitting_state.estimation_list)
    };
    return CollectLocalFittingObjectiveSamples(context, snapshot, sample_ref_list);
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
    if (!std::isfinite(scale_sample)) return std::nullopt;
    return scale_sample;
}

LocalFittingObjectiveStats CalculateLocalFittingObjectiveStats(
    const LocalFittingObjectiveSamples & objective_samples,
    const algorithm::ScaleReference & objective_reference,
    double residual_scale_sample)
{
    LocalFittingObjectiveStats stats;
    if (!objective_reference.has_reference || !std::isfinite(residual_scale_sample))
    {
        return stats;
    }

    double loss_sum{ 0.0 };
    for (const auto residual : objective_samples.residual_list)
    {
        const auto normalized_residual{ residual / objective_reference.scale };
        loss_sum += CalculateHuberLoss(normalized_residual, kHuberCutoffMultiplier);
    }
    const auto quality_objective{ loss_sum / static_cast<double>(objective_samples.residual_list.size()) };
    if (!std::isfinite(quality_objective)) return stats;

    stats.has_quality_objective = true;
    stats.quality_objective = quality_objective;
    return stats;
}

LocalFittingObjectiveStats CalculateLocalFittingObjectiveStats(
    const LocalFittingObjectiveSamples & objective_samples,
    const algorithm::ScaleReference & objective_reference)
{
    const auto residual_scale_sample{
        CalculateLocalFittingResidualScaleSample(objective_samples)
    };
    if (!residual_scale_sample.has_value())
    {
        return LocalFittingObjectiveStats{};
    }
    return CalculateLocalFittingObjectiveStats(
        objective_samples, objective_reference, *residual_scale_sample);
}

algorithm::ClusteredFittingQualityCandidateScore<LocalFittingObjectiveSamples>
ScoreLocalFittingClusterCandidate(
    const SecondStageLocalFittingContext & context,
    const GaussianFittingState & candidate_state,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const algorithm::ScaleReferenceTracker & objective_scale_tracker,
    algorithm::FittingQualityCandidateStats & previous_candidate_stats,
    const std::optional<LocalFittingObjectiveSamples> & previous_objective_samples,
    const std::optional<algorithm::ClusteredFittingQualityTrackedCandidate<LocalFittingObjectiveSamples>> & best_candidate,
    const algorithm::ParameterChangeStats & normalized_change_stats)
{
    algorithm::ClusteredFittingQualityCandidateScore<LocalFittingObjectiveSamples> score;
    auto objective_reference{ objective_scale_tracker.GetCommittedReference() };
    score.has_objective_reference = objective_reference.has_reference;
    if (best_candidate.has_value())
    {
        score.best_candidate_stats = best_candidate->candidate_stats;
    }

    if (objective_scale_tracker.HasReference())
    {
        score.objective_samples =
            CollectLocalFittingObjectiveSamples(context, candidate_state, sample_ref_list);
        if (score.objective_samples.has_value())
        {
            const auto scale_sample{
                CalculateLocalFittingResidualScaleSample(*score.objective_samples)
            };
            if (scale_sample.has_value())
            {
                objective_reference = objective_scale_tracker.GetProvisionalReference(*scale_sample);
                score.has_objective_reference = objective_reference.has_reference;
                const auto objective_stats{
                    CalculateLocalFittingObjectiveStats(
                        *score.objective_samples,
                        objective_reference,
                        *scale_sample)
                };
                score.objective_scale_sample = *scale_sample;
                if (objective_stats.has_quality_objective)
                {
                    if (previous_objective_samples.has_value())
                    {
                        const auto recalculated_previous_objective_stats{
                            CalculateLocalFittingObjectiveStats(
                                *previous_objective_samples,
                                objective_reference)
                        };
                        previous_candidate_stats.has_quality_objective =
                            recalculated_previous_objective_stats.has_quality_objective;
                        previous_candidate_stats.quality_objective =
                            recalculated_previous_objective_stats.quality_objective;
                    }
                    if (best_candidate.has_value() &&
                        best_candidate->objective_samples.has_value() &&
                        score.best_candidate_stats.has_value())
                    {
                        const auto recalculated_best_objective_stats{
                            CalculateLocalFittingObjectiveStats(
                                *best_candidate->objective_samples,
                                objective_reference)
                        };
                        score.best_candidate_stats->has_quality_objective = recalculated_best_objective_stats.has_quality_objective;
                        score.best_candidate_stats->quality_objective = recalculated_best_objective_stats.quality_objective;
                    }
                }
                score.candidate_stats.has_quality_objective = objective_stats.has_quality_objective;
                score.candidate_stats.quality_objective = objective_stats.quality_objective;
            }
        }
    }

    score.candidate_stats.parameter_change_stats = normalized_change_stats;
    return score;
}

ParameterSummaryStats SummarizeParameterValues(const std::vector<double> & value_list)
{
    if (value_list.empty()) return {};

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
    std::map<std::string, GaussianModelParameterSamples> spot_sample_map;
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
        Logger::Log(LogLevel::Info, "Group fitting prior summary by Spot: no atom groups available.");
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
        const auto & current_estimation{ current_estimation_list.at(i) };
        const auto & previous_estimation{ previous_estimation_list.at(i) };
        change_list.at(i).value_list = {
            std::abs(
                current_estimation(GaussianModel3D::AmplitudeIndex()) -
                previous_estimation(GaussianModel3D::AmplitudeIndex())),
            std::abs(
                current_estimation(GaussianModel3D::WidthIndex()) -
                previous_estimation(GaussianModel3D::WidthIndex())),
            std::abs(
                current_estimation(GaussianModel3D::OffsetIndex()) -
                previous_estimation(GaussianModel3D::OffsetIndex()))
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
        throw std::invalid_argument("Local fitting normalized parameter change input sizes are inconsistent.");
    }

    std::vector<algorithm::ParameterChange> change_list(current_estimation_list.size());
    for (size_t i = 0; i < current_estimation_list.size(); i++)
    {
        const auto & current_estimation{ current_estimation_list.at(i) };
        const auto & previous_estimation{ previous_estimation_list.at(i) };
        change_list.at(i).value_list = {
            algorithm::CalculateNormalizedChange(
                current_estimation(GaussianModel3D::AmplitudeIndex()),
                previous_estimation(GaussianModel3D::AmplitudeIndex()),
                kLocalFittingNormalizedChangeScaleFloor),
            algorithm::CalculateNormalizedChange(
                current_estimation(GaussianModel3D::WidthIndex()),
                previous_estimation(GaussianModel3D::WidthIndex()),
                kLocalFittingNormalizedChangeScaleFloor),
            algorithm::CalculateNormalizedChange(
                current_estimation(GaussianModel3D::OffsetIndex()),
                previous_estimation(GaussianModel3D::OffsetIndex()),
                kLocalFittingNormalizedChangeScaleFloor)
        };
    }
    return change_list;
}

bool IsLocalFittingNormalizedParameterChangeConverged(const algorithm::ParameterChangeStats & stats)
{
    for (std::size_t i = 0; i < stats.percentile_list.size(); i++)
    {
        if (stats.percentile_list.at(i) >= kLocalFittingNormalizedChangeTolerance) return false;
    }
    return true;
}

void ApplyLocalFittingDampedFixedPoint(
    GaussianFittingState & current_state,
    const GaussianFittingState & raw_state,
    const GaussianFittingState & previous_state,
    const std::vector<std::size_t> & active_index_list,
    double damping)
{
    if (current_state.estimation_list.size() != previous_state.estimation_list.size() ||
        raw_state.estimation_list.size() != previous_state.estimation_list.size() ||
        current_state.result_list.size() != previous_state.result_list.size() ||
        raw_state.result_list.size() != previous_state.result_list.size())
    {
        throw std::invalid_argument("Local fitting fixed-point damping input sizes are inconsistent.");
    }
    for (const auto active_index : active_index_list)
    {
        if (active_index >= current_state.estimation_list.size())
        {
            throw std::invalid_argument("Local fitting fixed-point damping active index is out of range.");
        }
        auto damped_estimation{
            (damping * raw_state.estimation_list.at(active_index) +
            (1.0 - damping) * previous_state.estimation_list.at(active_index)).eval()
        };
        const auto damped_model{ GaussianModel3D::FromVector(damped_estimation) };
        auto & result{ current_state.result_list.at(active_index) };
        result.mdpde = GaussianModel3DWithUncertainty{
            damped_model,
            raw_state.result_list.at(active_index).mdpde.GetStandardDeviationModel()
        };
        current_state.estimation_list.at(active_index) = damped_estimation;
    }
}

const char * GetLocalFittingAccelerationText(LocalFittingAccelerationKind kind)
{
    switch (kind)
    {
        case LocalFittingAccelerationKind::Anderson: return "aa";
        case LocalFittingAccelerationKind::DampedAnderson: return "damped-aa";
        case LocalFittingAccelerationKind::DampedFixedPoint: return "damped-fixed-point";
    }
    return "unknown";
}

LocalFittingAccelerationAttempt BuildLocalFittingAccelerationAttempt(
    bool use_anderson,
    double damping)
{
    return LocalFittingAccelerationAttempt{
        use_anderson && damping == 1.0 ?
            LocalFittingAccelerationKind::Anderson :
            (use_anderson ?
                LocalFittingAccelerationKind::DampedAnderson :
                LocalFittingAccelerationKind::DampedFixedPoint),
        damping
    };
}

bool IsLocalFittingActiveEstimationFinitePositive(
    const Eigen::VectorXd & estimation)
{
    return estimation.size() == GaussianModel3D::ParameterSize() &&
        estimation.allFinite() &&
        estimation(GaussianModel3D::WidthIndex()) > 0.0;
}

bool TryApplyLocalFittingAndersonCandidate(
    GaussianFittingState & current_state,
    const GaussianFittingState & previous_state,
    const std::vector<Eigen::VectorXd> & candidate_estimation_list,
    const std::vector<std::size_t> & active_index_list,
    double damping)
{
    if (!std::isfinite(damping) || damping <= 0.0 || damping > 1.0)
    {
        throw std::invalid_argument("Local fitting acceleration damping is out of range.");
    }
    if (current_state.estimation_list.size() != previous_state.estimation_list.size() ||
        current_state.result_list.size() != previous_state.result_list.size() ||
        candidate_estimation_list.size() != previous_state.estimation_list.size())
    {
        throw std::invalid_argument("Local fitting acceleration input sizes are inconsistent.");
    }

    std::vector<std::pair<std::size_t, Eigen::VectorXd>> accelerated_estimation_list;
    accelerated_estimation_list.reserve(active_index_list.size());
    for (const auto active_index : active_index_list)
    {
        if (active_index >= previous_state.estimation_list.size() ||
            active_index >= current_state.estimation_list.size())
        {
            throw std::invalid_argument("Local fitting acceleration active index is out of range.");
        }
        const auto accelerated_estimation{
            (previous_state.estimation_list.at(active_index) +
                damping * (
                    candidate_estimation_list.at(active_index) -
                    previous_state.estimation_list.at(active_index))).eval()
        };
        if (!IsLocalFittingActiveEstimationFinitePositive(accelerated_estimation))
        {
            return false;
        }
        accelerated_estimation_list.emplace_back(active_index, accelerated_estimation);
    }

    for (const auto & [active_index, accelerated_estimation] : accelerated_estimation_list)
    {
        const auto accelerated_model{ GaussianModel3D::FromVector(accelerated_estimation) };
        auto & result{ current_state.result_list.at(active_index) };
        result.mdpde = GaussianModel3DWithUncertainty{
            accelerated_model,
            result.mdpde.GetStandardDeviationModel()
        };
        current_state.estimation_list.at(active_index) = accelerated_estimation;
    }
    return true;
}

LocalFittingClusterLayout BuildLocalFittingClusterLayout(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list)
{
    LocalFittingClusterLayout layout;
    layout.cluster_list = BuildLocalFittingAndersonClusters(context, active_index_list);
    layout.key_list.reserve(layout.cluster_list.size());
    for (std::size_t cluster_index = 0; cluster_index < layout.cluster_list.size(); cluster_index++)
    {
        const auto & key{ layout.cluster_list.at(cluster_index).active_index_list };
        layout.key_list.emplace_back(key);
        const auto inserted{ layout.index_by_key.emplace(key, cluster_index) };
        if (!inserted.second)
        {
            throw std::invalid_argument("Local fitting cluster key is duplicated.");
        }
    }
    return layout;
}

template <typename Predicate>
std::vector<LocalFittingAndersonClusterKey> BuildLocalFittingClusterStatusKeyList(
    const LocalFittingClusterLayout & layout,
    const std::vector<LocalFittingClusterAttemptStatus> & status_list,
    Predicate predicate)
{
    if (layout.key_list.size() != status_list.size())
    {
        throw std::invalid_argument("Local fitting cluster status size is inconsistent.");
    }
    std::vector<LocalFittingAndersonClusterKey> key_list;
    for (std::size_t cluster_index = 0; cluster_index < status_list.size(); cluster_index++)
    {
        if (!predicate(status_list.at(cluster_index))) continue;
        key_list.emplace_back(layout.key_list.at(cluster_index));
    }
    return key_list;
}

template <typename Predicate>
bool HasLocalFittingClusterStatus(
    const std::vector<LocalFittingClusterAttemptStatus> & status_list,
    Predicate predicate)
{
    return std::any_of(status_list.begin(), status_list.end(), predicate);
}

algorithm::FittingQualityCandidateStats BuildInitialLocalFittingClusterCandidateStats(
    const std::optional<LocalFittingObjectiveSamples> & objective_samples,
    const algorithm::ScaleReference & objective_reference,
    std::optional<double> objective_scale_sample)
{
    LocalFittingObjectiveStats objective_stats;
    if (objective_samples.has_value() && objective_scale_sample.has_value())
    {
        objective_stats = CalculateLocalFittingObjectiveStats(
            *objective_samples,
            objective_reference,
            *objective_scale_sample);
    }
    return algorithm::FittingQualityCandidateStats{
        objective_stats.has_quality_objective,
        objective_stats.quality_objective,
        algorithm::ParameterChangeStats{
            std::vector<double>(
                static_cast<std::size_t>(GaussianModel3D::ParameterSize()),
                0.0)
        }
    };
}

algorithm::ClusteredFittingQualityInitialState<LocalFittingObjectiveSamples>
BuildInitialLocalFittingClusterQualityState(
    const SecondStageLocalFittingContext & context,
    const GaussianFittingState & previous_state,
    const LocalFittingAndersonCluster & cluster)
{
    auto initial_objective_samples{
        CollectLocalFittingObjectiveSamples(
            context,
            previous_state,
            cluster.objective_sample_ref_list)
    };
    std::optional<double> initial_objective_scale_sample;
    if (initial_objective_samples.has_value())
    {
        initial_objective_scale_sample =
            CalculateLocalFittingResidualScaleSample(*initial_objective_samples);
    }
    algorithm::ScaleReferenceTracker initial_tracker{
        kLocalFittingObjectiveScaleWarmupCount,
        initial_objective_scale_sample
    };
    auto initial_candidate_stats{
        BuildInitialLocalFittingClusterCandidateStats(
            initial_objective_samples,
            initial_tracker.GetCommittedReference(),
            initial_objective_scale_sample)
    };
    return algorithm::ClusteredFittingQualityInitialState<LocalFittingObjectiveSamples>{
        initial_objective_scale_sample,
        std::move(initial_candidate_stats),
        std::move(initial_objective_samples)
    };
}

LocalRefitResult FitAtomWithJointOffsetFallback(
    const SecondStageLocalFittingContext & context,
    std::size_t atom_index,
    const LocalGaussianResult & previous_result,
    const FittedGaussianSnapshot & offset_snapshot,
    const FitOptions & options)
{
    auto sample_entries{ BuildSecondStageAdjustedSamples(context, atom_index, offset_snapshot) };
    const auto & offset_model{ offset_snapshot.at(atom_index) };
    try
    {
        auto candidate_result{
            EstimateLocalGaussian(
                sample_entries,
                context.atom_context_list.at(atom_index).alpha_r,
                options,
                offset_model)
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
    result.ols = GaussianModel3DWithUncertainty{
        result.ols.GetModel().WithOffset(offset_model.GetOffset()),
        result.ols.GetStandardDeviationModel()
    };
    result.mdpde = GaussianModel3DWithUncertainty{
        result.mdpde.GetModel().WithOffset(offset_model.GetOffset()),
        result.mdpde.GetStandardDeviationModel()
    };
    if (IsSuspiciousJointOffset(
            sample_entries, previous_result.mdpde.GetModel(), result.mdpde.GetModel()))
    {
        return LocalRefitResult{ previous_result, true, true };
    }
    return LocalRefitResult{ result, true, false };
}

void ThawChangedActiveAtomNeighbors(
    const SecondStageLocalFittingContext & context,
    const std::vector<algorithm::ParameterChange> & change_list,
    const std::vector<std::size_t> & active_index_list,
    algorithm::ConvergenceFreezeTracker & freeze_tracker,
    algorithm::DependencyThawHysteresisTracker & thaw_hysteresis_tracker)
{
    if (change_list.size() != context.AtomSize())
    {
        throw std::invalid_argument("Local fitting dependency thaw input size is inconsistent.");
    }

    const double thaw_threshold{ std::sqrt(kLocalFittingParameterChangeTolerance) };
    for (const auto active_index : active_index_list)
    {
        if (active_index >= context.AtomSize())
        {
            throw std::invalid_argument("Local fitting dependency thaw active index is out of range.");
        }
        const auto active_change{ algorithm::GetMaximumParameterChange(change_list.at(active_index)) };
        for (const auto neighbor_index : context.atom_context_list.at(active_index).selected_neighbor_index_list)
        {
            if (!freeze_tracker.IsFrozen(neighbor_index)) continue;
            if (!thaw_hysteresis_tracker.ShouldThaw(neighbor_index, active_change, thaw_threshold)) continue;
            if (!thaw_hysteresis_tracker.CanDependencyThaw(neighbor_index)) continue;
            if (freeze_tracker.Thaw(neighbor_index))
            {
                thaw_hysteresis_tracker.RecordDependencyThaw(neighbor_index);
            }
        }
    }
}

void ExpandSuspiciousOffsetClusters(
    const ActiveCouplingGraph & active_coupling_graph,
    const std::vector<int> & suspicious_offset_seed_flag_list,
    std::vector<int> & suspicious_offset_flag_list,
    std::vector<int> & refit_fallback_flag_list)
{
    if (suspicious_offset_seed_flag_list.size() != suspicious_offset_flag_list.size() ||
        refit_fallback_flag_list.size() != suspicious_offset_flag_list.size())
    {
        throw std::invalid_argument("Suspicious offset cluster input sizes are inconsistent.");
    }
    if (!active_coupling_graph.empty() &&
        active_coupling_graph.size() != suspicious_offset_flag_list.size())
    {
        throw std::invalid_argument("Suspicious offset cluster graph size is inconsistent.");
    }

    for (std::size_t seed_index = 0; seed_index < suspicious_offset_seed_flag_list.size(); seed_index++)
    {
        if (suspicious_offset_seed_flag_list.at(seed_index) == 0) continue;

        std::vector<char> visited(suspicious_offset_seed_flag_list.size(), 0);
        std::vector<std::pair<std::size_t, std::size_t>> queue{ std::pair<std::size_t, std::size_t>{ seed_index, 0 } };
        visited.at(seed_index) = 1;
        for (std::size_t queue_index = 0; queue_index < queue.size(); queue_index++)
        {
            const auto [active_index, depth]{ queue.at(queue_index) };
            suspicious_offset_flag_list.at(active_index) = 1;
            refit_fallback_flag_list.at(active_index) = 1;
            if (active_coupling_graph.empty() || depth >= kSuspiciousOffsetClusterMaxDepth) continue;

            for (const auto & edge : active_coupling_graph.at(active_index))
            {
                if (edge.neighbor_index >= active_coupling_graph.size())
                {
                    throw std::invalid_argument("Suspicious offset cluster adjacency index is out of range.");
                }
                if (!std::isfinite(edge.overlap) || edge.overlap < kSuspiciousOffsetClusterMinimumOverlap) continue;
                if (visited.at(edge.neighbor_index) != 0) continue;
                visited.at(edge.neighbor_index) = 1;
                queue.emplace_back(edge.neighbor_index, depth + 1);
            }
        }
    }
}

void RollBackSuspiciousOffsetClusters(
    const SecondStageLocalFittingContext & context,
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
        if (state_index >= context.AtomSize() ||
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
        current_snapshot.at(state_index) = previous_model;
        iteration_state.estimation_list.at(state_index) = previous_state.estimation_list.at(state_index);
        iteration_state.result_list.at(state_index) = previous_state.result_list.at(state_index);
    }
}

void ExpandAndRollBackSuspiciousOffsetClusters(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const GaussianFittingState & previous_state,
    const ActiveCouplingGraph & active_coupling_graph,
    const std::vector<int> & suspicious_offset_seed_flag_list,
    std::vector<int> & suspicious_offset_flag_list,
    std::vector<int> & refit_fallback_flag_list,
    FittedGaussianSnapshot & current_snapshot,
    GaussianFittingState & iteration_state)
{
    ExpandSuspiciousOffsetClusters(
        active_coupling_graph,
        suspicious_offset_seed_flag_list,
        suspicious_offset_flag_list,
        refit_fallback_flag_list);
    RollBackSuspiciousOffsetClusters(
        context,
        active_index_list,
        previous_state,
        suspicious_offset_flag_list,
        current_snapshot,
        iteration_state);
}

LocalFittingIterationResult RunLocalFittingIteration(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const GaussianFittingState & previous_state,
    const FitOptions & options,
    double ridge_ratio,
    const std::vector<double> & ridge_multiplier_list)
{
    const auto selected_atom_size{ context.AtomSize() };
    if (previous_state.result_list.size() != selected_atom_size ||
        previous_state.estimation_list.size() != selected_atom_size ||
        ridge_multiplier_list.size() != selected_atom_size)
    {
        throw std::invalid_argument("Local fitting iteration input sizes are inconsistent.");
    }
    auto current_snapshot{
        BuildFittedGaussianSnapshot(context, previous_state.estimation_list)
    };
    const auto joint_offset_result{
        EstimateJointOffsets(
            context, active_index_list, current_snapshot, ridge_ratio, ridge_multiplier_list)
    };
    ApplyJointOffsetsToSnapshot(active_index_list, joint_offset_result.offset, current_snapshot);
    auto iteration_state{ previous_state };
    std::vector<int> refit_fallback_flag_list(active_index_list.size(), 0);
    std::vector<int> suspicious_offset_flag_list(active_index_list.size(), 0);
    std::vector<int> pre_refit_suspicious_seed_flag_list(active_index_list.size(), 0);

    for (std::size_t i = 0; i < active_index_list.size(); i++)
    {
        const auto state_index{ active_index_list.at(i) };
        const auto previous_model{
            GaussianModel3D::FromVector(previous_state.estimation_list.at(state_index))
        };
        const auto & offset_model{ current_snapshot.at(state_index) };
        if (!IsSuspiciousJointOffset(
                context.atom_context_list.at(state_index).sample_entries,
                previous_model,
                offset_model))
        {
            continue;
        }
        pre_refit_suspicious_seed_flag_list.at(i) = 1;
    }
    ExpandAndRollBackSuspiciousOffsetClusters(
        context,
        active_index_list,
        previous_state,
        joint_offset_result.active_coupling_graph,
        pre_refit_suspicious_seed_flag_list,
        suspicious_offset_flag_list,
        refit_fallback_flag_list,
        current_snapshot,
        iteration_state);

    std::vector<int> post_refit_suspicious_seed_flag_list(active_index_list.size(), 0);
    for (size_t i = 0; i < active_index_list.size(); i++)
    {
        if (suspicious_offset_flag_list.at(i) != 0) continue;

        const auto state_index{ active_index_list.at(i) };
        auto refit_result{
            FitAtomWithJointOffsetFallback(
                context,
                state_index,
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
            post_refit_suspicious_seed_flag_list.at(i) = 1;
        }
        auto result{ std::move(refit_result.result) };
        const auto fitted_model{ result.mdpde.GetModel() };
        iteration_state.estimation_list.at(state_index) = fitted_model.ToVector();
        iteration_state.result_list.at(state_index) = std::move(result);
    }
    ExpandAndRollBackSuspiciousOffsetClusters(
        context,
        active_index_list,
        previous_state,
        joint_offset_result.active_coupling_graph,
        post_refit_suspicious_seed_flag_list,
        suspicious_offset_flag_list,
        refit_fallback_flag_list,
        current_snapshot,
        iteration_state);

    LocalFittingIterationResult iteration_result;
    iteration_result.state = std::move(iteration_state);
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
        throw std::invalid_argument("local_editor_list and local fitting state sizes are inconsistent.");
    }

    for (std::size_t i = 0; i < local_editor_list.size(); i++)
    {
        local_editor_list.at(i).SetGaussianResult(iteration_state.result_list.at(i));
    }
}

void LogLocalFittingAllAtomsFrozen(const FitOptions & options, std::size_t accepted_iteration_count)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info,
        "Converged after " + std::to_string(accepted_iteration_count) +
        " iterations because all local atoms are frozen.");
}

void LogLocalFittingAndersonFallbackSwitch(
    const FitOptions & options,
    std::size_t iteration,
    bool has_invalid_anderson_candidate)
{
    if (options.quiet_mode) return;

    std::ostringstream progress_message;
    progress_message << "Local fitting iteration " << iteration + 1 << '/'
        << kLocalFittingMaximumIterations
        << " switching from Anderson acceleration to damped fixed-point fallback"
        << ", next acceleration = "
        << GetLocalFittingAccelerationText(LocalFittingAccelerationKind::DampedFixedPoint)
        << (has_invalid_anderson_candidate ? " after invalid Anderson candidate" : "");
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info, progress_message.str());
}

void LogLocalFittingBacktrackingRetry(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    double ridge_ratio,
    bool uses_cluster_local_objective_ridge)
{
    if (options.quiet_mode) return;

    std::ostringstream progress_message;
    progress_message
        << "Objective backtracking rejected all attempts; retrying after local fitting iteration "
        << accepted_iteration_count
        << std::fixed << std::setprecision(5)
        << "; acceleration history reset";
    if (uses_cluster_local_objective_ridge)
    {
        progress_message
            << ", next attempt uses increased cluster-local objective ridge"
            << ", global ridge ratio remains = " << ridge_ratio;
    }
    else
    {
        progress_message
            << ", next attempt uses increased global ridge ratio = " << ridge_ratio;
    }
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info, progress_message.str());
}

void LogLocalFittingBacktrackingStop(
    const FitOptions & options,
    double ridge_ratio)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message
        << "Stopped local fitting because objective backtracking rejected all "
        << "acceleration and fixed-point attempts "
        << (ridge_ratio >= kJointOffsetRidgeRatioMax ?
            "at the maximum joint-offset ridge ratio" : "at the maximum iteration limit")
        << "; applying previous state.";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogLocalFittingProgress(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const LocalFittingAccelerationAttempt & current_acceleration_attempt,
    const algorithm::ConvergenceFreezeTracker & freeze_tracker)
{
    if (options.quiet_mode) return;

    std::ostringstream progress_message;
    progress_message << "Iter. " << accepted_iteration_count
        << '/' << kLocalFittingMaximumIterations
        << std::fixed << std::setprecision(4)
        << ", acceleration = "<< GetLocalFittingAccelerationText(current_acceleration_attempt.kind)
        << ", damping = "<< current_acceleration_attempt.damping
        << ", active/frozen atoms = "<< freeze_tracker.GetActiveCount()
        << "/" << freeze_tracker.GetFrozenCount();
    Logger::ProgressLine(progress_message.str());
}

void LogLocalFittingConverged(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const algorithm::ParameterChangeStats & normalized_change_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info,
        "Converged after " + std::to_string(accepted_iteration_count) +
        " iterations with normalized percentile amplitude change = " +
        std::to_string(normalized_change_stats.percentile_list.at(GaussianModel3D::AmplitudeIndex())) +
        ", normalized percentile width change = " +
        std::to_string(normalized_change_stats.percentile_list.at(GaussianModel3D::WidthIndex())) +
        ", and normalized percentile offset change = " +
        std::to_string(normalized_change_stats.percentile_list.at(GaussianModel3D::OffsetIndex())) +
        ".");
}

void LogLocalFittingMaximumIterations(
    const FitOptions & options,
    const algorithm::ParameterChangeStats & normalized_change_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Warning,
        "Reached maximum iteration size; refitting at current accepted candidate "
        "with normalized percentile amplitude change = " +
        std::to_string(
            normalized_change_stats.percentile_list.at(GaussianModel3D::AmplitudeIndex())) +
        ", normalized percentile width change = " +
        std::to_string(
            normalized_change_stats.percentile_list.at(GaussianModel3D::WidthIndex())) +
        ", and normalized percentile offset change = " +
        std::to_string(
            normalized_change_stats.percentile_list.at(GaussianModel3D::OffsetIndex())));
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

void RunFixedOffsetLocalFittingPass(
    ModelObject & model_object,
    const FitOptions & options,
    LocalFittingPass pass)
{
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    const auto selected_atom_size{ atom_list.size() };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, atom_list) };
    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto stage_label{
        pass == LocalFittingPass::FirstStage ? "1st-stage" : "3rd-stage"
    };
    std::atomic<size_t> atom_count{ 0 };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info,
            "Run " + std::string{ stage_label } + " local atom fitting for " +
            std::to_string(selected_atom_size) + " atoms.");
    }

#ifdef USE_OPENMP
    #pragma omp parallel for num_threads(options.thread_size)
#endif
    for (size_t i = 0; i < selected_atom_size; i++)
    {
        auto & atom{ *atom_list[i] };
        const auto local_view{ AtomLocalPotentialView::RequireFor(atom) };
        LocalPotentialSampleList sample_entries;
        GaussianModel3D offset_model;
        switch (pass)
        {
        case LocalFittingPass::FirstStage:
            sample_entries = local_view.GetSamplingEntries();
            offset_model = local_view.GetGaussianResult().mdpde.GetModel();
            break;
        case LocalFittingPass::ThirdStage:
            sample_entries = local_view.GetSamplingEntries(false, true);
            offset_model = analysis_view.GetAtomGroupPrior(data_internal::GetGroupKey(&atom));
            break;
        }

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

} // namespace

double TrainAlphaR(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    const FitOptions & options)
{
    numeric_validation::RequireFiniteNonNegativeRange(
        options.distance_min, options.distance_max, "fit range");

    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (const auto & sample_entries : sample_entries_list)
    {
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(
                sample_entries,
                options.distance_min,
                options.distance_max));
    }
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
    const auto range_min{ options.distance_min };
    const auto range_max{ options.distance_max };
    std::vector<RHBMMemberDataset> dataset_list;
    dataset_list.reserve(sample_entries_list.size());
    for (std::size_t i = 0; i < sample_entries_list.size(); i++)
    {
        const auto sampling_entries{
            BuildSamplesForZeroOffsetGaussianFit(
                sample_entries_list.at(i),
                member_result_list.at(i).mdpde.GetModel())
        };
        dataset_list.emplace_back(
            rhbm_helper::BuildMemberDataset(sampling_entries, range_min, range_max));
    }
    if (dataset_list.size() != member_result_list.size())
    {
        throw std::invalid_argument("dataset_list and member_result_list sizes are inconsistent.");
    }

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

void RunLocalAlphaTraining(
    ModelObject & model_object,
    const FitOptions & options,
    bool apply_selection,
    bool use_updated_sample)
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
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(group_key) };
        std::vector<LocalPotentialSampleList> sample_entries_list;
        sample_entries_list.reserve(group_atom_list.size());
        for (auto * atom : group_atom_list)
        {
            analysis.EnsureAtomLocalPotential(*atom);
            const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
            const auto & sample_entries{ local_view.GetSamplingEntries(apply_selection, use_updated_sample) };
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
        const auto & group_atom_list{ analysis_view.GetAtomObjectList(group_key) };
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
    RunFixedOffsetLocalFittingPass(model_object, options, LocalFittingPass::FirstStage);
}

void RunSecondStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    const auto context{ BuildSecondStageLocalFittingContext(model_object) };
    auto local_editor_list{ BuildAtomLocalEditors(model_object, context.atom_list) };
    const auto atom_size{ context.AtomSize() };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run 2nd-stage local atom fitting with iterations...");
    }

    auto previous_state{ BuildInitialLocalFittingState(context) };
    algorithm::ClusteredAndersonAccelerationHistorySet acceleration_history{
        algorithm::AndersonAccelerationOptions{
            kLocalFittingAndersonHistoryDepth,
            kLocalFittingNormalizedChangeScaleFloor,
            kLocalFittingAndersonCoefficientL1Limit,
            kLocalFittingAndersonRegularization
        }
    };
    algorithm::ConvergenceFreezeTracker freeze_tracker{
        atom_size,
        kLocalFittingParameterChangeTolerance,
        kLocalFittingFreezeChangeRatio,
        kLocalFittingFreezeStableIterations
    };
    algorithm::DependencyThawHysteresisTracker thaw_hysteresis_tracker{
        atom_size,
        kLocalFittingDependencyThawHysteresisGrowth,
        kLocalFittingDependencyThawHysteresisMax,
        kLocalFittingDependencyThawHysteresisFrozenDecay,
        kLocalFittingDependencyThawMaximumCount
    };
    double ridge_ratio{ kJointOffsetRidgeRatio };
    std::vector<double> suspicious_joint_offset_ridge_multiplier_list(atom_size, 1.0);
    algorithm::ClusteredFittingQualityStateSet<LocalFittingObjectiveSamples> cluster_quality_state{
        algorithm::ClusteredFittingQualityOptions{
            kLocalFittingObjectiveScaleWarmupCount,
            kLocalFittingConvergenceObjectiveRelativeTolerance,
            kLocalFittingObjectiveTieRelativeTolerance,
            1.0,
            kSuspiciousJointOffsetRidgeMultiplier,
            kJointOffsetRidgeGrowth,
            kJointOffsetRidgeShrink
        }
    };
    std::size_t accepted_iteration_count{ 0 };
    for (size_t iter = 0; iter < kLocalFittingMaximumIterations; iter++)
    {
        const auto active_index_list{ freeze_tracker.BuildActiveIndexList() };
        if (active_index_list.empty())
        {
            ApplyLocalFittingState(previous_state, local_editor_list);
            LogLocalFittingAllAtomsFrozen(options, accepted_iteration_count);
            break;
        }

        const auto cluster_layout{
            BuildLocalFittingClusterLayout(context, active_index_list)
        };
        cluster_quality_state.Reconcile(
            cluster_layout.key_list,
            [&](const LocalFittingAndersonClusterKey & key)
            {
                const auto cluster_index{ cluster_layout.index_by_key.at(key) };
                return BuildInitialLocalFittingClusterQualityState(
                    context,
                    previous_state,
                    cluster_layout.cluster_list.at(cluster_index));
            });
        acceleration_history.Reconcile(cluster_layout.key_list);
        const auto joint_offset_ridge_multiplier_list{
            CombineLocalFittingRidgeMultiplierLists(
                suspicious_joint_offset_ridge_multiplier_list,
                cluster_quality_state.BuildObjectiveRidgeMultiplierList(atom_size))
        };
        auto iteration_result{
            RunLocalFittingIteration(
                context,
                active_index_list,
                previous_state,
                options,
                ridge_ratio,
                joint_offset_ridge_multiplier_list)
        };
        const auto has_suspicious_offset_fallback{ !iteration_result.suspicious_offset_state_index_list.empty() };
        std::vector<bool> suspicious_offset_atom_seen(atom_size, false);
        for (const auto state_index : iteration_result.suspicious_offset_state_index_list)
        {
            if (state_index >= suspicious_offset_atom_seen.size())
            {
                throw std::invalid_argument("Local fitting suspicious offset atom index is out of range.");
            }
            suspicious_offset_atom_seen.at(state_index) = true;
            suspicious_joint_offset_ridge_multiplier_list.at(state_index) =
                kSuspiciousJointOffsetRidgeMultiplier;
        }
        const auto raw_state{ std::move(iteration_result.state) };
        bool has_current_candidate{ false };
        bool has_objective_backtracking_rejection{ false };
        bool has_invalid_anderson_candidate{ false };
        const auto localized_anderson_candidate{
            acceleration_history.BuildCandidate(
                cluster_layout.key_list,
                previous_state.estimation_list,
                raw_state.estimation_list)
        };
        std::optional<std::vector<Eigen::VectorXd>> anderson_candidate_estimation_list;
        if (localized_anderson_candidate.has_value())
        {
            anderson_candidate_estimation_list = localized_anderson_candidate->state_list;
        }

        auto assembled_state{ previous_state };
        std::vector<LocalFittingClusterAttemptStatus> cluster_status_list(
            cluster_layout.cluster_list.size());
        std::vector<
            algorithm::ClusteredFittingQualityAcceptedScore<LocalFittingObjectiveSamples>>
            accepted_cluster_score_list;
        LocalFittingAccelerationAttempt accepted_acceleration_attempt;
        bool has_accepted_acceleration_attempt{ false };

        const auto copy_cluster_state = [](
            GaussianFittingState & target_state,
            const GaussianFittingState & source_state,
            const LocalFittingAndersonClusterKey & key)
        {
            for (const auto active_id : key)
            {
                target_state.result_list.at(active_id) = source_state.result_list.at(active_id);
                target_state.estimation_list.at(active_id) = source_state.estimation_list.at(active_id);
            }
        };

        const auto run_attempt_group = [&](bool use_anderson)
        {
            const auto attempt_size{ static_cast<int>(kLocalFittingAccelerationDampingList.size()) };
            for (auto & status : cluster_status_list)
            {
                if (!status.accepted)
                {
                    status.stopped = false;
                }
            }
            for (int attempt = 0; attempt < attempt_size; attempt++)
            {
                const auto acceleration_attempt{
                    BuildLocalFittingAccelerationAttempt(
                        use_anderson,
                        kLocalFittingAccelerationDampingList.at(static_cast<std::size_t>(attempt)))
                };
                for (std::size_t cluster_index = 0; cluster_index < cluster_layout.cluster_list.size(); cluster_index++)
                {
                    auto & status{ cluster_status_list.at(cluster_index) };
                    if (status.accepted || status.stopped)
                    {
                        continue;
                    }

                    const auto & cluster{ cluster_layout.cluster_list.at(cluster_index) };
                    const auto & key{ cluster.active_index_list };
                    const auto & key_list{ localized_anderson_candidate->used_cluster_key_list };
                    if (use_anderson &&
                        (!localized_anderson_candidate.has_value() ||
                         !(std::find(key_list.begin(), key_list.end(), key) != key_list.end())))
                    {
                        continue;
                    }

                    auto attempt_state{ previous_state };
                    bool valid_attempt{ true };
                    if (acceleration_attempt.kind == LocalFittingAccelerationKind::DampedFixedPoint)
                    {
                        ApplyLocalFittingDampedFixedPoint(
                            attempt_state,
                            raw_state,
                            previous_state,
                            key,
                            acceleration_attempt.damping);
                    }
                    else
                    {
                        valid_attempt =
                            anderson_candidate_estimation_list.has_value() &&
                            TryApplyLocalFittingAndersonCandidate(
                                attempt_state,
                                previous_state,
                                *anderson_candidate_estimation_list,
                                key,
                                acceleration_attempt.damping);
                    }
                    if (!valid_attempt)
                    {
                        has_invalid_anderson_candidate = true;
                        status.rejected_anderson = true;
                        continue;
                    }

                    const auto normalized_change_list{
                        CalculateLocalFittingNormalizedParameterChanges(
                            attempt_state.estimation_list,
                            previous_state.estimation_list)
                    };
                    auto normalized_change_stats{
                        algorithm::SummarizeParameterChangeStats(
                            normalized_change_list,
                            key,
                            kLocalFittingChangePercentile)
                    };
                    auto evaluated_attempt{
                        cluster_quality_state.EvaluateCandidate(
                            key,
                            attempt,
                            attempt_size,
                            [&](const algorithm::ScaleReferenceTracker & objective_scale_tracker,
                                algorithm::FittingQualityCandidateStats & previous_candidate_stats,
                                const std::optional<LocalFittingObjectiveSamples> & previous_objective_samples,
                                const std::optional<
                                    algorithm::ClusteredFittingQualityTrackedCandidate<LocalFittingObjectiveSamples>> & best_candidate)
                            {
                                return ScoreLocalFittingClusterCandidate(
                                    context,
                                    attempt_state,
                                    cluster.objective_sample_ref_list,
                                    objective_scale_tracker,
                                    previous_candidate_stats,
                                    previous_objective_samples,
                                    best_candidate,
                                    normalized_change_stats);
                            })
                    };
                    if (evaluated_attempt.outcome == algorithm::ClusteredFittingQualityAttemptOutcome::Accepted)
                    {
                        copy_cluster_state(assembled_state, attempt_state, key);
                        status.accepted = true;
                        status.rejected = false;
                        accepted_cluster_score_list.emplace_back(std::move(evaluated_attempt.accepted_score));
                        if (use_anderson)
                        {
                            status.accepted_by_anderson = true;
                        }
                        if (!use_anderson)
                        {
                            status.rejected_anderson = false;
                        }
                        if (!has_accepted_acceleration_attempt)
                        {
                            accepted_acceleration_attempt = acceleration_attempt;
                            has_accepted_acceleration_attempt = true;
                        }
                        continue;
                    }

                    has_objective_backtracking_rejection = true;
                    status.rejected = true;
                    if (use_anderson)
                    {
                        status.rejected_anderson = true;
                    }
                    if (evaluated_attempt.outcome == algorithm::ClusteredFittingQualityAttemptOutcome::ObjectiveStop)
                    {
                        status.stopped = true;
                    }
                }
            }
        };

        if (anderson_candidate_estimation_list.has_value())
        {
            run_attempt_group(true);
        }

        const auto has_pending_cluster = [&]()
        {
            return HasLocalFittingClusterStatus(
                cluster_status_list,
                [](const LocalFittingClusterAttemptStatus & status)
                {
                    return !status.accepted;
                });
        };

        if (has_pending_cluster())
        {
            const auto rejected_anderson_cluster_key_list{
                BuildLocalFittingClusterStatusKeyList(
                    cluster_layout,
                    cluster_status_list,
                    [](const LocalFittingClusterAttemptStatus & status)
                    {
                        return status.rejected_anderson;
                    })
            };
            if (localized_anderson_candidate.has_value() && !rejected_anderson_cluster_key_list.empty())
            {
                acceleration_history.ClearAndSuppress(rejected_anderson_cluster_key_list);
                LogLocalFittingAndersonFallbackSwitch(options, iter, has_invalid_anderson_candidate);
            }
            run_attempt_group(false);
        }

        has_current_candidate =
            HasLocalFittingClusterStatus(
                cluster_status_list,
                [](const LocalFittingClusterAttemptStatus & status)
                {
                    return status.accepted;
                });

        const auto accepted_cluster_key_list{
            BuildLocalFittingClusterStatusKeyList(
                cluster_layout,
                cluster_status_list,
                [](const LocalFittingClusterAttemptStatus & status)
                {
                    return status.accepted;
                })
        };
        const auto rejected_cluster_key_list{
            BuildLocalFittingClusterStatusKeyList(
                cluster_layout,
                cluster_status_list,
                [](const LocalFittingClusterAttemptStatus & status)
                {
                    return status.rejected;
                })
        };
        const auto rejected_anderson_cluster_key_list{
            BuildLocalFittingClusterStatusKeyList(
                cluster_layout,
                cluster_status_list,
                [](const LocalFittingClusterAttemptStatus & status)
                {
                    return status.rejected_anderson;
                })
        };

        if (!has_current_candidate)
        {
            bool increased_cluster_objective_ridge{ false };
            bool increased_global_ridge_ratio{ false };
            if (!rejected_anderson_cluster_key_list.empty())
            {
                acceleration_history.ClearAndSuppress(rejected_anderson_cluster_key_list);
            }
            if (!rejected_cluster_key_list.empty())
            {
                increased_cluster_objective_ridge = cluster_quality_state.IncreaseObjectiveRidge(rejected_cluster_key_list);
            }
            if (has_objective_backtracking_rejection && !increased_cluster_objective_ridge)
            {
                const auto previous_ridge_ratio{ ridge_ratio };
                ridge_ratio = std::min(kJointOffsetRidgeRatioMax, ridge_ratio * kJointOffsetRidgeGrowth);
                increased_global_ridge_ratio = ridge_ratio > previous_ridge_ratio;
            }
            if ((increased_cluster_objective_ridge ||
                (increased_global_ridge_ratio && ridge_ratio < kJointOffsetRidgeRatioMax)) &&
                iter + 1 < kLocalFittingMaximumIterations)
            {
                LogLocalFittingBacktrackingRetry(
                    options, accepted_iteration_count, ridge_ratio,
                    increased_cluster_objective_ridge);
                continue;
            }

            ApplyLocalFittingState(previous_state, local_editor_list);
            LogLocalFittingBacktrackingStop(options, ridge_ratio);
            break;
        }

        auto change_list{
            CalculateLocalFittingParameterChanges(
                assembled_state.estimation_list,
                previous_state.estimation_list)
        };
        const auto normalized_change_list{
            CalculateLocalFittingNormalizedParameterChanges(
                assembled_state.estimation_list,
                previous_state.estimation_list)
        };
        auto normalized_change_stats{
            algorithm::SummarizeParameterChangeStats(
                normalized_change_list,
                active_index_list,
                kLocalFittingChangePercentile)
        };
        accepted_iteration_count++;
        cluster_quality_state.CommitAccepted(accepted_cluster_score_list);
        cluster_quality_state.DecreaseObjectiveRidge(accepted_cluster_key_list);
        if (!rejected_cluster_key_list.empty())
        {
            cluster_quality_state.IncreaseObjectiveRidge(rejected_cluster_key_list);
        }

        const auto accepted_fixed_point_cluster_key_list{
            BuildLocalFittingClusterStatusKeyList(
                cluster_layout,
                cluster_status_list,
                [](const LocalFittingClusterAttemptStatus & status)
                {
                    return status.accepted && !status.accepted_by_anderson;
                })
        };
        acceleration_history.ClearAndSuppress(rejected_anderson_cluster_key_list);
        acceleration_history.ReleaseSuppression(accepted_fixed_point_cluster_key_list);

        for (const auto active_index : active_index_list)
        {
            if (!suspicious_offset_atom_seen.at(active_index))
            {
                suspicious_joint_offset_ridge_multiplier_list.at(active_index) = 1.0;
            }
        }
        if (!has_objective_backtracking_rejection)
        {
            ridge_ratio = std::max(kJointOffsetRidgeRatioMin, ridge_ratio * kJointOffsetRidgeShrink);
        }
        if (has_suspicious_offset_fallback)
        {
            acceleration_history.ClearAndSuppressContaining(iteration_result.suspicious_offset_state_index_list);
        }
        acceleration_history.Commit(
            accepted_cluster_key_list,
            previous_state.estimation_list,
            raw_state.estimation_list);
        std::vector<std::size_t> accepted_active_index_list;
        for (const auto & key : accepted_cluster_key_list)
        {
            accepted_active_index_list.insert(
                accepted_active_index_list.end(), key.begin(), key.end());
        }
        std::sort(accepted_active_index_list.begin(), accepted_active_index_list.end());
        accepted_active_index_list.erase(
            std::unique(accepted_active_index_list.begin(), accepted_active_index_list.end()),
            accepted_active_index_list.end());
        freeze_tracker.Update(change_list, accepted_active_index_list);
        ThawChangedActiveAtomNeighbors(
            context, change_list, accepted_active_index_list,
            freeze_tracker, thaw_hysteresis_tracker);
        for (const auto state_index : iteration_result.suspicious_offset_state_index_list)
        {
            freeze_tracker.Thaw(state_index);
        }
        for (std::size_t state_index = 0; state_index < atom_size; state_index++)
        {
            if (freeze_tracker.IsFrozen(state_index))
            {
                thaw_hysteresis_tracker.DecayFrozen(state_index);
            }
        }

        LogLocalFittingProgress(
            options, accepted_iteration_count, accepted_acceleration_attempt, freeze_tracker);

        if (freeze_tracker.GetActiveCount() == 0)
        {
            ApplyLocalFittingState(assembled_state, local_editor_list);
            LogLocalFittingAllAtomsFrozen(options, accepted_iteration_count);
            break;
        }

        const auto converged{
            !has_suspicious_offset_fallback &&
            !has_objective_backtracking_rejection &&
            cluster_quality_state.AllActiveReferencesLocked(cluster_layout.key_list) &&
            IsLocalFittingNormalizedParameterChangeConverged(normalized_change_stats)
        };
        if (converged)
        {
            ApplyLocalFittingState(assembled_state, local_editor_list);
            LogLocalFittingConverged(options, accepted_iteration_count, normalized_change_stats);
            break;
        }

        if (iter + 1 == kLocalFittingMaximumIterations)
        {
            ApplyLocalFittingState(assembled_state, local_editor_list);
            LogLocalFittingMaximumIterations(options, normalized_change_stats);
        }
        previous_state = std::move(assembled_state);
    }
}

void RunThirdStageLocalFitting(ModelObject & model_object, const FitOptions & options)
{
    RunFixedOffsetLocalFittingPass(model_object, options, LocalFittingPass::ThirdStage);
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
            sample_entries_list.emplace_back(local_view.GetSamplingEntries(false, true));
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
    RunLocalAlphaTraining(model_object, options, true, false);

    InitializeLocalFittingSeedModels(model_object);
    RunFirstStageLocalFitting(model_object, options);
    RunSecondStageLocalFitting(model_object, options);
    RunGroupAlphaTraining(model_object, options);
    SetUpdatedSamplingEntriesFromGroupMedianGaussian(model_object);
    RunGroupPotentialFitting(model_object, options);
    SetUpdatedSamplingEntriesFromFittedGroupGaussian(model_object);
    RunLocalAlphaTraining(model_object, options, false, true);
    RunThirdStageLocalFitting(model_object, options);

    RunGroupAlphaTraining(model_object, options);
    RunGroupPotentialFitting(model_object, options);
    if (!options.quiet_mode)
    {
        LogGroupPriorSpotSummary(model_object);
    }
}

} // namespace rhbm_gem::core
