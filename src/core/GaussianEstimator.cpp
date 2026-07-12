#include <cstddef>
#include <rhbm_gem/core/GaussianEstimator.hpp>

#include "core/detail/GaussianEstimatorStages.hpp"
#include "core/detail/LocalFittingAudit.hpp"
#include "core/detail/LocalFittingAndersonRegime.hpp"
#include "core/detail/LocalFittingHealth.hpp"
#include "core/detail/LocalFittingSeedRepair.hpp"
#include "core/detail/LocalFittingTransformedChange.hpp"
#include "core/detail/PostRefitRollback.hpp"
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
#include <rhbm_gem/utils/algorithm/NormalizedChange.hpp>
#include <rhbm_gem/utils/algorithm/ParameterChangeStats.hpp>
#include <rhbm_gem/utils/algorithm/RobustLoss.hpp>
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
constexpr double kLocalFittingTransformedChangeTolerance{ 1.0e-4 };
constexpr double kLocalFittingAndersonScaleFloor{ 1.0 };
constexpr std::size_t kMinimumAlphaRTrainingSampleCount{ 10 };
constexpr std::size_t kMinimumAlphaGTrainingMemberCount{ 10 };
constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr double kNeighborAtomSearchRange{ 2.0 * kNeighborContributionDistanceMax };
constexpr std::size_t kLocalFittingMaximumIterations{ 200 };
constexpr double kLocalFittingFreezeTrackerChangeTolerance{ 1.0e-6 };
constexpr double kLocalFittingChangePercentile{ 0.99 };
constexpr std::array<Spot, 4> kGroupPriorSummarySpotList{
    Spot::C,
    Spot::CA,
    Spot::N,
    Spot::O
};
constexpr algorithm::RobustLossKind kSecondStageRobustLossKind{ algorithm::RobustLossKind::Cauchy };
constexpr int kRobustLossMaximumIterations{ 50 };
constexpr double kRobustScaleMultiplier{ 1.4826 };
constexpr double kRobustScaleMin{ 1.0e-12 };
constexpr double kRobustLossCutoffMultiplier{ 1.345 };
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
constexpr double kLocalFittingAndersonRegularization{ 1.0e-4 };
constexpr double kLocalFittingAndersonCoefficientAbsLimit{ 5.0 };
constexpr std::array<double, 5> kLocalFittingAccelerationDampingList{ 1.0, 0.5, 0.25, 0.125, 0.0625 };
constexpr double kLocalFittingFreezeChangeRatio{ 0.1 };
constexpr int kLocalFittingFreezeStableIterations{ 3 };
constexpr double kLocalFittingDependencyThawChangeThreshold{ 1.0e-3 };
constexpr double kLocalFittingDependencyThawHysteresisGrowth{ 2.0 };
constexpr double kLocalFittingDependencyThawHysteresisMax{ 8.0 };
constexpr double kLocalFittingDependencyThawHysteresisFrozenDecay{ 0.9 };
constexpr int kLocalFittingDependencyThawMaximumCount{ 5 };
constexpr double kLocalFittingObjectiveTieRelativeTolerance{ 1.0e-8 };
constexpr double kLocalFittingConvergenceObjectiveRelativeTolerance{ 1.0e-3 };
constexpr std::size_t kLocalFittingObjectiveScaleWarmupCount{ 5 };
constexpr double kLocalFittingObjectiveResidualScaleFloorRatio{ 1.0e-6 };
constexpr double kLocalFittingWidthPriorPenaltyWeight{ 1.0e-2 };
constexpr double kLocalFittingOffsetPlausibilityPenaltyWeight{ 1.0e-2 };
constexpr double kLocalFittingMovementPenaltyWeight{ 1.0e-2 };
constexpr double kLocalFittingWidthPriorLogScale{ 0.35 };
constexpr double kLocalFittingMovementAmplitudeLogScale{ 0.50 };
constexpr double kLocalFittingMovementWidthLogScale{ 0.35 };
constexpr double kLocalFittingOffsetPeakRatioMax{ 1.0 };
constexpr std::size_t kSuspiciousOffsetClusterMaxDepth{ 2 };
constexpr double kSuspiciousOffsetClusterMinimumOverlap{ 0.05 };
constexpr std::size_t kSuspiciousProfileMinimumRadiusCount{ 3 };
constexpr double kSuspiciousProfileDistanceTolerance{ 1.0e-6 };
constexpr double kSuspiciousProfileCenterSignFlipRatio{ 0.25 };
constexpr double kSuspiciousProfileReboundCenterRatio{ 1.5 };
constexpr double kSuspiciousProfileReboundReferenceRatio{ 0.25 };
constexpr double kSuspiciousProfileUpwardExcursionReferenceRatio{ 0.20 };
constexpr int kSuspiciousProfileMaximumUpwardExcursions{ 1 };
constexpr double kSuspiciousWidthGrowthLimit{ 1.5 };
constexpr double kSuspiciousWidthRangeLimitRatio{ 1.5 };
constexpr double kSuspiciousCompensationResponseRatio{ 2.0 };
constexpr std::size_t kPersistentSuspiciousRollbackIterationLimit{ 5 };
constexpr std::size_t kPersistentJointOffsetFailureIterationLimit{ 5 };

struct GaussianFittingState
{
    std::vector<LocalGaussianResult> result_list{};
    std::vector<Eigen::VectorXd> estimation_list{};
};

enum class LocalFittingCandidateKind
{
    Anderson,
    FixedPoint
};

enum class LocalFittingBacktrackingStopReason
{
    MaximumGlobalRidge,
    MaximumIterationLimit
};

struct LocalFittingCandidateAttempt
{
    LocalFittingCandidateKind kind{ LocalFittingCandidateKind::FixedPoint };
    double damping{ 1.0 };
};

struct ActiveCouplingEdge
{
    std::size_t neighbor_index{ 0 };
    double overlap{ 0.0 };
};

using ActiveCouplingGraph = std::vector<std::vector<ActiveCouplingEdge>>;
using LocalFittingClusterKey = algorithm::ClusterKey;
using detail::JointOffsetSolveStatus;
using detail::IsJointOffsetSolveHardFailure;
using detail::IsJointOffsetSolveProgressEligible;
using detail::IsJointOffsetSolveStationarityEligible;

const char * GetJointOffsetSolveStatusText(JointOffsetSolveStatus status)
{
    switch (status)
    {
    case JointOffsetSolveStatus::Converged:
        return "converged";
    case JointOffsetSolveStatus::SystemBuildFailed:
        return "system-build-failed";
    case JointOffsetSolveStatus::EmptySystem:
        return "empty-system";
    case JointOffsetSolveStatus::InitialSolveFailed:
        return "initial-solve-failed";
    case JointOffsetSolveStatus::IrlsSolveFailed:
        return "irls-solve-failed";
    case JointOffsetSolveStatus::IrlsObjectiveDeteriorated:
        return "irls-objective-deteriorated";
    case JointOffsetSolveStatus::IrlsMaximumIterationsReached:
        return "irls-maximum-iterations-reached";
    }
    throw std::logic_error("Joint offset solve status is invalid.");
}

struct ZeroOffsetProfileDiagnostics
{
    double distance_min{ 0.0 };
    double distance_max{ 0.0 };
    double center_response{ 0.0 };
    double max_abs_response{ 0.0 };
    std::vector<double> radius_response_median_list{};
};

struct JointOffsetSolveResult
{
    JointOffsetSolveStatus status{ JointOffsetSolveStatus::SystemBuildFailed };
    Eigen::VectorXd offset{};
    ActiveCouplingGraph active_coupling_graph{};
    std::vector<double> effective_ridge_multiplier_list{};
};

struct JointOffsetBuildResult
{
    algorithm::WeightedRidgeSystem system{};
    ActiveCouplingGraph active_coupling_graph{};
    std::vector<double> effective_ridge_multiplier_list{};
};

struct LocalFittingClusterHealth
{
    JointOffsetSolveStatus joint_offset_status{ JointOffsetSolveStatus::SystemBuildFailed };
    std::size_t unhealthy_refit_atom_count{ 0 };
};

using LocalFittingClusterHealthMap =
    std::map<LocalFittingClusterKey, LocalFittingClusterHealth>;

struct LocalFittingClusterHealthSummary
{
    std::size_t unhealthy_cluster_count{ 0 };
    std::size_t unhealthy_atom_count{ 0 };
    std::map<JointOffsetSolveStatus, std::size_t> unhealthy_joint_status_count{};
    std::map<JointOffsetSolveStatus, std::size_t> joint_status_cluster_count{};
    std::map<JointOffsetSolveStatus, std::size_t> joint_status_atom_count{};
    std::size_t unhealthy_refit_cluster_count{ 0 };
    std::size_t unhealthy_refit_atom_count{ 0 };
};

bool IsLocalFittingClusterHealthy(const LocalFittingClusterHealth & health)
{
    return IsJointOffsetSolveStationarityEligible(health.joint_offset_status) &&
        health.unhealthy_refit_atom_count == 0;
}

void RecordUnhealthyLocalRefit(
    LocalFittingClusterHealthMap & health_by_key,
    std::size_t atom_index)
{
    for (auto & [key, health] : health_by_key)
    {
        if (!std::binary_search(key.begin(), key.end(), atom_index)) continue;
        health.unhealthy_refit_atom_count++;
        return;
    }
    throw std::invalid_argument("Unhealthy local refit atom is not in an active cluster.");
}

std::vector<LocalFittingClusterKey> CollectUnhealthyLocalFittingClusterKeys(
    const LocalFittingClusterHealthMap & health_by_key)
{
    std::vector<LocalFittingClusterKey> key_list;
    for (const auto & [key, health] : health_by_key)
    {
        if (!IsLocalFittingClusterHealthy(health))
        {
            key_list.emplace_back(key);
        }
    }
    return key_list;
}

std::vector<LocalFittingClusterKey> CollectProgressIneligibleLocalFittingClusterKeys(
    const LocalFittingClusterHealthMap & health_by_key)
{
    std::vector<LocalFittingClusterKey> key_list;
    for (const auto & [key, health] : health_by_key)
    {
        if (!IsJointOffsetSolveProgressEligible(health.joint_offset_status) ||
            health.unhealthy_refit_atom_count > 0)
        {
            key_list.emplace_back(key);
        }
    }
    return key_list;
}

std::vector<std::size_t> CollectLocalFittingClusterAtomIndexes(
    const std::vector<LocalFittingClusterKey> & key_list)
{
    std::vector<std::size_t> atom_index_list;
    for (const auto & key : key_list)
    {
        atom_index_list.insert(atom_index_list.end(), key.begin(), key.end());
    }
    return atom_index_list;
}

LocalFittingClusterHealthSummary SummarizeLocalFittingClusterHealth(
    const LocalFittingClusterHealthMap & health_by_key)
{
    LocalFittingClusterHealthSummary summary;
    for (const auto & [key, health] : health_by_key)
    {
        summary.joint_status_cluster_count[health.joint_offset_status]++;
        summary.joint_status_atom_count[health.joint_offset_status] += key.size();
        if (!IsLocalFittingClusterHealthy(health))
        {
            summary.unhealthy_cluster_count++;
            summary.unhealthy_atom_count += key.size();
        }
        if (!IsJointOffsetSolveStationarityEligible(health.joint_offset_status))
        {
            summary.unhealthy_joint_status_count[health.joint_offset_status]++;
        }
        if (health.unhealthy_refit_atom_count > 0)
        {
            summary.unhealthy_refit_cluster_count++;
            summary.unhealthy_refit_atom_count += health.unhealthy_refit_atom_count;
        }
    }
    return summary;
}

class LocalFittingJointOffsetStatusTracker
{
    std::map<LocalFittingClusterKey, JointOffsetSolveStatus> m_status_by_key{};

public:
    void Reconcile(const std::vector<LocalFittingClusterKey> & key_list)
    {
        std::map<LocalFittingClusterKey, JointOffsetSolveStatus> reconciled_status_by_key;
        for (const auto & key : key_list)
        {
            const auto iter{ m_status_by_key.find(key) };
            if (iter == m_status_by_key.end()) continue;
            reconciled_status_by_key.emplace(key, iter->second);
        }
        m_status_by_key = std::move(reconciled_status_by_key);
    }

    std::vector<LocalFittingClusterKey> FindIncompatible(
        const LocalFittingClusterHealthMap & health_by_key) const
    {
        std::vector<LocalFittingClusterKey> incompatible_key_list;
        for (const auto & [key, health] : health_by_key)
        {
            if (!IsJointOffsetSolveProgressEligible(health.joint_offset_status)) continue;
            const auto committed_iter{ m_status_by_key.find(key) };
            if (committed_iter != m_status_by_key.end() &&
                committed_iter->second != health.joint_offset_status)
            {
                incompatible_key_list.emplace_back(key);
            }
        }
        return incompatible_key_list;
    }

    void Invalidate(const std::vector<LocalFittingClusterKey> & key_list)
    {
        for (const auto & key : key_list)
        {
            m_status_by_key.erase(key);
        }
    }

    void InvalidateContaining(const std::vector<std::size_t> & atom_index_list)
    {
        for (auto iter = m_status_by_key.begin(); iter != m_status_by_key.end();)
        {
            const auto contains_affected_atom{
                std::any_of(
                    atom_index_list.begin(), atom_index_list.end(),
                    [&](std::size_t atom_index)
                    {
                        return std::binary_search(
                            iter->first.begin(), iter->first.end(), atom_index);
                    })
            };
            if (contains_affected_atom)
            {
                iter = m_status_by_key.erase(iter);
            }
            else
            {
                ++iter;
            }
        }
    }

    void Commit(
        const std::vector<LocalFittingClusterKey> & key_list,
        const LocalFittingClusterHealthMap & health_by_key)
    {
        for (const auto & key : key_list)
        {
            const auto health_iter{ health_by_key.find(key) };
            if (health_iter == health_by_key.end() ||
                !IsJointOffsetSolveProgressEligible(
                    health_iter->second.joint_offset_status))
            {
                throw std::invalid_argument(
                    "Joint offset status is missing or ineligible for commit.");
            }
            m_status_by_key[key] = health_iter->second.joint_offset_status;
        }
    }
};

struct ClusteredJointOffsetSolveResult
{
    Eigen::VectorXd offset{};
    ActiveCouplingGraph active_coupling_graph{};
    LocalFittingClusterHealthMap health_by_key{};
    detail::LocalFittingAndersonRegimeSignatureMap anderson_regime_signature_by_key{};
};

struct LocalFittingIterationResult
{
    GaussianFittingState state{};
    std::vector<std::size_t> suspicious_offset_state_index_list{};
    LocalFittingClusterHealthMap health_by_key{};
    detail::LocalFittingAndersonRegimeSignatureMap anderson_regime_signature_by_key{};
};

struct LocalAtomRefitResult
{
    LocalGaussianResult result{};
    bool is_health_eligible{ false };
};

struct LocalFittingObjectiveAtomModel
{
    std::size_t atom_index{ 0 };
    GaussianModel3D model{};
};

struct LocalFittingObjectiveSamples
{
    std::vector<double> residual_list{};
    std::vector<LocalFittingObjectiveAtomModel> active_model_list{};
    double scale_sample{ 0.0 };
};

struct LocalFittingObjectiveSampleRef
{
    std::size_t atom_index{ 0 };
    std::size_t sample_index{ 0 };
};

struct LocalFittingParameterPenaltyComponents
{
    double width_prior_penalty_sum{ 0.0 };
    double offset_plausibility_penalty_sum{ 0.0 };
    double movement_penalty_sum{ 0.0 };
};

struct LocalFittingBestAuditState
{
    std::vector<LocalFittingObjectiveSampleRef> sample_ref_list{};
    std::vector<std::size_t> atom_index_list{};
    std::optional<double> fixed_objective_scale{};
    std::optional<detail::LocalFittingObjectiveBreakdown> best_objective{};
    std::optional<GaussianFittingState> best_state{};
    std::optional<std::size_t> best_accepted_iteration{};
};

struct LocalFittingOffsetStats
{
    std::size_t atom_count{ 0 };
    std::size_t finite_count{ 0 };
    std::size_t exact_zero_count{ 0 };
    double median_absolute_offset{ 0.0 };
    double percentile_absolute_offset{ 0.0 };
    double maximum_absolute_offset{ 0.0 };
};

LocalFittingOffsetStats SummarizeLocalFittingOffsetValues(
    const std::vector<double> & offset_list);

void AppendLocalFittingOffsetSummary(
    std::ostringstream & stream,
    const LocalFittingOffsetStats & stats);

struct PersistentSuspiciousRollbackState
{
    std::vector<std::size_t> suspicious_atom_index_list{};
    std::size_t stable_iteration_count{ 0 };
};

using PersistentSuspiciousRollbackStateMap =
    std::map<LocalFittingClusterKey, PersistentSuspiciousRollbackState>;

struct PersistentJointOffsetFailureState
{
    JointOffsetSolveStatus status{ JointOffsetSolveStatus::SystemBuildFailed };
    std::size_t stable_iteration_count{ 0 };
};

using PersistentJointOffsetFailureStateMap =
    std::map<LocalFittingClusterKey, PersistentJointOffsetFailureState>;
using TerminalJointOffsetFailureMap =
    std::map<LocalFittingClusterKey, JointOffsetSolveStatus>;

struct LocalFittingTerminalSummary
{
    std::size_t suspicious_cluster_count{ 0 };
    std::size_t suspicious_atom_count{ 0 };
    std::size_t joint_offset_failure_cluster_count{ 0 };
    std::size_t joint_offset_failure_atom_count{ 0 };
    std::map<JointOffsetSolveStatus, std::size_t> joint_offset_failure_status_count{};

    std::size_t AtomCount() const
    {
        return suspicious_atom_count + joint_offset_failure_atom_count;
    }
};

enum class LocalFittingObjectiveAttemptDiagnosticStatus
{
    Scored,
    InvalidModel,
    ObjectiveUnavailable
};

enum class LocalFittingCandidateBuildFailureReason
{
    ParameterSize,
    NonFiniteParameter,
    NonPositiveWidth
};

struct LocalFittingCandidateBuildFailure
{
    std::size_t atom_index{ 0 };
    LocalFittingCandidateBuildFailureReason reason{
        LocalFittingCandidateBuildFailureReason::ParameterSize
    };
    Eigen::VectorXd estimation{};
};

struct LocalFittingObjectiveAttemptDiagnostic
{
    LocalFittingCandidateAttempt attempt{};
    LocalFittingObjectiveAttemptDiagnosticStatus status{
        LocalFittingObjectiveAttemptDiagnosticStatus::ObjectiveUnavailable
    };
    std::optional<double> objective_scale{};
    std::optional<detail::LocalFittingObjectiveBreakdown> candidate_objective{};
    std::optional<detail::LocalFittingObjectiveBreakdown> previous_objective{};
    std::optional<detail::LocalFittingObjectiveBreakdown> best_objective{};
    bool rejected_by_previous{ false };
    bool rejected_by_best{ false };
    std::optional<LocalFittingCandidateBuildFailure> build_failure{};
};

struct LocalFittingRejectedClusterDiagnostic
{
    LocalFittingClusterKey key{};
    std::vector<LocalFittingObjectiveAttemptDiagnostic> attempt_list{};
};

struct LocalFittingClusterWork
{
    std::vector<LocalFittingObjectiveSampleRef> objective_sample_ref_list{};
    std::optional<LocalFittingCandidateKind> accepted_kind{};
    bool rejected_anderson{ false };
    std::vector<LocalFittingObjectiveAttemptDiagnostic> objective_attempt_list{};
};

using LocalFittingClusterMap = std::map<LocalFittingClusterKey, LocalFittingClusterWork>;

struct LocalFittingCandidateSelection
{
    GaussianFittingState assembled_state{};
    std::vector<LocalFittingClusterKey> accepted_key_list{};
    std::vector<LocalFittingClusterKey> rejected_key_list{};
    std::vector<LocalFittingRejectedClusterDiagnostic> rejected_cluster_diagnostic_list{};
    std::optional<LocalFittingCandidateAttempt> accepted_candidate_attempt{};
    bool has_objective_backtracking_rejection{ false };
};

struct LocalFittingClusterSelectionSummary
{
    std::size_t accepted_cluster_count{ 0 };
    std::size_t rejected_cluster_count{ 0 };
    std::size_t accepted_atom_count{ 0 };
    std::size_t rejected_atom_count{ 0 };
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
    double prior_width{ 1.0 };
};

struct SecondStageLocalFittingContext
{
    std::vector<AtomObject *> atom_list{};
    std::vector<SecondStageAtomContext> atom_context_list{};

    std::size_t AtomSize() const { return atom_list.size(); }
};

struct GaussianModelParameterSamples
{
    std::vector<double> amplitude_list{};
    std::vector<double> width_list{};
    std::vector<double> offset_list{};
};

struct SecondStageSeedRepairRecord
{
    std::size_t atom_index{ 0 };
    detail::SecondStageSeedRepairSource source{
        detail::SecondStageSeedRepairSource::GlobalMedian
    };
    GaussianModel3D original_model{};
    GaussianModel3D repaired_model{};
};

struct SecondStageInitialStateBuildResult
{
    std::optional<GaussianFittingState> state{};
    std::vector<SecondStageSeedRepairRecord> repair_record_list{};
};

const char * GetSecondStageSeedRepairSourceText(
    detail::SecondStageSeedRepairSource source)
{
    switch (source)
    {
    case detail::SecondStageSeedRepairSource::GroupPosterior:
        return "group-posterior";
    case detail::SecondStageSeedRepairSource::GroupPrior:
        return "group-prior";
    case detail::SecondStageSeedRepairSource::LocalOls:
        return "local-ols";
    case detail::SecondStageSeedRepairSource::GroupMedian:
        return "group-median";
    case detail::SecondStageSeedRepairSource::GlobalMedian:
        return "global-median";
    }
    throw std::logic_error("Unknown second-stage seed repair source.");
}

std::size_t GetSecondStageSeedRepairSourceIndex(
    detail::SecondStageSeedRepairSource source)
{
    return static_cast<std::size_t>(source);
}

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

bool IsSameSuspiciousProfileRadius(double lhs, double rhs)
{
    const auto scale{ std::max({ std::abs(lhs), std::abs(rhs), 1.0 }) };
    return std::abs(lhs - rhs) <= kSuspiciousProfileDistanceTolerance * scale;
}

std::optional<ZeroOffsetProfileDiagnostics> BuildZeroOffsetProfileDiagnostics(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model,
    const FitOptions & options)
{
    std::vector<std::pair<double, double>> profile_samples;
    profile_samples.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        const auto distance{ static_cast<double>(sample.point.distance) };
        if (distance < options.distance_min || distance > options.distance_max) continue;
        const auto response{ CalculateZeroOffsetResponse(sample, model) };
        if (!std::isfinite(response)) continue;
        if (std::abs(response) > static_cast<double>(std::numeric_limits<float>::max())) continue;
        profile_samples.emplace_back(distance, response);
    }

    if (profile_samples.empty()) return std::nullopt;
    ZeroOffsetProfileDiagnostics diagnostics;
    std::sort(
        profile_samples.begin(),
        profile_samples.end(),
        [](const auto & lhs, const auto & rhs)
        {
            return lhs.first < rhs.first;
        });

    diagnostics.distance_min = profile_samples.front().first;
    diagnostics.distance_max = profile_samples.back().first;
    for (std::size_t i = 0; i < profile_samples.size();)
    {
        const auto radius{ profile_samples.at(i).first };
        std::vector<double> response_list;
        while (i < profile_samples.size() &&
            IsSameSuspiciousProfileRadius(profile_samples.at(i).first, radius))
        {
            const auto response{ profile_samples.at(i).second };
            diagnostics.max_abs_response = std::max(diagnostics.max_abs_response, std::abs(response));
            response_list.emplace_back(response);
            i++;
        }
        diagnostics.radius_response_median_list.emplace_back(array_helper::ComputeMedian(response_list));
    }
    diagnostics.center_response = diagnostics.radius_response_median_list.front();
    return diagnostics;
}

bool HasUsableSuspiciousProfileBaseline(
    const GaussianModel3D & previous_model,
    const ZeroOffsetProfileDiagnostics & previous_profile)
{
    if (previous_profile.radius_response_median_list.size() < kSuspiciousProfileMinimumRadiusCount)
    {
        return false;
    }
    if (!std::isfinite(previous_model.GetAmplitude()) ||
        !std::isfinite(previous_model.GetWidth()) ||
        !std::isfinite(previous_model.GetOffset()) ||
        previous_model.GetWidth() <= 0.0 ||
        previous_profile.max_abs_response <= kRobustScaleMin)
    {
        return false;
    }
    const auto center_scale{
        std::max(std::abs(previous_profile.center_response), kRobustScaleMin)
    };
    for (std::size_t i = 1; i < previous_profile.radius_response_median_list.size(); i++)
    {
        const auto current_scale{
            std::abs(previous_profile.radius_response_median_list.at(i))
        };
        if (current_scale > kSuspiciousProfileReboundCenterRatio * center_scale) return false;
    }
    return true;
}

bool HasSuspiciousOffsetMagnitude(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    double previous_profile_max_abs_response)
{
    const auto previous_offset_response{
        previous_model.GetOffset() * previous_model.OffsetBasisAtDistance(0.0)
    };
    const auto candidate_offset_response{
        candidate_model.GetOffset() * candidate_model.OffsetBasisAtDistance(0.0)
    };
    if (!std::isfinite(previous_offset_response) || !std::isfinite(candidate_offset_response))
    {
        return true;
    }
    const auto reference_scale{
        std::max({
            std::abs(previous_model.SignalAtDistance(0.0)),
            std::abs(previous_offset_response),
            previous_profile_max_abs_response,
            kRobustScaleMin
        })
    };
    return std::abs(candidate_offset_response) >
        kSuspiciousCompensationResponseRatio * reference_scale;
}

bool HasSuspiciousCenterSignFlip(
    const ZeroOffsetProfileDiagnostics & previous_profile,
    const ZeroOffsetProfileDiagnostics & candidate_profile)
{
    const auto previous_center{ previous_profile.center_response };
    const auto candidate_center{ candidate_profile.center_response };
    const auto reference_scale{ std::max(std::abs(previous_center), kRobustScaleMin) };
    return previous_center * candidate_center < 0.0 &&
        std::abs(candidate_center) > kSuspiciousProfileCenterSignFlipRatio * reference_scale;
}

bool HasSuspiciousRadialRebound(
    const ZeroOffsetProfileDiagnostics & previous_profile,
    const ZeroOffsetProfileDiagnostics & candidate_profile)
{
    if (candidate_profile.radius_response_median_list.size() < kSuspiciousProfileMinimumRadiusCount)
    {
        return false;
    }
    const auto reference_center_scale{
        std::max(std::abs(previous_profile.center_response), kRobustScaleMin)
    };
    const auto candidate_center_scale{
        std::max(std::abs(candidate_profile.center_response), kRobustScaleMin)
    };
    int upward_excursion_count{ 0 };
    auto previous_abs_response{ std::abs(candidate_profile.radius_response_median_list.front()) };
    for (std::size_t i = 1; i < candidate_profile.radius_response_median_list.size(); i++)
    {
        const auto current_abs_response{
            std::abs(candidate_profile.radius_response_median_list.at(i))
        };
        if (current_abs_response > kSuspiciousProfileReboundCenterRatio * candidate_center_scale &&
            current_abs_response > kSuspiciousProfileReboundReferenceRatio * reference_center_scale)
        {
            return true;
        }
        if (current_abs_response >
            previous_abs_response + kSuspiciousProfileUpwardExcursionReferenceRatio * reference_center_scale)
        {
            upward_excursion_count++;
        }
        previous_abs_response = current_abs_response;
    }
    return upward_excursion_count > kSuspiciousProfileMaximumUpwardExcursions;
}

bool HasSuspiciousWidthGrowth(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    const ZeroOffsetProfileDiagnostics & previous_profile)
{
    if (!std::isfinite(candidate_model.GetWidth()) || candidate_model.GetWidth() <= 0.0) return true;
    if (candidate_model.GetWidth() > kSuspiciousWidthGrowthLimit * previous_model.GetWidth()) return true;
    const auto distance_range{ previous_profile.distance_max - previous_profile.distance_min };
    return distance_range > 0.0 && candidate_model.GetWidth() > kSuspiciousWidthRangeLimitRatio * distance_range;
}

bool HasSuspiciousAmplitudeOffsetCompensation(
    const GaussianModel3D & previous_model,
    const GaussianModel3D & candidate_model,
    const ZeroOffsetProfileDiagnostics & previous_profile)
{
    const auto signal_delta{
        candidate_model.SignalAtDistance(0.0) - previous_model.SignalAtDistance(0.0)
    };
    const auto offset_delta_response{
        (candidate_model.GetOffset() - previous_model.GetOffset()) * previous_model.OffsetBasisAtDistance(0.0)
    };
    const auto reference_scale{
        std::max({
            std::abs(previous_profile.center_response),
            std::abs(previous_model.SignalAtDistance(0.0)),
            kRobustScaleMin
        })
    };
    return signal_delta * offset_delta_response < 0.0 &&
        std::abs(signal_delta) > kSuspiciousCompensationResponseRatio * reference_scale &&
        std::abs(offset_delta_response) > kSuspiciousCompensationResponseRatio * reference_scale;
}

bool IsSuspiciousJointOffset(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & previous_model,
    const GaussianModel3D & offset_model,
    const FitOptions & options)
{
    if (!CanBuildFiniteZeroOffsetSamples(sample_entries, previous_model)) return false;
    if (!CanBuildFiniteZeroOffsetSamples(sample_entries, offset_model)) return true;
    const auto previous_profile{ BuildZeroOffsetProfileDiagnostics(sample_entries, previous_model, options) };
    if (HasSuspiciousOffsetMagnitude(
            previous_model,
            offset_model,
            previous_profile.has_value() ? previous_profile->max_abs_response : 0.0))
    {
        return true;
    }
    if (!previous_profile.has_value() ||
        !HasUsableSuspiciousProfileBaseline(previous_model, *previous_profile))
    {
        return false;
    }
    const auto candidate_profile{ BuildZeroOffsetProfileDiagnostics(sample_entries, offset_model, options) };
    if (!candidate_profile.has_value()) return true;
    return HasSuspiciousCenterSignFlip(*previous_profile, *candidate_profile) ||
        HasSuspiciousRadialRebound(*previous_profile, *candidate_profile) ||
        HasSuspiciousWidthGrowth(previous_model, offset_model, *previous_profile) ||
        HasSuspiciousAmplitudeOffsetCompensation(previous_model, offset_model, *previous_profile);
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
        linearization_service::DecodeParameterVector(fit_result.beta_ols).WithOffset(offset)
    };
    const auto mdpde_model{
        linearization_service::DecodeParameterVector(fit_result.beta_mdpde).WithOffset(offset)
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

SecondStageLocalFittingContext BuildSecondStageLocalFittingContext(ModelObject & model_object)
{
    SecondStageLocalFittingContext context;
    context.atom_list = model_object.GetSelectedAtoms();
    context.atom_context_list.resize(context.atom_list.size());
    std::unordered_map<const AtomObject *, std::size_t> atom_index_map;
    atom_index_map.reserve(context.atom_list.size());
    for (std::size_t i = 0; i < context.atom_list.size(); i++)
    {
        atom_index_map.emplace(context.atom_list.at(i), i);
    }
    const auto analysis_view{ model_object.GetAnalysisView() };

    std::unordered_map<GroupKey, std::vector<double>> width_samples_by_group;
    width_samples_by_group.reserve(context.atom_list.size());
    for (const auto * atom : context.atom_list)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        const auto width{ local_view.GetEstimateMDPDE().GetWidth() };
        if (numeric_validation::IsFinitePositive(width))
        {
            width_samples_by_group[data_internal::GetGroupKey(atom)].emplace_back(width);
        }
    }

    std::unordered_map<GroupKey, double> median_width_by_group;
    median_width_by_group.reserve(width_samples_by_group.size());
    for (const auto & [group_key, width_samples] : width_samples_by_group)
    {
        if (width_samples.empty()) continue;
        const auto median_width{ array_helper::ComputeMedian(width_samples) };
        if (numeric_validation::IsFinitePositive(median_width))
        {
            median_width_by_group.emplace(group_key, median_width);
        }
    }

    for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
    {
        const auto * atom{ context.atom_list.at(atom_index) };
        auto & atom_context{ context.atom_context_list.at(atom_index) };
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        atom_context.sample_entries = local_view.GetSamplingEntries(false);
        atom_context.alpha_r = local_view.GetAlphaR();
        const auto group_key{ data_internal::GetGroupKey(atom) };
        const auto current_width{ local_view.GetEstimateMDPDE().GetWidth() };
        atom_context.prior_width = numeric_validation::IsFinitePositive(current_width) ?
            current_width : 1.0;
        const auto median_width_iter{ median_width_by_group.find(group_key) };
        if (median_width_iter != median_width_by_group.end())
        {
            atom_context.prior_width = median_width_iter->second;
        }
        if (analysis_view.HasAtomGroup(group_key))
        {
            const auto group_prior_width{ analysis_view.GetAtomGroupPrior(group_key).GetWidth() };
            if (numeric_validation::IsFinitePositive(group_prior_width))
            {
                atom_context.prior_width = group_prior_width;
            }
        }
    }

    for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
    {
        const auto * atom{ context.atom_list.at(atom_index) };
        auto & atom_context{ context.atom_context_list.at(atom_index) };
        const auto neighbor_atom_list{ atom->FindNeighborAtoms(kNeighborAtomSearchRange) };
        atom_context.selected_neighbor_index_list.reserve(neighbor_atom_list.size());
        for (const auto * neighbor_atom : neighbor_atom_list)
        {
            const auto neighbor_iter{ atom_index_map.find(neighbor_atom) };
            if (neighbor_iter == atom_index_map.end()) continue;

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

std::size_t FindLocalFittingClusterRoot(std::vector<std::size_t> & parent_list, std::size_t index)
{
    if (index >= parent_list.size())
    {
        throw std::invalid_argument("Local fitting cluster index is out of range.");
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

void MergeLocalFittingClusters(
    std::vector<std::size_t> & parent_list,
    std::size_t left_index,
    std::size_t right_index)
{
    const auto left_root{ FindLocalFittingClusterRoot(parent_list, left_index) };
    const auto right_root{ FindLocalFittingClusterRoot(parent_list, right_index) };
    if (left_root == right_root) return;

    parent_list.at(right_root) = left_root;
}

LocalFittingClusterMap BuildLocalFittingClusters(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list)
{
    std::vector<std::optional<std::size_t>> active_position_by_atom_index(
        context.AtomSize());
    for (std::size_t active_position = 0; active_position < active_index_list.size(); active_position++)
    {
        const auto active_index{ active_index_list.at(active_position) };
        if (active_index >= context.AtomSize())
        {
            throw std::invalid_argument("Local fitting active index is out of range.");
        }
        active_position_by_atom_index.at(active_index) = active_position;
    }

    std::vector<std::size_t> parent_list(active_index_list.size());
    for (std::size_t i = 0; i < parent_list.size(); i++)
    {
        parent_list.at(i) = i;
    }

    std::vector<std::pair<LocalFittingObjectiveSampleRef, std::size_t>>
        sample_representative_position_list;
    for (std::size_t atom_index = 0; atom_index < context.AtomSize(); atom_index++)
    {
        const auto & atom_context{ context.atom_context_list.at(atom_index) };
        for (std::size_t sample_index = 0; sample_index < atom_context.sample_entries.size(); sample_index++)
        {
            std::optional<std::size_t> representative_position;
            const auto target_active_position{ active_position_by_atom_index.at(atom_index) };
            if (target_active_position.has_value())
            {
                representative_position = target_active_position;
            }
            for (const auto & neighbor_sample : atom_context.sample_neighbor_list.at(sample_index))
            {
                if (neighbor_sample.atom_index >= active_position_by_atom_index.size())
                {
                    throw std::invalid_argument("Local fitting neighbor index is out of range.");
                }
                const auto neighbor_active_position{
                    active_position_by_atom_index.at(neighbor_sample.atom_index)
                };
                if (!neighbor_active_position.has_value()) continue;
                if (!representative_position.has_value())
                {
                    representative_position = neighbor_active_position;
                    continue;
                }
                MergeLocalFittingClusters(
                    parent_list,
                    *representative_position,
                    *neighbor_active_position);
            }

            if (representative_position.has_value())
            {
                sample_representative_position_list.emplace_back(
                    LocalFittingObjectiveSampleRef{ atom_index, sample_index },
                    *representative_position);
            }
        }
    }

    std::map<std::size_t, LocalFittingClusterKey> active_index_list_by_root;
    for (std::size_t active_position = 0; active_position < active_index_list.size(); active_position++)
    {
        const auto root{ FindLocalFittingClusterRoot(parent_list, active_position) };
        active_index_list_by_root[root].emplace_back(active_index_list.at(active_position));
    }

    std::map<std::size_t, std::vector<LocalFittingObjectiveSampleRef>> sample_ref_list_by_root;
    for (const auto & [sample_ref, representative_position] :
        sample_representative_position_list)
    {
        const auto root{
            FindLocalFittingClusterRoot(parent_list, representative_position)
        };
        sample_ref_list_by_root[root].emplace_back(sample_ref);
    }

    LocalFittingClusterMap cluster_map;
    for (auto & [root, cluster_active_index_list] : active_index_list_by_root)
    {
        std::sort(cluster_active_index_list.begin(), cluster_active_index_list.end());
        const auto inserted{
            cluster_map.emplace(
                std::move(cluster_active_index_list),
                LocalFittingClusterWork{ std::move(sample_ref_list_by_root[root]) })
        };
        if (!inserted.second)
        {
            throw std::logic_error("Local fitting cluster key is duplicated.");
        }
    }
    return cluster_map;
}

void AddGaussianModelParameterSample(
    GaussianModelParameterSamples & samples,
    const GaussianModel3D & model)
{
    if (!detail::IsValidSecondStageGaussianModel(model)) return;
    samples.amplitude_list.emplace_back(model.GetAmplitude());
    samples.width_list.emplace_back(model.GetWidth());
    samples.offset_list.emplace_back(model.GetOffset());
}

std::optional<GaussianModel3DWithUncertainty> BuildValidGaussianParameterMedian(
    const GaussianModelParameterSamples & samples)
{
    if (samples.amplitude_list.empty() ||
        samples.width_list.empty() ||
        samples.offset_list.empty())
    {
        return std::nullopt;
    }
    const GaussianModel3D median_model{
        array_helper::ComputeMedian(samples.amplitude_list),
        array_helper::ComputeMedian(samples.width_list),
        array_helper::ComputeMedian(samples.offset_list)
    };
    if (!detail::IsValidSecondStageGaussianModel(median_model))
    {
        return std::nullopt;
    }
    return GaussianModel3DWithUncertainty{
        median_model,
        GaussianModel3DUncertainty{}
    };
}

SecondStageInitialStateBuildResult BuildInitialLocalFittingState(
    const SecondStageLocalFittingContext & context,
    const ModelAnalysisView & analysis_view)
{
    SecondStageInitialStateBuildResult build_result;
    GaussianFittingState state{
        std::vector<LocalGaussianResult>(context.AtomSize()),
        std::vector<Eigen::VectorXd>(context.AtomSize())
    };
    std::vector<std::optional<GaussianModel3DWithUncertainty>> group_prior_list(
        context.AtomSize());
    std::unordered_map<GroupKey, GaussianModelParameterSamples> samples_by_group;
    GaussianModelParameterSamples global_samples;

    for (std::size_t i = 0; i < context.AtomSize(); i++)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*context.atom_list.at(i)) };
        state.result_list.at(i) = local_view.GetGaussianResult();
        const auto group_key{ data_internal::GetGroupKey(context.atom_list.at(i)) };
        if (analysis_view.HasAtomGroup(group_key))
        {
            group_prior_list.at(i) =
                analysis_view.GetAtomGroupPriorWithUncertainty(group_key);
        }

        const auto & result{ state.result_list.at(i) };
        const auto & mdpde_model{ result.mdpde.GetModel() };
        std::optional<GaussianModel3DWithUncertainty> preferred_model;
        if (detail::IsValidSecondStageGaussianModel(mdpde_model))
        {
            preferred_model = result.mdpde;
        }
        else
        {
            const auto direct_selection{
                detail::SelectSecondStageSeedRepair(
                    detail::SecondStageSeedRepairCandidates{
                        result.posterior,
                        group_prior_list.at(i),
                        result.ols,
                        std::nullopt,
                        std::nullopt
                    })
            };
            if (direct_selection.has_value())
            {
                preferred_model = direct_selection->model;
            }
        }
        if (!preferred_model.has_value()) continue;

        AddGaussianModelParameterSample(
            samples_by_group[group_key],
            preferred_model->GetModel());
        AddGaussianModelParameterSample(
            global_samples,
            preferred_model->GetModel());
    }

    std::unordered_map<GroupKey, GaussianModel3DWithUncertainty> median_by_group;
    median_by_group.reserve(samples_by_group.size());
    for (const auto & [group_key, samples] : samples_by_group)
    {
        const auto median_model{ BuildValidGaussianParameterMedian(samples) };
        if (median_model.has_value())
        {
            median_by_group.emplace(group_key, *median_model);
        }
    }
    const auto global_median{ BuildValidGaussianParameterMedian(global_samples) };

    for (std::size_t i = 0; i < context.AtomSize(); i++)
    {
        auto & result{ state.result_list.at(i) };
        const auto original_model{ result.mdpde.GetModel() };
        if (!detail::IsValidSecondStageGaussianModel(original_model))
        {
            const auto group_key{ data_internal::GetGroupKey(context.atom_list.at(i)) };
            std::optional<GaussianModel3DWithUncertainty> group_median;
            const auto group_median_iter{ median_by_group.find(group_key) };
            if (group_median_iter != median_by_group.end())
            {
                group_median = group_median_iter->second;
            }
            const auto selection{
                detail::SelectSecondStageSeedRepair(
                    detail::SecondStageSeedRepairCandidates{
                        result.posterior,
                        group_prior_list.at(i),
                        result.ols,
                        group_median,
                        global_median
                    })
            };
            if (!selection.has_value())
            {
                return build_result;
            }

            const auto repaired_seed{
                detail::BuildRepairedSecondStageSeed(original_model, *selection)
            };
            const auto repaired_model{ repaired_seed.GetModel() };
            if (!detail::IsValidSecondStageGaussianModel(repaired_model))
            {
                return build_result;
            }
            result.mdpde = repaired_seed;
            build_result.repair_record_list.emplace_back(
                SecondStageSeedRepairRecord{
                    i,
                    selection->source,
                    original_model,
                    repaired_model
                });
        }
        state.estimation_list.at(i) = result.mdpde.GetModel().ToVector();
    }
    build_result.state = std::move(state);
    return build_result;
}

void LogSecondStageSeedRepairs(
    const std::vector<SecondStageSeedRepairRecord> & repair_record_list,
    const FitOptions & options)
{
    if (options.quiet_mode || repair_record_list.empty()) return;

    constexpr std::array<detail::SecondStageSeedRepairSource, 5> source_list{
        detail::SecondStageSeedRepairSource::GroupPosterior,
        detail::SecondStageSeedRepairSource::GroupPrior,
        detail::SecondStageSeedRepairSource::LocalOls,
        detail::SecondStageSeedRepairSource::GroupMedian,
        detail::SecondStageSeedRepairSource::GlobalMedian
    };
    std::array<std::size_t, source_list.size()> source_count{};
    for (const auto & record : repair_record_list)
    {
        source_count.at(GetSecondStageSeedRepairSourceIndex(record.source))++;
    }

    std::ostringstream summary;
    summary << "Repaired invalid second-stage seed atoms = "
        << repair_record_list.size() << ", sources = ";
    for (std::size_t i = 0; i < source_list.size(); i++)
    {
        if (i != 0) summary << ", ";
        summary << GetSecondStageSeedRepairSourceText(source_list.at(i))
            << ":" << source_count.at(i);
    }
    summary << ".";
    Logger::Log(LogLevel::Info, summary.str());

    for (const auto & record : repair_record_list)
    {
        std::ostringstream detail_message;
        detail_message << "Second-stage seed repair: atom index = "
            << record.atom_index
            << ", source = " << GetSecondStageSeedRepairSourceText(record.source)
            << std::scientific
            << ", original A/B/C = "
            << record.original_model.GetAmplitude() << "/"
            << record.original_model.GetWidth() << "/"
            << record.original_model.GetOffset()
            << ", repaired A/B/C = "
            << record.repaired_model.GetAmplitude() << "/"
            << record.repaired_model.GetWidth() << "/"
            << record.repaired_model.GetOffset() << ".";
        Logger::Log(LogLevel::Debug, detail_message.str());
    }
}

FittedGaussianSnapshot BuildFittedGaussianSnapshot(
    const std::vector<Eigen::VectorXd> & estimation_list)
{
    FittedGaussianSnapshot snapshot;
    snapshot.reserve(estimation_list.size());
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
    std::vector<double> effective_ridge_multiplier_list(
        static_cast<std::size_t>(column_count), 1.0);
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
        effective_ridge_multiplier_list.at(static_cast<std::size_t>(column_index)) =
            combined_multiplier;
        const auto base_ridge{
            square_sum > std::numeric_limits<double>::epsilon()
                ? ridge_ratio * square_sum
                : ridge_ratio / kJointOffsetRidgeRatio
        };
        system.ridge_diagonal(column_index) = combined_multiplier * base_ridge;
    }
    return JointOffsetBuildResult{
        std::move(system),
        std::move(active_coupling_graph),
        std::move(effective_ridge_multiplier_list)
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
    return kRobustScaleMultiplier * array_helper::ComputeMedian(deviation_list);
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
    catch (const std::runtime_error &)
    {
        return JointOffsetSolveResult{
            JointOffsetSolveStatus::SystemBuildFailed,
            previous_offset,
            {},
            {}
        };
    }
    auto system{ std::move(build_result.system) };
    auto active_coupling_graph{ std::move(build_result.active_coupling_graph) };
    auto effective_ridge_multiplier_list{
        std::move(build_result.effective_ridge_multiplier_list)
    };
    if (system.response.size() == 0 || system.previous_parameter.size() == 0)
    {
        return JointOffsetSolveResult{
            JointOffsetSolveStatus::EmptySystem,
            previous_offset,
            std::move(active_coupling_graph),
            std::move(effective_ridge_multiplier_list)
        };
    }

    Eigen::VectorXd weight{ Eigen::VectorXd::Ones(system.response.size()) };
    algorithm::WeightedRidgeSolver solver{ system };
    Eigen::VectorXd offset;
    if (!solver.Solve(system, weight, offset))
    {
        return JointOffsetSolveResult{
            JointOffsetSolveStatus::InitialSolveFailed,
            previous_offset,
            std::move(active_coupling_graph),
            std::move(effective_ridge_multiplier_list)
        };
    }

    for (int iteration = 0; iteration < kRobustLossMaximumIterations; iteration++)
    {
        const Eigen::VectorXd residual{ system.response - system.design_matrix * offset };
        std::vector<double> residual_list(residual.data(), residual.data() + residual.size());
        const auto residual_scale{
            std::max(CalculateMedianAbsoluteDeviationScale(residual_list), kRobustScaleMin)
        };
        for (Eigen::Index i = 0; i < residual.size(); i++)
        {
            weight(i) = algorithm::CalculateRobustWeight(
                kSecondStageRobustLossKind,
                residual(i),
                residual_scale,
                kRobustLossCutoffMultiplier);
        }

        Eigen::VectorXd updated_offset;
        if (!solver.Solve(system, weight, updated_offset))
        {
            return JointOffsetSolveResult{
                JointOffsetSolveStatus::IrlsSolveFailed,
                previous_offset,
                std::move(active_coupling_graph),
                std::move(effective_ridge_multiplier_list)
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
            return JointOffsetSolveResult{
                JointOffsetSolveStatus::IrlsObjectiveDeteriorated,
                offset,
                std::move(active_coupling_graph),
                std::move(effective_ridge_multiplier_list)
            };
        }
        const auto maximum_change{
            algorithm::CalculateMaximumNormalizedVectorChange(updated_offset, offset, kJointOffsetIrlsScaleFloor)
        };
        offset = std::move(updated_offset);
        if (maximum_change < kJointOffsetIrlsNormalizedChangeTolerance)
        {
            return JointOffsetSolveResult{
                JointOffsetSolveStatus::Converged,
                offset,
                std::move(active_coupling_graph),
                std::move(effective_ridge_multiplier_list)
            };
        }
    }

    return JointOffsetSolveResult{
        JointOffsetSolveStatus::IrlsMaximumIterationsReached,
        offset,
        std::move(active_coupling_graph),
        std::move(effective_ridge_multiplier_list)
    };
}

ClusteredJointOffsetSolveResult EstimateClusteredJointOffsets(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const std::vector<LocalFittingClusterKey> & cluster_key_list,
    const FittedGaussianSnapshot & snapshot,
    double ridge_ratio,
    const std::vector<double> & ridge_multiplier_list)
{
    std::unordered_map<std::size_t, std::size_t> active_position_by_atom_index;
    active_position_by_atom_index.reserve(active_index_list.size());
    for (std::size_t active_position = 0;
        active_position < active_index_list.size();
        active_position++)
    {
        if (!active_position_by_atom_index.emplace(
                active_index_list.at(active_position), active_position).second)
        {
            throw std::invalid_argument(
                "Clustered joint offset active atom indexes must be unique.");
        }
    }

    ClusteredJointOffsetSolveResult clustered_result;
    clustered_result.offset = Eigen::VectorXd::Zero(
        static_cast<Eigen::Index>(active_index_list.size()));
    clustered_result.active_coupling_graph.resize(active_index_list.size());
    std::vector<char> covered_active_position_list(active_index_list.size(), 0);
    for (const auto & key : cluster_key_list)
    {
        if (key.empty() || !std::is_sorted(key.begin(), key.end()) ||
            std::adjacent_find(key.begin(), key.end()) != key.end())
        {
            throw std::invalid_argument(
                "Clustered joint offset key must be non-empty and canonical.");
        }

        const auto cluster_result{
            EstimateJointOffsets(
                context,
                key,
                snapshot,
                ridge_ratio,
                ridge_multiplier_list)
        };
        if (cluster_result.offset.size() != static_cast<Eigen::Index>(key.size()))
        {
            throw std::logic_error("Clustered joint offset result size is inconsistent.");
        }
        if (!cluster_result.active_coupling_graph.empty() &&
            cluster_result.active_coupling_graph.size() != key.size())
        {
            throw std::logic_error("Clustered joint offset graph size is inconsistent.");
        }
        if (!clustered_result.health_by_key.emplace(
                key,
                LocalFittingClusterHealth{ cluster_result.status, 0 }).second)
        {
            throw std::invalid_argument("Clustered joint offset keys must be unique.");
        }

        for (std::size_t local_position = 0; local_position < key.size(); local_position++)
        {
            const auto atom_index{ key.at(local_position) };
            const auto active_iter{ active_position_by_atom_index.find(atom_index) };
            if (active_iter == active_position_by_atom_index.end())
            {
                throw std::invalid_argument("Clustered joint offset atom is not active.");
            }
            const auto active_position{ active_iter->second };
            if (covered_active_position_list.at(active_position) != 0)
            {
                throw std::invalid_argument(
                    "Clustered joint offset keys must not share active atoms.");
            }
            covered_active_position_list.at(active_position) = 1;
            clustered_result.offset(static_cast<Eigen::Index>(active_position)) =
                cluster_result.offset(static_cast<Eigen::Index>(local_position));

            if (cluster_result.active_coupling_graph.empty()) continue;
            for (const auto & edge : cluster_result.active_coupling_graph.at(local_position))
            {
                if (edge.neighbor_index >= key.size())
                {
                    throw std::logic_error(
                        "Clustered joint offset graph neighbor is out of range.");
                }
                const auto neighbor_iter{
                    active_position_by_atom_index.find(key.at(edge.neighbor_index))
                };
                if (neighbor_iter == active_position_by_atom_index.end())
                {
                    throw std::logic_error(
                        "Clustered joint offset graph neighbor is not active.");
                }
                clustered_result.active_coupling_graph.at(active_position).emplace_back(
                    ActiveCouplingEdge{ neighbor_iter->second, edge.overlap });
            }
        }

        if (!cluster_result.effective_ridge_multiplier_list.empty())
        {
            if (cluster_result.effective_ridge_multiplier_list.size() != key.size())
            {
                throw std::logic_error(
                    "Clustered joint offset ridge multiplier size is inconsistent.");
            }
            detail::LocalFittingAndersonRegimeSignature signature{
                ridge_ratio,
                cluster_result.effective_ridge_multiplier_list
            };
            detail::ValidateLocalFittingAndersonRegimeSignature(signature);
            clustered_result.anderson_regime_signature_by_key.emplace(
                key,
                std::move(signature));
        }
    }

    if (std::find(
            covered_active_position_list.begin(),
            covered_active_position_list.end(),
            0) != covered_active_position_list.end())
    {
        throw std::invalid_argument(
            "Clustered joint offset keys must cover active atoms exactly once.");
    }
    return clustered_result;
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

std::optional<LocalFittingObjectiveSamples> CollectLocalFittingObjectiveSamples(
    const SecondStageLocalFittingContext & context,
    const std::vector<Eigen::VectorXd> & estimation_list,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const std::vector<std::size_t> & active_index_list)
{
    if (estimation_list.size() != context.AtomSize())
    {
        throw std::invalid_argument(
            "Local fitting objective estimation size is inconsistent.");
    }
    const auto snapshot{ BuildFittedGaussianSnapshot(estimation_list) };

    LocalFittingObjectiveSamples objective_samples;
    objective_samples.residual_list.reserve(sample_ref_list.size());
    objective_samples.active_model_list.reserve(active_index_list.size());
    for (const auto active_index : active_index_list)
    {
        if (active_index >= snapshot.size())
        {
            throw std::invalid_argument("Local fitting objective active index is out of range.");
        }
        objective_samples.active_model_list.emplace_back(LocalFittingObjectiveAtomModel{
            active_index,
            snapshot.at(active_index)
        });
    }
    std::vector<double> response_list;
    response_list.reserve(sample_ref_list.size());
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
        response_list.emplace_back(response);
    }
    if (objective_samples.residual_list.empty())
    {
        return std::nullopt;
    }

    const auto residual_scale{
        CalculateMedianAbsoluteDeviationScale(objective_samples.residual_list)
    };
    const auto response_scale_floor{
        kLocalFittingObjectiveResidualScaleFloorRatio *
        CalculateMedianAbsoluteDeviationScale(response_list)
    };
    objective_samples.scale_sample =
        std::max({ residual_scale, response_scale_floor, kRobustScaleMin });
    if (!std::isfinite(objective_samples.scale_sample)) return std::nullopt;
    return std::move(objective_samples);
}

double CalculateSquaredValue(double value)
{
    return value * value;
}

std::optional<LocalFittingParameterPenaltyComponents>
CalculateLocalFittingParameterPenaltyComponents(
    const SecondStageLocalFittingContext & context,
    const LocalFittingObjectiveSamples & objective_samples,
    double objective_scale,
    const std::vector<LocalFittingObjectiveAtomModel> * movement_reference_model_list)
{
    if (!numeric_validation::IsFinitePositive(objective_scale))
    {
        return std::nullopt;
    }
    if (movement_reference_model_list != nullptr &&
        movement_reference_model_list->size() != objective_samples.active_model_list.size())
    {
        return std::nullopt;
    }

    LocalFittingParameterPenaltyComponents components;
    const auto residual_scale_floor{ std::max(objective_scale, kRobustScaleMin) };
    for (std::size_t model_position = 0;
        model_position < objective_samples.active_model_list.size();
        model_position++)
    {
        const auto & atom_model{
            objective_samples.active_model_list.at(model_position)
        };
        const auto active_index{ atom_model.atom_index };
        if (active_index >= context.AtomSize()) return std::nullopt;

        const auto & model{ atom_model.model };
        if (!detail::IsValidSecondStageGaussianModel(model)) return std::nullopt;

        const auto prior_width{ context.atom_context_list.at(active_index).prior_width };
        if (!numeric_validation::IsFinitePositive(prior_width)) return std::nullopt;
        components.width_prior_penalty_sum += CalculateSquaredValue(
            (std::log(model.GetWidth()) - std::log(prior_width)) / kLocalFittingWidthPriorLogScale);

        const auto peak_signal{ model.SignalAtDistance(0.0) };
        const auto offset_peak{ model.GetOffset() * model.OffsetBasisAtDistance(0.0) };
        if (!std::isfinite(peak_signal) || !std::isfinite(offset_peak))
        {
            return std::nullopt;
        }
        const auto offset_ratio{
            std::abs(offset_peak) /
            std::max({ std::abs(peak_signal), residual_scale_floor, kRobustScaleMin })
        };
        if (!std::isfinite(offset_ratio)) return std::nullopt;
        const auto offset_excess{
            std::max(0.0, offset_ratio - kLocalFittingOffsetPeakRatioMax)
        };
        components.offset_plausibility_penalty_sum += CalculateSquaredValue(offset_excess);

        if (movement_reference_model_list != nullptr)
        {
            const auto & previous_atom_model{
                movement_reference_model_list->at(model_position)
            };
            if (previous_atom_model.atom_index != active_index) return std::nullopt;
            const auto & previous_model{ previous_atom_model.model };
            if (!detail::IsValidSecondStageGaussianModel(previous_model))
            {
                return std::nullopt;
            }

            const auto previous_offset_basis{ previous_model.OffsetBasisAtDistance(0.0) };
            const auto previous_peak_signal{ previous_model.SignalAtDistance(0.0) };
            if (!std::isfinite(previous_offset_basis) ||
                previous_offset_basis == 0.0 ||
                !std::isfinite(previous_peak_signal))
            {
                return std::nullopt;
            }
            const auto offset_scale{
                std::max({
                    std::abs(previous_model.GetOffset()),
                    std::abs(previous_peak_signal / previous_offset_basis),
                    residual_scale_floor,
                    kRobustScaleMin
                })
            };
            components.movement_penalty_sum += CalculateSquaredValue(
                    (std::log(model.GetAmplitude()) - std::log(previous_model.GetAmplitude())) /
                    kLocalFittingMovementAmplitudeLogScale) +
                CalculateSquaredValue(
                    (std::log(model.GetWidth()) - std::log(previous_model.GetWidth())) /
                    kLocalFittingMovementWidthLogScale) +
                CalculateSquaredValue(
                    (model.GetOffset() - previous_model.GetOffset()) / offset_scale);
        }
    }

    if (!std::isfinite(components.width_prior_penalty_sum) ||
        !std::isfinite(components.offset_plausibility_penalty_sum) ||
        !std::isfinite(components.movement_penalty_sum))
    {
        return std::nullopt;
    }
    return components;
}

std::optional<double> CalculateLocalFittingResidualObjective(
    const LocalFittingObjectiveSamples & objective_samples,
    double objective_scale)
{
    if (!numeric_validation::IsFinitePositive(objective_scale) ||
        objective_samples.residual_list.empty())
    {
        return std::nullopt;
    }

    double loss_sum{ 0.0 };
    for (const auto residual : objective_samples.residual_list)
    {
        const auto normalized_residual{ residual / objective_scale };
        loss_sum += algorithm::CalculateRobustLoss(
            kSecondStageRobustLossKind,
            normalized_residual,
            kRobustLossCutoffMultiplier);
    }
    const auto residual_objective{
        loss_sum / static_cast<double>(objective_samples.residual_list.size())
    };
    return std::isfinite(residual_objective) ?
        std::optional<double>{ residual_objective } :
        std::nullopt;
}

std::optional<detail::LocalFittingObjectiveBreakdown>
CalculateLocalFittingObjectiveBreakdown(
    const SecondStageLocalFittingContext & context,
    const LocalFittingObjectiveSamples & objective_samples,
    double objective_scale,
    const std::vector<LocalFittingObjectiveAtomModel> * movement_reference_model_list)
{
    const auto residual_objective{
        CalculateLocalFittingResidualObjective(
            objective_samples,
            objective_scale)
    };
    if (!residual_objective.has_value()) return std::nullopt;

    const auto parameter_penalty_components{
        CalculateLocalFittingParameterPenaltyComponents(
            context,
            objective_samples,
            objective_scale,
            movement_reference_model_list)
    };
    if (!parameter_penalty_components.has_value()) return std::nullopt;

    return detail::BuildLocalFittingMeanObjectiveBreakdown(
            *residual_objective,
            parameter_penalty_components->width_prior_penalty_sum,
            parameter_penalty_components->offset_plausibility_penalty_sum,
            parameter_penalty_components->movement_penalty_sum,
            objective_samples.active_model_list.size(),
            kLocalFittingWidthPriorPenaltyWeight,
            kLocalFittingOffsetPlausibilityPenaltyWeight,
            kLocalFittingMovementPenaltyWeight);
}

std::optional<double> CalculateLocalFittingObjective(
    const SecondStageLocalFittingContext & context,
    const LocalFittingObjectiveSamples & objective_samples,
    double objective_scale,
    const std::vector<LocalFittingObjectiveAtomModel> * movement_reference_model_list)
{
    const auto objective{
        CalculateLocalFittingObjectiveBreakdown(
            context,
            objective_samples,
            objective_scale,
            movement_reference_model_list)
    };
    return objective.has_value() ?
        std::optional<double>{ objective->total_objective } :
        std::nullopt;
}

std::optional<detail::LocalFittingObjectiveBreakdown>
EvaluateLocalFittingAuditObjective(
    const SecondStageLocalFittingContext & context,
    const GaussianFittingState & state,
    LocalFittingBestAuditState & audit_state)
{
    auto objective_samples{
        CollectLocalFittingObjectiveSamples(
            context,
            state.estimation_list,
            audit_state.sample_ref_list,
            audit_state.atom_index_list)
    };
    if (!objective_samples.has_value()) return std::nullopt;

    const auto objective_scale{
        audit_state.fixed_objective_scale.value_or(
            objective_samples->scale_sample)
    };
    const auto residual_objective{
        CalculateLocalFittingResidualObjective(
            *objective_samples,
            objective_scale)
    };
    const auto penalty_components{
        CalculateLocalFittingParameterPenaltyComponents(
            context,
            *objective_samples,
            objective_scale,
            nullptr)
    };
    if (!residual_objective.has_value() || !penalty_components.has_value())
    {
        return std::nullopt;
    }
    const auto objective{
        detail::BuildLocalFittingMeanObjectiveBreakdown(
            *residual_objective,
            penalty_components->width_prior_penalty_sum,
            penalty_components->offset_plausibility_penalty_sum,
            penalty_components->movement_penalty_sum,
            objective_samples->active_model_list.size(),
            kLocalFittingWidthPriorPenaltyWeight,
            kLocalFittingOffsetPlausibilityPenaltyWeight,
            0.0)
    };
    if (!objective.has_value()) return std::nullopt;
    if (!audit_state.fixed_objective_scale.has_value())
    {
        audit_state.fixed_objective_scale = objective_scale;
    }
    return objective;
}

void TryUpdateLocalFittingBestAuditState(
    const SecondStageLocalFittingContext & context,
    const GaussianFittingState & candidate_state,
    const std::optional<std::size_t> & accepted_iteration,
    LocalFittingBestAuditState & audit_state)
{
    const auto candidate_objective{
        EvaluateLocalFittingAuditObjective(
            context,
            candidate_state,
            audit_state)
    };
    if (!candidate_objective.has_value()) return;
    if (audit_state.best_objective.has_value() &&
        !detail::IsBetterLocalFittingAuditObjective(
            candidate_objective->total_objective,
            audit_state.best_objective->total_objective,
            kLocalFittingObjectiveTieRelativeTolerance))
    {
        return;
    }
    audit_state.best_objective = *candidate_objective;
    audit_state.best_state = candidate_state;
    audit_state.best_accepted_iteration = accepted_iteration;
}

LocalFittingBestAuditState BuildInitialLocalFittingBestAuditState(
    const SecondStageLocalFittingContext & context,
    const GaussianFittingState & initial_state)
{
    LocalFittingBestAuditState audit_state;
    audit_state.atom_index_list.reserve(context.AtomSize());
    for (std::size_t atom_index = 0; atom_index < context.AtomSize(); atom_index++)
    {
        audit_state.atom_index_list.emplace_back(atom_index);
        const auto sample_size{
            context.atom_context_list.at(atom_index).sample_entries.size()
        };
        for (std::size_t sample_index = 0; sample_index < sample_size; sample_index++)
        {
            audit_state.sample_ref_list.emplace_back(
                LocalFittingObjectiveSampleRef{ atom_index, sample_index });
        }
    }
    TryUpdateLocalFittingBestAuditState(
        context,
        initial_state,
        std::nullopt,
        audit_state);
    return audit_state;
}

void ReconcileLocalFittingBestAuditTerminalFallback(
    const SecondStageLocalFittingContext & context,
    const std::vector<LocalFittingClusterKey> & terminal_key_list,
    const GaussianFittingState & terminal_fallback_state,
    LocalFittingBestAuditState & audit_state)
{
    if (terminal_key_list.empty() || !audit_state.best_state.has_value()) return;
    if (terminal_fallback_state.estimation_list.size() != context.AtomSize() ||
        terminal_fallback_state.result_list.size() != context.AtomSize())
    {
        throw std::invalid_argument(
            "Local fitting audit terminal fallback state sizes are inconsistent.");
    }

    auto reconciled_best_state{ *audit_state.best_state };
    for (const auto & key : terminal_key_list)
    {
        for (const auto atom_index : key)
        {
            if (atom_index >= context.AtomSize())
            {
                throw std::invalid_argument(
                    "Local fitting audit terminal atom index is out of range.");
            }
            reconciled_best_state.estimation_list.at(atom_index) =
                terminal_fallback_state.estimation_list.at(atom_index);
            reconciled_best_state.result_list.at(atom_index) =
                terminal_fallback_state.result_list.at(atom_index);
        }
    }

    const auto reconciled_objective{
        EvaluateLocalFittingAuditObjective(
            context,
            reconciled_best_state,
            audit_state)
    };
    if (!reconciled_objective.has_value())
    {
        audit_state.best_objective.reset();
        audit_state.best_state.reset();
        audit_state.best_accepted_iteration.reset();
        return;
    }
    audit_state.best_objective = *reconciled_objective;
    audit_state.best_state = std::move(reconciled_best_state);
}

algorithm::ClusteredFittingQualityCandidateScore<LocalFittingObjectiveSamples>
ScoreLocalFittingClusterCandidate(
    const SecondStageLocalFittingContext & context,
    const std::vector<Eigen::VectorXd> & candidate_estimation_list,
    const std::vector<std::size_t> & active_index_list,
    const std::vector<LocalFittingObjectiveSampleRef> & sample_ref_list,
    const algorithm::ScaleReferenceTracker & objective_scale_tracker,
    algorithm::FittingQualityCandidateStats & previous_candidate_stats,
    const std::optional<LocalFittingObjectiveSamples> & previous_objective_samples,
    const std::optional<algorithm::ClusteredFittingQualityTrackedCandidate<LocalFittingObjectiveSamples>> & best_candidate,
    const algorithm::ParameterChangeStats & transformed_change_stats,
    LocalFittingObjectiveAttemptDiagnostic & diagnostic)
{
    algorithm::ClusteredFittingQualityCandidateScore<LocalFittingObjectiveSamples> score;
    score.candidate_stats.parameter_change_stats = transformed_change_stats;
    auto objective_scale{ objective_scale_tracker.GetCommittedReference() };
    score.has_objective_reference = objective_scale.has_value();
    if (best_candidate.has_value())
    {
        score.best_candidate_stats = best_candidate->candidate_stats;
    }

    if (!objective_scale.has_value()) return score;

    score.objective_samples =
        CollectLocalFittingObjectiveSamples(
            context,
            candidate_estimation_list,
            sample_ref_list,
            active_index_list);
    if (!score.objective_samples.has_value()) return score;

    const auto scale_sample{ score.objective_samples->scale_sample };
    objective_scale = objective_scale_tracker.GetProvisionalReference(scale_sample);
    score.has_objective_reference = objective_scale.has_value();
    if (!objective_scale.has_value()) return score;
    diagnostic.objective_scale = *objective_scale;
    score.objective_scale_sample = scale_sample;
    if (!previous_objective_samples.has_value()) return score;

    diagnostic.candidate_objective =
        CalculateLocalFittingObjectiveBreakdown(
            context,
            *score.objective_samples,
            *objective_scale,
            &previous_objective_samples->active_model_list);
    score.candidate_stats.quality_objective =
        diagnostic.candidate_objective.has_value() ?
            std::optional<double>{ diagnostic.candidate_objective->total_objective } :
            std::nullopt;
    score.committed_quality_objective = CalculateLocalFittingObjective(
        context,
        *score.objective_samples,
        *objective_scale,
        nullptr);
    if (diagnostic.candidate_objective.has_value())
    {
        diagnostic.status = LocalFittingObjectiveAttemptDiagnosticStatus::Scored;
        diagnostic.previous_objective =
            CalculateLocalFittingObjectiveBreakdown(
                context,
                *previous_objective_samples,
                *objective_scale,
                nullptr);
        previous_candidate_stats.quality_objective =
            diagnostic.previous_objective.has_value() ?
                std::optional<double>{ diagnostic.previous_objective->total_objective } :
                std::nullopt;
        if (best_candidate.has_value() &&
            best_candidate->objective_samples.has_value() &&
            score.best_candidate_stats.has_value())
        {
            diagnostic.best_objective =
                CalculateLocalFittingObjectiveBreakdown(
                    context,
                    *best_candidate->objective_samples,
                    *objective_scale,
                    nullptr);
            score.best_candidate_stats->quality_objective =
                diagnostic.best_objective.has_value() ?
                    std::optional<double>{ diagnostic.best_objective->total_objective } :
                    std::nullopt;
        }
    }
    if (score.has_objective_reference)
    {
        diagnostic.rejected_by_previous =
            algorithm::IsFittingQualityObjectiveDeteriorated(
                score.candidate_stats,
                previous_candidate_stats,
                kLocalFittingConvergenceObjectiveRelativeTolerance);
        if (score.best_candidate_stats.has_value())
        {
            diagnostic.rejected_by_best =
                algorithm::IsFittingQualityObjectiveDeteriorated(
                    score.candidate_stats,
                    *score.best_candidate_stats,
                    kLocalFittingConvergenceObjectiveRelativeTolerance);
        }
    }
    return score;
}

void LogSelectedAtomOffsetSummary(
    const ModelObject & model_object,
    const FitOptions & options,
    const std::string & stage_label)
{
    if (options.quiet_mode) return;

    std::vector<double> offset_list;
    const auto & atom_list{ model_object.GetSelectedAtoms() };
    offset_list.reserve(atom_list.size());
    for (const auto * atom : atom_list)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        offset_list.emplace_back(
            local_view.GetEstimateMDPDE().GetOffset());
    }

    std::ostringstream message;
    message << stage_label;
    AppendLocalFittingOffsetSummary(
        message,
        SummarizeLocalFittingOffsetValues(offset_list));
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
}

void LogAtomGroupPriorOffsetSummary(
    const ModelObject & model_object,
    const FitOptions & options,
    const std::string & stage_label)
{
    if (options.quiet_mode) return;

    const auto analysis_view{ model_object.GetAnalysisView() };
    const auto group_key_list{ analysis_view.CollectAtomGroupKeys() };
    std::vector<double> offset_list;
    offset_list.reserve(group_key_list.size());
    for (const auto group_key : group_key_list)
    {
        offset_list.emplace_back(
            analysis_view.GetAtomGroupPrior(group_key).GetOffset());
    }

    std::ostringstream message;
    message << stage_label;
    AppendLocalFittingOffsetSummary(
        message,
        SummarizeLocalFittingOffsetValues(offset_list));
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
}

std::vector<std::string> BuildGroupPriorSpotSummaryLines(const ModelObject & model_object)
{
    const auto analysis_view{ model_object.GetAnalysisView() };
    std::map<Spot, GaussianModelParameterSamples> spot_sample_map;
    for (const auto group_key : analysis_view.CollectAtomGroupKeys())
    {
        const auto & atom_list{ analysis_view.GetAtomObjectList(group_key) };
        if (atom_list.empty()) continue;

        const auto spot{ atom_list.front()->GetSpot() };
        if (std::find(
                kGroupPriorSummarySpotList.begin(),
                kGroupPriorSummarySpotList.end(),
                spot) == kGroupPriorSummarySpotList.end())
        {
            continue;
        }
        const auto & prior{ analysis_view.GetAtomGroupPrior(group_key) };
        auto & sample_list{ spot_sample_map[spot] };
        sample_list.amplitude_list.emplace_back(prior.GetAmplitude());
        sample_list.width_list.emplace_back(prior.GetWidth());
        sample_list.offset_list.emplace_back(prior.GetOffset());
    }

    std::vector<std::string> summary_lines;
    summary_lines.reserve(spot_sample_map.size());
    for (const auto spot : kGroupPriorSummarySpotList)
    {
        const auto sample_iter{ spot_sample_map.find(spot) };
        if (sample_iter == spot_sample_map.end()) continue;
        const auto & sample_list{ sample_iter->second };

        const auto amplitude_mean{
            array_helper::ComputeMean(
                sample_list.amplitude_list.data(), sample_list.amplitude_list.size())
        };
        const auto width_mean{
            array_helper::ComputeMean(
                sample_list.width_list.data(), sample_list.width_list.size())
        };
        const auto offset_mean{
            array_helper::ComputeMean(
                sample_list.offset_list.data(), sample_list.offset_list.size())
        };

        std::ostringstream stream;
        stream << "Spot::" << ChemicalDataHelper::GetLabel(spot)
            << std::fixed << std::setprecision(2)
            << " , amplitude mean = " << amplitude_mean
            << ", amplitude s.d. = "
            << array_helper::ComputeStandardDeviation(
                sample_list.amplitude_list.data(),
                sample_list.amplitude_list.size(),
                amplitude_mean)
            << ", width mean = " << width_mean
            << ", width s.d. = "
            << array_helper::ComputeStandardDeviation(
                sample_list.width_list.data(),
                sample_list.width_list.size(),
                width_mean)
            << ", offset mean = " << offset_mean
            << ", offset s.d. = "
            << array_helper::ComputeStandardDeviation(
                sample_list.offset_list.data(),
                sample_list.offset_list.size(),
                offset_mean);
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

algorithm::ParameterChange CalculateLocalFittingTransformedChange(
    const Eigen::VectorXd & current_estimation,
    const Eigen::VectorXd & previous_estimation)
{
    if (current_estimation.size() != GaussianModel3D::ParameterSize() ||
        previous_estimation.size() != GaussianModel3D::ParameterSize())
    {
        throw std::invalid_argument(
            "Local fitting transformed change parameter sizes are inconsistent.");
    }
    return detail::CalculateLocalFittingTransformedChange(
        GaussianModel3D{
            current_estimation(GaussianModel3D::AmplitudeIndex()),
            current_estimation(GaussianModel3D::WidthIndex()),
            current_estimation(GaussianModel3D::OffsetIndex())
        },
        GaussianModel3D{
            previous_estimation(GaussianModel3D::AmplitudeIndex()),
            previous_estimation(GaussianModel3D::WidthIndex()),
            previous_estimation(GaussianModel3D::OffsetIndex())
        });
}

algorithm::ParameterChange CalculateLocalFittingFreezeEvidenceChange(
    const Eigen::VectorXd & accepted_estimation,
    const Eigen::VectorXd & raw_fixed_point_estimation,
    const Eigen::VectorXd & previous_estimation)
{
    if (accepted_estimation.size() != GaussianModel3D::ParameterSize() ||
        raw_fixed_point_estimation.size() != GaussianModel3D::ParameterSize() ||
        previous_estimation.size() != GaussianModel3D::ParameterSize())
    {
        throw std::invalid_argument(
            "Local fitting freeze evidence parameter sizes are inconsistent.");
    }
    const auto build_model = [](const Eigen::VectorXd & estimation)
    {
        return GaussianModel3D{
            estimation(GaussianModel3D::AmplitudeIndex()),
            estimation(GaussianModel3D::WidthIndex()),
            estimation(GaussianModel3D::OffsetIndex())
        };
    };
    return detail::CalculateLocalFittingFreezeEvidenceChange(
        build_model(accepted_estimation),
        build_model(raw_fixed_point_estimation),
        build_model(previous_estimation));
}

std::vector<algorithm::ParameterChange> CalculateLocalFittingTransformedChanges(
    const std::vector<Eigen::VectorXd> & current_estimation_list,
    const std::vector<Eigen::VectorXd> & previous_estimation_list)
{
    if (current_estimation_list.size() != previous_estimation_list.size())
    {
        throw std::invalid_argument(
            "Local fitting transformed change input sizes are inconsistent.");
    }

    std::vector<algorithm::ParameterChange> change_list(current_estimation_list.size());
    for (size_t i = 0; i < current_estimation_list.size(); i++)
    {
        change_list.at(i) = CalculateLocalFittingTransformedChange(
            current_estimation_list.at(i),
            previous_estimation_list.at(i));
    }
    return change_list;
}

std::vector<algorithm::ParameterChange> CalculateLocalFittingFreezeEvidenceChanges(
    const std::vector<Eigen::VectorXd> & accepted_estimation_list,
    const std::vector<Eigen::VectorXd> & raw_fixed_point_estimation_list,
    const std::vector<Eigen::VectorXd> & previous_estimation_list)
{
    if (accepted_estimation_list.size() != previous_estimation_list.size() ||
        raw_fixed_point_estimation_list.size() != previous_estimation_list.size())
    {
        throw std::invalid_argument(
            "Local fitting freeze evidence input sizes are inconsistent.");
    }

    std::vector<algorithm::ParameterChange> change_list(
        previous_estimation_list.size());
    for (std::size_t i = 0; i < previous_estimation_list.size(); i++)
    {
        change_list.at(i) = CalculateLocalFittingFreezeEvidenceChange(
            accepted_estimation_list.at(i),
            raw_fixed_point_estimation_list.at(i),
            previous_estimation_list.at(i));
    }
    return change_list;
}

algorithm::ParameterChangeStats SummarizeLocalFittingTransformedChanges(
    const std::vector<Eigen::VectorXd> & current_estimation_list,
    const std::vector<Eigen::VectorXd> & previous_estimation_list,
    const std::vector<std::size_t> & index_list)
{
    if (current_estimation_list.size() != previous_estimation_list.size())
    {
        throw std::invalid_argument(
            "Local fitting transformed change input sizes are inconsistent.");
    }

    std::vector<algorithm::ParameterChange> change_list;
    change_list.reserve(index_list.size());
    for (const auto i : index_list)
    {
        if (i >= current_estimation_list.size())
        {
            throw std::invalid_argument(
                "Local fitting transformed change index is out of range.");
        }
        change_list.emplace_back(CalculateLocalFittingTransformedChange(
            current_estimation_list.at(i),
            previous_estimation_list.at(i)));
    }

    std::vector<std::size_t> local_index_list(change_list.size());
    for (std::size_t i = 0; i < local_index_list.size(); i++)
    {
        local_index_list.at(i) = i;
    }
    return algorithm::SummarizeParameterChangeStats(
        change_list,
        local_index_list,
        kLocalFittingChangePercentile);
}

bool IsLocalFittingTransformedChangeConverged(const algorithm::ParameterChangeStats & stats)
{
    for (std::size_t i = 0; i < stats.percentile_list.size(); i++)
    {
        if (stats.percentile_list.at(i) >= kLocalFittingTransformedChangeTolerance) return false;
    }
    return true;
}

std::vector<std::size_t> CollectClusterSuspiciousAtomIndexes(
    const LocalFittingClusterKey & key,
    const std::vector<std::size_t> & suspicious_atom_index_list)
{
    std::vector<std::size_t> cluster_suspicious_atom_index_list;
    for (const auto atom_index : key)
    {
        if (std::binary_search(
                suspicious_atom_index_list.begin(),
                suspicious_atom_index_list.end(),
                atom_index))
        {
            cluster_suspicious_atom_index_list.emplace_back(atom_index);
        }
    }
    return cluster_suspicious_atom_index_list;
}

std::vector<LocalFittingClusterKey> UpdatePersistentSuspiciousRollbackState(
    const std::vector<LocalFittingClusterKey> & cluster_key_list,
    const std::vector<LocalFittingClusterKey> & accepted_key_list,
    const std::vector<std::size_t> & suspicious_atom_index_list,
    const GaussianFittingState & assembled_state,
    const GaussianFittingState & previous_state,
    PersistentSuspiciousRollbackStateMap & state_by_key)
{
    PersistentSuspiciousRollbackStateMap next_state_by_key;
    std::vector<LocalFittingClusterKey> terminal_key_list;
    for (const auto & key : cluster_key_list)
    {
        if (std::find(accepted_key_list.begin(), accepted_key_list.end(), key) ==
            accepted_key_list.end())
        {
            continue;
        }

        auto cluster_suspicious_atom_index_list{
            CollectClusterSuspiciousAtomIndexes(key, suspicious_atom_index_list)
        };
        if (cluster_suspicious_atom_index_list.empty()) continue;

        const auto transformed_change_stats{
            SummarizeLocalFittingTransformedChanges(
                assembled_state.estimation_list,
                previous_state.estimation_list,
                key)
        };
        if (!IsLocalFittingTransformedChangeConverged(transformed_change_stats))
        {
            continue;
        }

        PersistentSuspiciousRollbackState next_state{
            std::move(cluster_suspicious_atom_index_list),
            1
        };
        const auto previous_iter{ state_by_key.find(key) };
        if (previous_iter != state_by_key.end() &&
            previous_iter->second.suspicious_atom_index_list ==
                next_state.suspicious_atom_index_list)
        {
            next_state.stable_iteration_count =
                previous_iter->second.stable_iteration_count + 1;
        }

        if (next_state.stable_iteration_count >=
            kPersistentSuspiciousRollbackIterationLimit)
        {
            terminal_key_list.emplace_back(key);
            continue;
        }
        next_state_by_key.emplace(key, std::move(next_state));
    }
    state_by_key = std::move(next_state_by_key);
    return terminal_key_list;
}

TerminalJointOffsetFailureMap UpdatePersistentJointOffsetFailureState(
    const std::vector<LocalFittingClusterKey> & cluster_key_list,
    const std::vector<LocalFittingClusterKey> & accepted_key_list,
    const std::vector<std::size_t> & suspicious_atom_index_list,
    const LocalFittingClusterHealthMap & health_by_key,
    const GaussianFittingState & assembled_state,
    const GaussianFittingState & previous_state,
    PersistentJointOffsetFailureStateMap & state_by_key)
{
    PersistentJointOffsetFailureStateMap next_state_by_key;
    TerminalJointOffsetFailureMap terminal_failure_by_key;
    for (const auto & key : cluster_key_list)
    {
        if (std::find(accepted_key_list.begin(), accepted_key_list.end(), key) ==
            accepted_key_list.end())
        {
            continue;
        }
        if (!CollectClusterSuspiciousAtomIndexes(
                key,
                suspicious_atom_index_list).empty())
        {
            continue;
        }

        const auto health_iter{ health_by_key.find(key) };
        if (health_iter == health_by_key.end())
        {
            throw std::invalid_argument(
                "Persistent joint-offset failure cluster health is missing.");
        }
        const auto status{ health_iter->second.joint_offset_status };
        if (!IsJointOffsetSolveHardFailure(status)) continue;

        const auto transformed_change_stats{
            SummarizeLocalFittingTransformedChanges(
                assembled_state.estimation_list,
                previous_state.estimation_list,
                key)
        };
        if (!IsLocalFittingTransformedChangeConverged(transformed_change_stats))
        {
            continue;
        }

        PersistentJointOffsetFailureState next_state{ status, 1 };
        const auto previous_iter{ state_by_key.find(key) };
        if (previous_iter != state_by_key.end() &&
            previous_iter->second.status == status)
        {
            next_state.stable_iteration_count =
                previous_iter->second.stable_iteration_count + 1;
        }

        if (next_state.stable_iteration_count >=
            kPersistentJointOffsetFailureIterationLimit)
        {
            terminal_failure_by_key.emplace(key, status);
            continue;
        }
        next_state_by_key.emplace(key, next_state);
    }
    state_by_key = std::move(next_state_by_key);
    return terminal_failure_by_key;
}

bool ContainsLocalFittingClusterKey(
    const std::vector<LocalFittingClusterKey> & key_list,
    const LocalFittingClusterKey & key)
{
    return std::find(key_list.begin(), key_list.end(), key) != key_list.end();
}

std::vector<LocalFittingClusterKey> ExcludeLocalFittingClusterKeys(
    const std::vector<LocalFittingClusterKey> & key_list,
    const std::vector<LocalFittingClusterKey> & excluded_key_list)
{
    std::vector<LocalFittingClusterKey> result;
    result.reserve(key_list.size());
    for (const auto & key : key_list)
    {
        if (!ContainsLocalFittingClusterKey(excluded_key_list, key))
        {
            result.emplace_back(key);
        }
    }
    return result;
}

std::vector<LocalFittingClusterKey> ExcludeLocalFittingClusterKeysContainingAtoms(
    const std::vector<LocalFittingClusterKey> & key_list,
    const std::vector<std::size_t> & excluded_atom_index_list)
{
    std::vector<LocalFittingClusterKey> result;
    result.reserve(key_list.size());
    for (const auto & key : key_list)
    {
        const auto contains_excluded_atom{
            std::any_of(
                excluded_atom_index_list.begin(), excluded_atom_index_list.end(),
                [&](std::size_t atom_index)
                {
                    return std::binary_search(key.begin(), key.end(), atom_index);
                })
        };
        if (!contains_excluded_atom)
        {
            result.emplace_back(key);
        }
    }
    return result;
}

std::vector<LocalFittingClusterKey> CollectLocalFittingClusterKeysWithRegimeSignatures(
    const std::vector<LocalFittingClusterKey> & key_list,
    const detail::LocalFittingAndersonRegimeSignatureMap & signature_by_key)
{
    std::vector<LocalFittingClusterKey> result;
    result.reserve(key_list.size());
    for (const auto & key : key_list)
    {
        if (signature_by_key.find(key) != signature_by_key.end())
        {
            result.emplace_back(key);
        }
    }
    return result;
}

void ApplyTerminalFallbackClusters(
    const std::vector<LocalFittingClusterKey> & terminal_key_list,
    const GaussianFittingState & previous_state,
    std::vector<char> & terminal_atom_mask,
    GaussianFittingState & assembled_state)
{
    if (previous_state.estimation_list.size() != terminal_atom_mask.size() ||
        previous_state.result_list.size() != terminal_atom_mask.size() ||
        assembled_state.estimation_list.size() != terminal_atom_mask.size() ||
        assembled_state.result_list.size() != terminal_atom_mask.size())
    {
        throw std::invalid_argument(
            "Local fitting terminal fallback state sizes are inconsistent.");
    }
    for (const auto & key : terminal_key_list)
    {
        for (const auto atom_index : key)
        {
            if (atom_index >= terminal_atom_mask.size())
            {
                throw std::invalid_argument(
                    "Local fitting terminal fallback atom index is out of range.");
            }
            terminal_atom_mask.at(atom_index) = 1;
            assembled_state.estimation_list.at(atom_index) =
                previous_state.estimation_list.at(atom_index);
            assembled_state.result_list.at(atom_index) =
                previous_state.result_list.at(atom_index);
        }
    }
}

std::vector<std::size_t> BuildEligibleLocalFittingActiveIndexList(
    const algorithm::ConvergenceFreezeTracker & freeze_tracker,
    const std::vector<char> & terminal_atom_mask)
{
    auto active_index_list{ freeze_tracker.BuildActiveIndexList() };
    active_index_list.erase(
        std::remove_if(
            active_index_list.begin(),
            active_index_list.end(),
            [&](std::size_t atom_index)
            {
                if (atom_index >= terminal_atom_mask.size())
                {
                    throw std::invalid_argument(
                        "Local fitting terminal fallback active index is out of range.");
                }
                return terminal_atom_mask.at(atom_index) != 0;
            }),
        active_index_list.end());
    return active_index_list;
}

std::optional<GaussianFittingState> BuildLocalFittingCandidateState(
    const GaussianFittingState & previous_state,
    const std::vector<Eigen::VectorXd> & candidate_estimation_list,
    const std::vector<LocalGaussianResult> & uncertainty_result_list,
    const std::vector<std::size_t> & active_index_list,
    double damping,
    LocalFittingCandidateBuildFailure * build_failure)
{
    if (!std::isfinite(damping) || damping <= 0.0 || damping > 1.0)
    {
        throw std::invalid_argument("Local fitting candidate damping is out of range.");
    }
    if (previous_state.estimation_list.size() != previous_state.result_list.size() ||
        candidate_estimation_list.size() != previous_state.estimation_list.size() ||
        uncertainty_result_list.size() != previous_state.result_list.size())
    {
        throw std::invalid_argument("Local fitting candidate input sizes are inconsistent.");
    }

    auto candidate_state{ previous_state };
    for (const auto active_index : active_index_list)
    {
        if (active_index >= candidate_state.estimation_list.size())
        {
            throw std::invalid_argument("Local fitting candidate active index is out of range.");
        }
        const auto damped_estimation{
            (previous_state.estimation_list.at(active_index) +
                damping * (
                    candidate_estimation_list.at(active_index) -
                    previous_state.estimation_list.at(active_index))).eval()
        };
        const auto set_build_failure = [&](LocalFittingCandidateBuildFailureReason reason)
        {
            if (build_failure == nullptr) return;
            *build_failure = LocalFittingCandidateBuildFailure{
                active_index,
                reason,
                damped_estimation
            };
        };
        if (damped_estimation.size() != GaussianModel3D::ParameterSize())
        {
            set_build_failure(LocalFittingCandidateBuildFailureReason::ParameterSize);
            return std::nullopt;
        }
        if (!damped_estimation.allFinite())
        {
            set_build_failure(LocalFittingCandidateBuildFailureReason::NonFiniteParameter);
            return std::nullopt;
        }
        if (damped_estimation(GaussianModel3D::WidthIndex()) <= 0.0)
        {
            set_build_failure(LocalFittingCandidateBuildFailureReason::NonPositiveWidth);
            return std::nullopt;
        }

        const auto damped_model{ GaussianModel3D::FromVector(damped_estimation) };
        auto & result{ candidate_state.result_list.at(active_index) };
        result.mdpde = GaussianModel3DWithUncertainty{
            damped_model,
            uncertainty_result_list.at(active_index).mdpde.GetStandardDeviationModel()
        };
        candidate_state.estimation_list.at(active_index) = damped_estimation;
    }
    return std::move(candidate_state);
}

const char * GetLocalFittingCandidateText(const LocalFittingCandidateAttempt & attempt)
{
    if (attempt.kind == LocalFittingCandidateKind::FixedPoint) return "damped-fixed-point";
    return attempt.damping == 1.0 ? "aa" : "damped-aa";
}

algorithm::ClusteredFittingQualityInitialState<LocalFittingObjectiveSamples>
BuildInitialLocalFittingClusterQualityState(
    const SecondStageLocalFittingContext & context,
    const GaussianFittingState & previous_state,
    const LocalFittingClusterKey & key,
    const LocalFittingClusterWork & cluster)
{
    auto initial_objective_samples{
        CollectLocalFittingObjectiveSamples(
            context,
            previous_state.estimation_list,
            cluster.objective_sample_ref_list,
            key)
    };
    std::optional<double> initial_objective;
    if (initial_objective_samples.has_value())
    {
        initial_objective = CalculateLocalFittingObjective(
            context,
            *initial_objective_samples,
            initial_objective_samples->scale_sample,
            nullptr);
    }
    algorithm::FittingQualityCandidateStats initial_candidate_stats{
        initial_objective,
        algorithm::ParameterChangeStats{
            std::vector<double>(
                static_cast<std::size_t>(GaussianModel3D::ParameterSize()),
                0.0)
        }
    };
    return algorithm::ClusteredFittingQualityInitialState<LocalFittingObjectiveSamples>{
        initial_objective_samples.has_value() ?
            std::optional<double>{ initial_objective_samples->scale_sample } :
            std::nullopt,
        algorithm::ClusteredFittingQualityTrackedCandidate<LocalFittingObjectiveSamples>{
            std::move(initial_candidate_stats),
            std::move(initial_objective_samples)
        }
    };
}

LocalFittingCandidateSelection SelectLocalFittingClusterCandidates(
    const SecondStageLocalFittingContext & context,
    LocalFittingClusterMap cluster_map,
    const std::vector<LocalFittingClusterKey> & key_list,
    const GaussianFittingState & previous_state,
    const GaussianFittingState & raw_state,
    algorithm::ClusteredAndersonAccelerationHistorySet & acceleration_history,
    algorithm::ClusteredFittingQualityStateSet<LocalFittingObjectiveSamples> & cluster_quality_state)
{
    const auto localized_anderson_candidate{
        acceleration_history.BuildCandidate(
            key_list,
            previous_state.estimation_list,
            raw_state.estimation_list)
    };
    LocalFittingCandidateSelection selection;
    selection.assembled_state = previous_state;

    const auto run_attempt_group = [&](
        LocalFittingCandidateKind candidate_kind,
        const std::vector<LocalFittingClusterKey> & candidate_key_list,
        const std::vector<Eigen::VectorXd> & candidate_estimation_list,
        const std::vector<LocalGaussianResult> & uncertainty_result_list)
    {
        for (const auto damping : kLocalFittingAccelerationDampingList)
        {
            const auto candidate_attempt{
                LocalFittingCandidateAttempt{
                    candidate_kind,
                    damping
                }
            };
            for (const auto & key : candidate_key_list)
            {
                auto & cluster{ cluster_map.at(key) };
                if (cluster.accepted_kind.has_value()) continue;
                LocalFittingObjectiveAttemptDiagnostic attempt_diagnostic;
                attempt_diagnostic.attempt = candidate_attempt;
                LocalFittingCandidateBuildFailure build_failure;

                auto attempt_state{
                    BuildLocalFittingCandidateState(
                        previous_state,
                        candidate_estimation_list,
                        uncertainty_result_list,
                        key,
                        damping,
                        &build_failure)
                };
                if (!attempt_state.has_value())
                {
                    attempt_diagnostic.status =
                        LocalFittingObjectiveAttemptDiagnosticStatus::InvalidModel;
                    attempt_diagnostic.build_failure = std::move(build_failure);
                    cluster.objective_attempt_list.emplace_back(
                        std::move(attempt_diagnostic));
                    if (candidate_kind == LocalFittingCandidateKind::Anderson)
                    {
                        cluster.rejected_anderson = true;
                    }
                    continue;
                }

                auto transformed_change_stats{
                    SummarizeLocalFittingTransformedChanges(
                        attempt_state->estimation_list,
                        previous_state.estimation_list,
                        key)
                };
                const auto accepted{
                    cluster_quality_state.TryCommitCandidate(
                        key,
                        [&](const algorithm::ScaleReferenceTracker & objective_scale_tracker,
                            algorithm::FittingQualityCandidateStats & previous_candidate_stats,
                            const std::optional<LocalFittingObjectiveSamples> & previous_objective_samples,
                            const std::optional<
                                algorithm::ClusteredFittingQualityTrackedCandidate<LocalFittingObjectiveSamples>> & best_candidate)
                        {
                            return ScoreLocalFittingClusterCandidate(
                                context,
                                attempt_state->estimation_list,
                                key,
                                cluster.objective_sample_ref_list,
                                objective_scale_tracker,
                                previous_candidate_stats,
                                previous_objective_samples,
                                best_candidate,
                                transformed_change_stats,
                                attempt_diagnostic);
                        })
                };
                cluster.objective_attempt_list.emplace_back(
                    std::move(attempt_diagnostic));
                if (accepted)
                {
                    for (const auto active_index : key)
                    {
                        selection.assembled_state.result_list.at(active_index) =
                            attempt_state->result_list.at(active_index);
                        selection.assembled_state.estimation_list.at(active_index) =
                            attempt_state->estimation_list.at(active_index);
                    }
                    cluster.accepted_kind = candidate_kind;
                    if (!selection.accepted_candidate_attempt.has_value())
                    {
                        selection.accepted_candidate_attempt = candidate_attempt;
                    }
                    continue;
                }

                selection.has_objective_backtracking_rejection = true;
                if (candidate_kind == LocalFittingCandidateKind::Anderson)
                {
                    cluster.rejected_anderson = true;
                }
            }
        }
    };

    if (localized_anderson_candidate.has_value())
    {
        run_attempt_group(
            LocalFittingCandidateKind::Anderson,
            localized_anderson_candidate->used_cluster_key_list,
            localized_anderson_candidate->state_list,
            previous_state.result_list);
    }
    std::vector<LocalFittingClusterKey> anderson_failure_keys;
    for (auto & [key, cluster] : cluster_map)
    {
        if (cluster.rejected_anderson)
        {
            anderson_failure_keys.emplace_back(key);
        }
    }
    if (!anderson_failure_keys.empty())
    {
        acceleration_history.ClearAndSuppress(anderson_failure_keys);
    }
    run_attempt_group(
        LocalFittingCandidateKind::FixedPoint,
        key_list,
        raw_state.estimation_list,
        raw_state.result_list);

    std::vector<LocalFittingClusterKey> fixed_point_progress_keys;
    for (const auto & [key, cluster] : cluster_map)
    {
        if (cluster.accepted_kind.has_value())
        {
            selection.accepted_key_list.emplace_back(key);
            if (*cluster.accepted_kind == LocalFittingCandidateKind::FixedPoint)
            {
                fixed_point_progress_keys.emplace_back(key);
            }
        }
        else
        {
            selection.rejected_key_list.emplace_back(key);
            selection.rejected_cluster_diagnostic_list.emplace_back(
                LocalFittingRejectedClusterDiagnostic{
                    key,
                    cluster.objective_attempt_list
                });
        }
    }
    acceleration_history.ReleaseSuppression(fixed_point_progress_keys);
    return selection;
}

LocalFittingClusterSelectionSummary SummarizeLocalFittingClusterSelection(
    const LocalFittingCandidateSelection & selection)
{
    LocalFittingClusterSelectionSummary summary;
    summary.accepted_cluster_count = selection.accepted_key_list.size();
    summary.rejected_cluster_count = selection.rejected_key_list.size();
    for (const auto & key : selection.accepted_key_list)
    {
        summary.accepted_atom_count += key.size();
    }
    for (const auto & key : selection.rejected_key_list)
    {
        summary.rejected_atom_count += key.size();
    }
    return summary;
}

bool IsLocalGaussianRefitHealthEligible(const LocalGaussianResult & result)
{
    if (!result.fit_result.has_value()) return false;
    const auto & model{ result.mdpde.GetModel() };
    return detail::IsLocalGaussianRefitStatusHealthEligible(
            result.fit_result->status) &&
        detail::IsValidSecondStageGaussianModel(model);
}

std::optional<LocalAtomRefitResult> FitAtomWithJointOffsetFallback(
    const SecondStageLocalFittingContext & context,
    std::size_t atom_index,
    const LocalGaussianResult & previous_result,
    const FittedGaussianSnapshot & offset_snapshot,
    const FitOptions & options)
{
    auto sample_entries{ BuildSecondStageAdjustedSamples(context, atom_index, offset_snapshot) };
    const auto & offset_model{ offset_snapshot.at(atom_index) };
    const auto & previous_model{ previous_result.mdpde.GetModel() };
    const auto is_acceptable = [&](const GaussianModel3D & model)
    {
        return detail::IsValidSecondStageGaussianModel(model) &&
            CanBuildFiniteZeroOffsetSamples(sample_entries, model) &&
            !IsSuspiciousJointOffset(
                sample_entries,
                previous_model,
                model,
                options);
    };
    try
    {
        auto candidate_result{
            EstimateLocalGaussian(
                sample_entries,
                context.atom_context_list.at(atom_index).alpha_r,
                options,
                offset_model)
        };
        const auto candidate_model{ candidate_result.mdpde.GetModel() };
        if (is_acceptable(candidate_model))
        {
            const auto is_health_eligible{
                IsLocalGaussianRefitHealthEligible(candidate_result)
            };
            return LocalAtomRefitResult{
                std::move(candidate_result),
                is_health_eligible
            };
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
    if (!is_acceptable(result.mdpde.GetModel()))
    {
        return std::nullopt;
    }
    return LocalAtomRefitResult{ std::move(result), false };
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
            if (!thaw_hysteresis_tracker.ShouldThaw(
                    neighbor_index,
                    active_change,
                    kLocalFittingDependencyThawChangeThreshold))
            {
                continue;
            }
            if (!thaw_hysteresis_tracker.CanDependencyThaw(neighbor_index)) continue;
            if (freeze_tracker.Thaw(neighbor_index))
            {
                thaw_hysteresis_tracker.RecordDependencyThaw(neighbor_index);
            }
        }
    }
}

std::vector<std::size_t> ExpandSuspiciousOffsetClusters(
    const ActiveCouplingGraph & active_coupling_graph,
    const std::vector<std::size_t> & suspicious_offset_seed_position_list,
    std::vector<char> & suspicious_offset_mask)
{
    if (!active_coupling_graph.empty() &&
        active_coupling_graph.size() != suspicious_offset_mask.size())
    {
        throw std::invalid_argument("Suspicious offset cluster graph size is inconsistent.");
    }

    std::vector<char> visited(suspicious_offset_mask.size(), 0);
    std::vector<std::pair<std::size_t, std::size_t>> queue;
    queue.reserve(suspicious_offset_mask.size());
    for (const auto seed_position : suspicious_offset_seed_position_list)
    {
        if (seed_position >= suspicious_offset_mask.size())
        {
            throw std::invalid_argument("Suspicious offset seed position is out of range.");
        }
        if (visited.at(seed_position) != 0) continue;
        visited.at(seed_position) = 1;
        queue.emplace_back(seed_position, 0);
    }

    std::vector<std::size_t> added_position_list;
    for (std::size_t queue_index = 0; queue_index < queue.size(); queue_index++)
    {
        const auto [active_position, depth]{ queue.at(queue_index) };
        if (suspicious_offset_mask.at(active_position) == 0)
        {
            suspicious_offset_mask.at(active_position) = 1;
            added_position_list.emplace_back(active_position);
        }
        if (active_coupling_graph.empty() || depth >= kSuspiciousOffsetClusterMaxDepth) continue;

        for (const auto & edge : active_coupling_graph.at(active_position))
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
    return added_position_list;
}

void RollBackSuspiciousOffsetClusters(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const GaussianFittingState & previous_state,
    const std::vector<std::size_t> & suspicious_active_position_list,
    FittedGaussianSnapshot & current_snapshot,
    GaussianFittingState & iteration_state)
{
    for (const auto active_position : suspicious_active_position_list)
    {
        if (active_position >= active_index_list.size())
        {
            throw std::invalid_argument("Suspicious offset rollback position is out of range.");
        }

        const auto state_index{ active_index_list.at(active_position) };
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

LocalFittingIterationResult RunLocalFittingIteration(
    const SecondStageLocalFittingContext & context,
    const std::vector<std::size_t> & active_index_list,
    const std::vector<LocalFittingClusterKey> & cluster_key_list,
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
        BuildFittedGaussianSnapshot(previous_state.estimation_list)
    };
    auto clustered_joint_offset_result{
        EstimateClusteredJointOffsets(
            context,
            active_index_list,
            cluster_key_list,
            current_snapshot,
            ridge_ratio,
            ridge_multiplier_list)
    };
    ApplyJointOffsetsToSnapshot(
        active_index_list,
        clustered_joint_offset_result.offset,
        current_snapshot);
    auto iteration_state{ previous_state };
    std::vector<char> suspicious_offset_mask(active_index_list.size(), 0);
    std::vector<std::size_t> pre_refit_suspicious_seed_position_list;

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
                offset_model,
                options))
        {
            continue;
        }
        pre_refit_suspicious_seed_position_list.emplace_back(i);
    }
    const auto pre_refit_suspicious_position_list{
        ExpandSuspiciousOffsetClusters(
            clustered_joint_offset_result.active_coupling_graph,
            pre_refit_suspicious_seed_position_list,
            suspicious_offset_mask)
    };
    RollBackSuspiciousOffsetClusters(
        context,
        active_index_list,
        previous_state,
        pre_refit_suspicious_position_list,
        current_snapshot,
        iteration_state);

    const auto refit_snapshot{ current_snapshot };
    std::vector<std::size_t> post_refit_suspicious_seed_position_list;
    for (size_t i = 0; i < active_index_list.size(); i++)
    {
        if (suspicious_offset_mask.at(i) != 0) continue;

        const auto state_index{ active_index_list.at(i) };
        auto refit_result{
            FitAtomWithJointOffsetFallback(
                context,
                state_index,
                previous_state.result_list.at(state_index),
                refit_snapshot,
                options)
        };
        if (!refit_result.has_value())
        {
            RecordUnhealthyLocalRefit(
                clustered_joint_offset_result.health_by_key,
                state_index);
            post_refit_suspicious_seed_position_list.emplace_back(i);
            continue;
        }
        if (!refit_result->is_health_eligible)
        {
            RecordUnhealthyLocalRefit(
                clustered_joint_offset_result.health_by_key,
                state_index);
        }
        auto result{ std::move(refit_result->result) };
        const auto fitted_model{ result.mdpde.GetModel() };
        iteration_state.estimation_list.at(state_index) = fitted_model.ToVector();
        iteration_state.result_list.at(state_index) = std::move(result);
    }
    const auto post_refit_suspicious_position_list{
        detail::ExpandPostRefitRollbackClusters(
            active_index_list,
            cluster_key_list,
            post_refit_suspicious_seed_position_list,
            suspicious_offset_mask)
    };
    RollBackSuspiciousOffsetClusters(
        context,
        active_index_list,
        previous_state,
        post_refit_suspicious_position_list,
        current_snapshot,
        iteration_state);

    LocalFittingIterationResult iteration_result;
    iteration_result.state = std::move(iteration_state);
    iteration_result.health_by_key =
        std::move(clustered_joint_offset_result.health_by_key);
    iteration_result.anderson_regime_signature_by_key =
        std::move(clustered_joint_offset_result.anderson_regime_signature_by_key);
    for (std::size_t i = 0; i < suspicious_offset_mask.size(); i++)
    {
        if (suspicious_offset_mask.at(i) == 0) continue;
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

LocalFittingOffsetStats SummarizeLocalFittingOffsetValues(
    const std::vector<double> & offset_list)
{
    LocalFittingOffsetStats stats;
    stats.atom_count = offset_list.size();
    std::vector<double> absolute_offset_list;
    absolute_offset_list.reserve(offset_list.size());
    for (const auto offset : offset_list)
    {
        if (!std::isfinite(offset)) continue;
        stats.finite_count++;
        if (offset == 0.0) stats.exact_zero_count++;
        absolute_offset_list.emplace_back(std::abs(offset));
    }
    if (absolute_offset_list.empty()) return stats;

    stats.median_absolute_offset =
        array_helper::ComputeMedian(absolute_offset_list);
    stats.percentile_absolute_offset =
        array_helper::ComputePercentile(
            absolute_offset_list,
            kLocalFittingChangePercentile);
    stats.maximum_absolute_offset =
        *std::max_element(
            absolute_offset_list.begin(),
            absolute_offset_list.end());
    return stats;
}

LocalFittingOffsetStats SummarizeLocalFittingOffsets(
    const std::vector<Eigen::VectorXd> & estimation_list)
{
    std::vector<double> offset_list;
    offset_list.reserve(estimation_list.size());
    for (const auto & estimation : estimation_list)
    {
        GaussianModel3D::RequireParameterVector(
            estimation,
            "Local fitting offset summary model");
        offset_list.emplace_back(
            estimation(GaussianModel3D::OffsetIndex()));
    }
    return SummarizeLocalFittingOffsetValues(offset_list);
}

void AppendLocalFittingOffsetSummary(
    std::ostringstream & stream,
    const LocalFittingOffsetStats & stats)
{
    stream
        << std::scientific << std::setprecision(4)
        << "; offsets finite/exact-zero = "
        << stats.finite_count << "/" << stats.exact_zero_count
        << " of " << stats.atom_count
        << ", |C| median/p99/max = "
        << stats.median_absolute_offset << "/"
        << stats.percentile_absolute_offset << "/"
        << stats.maximum_absolute_offset;
}

void AppendLocalFittingAuditSummary(
    std::ostringstream & stream,
    const LocalFittingBestAuditState & audit_state)
{
    if (!audit_state.best_objective.has_value()) return;
    const auto & objective{ *audit_state.best_objective };
    stream
        << "; audit best source = ";
    if (audit_state.best_accepted_iteration.has_value())
    {
        stream << "accepted iteration " << *audit_state.best_accepted_iteration;
    }
    else
    {
        stream << "initial";
    }
    stream
        << std::scientific << std::setprecision(4)
        << ", fixed audit objective residual/width/offset/total = "
        << objective.residual_objective << "/"
        << objective.width_prior_penalty << "/"
        << objective.offset_plausibility_penalty << "/"
        << objective.total_objective;
}

void LogLocalFittingAllAtomsFrozen(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const LocalFittingOffsetStats & offset_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "Converged after " << accepted_iteration_count
        << " iterations because all local atoms are frozen";
    AppendLocalFittingOffsetSummary(message, offset_stats);
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
}

void AppendLocalFittingTerminalSummary(
    std::ostringstream & stream,
    const LocalFittingTerminalSummary & summary)
{
    if (summary.suspicious_atom_count > 0)
    {
        stream << "; terminal suspicious rollback fallback clusters/atoms = "
            << summary.suspicious_cluster_count
            << "/" << summary.suspicious_atom_count;
    }
    if (summary.joint_offset_failure_atom_count > 0)
    {
        stream << "; terminal joint-offset failure fallback clusters/atoms = "
            << summary.joint_offset_failure_cluster_count
            << "/" << summary.joint_offset_failure_atom_count;
        if (!summary.joint_offset_failure_status_count.empty())
        {
            stream << ", statuses = ";
            bool is_first_status{ true };
            for (const auto & [status, count] :
                summary.joint_offset_failure_status_count)
            {
                if (!is_first_status) stream << ",";
                stream << GetJointOffsetSolveStatusText(status) << ":" << count;
                is_first_status = false;
            }
        }
    }
}

void LogLocalFittingTerminalFallback(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const LocalFittingTerminalSummary & terminal_summary,
    const LocalFittingOffsetStats & offset_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message
        << "Completed local fitting after " << accepted_iteration_count
        << " accepted iterations with last validated states retained";
    AppendLocalFittingTerminalSummary(warning_message, terminal_summary);
    AppendLocalFittingOffsetSummary(warning_message, offset_stats);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void AppendLocalFittingClusterHealthSummary(
    std::ostringstream & stream,
    const LocalFittingClusterHealthSummary & summary)
{
    if (!summary.joint_status_cluster_count.empty())
    {
        stream << ", joint-offset statuses clusters/atoms = ";
        bool is_first_status{ true };
        for (const auto & [status, cluster_count] :
            summary.joint_status_cluster_count)
        {
            const auto atom_iter{ summary.joint_status_atom_count.find(status) };
            if (atom_iter == summary.joint_status_atom_count.end())
            {
                throw std::logic_error("Joint offset status atom count is missing.");
            }
            if (!is_first_status) stream << ",";
            stream
                << GetJointOffsetSolveStatusText(status) << ":"
                << cluster_count << "/" << atom_iter->second;
            is_first_status = false;
        }
    }

    if (summary.unhealthy_cluster_count == 0) return;

    stream
        << ", health-unhealthy clusters/atoms = "
        << summary.unhealthy_cluster_count << "/" << summary.unhealthy_atom_count;
    if (!summary.unhealthy_joint_status_count.empty())
    {
        stream << ", joint-offset = ";
        bool is_first_status{ true };
        for (const auto & [status, count] : summary.unhealthy_joint_status_count)
        {
            if (!is_first_status) stream << ",";
            stream << GetJointOffsetSolveStatusText(status) << ":" << count;
            is_first_status = false;
        }
    }
    if (summary.unhealthy_refit_cluster_count > 0)
    {
        stream
            << ", local-refit-fallback clusters/atoms = "
            << summary.unhealthy_refit_cluster_count
            << "/" << summary.unhealthy_refit_atom_count;
    }
}

void AppendLocalFittingClusterSelectionSummary(
    std::ostringstream & stream,
    const LocalFittingClusterSelectionSummary & summary)
{
    stream
        << ", objective accepted/rejected clusters = "
        << summary.accepted_cluster_count << "/"
        << summary.rejected_cluster_count
        << ", atoms = "
        << summary.accepted_atom_count << "/"
        << summary.rejected_atom_count;
}

void AppendLocalFittingObjectiveBreakdown(
    std::ostringstream & stream,
    const std::optional<detail::LocalFittingObjectiveBreakdown> & breakdown)
{
    if (!breakdown.has_value())
    {
        stream << "unavailable";
        return;
    }
    stream
        << breakdown->residual_objective << "/"
        << breakdown->width_prior_penalty << "/"
        << breakdown->offset_plausibility_penalty << "/"
        << breakdown->movement_penalty << "/"
        << breakdown->total_objective;
}

void LogRejectedLocalFittingClusterDiagnostics(
    const FitOptions & options,
    const std::vector<LocalFittingRejectedClusterDiagnostic> & diagnostic_list)
{
    if (options.quiet_mode || Logger::GetLogLevel() < LogLevel::Debug ||
        diagnostic_list.empty())
    {
        return;
    }

    Logger::FinishProgressLine();
    for (const auto & cluster_diagnostic : diagnostic_list)
    {
        if (cluster_diagnostic.key.empty())
        {
            throw std::logic_error(
                "Rejected local fitting cluster diagnostic key is empty.");
        }
        std::ostringstream header;
        header
            << "Rejected local fitting cluster objective diagnostics: atoms = "
            << cluster_diagnostic.key.size()
            << ", key first/last = "
            << cluster_diagnostic.key.front() << "/"
            << cluster_diagnostic.key.back()
            << ", breakdown order = residual/width/offset/movement/total";
        Logger::Log(LogLevel::Debug, header.str());

        for (const auto & diagnostic : cluster_diagnostic.attempt_list)
        {
            std::ostringstream message;
            message
                << std::scientific << std::setprecision(6)
                << "  kind = "
                << (diagnostic.attempt.kind == LocalFittingCandidateKind::Anderson ?
                    "anderson" : "fixed-point")
                << ", damping = " << diagnostic.attempt.damping;
            if (diagnostic.attempt.kind == LocalFittingCandidateKind::FixedPoint &&
                diagnostic.attempt.damping == 1.0)
            {
                message << " (raw)";
            }

            if (diagnostic.status ==
                LocalFittingObjectiveAttemptDiagnosticStatus::InvalidModel)
            {
                message << ", status = invalid-model";
                if (!diagnostic.build_failure.has_value())
                {
                    throw std::logic_error(
                        "Invalid local fitting candidate diagnostic has no build failure.");
                }
                const auto & failure{ *diagnostic.build_failure };
                message << ", first invalid atom = " << failure.atom_index
                    << ", reason = ";
                switch (failure.reason)
                {
                case LocalFittingCandidateBuildFailureReason::ParameterSize:
                    message << "parameter-size";
                    break;
                case LocalFittingCandidateBuildFailureReason::NonFiniteParameter:
                    message << "non-finite-parameter";
                    break;
                case LocalFittingCandidateBuildFailureReason::NonPositiveWidth:
                    message << "non-positive-width";
                    break;
                }
                if (failure.estimation.size() == GaussianModel3D::ParameterSize())
                {
                    message
                        << ", A/B/C = "
                        << failure.estimation(GaussianModel3D::AmplitudeIndex()) << "/"
                        << failure.estimation(GaussianModel3D::WidthIndex()) << "/"
                        << failure.estimation(GaussianModel3D::OffsetIndex());
                }
                else
                {
                    message << ", parameter size = " << failure.estimation.size();
                }
                Logger::Log(LogLevel::Debug, message.str());
                continue;
            }

            message << ", objective scale = ";
            if (diagnostic.objective_scale.has_value())
            {
                message << *diagnostic.objective_scale;
            }
            else
            {
                message << "unavailable";
            }
            message << ", candidate = ";
            AppendLocalFittingObjectiveBreakdown(
                message,
                diagnostic.candidate_objective);
            message << ", previous = ";
            AppendLocalFittingObjectiveBreakdown(
                message,
                diagnostic.previous_objective);
            message << ", best = ";
            AppendLocalFittingObjectiveBreakdown(
                message,
                diagnostic.best_objective);
            message << ", rejected-by = ";
            if (diagnostic.status ==
                LocalFittingObjectiveAttemptDiagnosticStatus::ObjectiveUnavailable)
            {
                message << "objective-unavailable";
            }
            else if (diagnostic.rejected_by_previous && diagnostic.rejected_by_best)
            {
                message << "previous+best";
            }
            else if (diagnostic.rejected_by_previous)
            {
                message << "previous";
            }
            else if (diagnostic.rejected_by_best)
            {
                message << "best";
            }
            else
            {
                message << "none";
            }
            Logger::Log(LogLevel::Debug, message.str());
        }
    }
}

void LogLocalFittingBacktrackingRetry(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    double ridge_ratio,
    bool uses_cluster_local_objective_ridge,
    const LocalFittingClusterHealthSummary & health_summary,
    const LocalFittingClusterSelectionSummary & selection_summary,
    double raw_offset_change_percentile,
    const LocalFittingOffsetStats & raw_offset_stats)
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
    progress_message
        << std::scientific << std::setprecision(4)
        << ", offset dQ_C p99 raw = " << raw_offset_change_percentile
        << ", exact-zero offsets raw = " << raw_offset_stats.exact_zero_count;
    AppendLocalFittingClusterSelectionSummary(
        progress_message,
        selection_summary);
    AppendLocalFittingClusterHealthSummary(progress_message, health_summary);
    Logger::ProgressLine(progress_message.str());
}

void LogLocalFittingBacktrackingStop(
    const FitOptions & options,
    LocalFittingBacktrackingStopReason reason,
    const LocalFittingBestAuditState * applied_audit_state,
    const LocalFittingTerminalSummary & terminal_summary,
    const LocalFittingOffsetStats & applied_offset_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message
        << "Stopped local fitting because objective backtracking rejected all "
        << "acceleration and fixed-point attempts "
        << (reason == LocalFittingBacktrackingStopReason::MaximumGlobalRidge ?
            "at the maximum joint-offset ridge ratio" : "at the maximum iteration limit");
    AppendLocalFittingTerminalSummary(warning_message, terminal_summary);
    if (applied_audit_state != nullptr)
    {
        warning_message << "; applying best validated audit state";
        AppendLocalFittingAuditSummary(
            warning_message,
            *applied_audit_state);
    }
    else
    {
        warning_message << "; applying previous state";
        if (reason == LocalFittingBacktrackingStopReason::MaximumIterationLimit)
        {
            warning_message
                << " because no finite fixed audit state is available";
        }
    }
    AppendLocalFittingOffsetSummary(warning_message, applied_offset_stats);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
}

void LogLocalFittingProgress(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const LocalFittingCandidateAttempt & current_candidate_attempt,
    const algorithm::ConvergenceFreezeTracker & freeze_tracker,
    const LocalFittingClusterHealthSummary & health_summary,
    const LocalFittingTerminalSummary & terminal_summary,
    const LocalFittingClusterSelectionSummary & selection_summary,
    double raw_offset_change_percentile,
    double accepted_offset_change_percentile,
    const LocalFittingOffsetStats & raw_offset_stats,
    const LocalFittingOffsetStats & accepted_offset_stats)
{
    if (options.quiet_mode) return;

    if (terminal_summary.AtomCount() > freeze_tracker.GetActiveCount())
    {
        throw std::invalid_argument(
            "Local fitting terminal atom count exceeds the active atom count.");
    }
    const auto effective_active_count{
        freeze_tracker.GetActiveCount() - terminal_summary.AtomCount()
    };
    std::ostringstream progress_message;
    progress_message << "Iter. " << accepted_iteration_count
        << '/' << kLocalFittingMaximumIterations
        << std::fixed << std::setprecision(4)
        << ", acceleration = "<< GetLocalFittingCandidateText(current_candidate_attempt)
        << ", damping = "<< current_candidate_attempt.damping
        << ", active/frozen atoms = "<< effective_active_count
        << "/" << freeze_tracker.GetFrozenCount();
    progress_message
        << std::scientific << std::setprecision(4)
        << ", offset dQ_C p99 raw/accepted = "
        << raw_offset_change_percentile << "/"
        << accepted_offset_change_percentile
        << ", exact-zero offsets raw/accepted = "
        << raw_offset_stats.exact_zero_count << "/"
        << accepted_offset_stats.exact_zero_count;
    AppendLocalFittingClusterSelectionSummary(
        progress_message,
        selection_summary);
    if (terminal_summary.suspicious_atom_count > 0)
    {
        progress_message
            << ", terminal-suspicious atoms = "
            << terminal_summary.suspicious_atom_count;
    }
    if (terminal_summary.joint_offset_failure_atom_count > 0)
    {
        progress_message
            << ", terminal-joint-offset-failure atoms = "
            << terminal_summary.joint_offset_failure_atom_count;
    }
    AppendLocalFittingClusterHealthSummary(progress_message, health_summary);
    Logger::ProgressLine(progress_message.str());
}

void LogLocalFittingConverged(
    const FitOptions & options,
    std::size_t accepted_iteration_count,
    const algorithm::ParameterChangeStats & transformed_change_stats,
    const LocalFittingOffsetStats & offset_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream message;
    message
        << "Converged after " << accepted_iteration_count
        << " iterations with percentile log-peak-height change = "
        << transformed_change_stats.percentile_list.at(
            detail::kLogPeakHeightChangeIndex)
        << ", percentile log-width change = "
        << transformed_change_stats.percentile_list.at(
            detail::kLogWidthChangeIndex)
        << ", and percentile offset-to-peak-ratio change = "
        << transformed_change_stats.percentile_list.at(
            detail::kOffsetToPeakRatioChangeIndex);
    AppendLocalFittingOffsetSummary(message, offset_stats);
    message << ".";
    Logger::Log(LogLevel::Info, message.str());
}

void LogLocalFittingMaximumIterations(
    const FitOptions & options,
    const LocalFittingBestAuditState * applied_audit_state,
    const LocalFittingTerminalSummary & terminal_summary,
    const LocalFittingOffsetStats & applied_offset_stats)
{
    if (options.quiet_mode) return;

    Logger::FinishProgressLine();
    std::ostringstream warning_message;
    warning_message << "Reached maximum iteration size";
    AppendLocalFittingTerminalSummary(warning_message, terminal_summary);
    if (applied_audit_state != nullptr)
    {
        warning_message << "; applying best validated audit state";
        AppendLocalFittingAuditSummary(
            warning_message,
            *applied_audit_state);
    }
    else
    {
        warning_message
            << "; applying current accepted candidate because no finite fixed audit state is available";
    }
    AppendLocalFittingOffsetSummary(warning_message, applied_offset_stats);
    warning_message << ".";
    Logger::Log(LogLevel::Warning, warning_message.str());
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

} // namespace

void RunFixedOffsetLocalFitting(
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
    size_t atom_count{ 0 };
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

        if (!options.quiet_mode)
        {
#ifdef USE_OPENMP
            #pragma omp critical
#endif
            {
                atom_count++;
                Logger::ProgressPercent(atom_count, selected_atom_size);
            }
        }
    }
}

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
    LocalFittingPass pass)
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
            auto sample_entries{
                local_view.GetSamplingEntries(
                    pass == LocalFittingPass::FirstStage,
                    pass == LocalFittingPass::ThirdStage)
            };
            if (!HasEnoughSamplesInFitRange(
                    sample_entries,
                    options.distance_min,
                    options.distance_max,
                    kMinimumAlphaRTrainingSampleCount)) continue;
            sample_entries_list.emplace_back(std::move(sample_entries));
        }
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

namespace {

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
    for (const auto group_key : group_key_list)
    {
        analysis.SetAtomGroupAlphaG(group_key, alpha_g);
    }
}

} // namespace

void RunSecondStageLocalFitting(
    ModelObject & model_object,
    const FitOptions & options,
    const SecondStageLocalFittingInternalOptions & internal_options)
{
    const auto context{ BuildSecondStageLocalFittingContext(model_object) };
    const auto atom_size{ context.AtomSize() };
    if (!options.quiet_mode)
    {
        Logger::Log(LogLevel::Info, "Run 2nd-stage local atom fitting with iterations...");
    }

    auto initial_state_build_result{
        BuildInitialLocalFittingState(
            context,
            model_object.GetAnalysisView())
    };
    if (!initial_state_build_result.state.has_value())
    {
        if (!options.quiet_mode)
        {
            Logger::Log(
                LogLevel::Warning,
                "Skip 2nd-stage local atom fitting because no valid Gaussian seed "
                "is available for every selected atom.");
        }
        return;
    }
    LogSecondStageSeedRepairs(
        initial_state_build_result.repair_record_list,
        options);
    auto local_editor_list{ BuildAtomLocalEditors(model_object, context.atom_list) };
    auto previous_state{ std::move(*initial_state_build_result.state) };
    auto best_audit_state{
        BuildInitialLocalFittingBestAuditState(context, previous_state)
    };
    algorithm::ClusteredAndersonAccelerationHistorySet acceleration_history{
        algorithm::AndersonAccelerationOptions{
            kLocalFittingAndersonHistoryDepth,
            kLocalFittingAndersonScaleFloor,
            kLocalFittingAndersonCoefficientL1Limit,
            kLocalFittingAndersonRegularization,
            kLocalFittingAndersonCoefficientAbsLimit
        }
    };
    detail::LocalFittingAndersonRegimeTracker anderson_regime_tracker;
    LocalFittingJointOffsetStatusTracker joint_offset_status_tracker;
    algorithm::ConvergenceFreezeTracker freeze_tracker{
        atom_size,
        kLocalFittingFreezeTrackerChangeTolerance,
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
    std::vector<std::size_t> suspicious_offset_state_index_list;
    PersistentSuspiciousRollbackStateMap persistent_suspicious_rollback_state_by_key;
    PersistentJointOffsetFailureStateMap persistent_joint_offset_failure_state_by_key;
    std::vector<char> terminal_fallback_atom_mask(atom_size, 0);
    LocalFittingTerminalSummary terminal_summary;
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
        const auto active_index_list{
            BuildEligibleLocalFittingActiveIndexList(
                freeze_tracker,
                terminal_fallback_atom_mask)
        };
        if (active_index_list.empty())
        {
            ApplyLocalFittingState(previous_state, local_editor_list);
            if (terminal_summary.AtomCount() > 0)
            {
                LogLocalFittingTerminalFallback(
                    options,
                    accepted_iteration_count,
                    terminal_summary,
                    SummarizeLocalFittingOffsets(
                        previous_state.estimation_list));
            }
            else
            {
                LogLocalFittingAllAtomsFrozen(
                    options,
                    accepted_iteration_count,
                    SummarizeLocalFittingOffsets(
                        previous_state.estimation_list));
            }
            break;
        }

        auto cluster_map{ BuildLocalFittingClusters(context, active_index_list) };
        std::vector<LocalFittingClusterKey> cluster_key_list;
        cluster_key_list.reserve(cluster_map.size());
        for (const auto & [key, cluster] : cluster_map)
        {
            static_cast<void>(cluster);
            cluster_key_list.emplace_back(key);
        }
        cluster_quality_state.Reconcile(
            cluster_key_list,
            [&](const LocalFittingClusterKey & key)
            {
                return BuildInitialLocalFittingClusterQualityState(
                    context,
                    previous_state,
                    key,
                    cluster_map.at(key));
            });
        acceleration_history.Reconcile(cluster_key_list);
        anderson_regime_tracker.Reconcile(cluster_key_list);
        joint_offset_status_tracker.Reconcile(cluster_key_list);
        auto joint_offset_ridge_multiplier_list{
            cluster_quality_state.BuildObjectiveRidgeMultiplierList(atom_size)
        };
        for (const auto state_index : suspicious_offset_state_index_list)
        {
            if (state_index >= atom_size)
            {
                throw std::invalid_argument(
                    "Local fitting suspicious offset atom index is out of range.");
            }
            joint_offset_ridge_multiplier_list.at(state_index) = std::max(
                joint_offset_ridge_multiplier_list.at(state_index),
                kSuspiciousJointOffsetRidgeMultiplier);
        }
        auto iteration_result{
            RunLocalFittingIteration(
                context,
                active_index_list,
                cluster_key_list,
                previous_state,
                options,
                ridge_ratio,
                joint_offset_ridge_multiplier_list)
        };
        const auto current_health_by_key{
            std::move(iteration_result.health_by_key)
        };
        const auto actual_unhealthy_key_list{
            CollectUnhealthyLocalFittingClusterKeys(current_health_by_key)
        };
        const auto actual_progress_ineligible_key_list{
            CollectProgressIneligibleLocalFittingClusterKeys(current_health_by_key)
        };
        const auto policy_stationarity_ineligible_key_list{
            internal_options.enable_end_to_end_health_status ?
                actual_unhealthy_key_list :
                std::vector<LocalFittingClusterKey>{}
        };
        const auto policy_progress_ineligible_key_list{
            internal_options.enable_end_to_end_health_status ?
                actual_progress_ineligible_key_list :
                std::vector<LocalFittingClusterKey>{}
        };
        const auto health_summary{
            internal_options.enable_end_to_end_health_status ?
                SummarizeLocalFittingClusterHealth(current_health_by_key) :
                LocalFittingClusterHealthSummary{}
        };
        const auto current_anderson_regime_signature_by_key{
            std::move(iteration_result.anderson_regime_signature_by_key)
        };
        suspicious_offset_state_index_list =
            std::move(iteration_result.suspicious_offset_state_index_list);
        const auto has_suspicious_offset_fallback{
            !suspicious_offset_state_index_list.empty()
        };
        acceleration_history.ClearAndSuppress(policy_progress_ineligible_key_list);
        anderson_regime_tracker.Invalidate(policy_progress_ineligible_key_list);
        joint_offset_status_tracker.Invalidate(policy_progress_ineligible_key_list);
        if (internal_options.enable_end_to_end_health_status)
        {
            const auto incompatible_status_key_list{
                joint_offset_status_tracker.FindIncompatible(current_health_by_key)
            };
            acceleration_history.ClearAndSuppress(incompatible_status_key_list);
            anderson_regime_tracker.Invalidate(incompatible_status_key_list);
            joint_offset_status_tracker.Invalidate(incompatible_status_key_list);
        }
        if (!current_anderson_regime_signature_by_key.empty())
        {
            const auto incompatible_regime_key_list{
                anderson_regime_tracker.FindIncompatible(
                    current_anderson_regime_signature_by_key)
            };
            acceleration_history.ClearAndSuppress(incompatible_regime_key_list);
            anderson_regime_tracker.Invalidate(incompatible_regime_key_list);
            joint_offset_status_tracker.Invalidate(incompatible_regime_key_list);
        }
        if (has_suspicious_offset_fallback)
        {
            acceleration_history.ClearAndSuppressContaining(
                suspicious_offset_state_index_list);
            anderson_regime_tracker.InvalidateContaining(
                suspicious_offset_state_index_list);
            joint_offset_status_tracker.InvalidateContaining(
                suspicious_offset_state_index_list);
        }
        const auto raw_state{ std::move(iteration_result.state) };
        const auto raw_fixed_point_residual_stats{
            SummarizeLocalFittingTransformedChanges(
                raw_state.estimation_list,
                previous_state.estimation_list,
                active_index_list)
        };
        const auto raw_offset_stats{
            SummarizeLocalFittingOffsets(raw_state.estimation_list)
        };
        auto selection{
            SelectLocalFittingClusterCandidates(
                context,
                std::move(cluster_map),
                cluster_key_list,
                previous_state,
                raw_state,
                acceleration_history,
                cluster_quality_state)
        };
        const auto objective_selection_summary{
            SummarizeLocalFittingClusterSelection(selection)
        };
        auto assembled_state{ std::move(selection.assembled_state) };
        const auto policy_stationarity_ineligible_atom_index_list{
            CollectLocalFittingClusterAtomIndexes(
                policy_stationarity_ineligible_key_list)
        };
        freeze_tracker.ResetStability(
            policy_stationarity_ineligible_atom_index_list);
        freeze_tracker.ResetStability(
            CollectLocalFittingClusterAtomIndexes(selection.rejected_key_list));
        acceleration_history.ClearAndSuppress(policy_progress_ineligible_key_list);
        anderson_regime_tracker.Invalidate(policy_progress_ineligible_key_list);
        joint_offset_status_tracker.Invalidate(policy_progress_ineligible_key_list);
        const auto terminal_suspicious_key_list{
            UpdatePersistentSuspiciousRollbackState(
                cluster_key_list,
                selection.accepted_key_list,
                suspicious_offset_state_index_list,
                assembled_state,
                previous_state,
                persistent_suspicious_rollback_state_by_key)
        };
        TerminalJointOffsetFailureMap terminal_joint_offset_failure_by_key;
        if (internal_options.enable_end_to_end_health_status)
        {
            terminal_joint_offset_failure_by_key =
                UpdatePersistentJointOffsetFailureState(
                    cluster_key_list,
                    selection.accepted_key_list,
                    suspicious_offset_state_index_list,
                    current_health_by_key,
                    assembled_state,
                    previous_state,
                    persistent_joint_offset_failure_state_by_key);
        }
        else
        {
            persistent_joint_offset_failure_state_by_key.clear();
        }

        std::vector<LocalFittingClusterKey> terminal_key_list{
            terminal_suspicious_key_list
        };
        terminal_summary.suspicious_cluster_count +=
            terminal_suspicious_key_list.size();
        for (const auto & key : terminal_suspicious_key_list)
        {
            terminal_summary.suspicious_atom_count += key.size();
        }
        for (const auto & [key, status] : terminal_joint_offset_failure_by_key)
        {
            terminal_key_list.emplace_back(key);
            terminal_summary.joint_offset_failure_cluster_count++;
            terminal_summary.joint_offset_failure_atom_count += key.size();
            terminal_summary.joint_offset_failure_status_count[status]++;
        }
        ApplyTerminalFallbackClusters(
            terminal_key_list,
            previous_state,
            terminal_fallback_atom_mask,
            assembled_state);
        ReconcileLocalFittingBestAuditTerminalFallback(
            context,
            terminal_key_list,
            assembled_state,
            best_audit_state);
        if (!terminal_key_list.empty())
        {
            acceleration_history.ClearAndSuppress(terminal_key_list);
            anderson_regime_tracker.Invalidate(terminal_key_list);
            joint_offset_status_tracker.Invalidate(terminal_key_list);
            for (const auto & key : terminal_key_list)
            {
                freeze_tracker.ResetStability(key);
            }
            suspicious_offset_state_index_list.erase(
                std::remove_if(
                    suspicious_offset_state_index_list.begin(),
                    suspicious_offset_state_index_list.end(),
                    [&](std::size_t atom_index)
                    {
                        return terminal_fallback_atom_mask.at(atom_index) != 0;
                    }),
                suspicious_offset_state_index_list.end());
        }
        const auto continuing_accepted_key_list{
            ExcludeLocalFittingClusterKeys(
                selection.accepted_key_list,
                terminal_key_list)
        };
        const auto progress_eligible_accepted_key_list{
            ExcludeLocalFittingClusterKeys(
                continuing_accepted_key_list,
                policy_progress_ineligible_key_list)
        };
        const auto progress_eligible_key_list{
            ExcludeLocalFittingClusterKeysContainingAtoms(
                progress_eligible_accepted_key_list,
                suspicious_offset_state_index_list)
        };
        const auto stationarity_eligible_accepted_key_list{
            ExcludeLocalFittingClusterKeys(
                continuing_accepted_key_list,
                policy_stationarity_ineligible_key_list)
        };
        const auto stationarity_eligible_key_list{
            ExcludeLocalFittingClusterKeysContainingAtoms(
                stationarity_eligible_accepted_key_list,
                suspicious_offset_state_index_list)
        };
        if (selection.accepted_key_list.empty())
        {
            bool increased_cluster_objective_ridge{ false };
            bool increased_global_ridge_ratio{ false };
            if (!selection.rejected_key_list.empty())
            {
                increased_cluster_objective_ridge =
                    cluster_quality_state.IncreaseObjectiveRidge(selection.rejected_key_list);
            }
            if (selection.has_objective_backtracking_rejection && !increased_cluster_objective_ridge)
            {
                const auto previous_ridge_ratio{ ridge_ratio };
                ridge_ratio = std::min(kJointOffsetRidgeRatioMax, ridge_ratio * kJointOffsetRidgeGrowth);
                increased_global_ridge_ratio = ridge_ratio > previous_ridge_ratio;
            }
            if ((increased_cluster_objective_ridge || increased_global_ridge_ratio) &&
                iter + 1 < kLocalFittingMaximumIterations)
            {
                LogLocalFittingBacktrackingRetry(
                    options, accepted_iteration_count, ridge_ratio,
                    increased_cluster_objective_ridge,
                    health_summary,
                    objective_selection_summary,
                    raw_fixed_point_residual_stats.percentile_list.at(
                        detail::kOffsetToPeakRatioChangeIndex),
                    raw_offset_stats);
                LogRejectedLocalFittingClusterDiagnostics(
                    options,
                    selection.rejected_cluster_diagnostic_list);
                continue;
            }

            const auto stop_reason{
                iter + 1 >= kLocalFittingMaximumIterations ?
                    LocalFittingBacktrackingStopReason::MaximumIterationLimit :
                    LocalFittingBacktrackingStopReason::MaximumGlobalRidge
            };
            const LocalFittingBestAuditState * applied_audit_state{ nullptr };
            const GaussianFittingState * applied_state{ &previous_state };
            if (stop_reason == LocalFittingBacktrackingStopReason::MaximumIterationLimit &&
                best_audit_state.best_state.has_value() &&
                best_audit_state.best_objective.has_value())
            {
                applied_audit_state = &best_audit_state;
                applied_state = &*best_audit_state.best_state;
            }
            ApplyLocalFittingState(*applied_state, local_editor_list);
            LogRejectedLocalFittingClusterDiagnostics(
                options,
                selection.rejected_cluster_diagnostic_list);
            LogLocalFittingBacktrackingStop(
                options,
                stop_reason,
                applied_audit_state,
                terminal_summary,
                SummarizeLocalFittingOffsets(
                    applied_state->estimation_list));
            break;
        }

        auto change_list{
            CalculateLocalFittingTransformedChanges(
                assembled_state.estimation_list,
                previous_state.estimation_list)
        };
        auto freeze_evidence_change_list{
            CalculateLocalFittingFreezeEvidenceChanges(
                assembled_state.estimation_list,
                raw_state.estimation_list,
                previous_state.estimation_list)
        };
        auto transformed_change_stats{
            SummarizeLocalFittingTransformedChanges(
                assembled_state.estimation_list,
                previous_state.estimation_list,
                active_index_list)
        };
        const auto accepted_offset_stats{
            SummarizeLocalFittingOffsets(assembled_state.estimation_list)
        };
        accepted_iteration_count++;
        TryUpdateLocalFittingBestAuditState(
            context,
            assembled_state,
            accepted_iteration_count,
            best_audit_state);
        cluster_quality_state.DecreaseObjectiveRidge(progress_eligible_key_list);
        if (!selection.rejected_key_list.empty())
        {
            cluster_quality_state.IncreaseObjectiveRidge(selection.rejected_key_list);
        }

        if (policy_progress_ineligible_key_list.empty() &&
            !selection.has_objective_backtracking_rejection)
        {
            ridge_ratio = std::max(kJointOffsetRidgeRatioMin, ridge_ratio * kJointOffsetRidgeShrink);
        }
        if (has_suspicious_offset_fallback)
        {
            acceleration_history.ClearAndSuppressContaining(
                suspicious_offset_state_index_list);
            anderson_regime_tracker.InvalidateContaining(
                suspicious_offset_state_index_list);
            joint_offset_status_tracker.InvalidateContaining(
                suspicious_offset_state_index_list);
        }
        acceleration_history.Commit(
            progress_eligible_key_list,
            previous_state.estimation_list,
            raw_state.estimation_list);
        const auto regime_signature_commit_key_list{
            CollectLocalFittingClusterKeysWithRegimeSignatures(
                progress_eligible_key_list,
                current_anderson_regime_signature_by_key)
        };
        anderson_regime_tracker.Commit(
            regime_signature_commit_key_list,
            current_anderson_regime_signature_by_key);
        if (internal_options.enable_end_to_end_health_status)
        {
            joint_offset_status_tracker.Commit(
                progress_eligible_key_list,
                current_health_by_key);
        }
        std::vector<std::size_t> accepted_active_index_list;
        for (const auto & key : continuing_accepted_key_list)
        {
            accepted_active_index_list.insert(
                accepted_active_index_list.end(), key.begin(), key.end());
        }
        std::vector<std::size_t> stationarity_eligible_active_index_list;
        for (const auto & key : stationarity_eligible_key_list)
        {
            stationarity_eligible_active_index_list.insert(
                stationarity_eligible_active_index_list.end(), key.begin(), key.end());
        }
        freeze_tracker.Update(
            freeze_evidence_change_list,
            stationarity_eligible_active_index_list);
        ThawChangedActiveAtomNeighbors(
            context, change_list, accepted_active_index_list,
            freeze_tracker, thaw_hysteresis_tracker);
        for (const auto state_index : suspicious_offset_state_index_list)
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
            options,
            accepted_iteration_count,
            selection.accepted_candidate_attempt.value(),
            freeze_tracker,
            health_summary,
            terminal_summary,
            objective_selection_summary,
            raw_fixed_point_residual_stats.percentile_list.at(
                detail::kOffsetToPeakRatioChangeIndex),
            transformed_change_stats.percentile_list.at(
                detail::kOffsetToPeakRatioChangeIndex),
            raw_offset_stats,
            accepted_offset_stats);
        LogRejectedLocalFittingClusterDiagnostics(
            options,
            selection.rejected_cluster_diagnostic_list);

        if (BuildEligibleLocalFittingActiveIndexList(
                freeze_tracker,
                terminal_fallback_atom_mask).empty())
        {
            ApplyLocalFittingState(assembled_state, local_editor_list);
            if (terminal_summary.AtomCount() > 0)
            {
                LogLocalFittingTerminalFallback(
                    options,
                    accepted_iteration_count,
                    terminal_summary,
                    accepted_offset_stats);
            }
            else
            {
                LogLocalFittingAllAtomsFrozen(
                    options,
                    accepted_iteration_count,
                    accepted_offset_stats);
            }
            break;
        }

        const auto converged{
            policy_stationarity_ineligible_key_list.empty() &&
            !has_suspicious_offset_fallback &&
            !selection.has_objective_backtracking_rejection &&
            cluster_quality_state.AllActiveReferencesLocked(cluster_key_list) &&
            IsLocalFittingTransformedChangeConverged(transformed_change_stats) &&
            IsLocalFittingTransformedChangeConverged(raw_fixed_point_residual_stats)
        };
        if (converged)
        {
            ApplyLocalFittingState(assembled_state, local_editor_list);
            if (terminal_summary.AtomCount() > 0)
            {
                LogLocalFittingTerminalFallback(
                    options,
                    accepted_iteration_count,
                    terminal_summary,
                    accepted_offset_stats);
            }
            else
            {
                LogLocalFittingConverged(
                    options,
                    accepted_iteration_count,
                    transformed_change_stats,
                    accepted_offset_stats);
            }
            break;
        }

        if (iter + 1 == kLocalFittingMaximumIterations)
        {
            const LocalFittingBestAuditState * applied_audit_state{ nullptr };
            const GaussianFittingState * applied_state{ &assembled_state };
            if (best_audit_state.best_state.has_value() &&
                best_audit_state.best_objective.has_value())
            {
                applied_audit_state = &best_audit_state;
                applied_state = &*best_audit_state.best_state;
            }
            ApplyLocalFittingState(*applied_state, local_editor_list);
            LogLocalFittingMaximumIterations(
                options,
                applied_audit_state,
                terminal_summary,
                SummarizeLocalFittingOffsets(
                    applied_state->estimation_list));
            break;
        }
        previous_state = std::move(assembled_state);
    }
}

namespace {

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
    size_t key_count{ 0 };

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

} // namespace

void RunPotentialFittingWorkflow(ModelObject & model_object, const FitOptions & options)
{
    RunLocalAlphaTraining(model_object, options, LocalFittingPass::FirstStage);

    InitializeLocalFittingSeedModels(model_object);
    RunFixedOffsetLocalFitting(model_object, options, LocalFittingPass::FirstStage);
    RunGroupAlphaTraining(model_object, options);
    RunGroupPotentialFitting(model_object, options);

    RunSecondStageLocalFitting(model_object, options);
    LogSelectedAtomOffsetSummary(
        model_object,
        options,
        "Offset summary after 2nd-stage local fitting");
    RunGroupAlphaTraining(model_object, options);
    SetUpdatedSamplingEntriesFromGroupMedianGaussian(model_object);
    RunGroupPotentialFitting(model_object, options);
    LogAtomGroupPriorOffsetSummary(
        model_object,
        options,
        "Offset summary for group priors before 3rd-stage local fitting");
    SetUpdatedSamplingEntriesFromFittedGroupGaussian(model_object);
    RunLocalAlphaTraining(model_object, options, LocalFittingPass::ThirdStage);
    RunFixedOffsetLocalFitting(model_object, options, LocalFittingPass::ThirdStage);
    LogSelectedAtomOffsetSummary(
        model_object,
        options,
        "Offset summary after 3rd-stage local fitting");

    RunGroupAlphaTraining(model_object, options);
    RunGroupPotentialFitting(model_object, options);
    if (!options.quiet_mode)
    {
        LogGroupPriorSpotSummary(model_object);
    }
}

} // namespace rhbm_gem::core
