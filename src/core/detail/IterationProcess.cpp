#include "core/detail/IterationProcess.hpp"

#include "core/detail/Diagnosis.hpp"
#include "core/detail/IterationProposal.hpp"
#include "core/detail/Quarantine.hpp"
#include "core/detail/BoundaryReconciliation.hpp"
#include "core/detail/DependencyPolish.hpp"
#include "core/detail/PreparedLocalGaussianFit.hpp"
#include "core/detail/CandidateSelection.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <ranges>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rhbm_gem::core::detail {

namespace {

constexpr std::size_t kAuditPatience{ 3 };
constexpr double kNeighborContributionDistanceMax{ 2.5 };
constexpr double kNeighborAtomSearchRange{ 2.0 * kNeighborContributionDistanceMax };
constexpr double kSuspiciousJointOffsetRidgeMultiplier{ 10.0 };
constexpr double kConvergencePercentile{ 0.99 };

bool UsesPolish(const PolishProvenance & provenance)
{
    return std::ranges::any_of(
        provenance,
        [](char value) { return value != 0; });
}

std::vector<double> BuildSuspiciousJointOffsetRidgeMultiplierList(
    const SuspiciousUpdateMask & rollback_atom_mask,
    const SuspiciousBlockActivity & block_activity,
    const std::set<std::size_t> & probation_atom_index_set)
{
    std::vector<double> ridge_multiplier_list(rollback_atom_mask.size(), 1.0);
    for (std::size_t atom_index = 0; atom_index < rollback_atom_mask.size(); atom_index++)
    {
        if (rollback_atom_mask.at(atom_index) != 0 ||
            !block_activity.HasActiveShape(atom_index) ||
            !block_activity.HasActiveOffset(atom_index) ||
            probation_atom_index_set.contains(atom_index))
        {
            ridge_multiplier_list.at(atom_index) =
                kSuspiciousJointOffsetRidgeMultiplier;
        }
    }
    return ridge_multiplier_list;
}

void ResetClusterSolverWorkspace(
    const std::vector<ClusterKey> & cluster_key_list,
    ClusterSolverWorkspaceMap & workspace_by_key)
{
    workspace_by_key.clear();
    for (const auto & key : cluster_key_list)
    {
        workspace_by_key.try_emplace(key);
    }
}

struct SecondStageInitializationResult
{
    SecondStageContext context{};
    FitState state{};
    std::vector<int> neighbor_count_list{};
    std::vector<SecondStageSeedSelectionRecord> selection_record_list{};
};

struct PendingTopology
{
    GraphTopology topology{};
    CouplingGraphPartition partition{};
};

struct IterationState
{
    FitState accepted_state{};
    FittedGaussianSnapshot topology_reference_state{};
    PolishProvenance previous_polish_provenance{};
    SuspiciousUpdateMask rollback_atom_mask{};
    std::vector<std::size_t> selected_atom_index_list{};
    CouplingGraphPartition graph_partition{};
    std::optional<PendingTopology> pending_topology{};
    ClusterSolverWorkspaceMap solver_workspace_by_key{};
    BoundaryJointCorrectionWorkspaceMap boundary_joint_correction_workspace_by_key{};
    ObjectiveDomain objective_domain{};
    BestAuditState best_audit_state{};
    QuarantineState quarantine_state{};
    ClusterObjectiveStateMap cluster_objective_state{};
    TrustRegionStateSet trust_region_state{};
    std::size_t accepted_iteration_count{ 0 };
    std::size_t accepted_iterations_since_topology_rebuild{ 0 };
    std::size_t audit_patience_count{ 0 };
};

static void ValidateBlockActivitySize(
    std::size_t atom_count,
    const SuspiciousBlockActivity & block_activity)
{
    if (block_activity.shape_fixed_atom_mask.size() != atom_count ||
        block_activity.offset_fixed_atom_mask.size() != atom_count ||
        block_activity.hard_failure_atom_mask.size() != atom_count)
    {
        throw std::invalid_argument("Active-coordinate convergence inputs are inconsistent.");
    }
}

static ConvergenceCertificate SummarizeFixedPointOperator(
    const FixedPointOperatorEvidence & evidence,
    const FitState & previous_state,
    const std::vector<std::size_t> & atom_index_list,
    std::string_view diagnostic_phase = "outer-operator")
{
    ConvergenceCertificate result;
    const ActiveCoordinatePopulation operator_nominal_population{
        atom_index_list,
        atom_index_list
    };

    std::vector<TransformedChange> change_list;
    change_list.reserve(previous_state.size());
    for (std::size_t atom_index = 0;
        atom_index < previous_state.size(); atom_index++)
    {
        auto change{ CalculateTransformedChange(
            GetFitModel(evidence.state, atom_index),
            GetFitModel(previous_state, atom_index)) };
        if (evidence.shape_available_atom_mask.at(atom_index) == 0)
        {
            change.at(GaussianModel3D::LogPeakHeightCoordinateIndex()) = std::numeric_limits<double>::infinity();
            change.at(GaussianModel3D::LogWidthCoordinateIndex()) = std::numeric_limits<double>::infinity();
        }
        if (evidence.offset_available_atom_mask.at(atom_index) == 0)
        {
            change.at(GaussianModel3D::OffsetToPeakRatioCoordinateIndex()) = std::numeric_limits<double>::infinity();
        }
        change_list.emplace_back(std::move(change));
    }
    result.operator_nominal_residual = SummarizeActiveDofChanges(change_list, operator_nominal_population);
    result.operator_complete = std::ranges::all_of(
        atom_index_list,
        [&](const auto atom_index)
        {
            return evidence.shape_available_atom_mask.at(atom_index) != 0 &&
                evidence.offset_available_atom_mask.at(atom_index) != 0;
        });
    if (Logger::GetLogLevel() >= LogLevel::Debug)
    {
        std::size_t shape_unavailable{ 0 };
        std::size_t offset_unavailable{ 0 };
        for (const auto atom_index : atom_index_list)
        {
            shape_unavailable += evidence.shape_available_atom_mask.at(atom_index) == 0;
            offset_unavailable += evidence.offset_available_atom_mask.at(atom_index) == 0;
        }
        std::ostringstream message;
        message << "Second-stage availability: schema=1, phase=" << diagnostic_phase
            << ", nominal-atoms=" << atom_index_list.size()
            << ", shape-unavailable=" << shape_unavailable
            << ", offset-unavailable=" << offset_unavailable << ".";
        Logger::Log(LogLevel::Debug, message.str());
    }
    return result;
}

static std::optional<SecondStageInitializationResult> BuildSecondStageInitialization(
    const ModelObject & model_object,
    const FitOptions & options)
{
    SecondStageInitializationResult build_result;
    auto & context{ build_result.context };
    auto & state{ build_result.state };

    const auto & atom_list{ model_object.GetSelectedAtoms() };
    context.atom_list.reserve(atom_list.size());
    for (auto * atom : atom_list)
    {
        context.atom_list.emplace_back(AtomContext{ atom });
    }
    state.resize(context.atom_list.size());

    std::unordered_map<const AtomObject *, std::size_t> atom_index_map;
    atom_index_map.reserve(context.atom_list.size());
    for (std::size_t i = 0; i < context.atom_list.size(); i++)
    {
        atom_index_map.emplace(context.atom_list.at(i).atom, i);
    }
    for (std::size_t atom_index = 0;
        atom_index < context.atom_list.size();
        atom_index++)
    {
        auto & atom_context{ context.atom_list.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        atom_context.raw_sampling_entries = local_view.GetRawSamplingEntries(false);
        state.at(atom_index) = local_view.GetGaussianResult(FittingStage::Second);
        atom_context.alpha_r = local_view.GetAlphaR(FittingStage::Second);
        atom_context.refit_design = PreparedLocalGaussianDesign{
            atom_context.raw_sampling_entries,
            options.distance_min,
            options.distance_max
        };
    }

    std::vector<GaussianModel3D> global_models;
    global_models.reserve(context.atom_list.size());
    for (std::size_t atom_index = 0;
        atom_index < context.atom_list.size();
        atom_index++)
    {
        global_models.emplace_back(state.at(atom_index).mdpde.GetModel());
    }
    const auto global_median{ BuildGaussianParameterMedian(global_models) };

    for (std::size_t i = 0; i < context.atom_list.size(); i++)
    {
        auto & result{ state.at(i) };
        const auto original_model{ result.mdpde.GetModel() };
        const auto selection{ SelectSecondStageSeed(result.mdpde, global_median) };
        if (!selection.has_value())
        {
            return std::nullopt;
        }

        result.mdpde = selection->model;
        build_result.selection_record_list.emplace_back(
            SecondStageSeedSelectionRecord{
                selection->source,
                original_model,
                selection->model.GetModel()
            });
    }

    build_result.neighbor_count_list.reserve(context.atom_list.size());
    for (std::size_t atom_index = 0;
        atom_index < context.atom_list.size();
        atom_index++)
    {
        auto & atom_context{ context.atom_list.at(atom_index) };
        const auto * atom{ atom_context.atom };
        const auto neighbor_atom_list{
            atom->FindNeighborAtoms(kNeighborAtomSearchRange)
        };
        std::unordered_set<const AtomObject *> neighbor_atom_set;
        atom_context.unselected_distance_list_by_sample.resize(atom_context.raw_sampling_entries.size());

        atom_context.neighbor_atom_sample_offset_list.reserve(atom_context.raw_sampling_entries.size() + 1);
        atom_context.neighbor_atom_sample_offset_list.emplace_back(0);
        atom_context.neighbor_atom_sample_list.reserve(
            atom_context.raw_sampling_entries.size() * neighbor_atom_list.size());
        for (std::size_t sample_index = 0;
            sample_index < atom_context.raw_sampling_entries.size();
            sample_index++)
        {
            const auto & sample{
                atom_context.raw_sampling_entries.at(sample_index)
            };
            std::unordered_set<const AtomObject *> sample_neighbor_set;
            for (auto * neighbor_atom : neighbor_atom_list)
            {
                if (options.exclude_hydrogen && neighbor_atom->GetElement() == Element::HYDROGEN)
                {
                    continue;
                }
                const auto distance{
                    array_helper::ComputeNorm(sample.point.position, neighbor_atom->GetPositionRef())
                };
                if (distance > kNeighborContributionDistanceMax ||
                    !sample_neighbor_set.emplace(neighbor_atom).second) continue;
                neighbor_atom_set.emplace(neighbor_atom);

                const auto selected_iter{
                    atom_index_map.find(neighbor_atom)
                };
                if (selected_iter != atom_index_map.end())
                {
                    atom_context.neighbor_atom_sample_list.emplace_back(
                        NeighborAtomSample{
                            selected_iter->second,
                            distance
                        });
                    continue;
                }

                atom_context.unselected_distance_list_by_sample.at(sample_index).emplace_back(distance);
            }
            atom_context.neighbor_atom_sample_offset_list.emplace_back(
                atom_context.neighbor_atom_sample_list.size());
        }
        build_result.neighbor_count_list.emplace_back(static_cast<int>(neighbor_atom_set.size()));
    }

    return build_result;
}

static void StoreSecondStageNeighborCounts(
    ModelObject & model_object,
    const SecondStageContext & context,
    const std::vector<int> & neighbor_count_list)
{
    auto analysis{ model_object.EditAnalysis() };
    for (std::size_t atom_index = 0;
        atom_index < context.atom_list.size();
        atom_index++)
    {
        analysis.SetAtomLocalNeighborCountForPeeling(
            *context.atom_list.at(atom_index).atom,
            neighbor_count_list.at(atom_index));
    }
}

static void RefreshBestAuditState(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot,
    IterationState & iteration_state)
{
    ReevaluateBestAuditState(context, iteration_state.objective_domain, iteration_state.best_audit_state);
    const auto audit_objective{
        EvaluateAuditObjective(
            iteration_state.objective_domain,
            context, model_snapshot)
    };
    if (audit_objective.has_value())
    {
        TryUpdateBestAuditState(
            iteration_state.accepted_state,
            UsesPolish(iteration_state.previous_polish_provenance),
            iteration_state.accepted_iteration_count,
            *audit_objective,
            iteration_state.best_audit_state);
    }
}

static void ResetIterationStateForPartition(
    const SecondStageContext & context,
    const FitOptions & options,
    CouplingGraphPartition partition,
    IterationState & iteration_state,
    PerformanceCounters & performance_counters)
{
    const auto cluster_key_list{ BuildGraphClusterKeyList(partition) };
    const auto model_snapshot{
        BuildSecondStageModelSnapshot(context, iteration_state.accepted_state)
    };
    iteration_state.objective_domain = BuildObjectiveDomain(
        context,
        model_snapshot,
        cluster_key_list,
        options.distance_min,
        options.distance_max);
    iteration_state.cluster_objective_state.clear();
    const auto objective_by_key{
        BuildObjectiveByKey(
            partition,
            iteration_state.objective_domain,
            context, model_snapshot)
    };
    ReconcileClusterObjectiveState(
        objective_by_key,
        iteration_state.cluster_objective_state);
    RefreshBestAuditState(context, model_snapshot, iteration_state);
    iteration_state.trust_region_state.Reconcile(cluster_key_list);
    performance_counters.RecordSolverWorkspaceReset();
    ResetClusterSolverWorkspace(
        cluster_key_list,
        iteration_state.solver_workspace_by_key);
    iteration_state.boundary_joint_correction_workspace_by_key.clear();
    iteration_state.graph_partition = std::move(partition);
    iteration_state.audit_patience_count = 0;
}

static bool TryRebuildAdaptiveTopology(
    const SecondStageContext & context,
    const FitOptions & options,
    const FitState & accepted_state,
    GraphTopology & graph_topology,
    IterationState & iteration_state,
    PerformanceCounters & performance_counters)
{
    const auto decision{
        EvaluateAdaptiveTopologyRebuildTrigger(
            accepted_state,
            iteration_state.topology_reference_state,
            iteration_state.selected_atom_index_list,
            iteration_state.accepted_iterations_since_topology_rebuild)
    };
    if (decision.trigger == AdaptiveTopologyRebuildTrigger::None) return false;

    FinishProgressLine(options.quiet_mode);
    const auto rebuild_start{ std::chrono::steady_clock::now() };
    auto rebuilt_topology{
        BuildSecondStageGraphTopology(
            context,
            accepted_state,
            options.quiet_mode,
            &graph_topology)
    };
    auto rebuilt_partition{
        BuildGraphPartition(
            rebuilt_topology,
            iteration_state.selected_atom_index_list)
    };
    const auto partition_changed{
        iteration_state.graph_partition != rebuilt_partition
    };
    const auto elapsed_milliseconds{
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - rebuild_start).count()
    };
    performance_counters.RecordTopologyRebuild(elapsed_milliseconds, partition_changed);

    LogAdaptiveTopologyRebuild(
        options.quiet_mode,
        iteration_state.accepted_iteration_count,
        decision,
        graph_topology,
        rebuilt_topology,
        iteration_state.graph_partition,
        rebuilt_partition,
        partition_changed);

    iteration_state.topology_reference_state = BuildSecondStageModelSnapshot(context, accepted_state).node;
    iteration_state.accepted_iterations_since_topology_rebuild = 0;
    if (partition_changed)
    {
        iteration_state.pending_topology = PendingTopology{
            std::move(rebuilt_topology), std::move(rebuilt_partition) };
    }
    else graph_topology = std::move(rebuilt_topology);
    return partition_changed;
}

static IterationState BuildIterationState(
    SecondStageContext & context,
    const GraphTopology & graph_topology,
    FitState initial_state,
    const FitOptions & options)
{
    IterationState iteration_state;
    iteration_state.accepted_state = std::move(initial_state);
    iteration_state.topology_reference_state =
        BuildSecondStageModelSnapshot(context, iteration_state.accepted_state).node;
    iteration_state.previous_polish_provenance.assign(context.atom_list.size(), 0);
    iteration_state.rollback_atom_mask.assign(context.atom_list.size(), 0);
    iteration_state.quarantine_state = QuarantineState(context.atom_list.size());
    iteration_state.selected_atom_index_list.resize(context.atom_list.size());
    std::iota(
        iteration_state.selected_atom_index_list.begin(),
        iteration_state.selected_atom_index_list.end(),
        0);
    iteration_state.graph_partition = BuildGraphPartition(graph_topology, iteration_state.selected_atom_index_list);
    const auto cluster_key_list{ BuildGraphClusterKeyList(iteration_state.graph_partition) };
    context.frozen_background = BuildFrozenBackground(context, iteration_state.accepted_state, cluster_key_list);
    if (!context.frozen_background)
        throw std::runtime_error("Second-stage initial fixed background is unavailable.");
    LogFrozenBackground(context, options.quiet_mode);
    ResetClusterSolverWorkspace(
        cluster_key_list,
        iteration_state.solver_workspace_by_key);
    const auto initial_model_snapshot{
        BuildSecondStageModelSnapshot(context, iteration_state.accepted_state)
    };
    iteration_state.objective_domain = BuildObjectiveDomain(
        context,
        initial_model_snapshot,
        cluster_key_list,
        options.distance_min,
        options.distance_max);
    const auto initial_audit_objective{
        EvaluateAuditObjective(
            iteration_state.objective_domain,
            context, initial_model_snapshot)
    };
    if (initial_audit_objective.has_value())
    {
        TryUpdateBestAuditState(
            iteration_state.accepted_state,
            UsesPolish(
                iteration_state.previous_polish_provenance),
            0,
            *initial_audit_objective,
            iteration_state.best_audit_state);
    }
    return iteration_state;
}

static bool BeginFrozenBackgroundIteration(
    SecondStageContext & context,
    GraphTopology & graph_topology,
    const FitOptions & options,
    IterationState & iteration_state,
    PerformanceCounters & performance_counters)
{
    const bool partition_changed{ iteration_state.pending_topology.has_value() };
    const auto & partition{ partition_changed ?
        iteration_state.pending_topology->partition : iteration_state.graph_partition };
    const auto background{ BuildFrozenBackground(
        context, iteration_state.accepted_state, BuildGraphClusterKeyList(partition)) };
    if (!background)
        throw std::runtime_error("Second-stage fixed background refresh is unavailable.");

    const auto previous_background{ context.frozen_background };
    context.frozen_background = background;
    LogFrozenBackground(context, options.quiet_mode);
    if (partition_changed)
    {
        auto pending{ std::move(*iteration_state.pending_topology) };
        iteration_state.pending_topology.reset();
        ResetIterationStateForPartition(context, options, std::move(pending.partition),
            iteration_state, performance_counters);
        graph_topology = std::move(pending.topology);
        LogObjectiveDomain(iteration_state.objective_domain, options.quiet_mode, true);
        return true;
    }
    if (previous_background && previous_background->response_by_atom == background->response_by_atom)
        return false;

    const auto previous_snapshot{ BuildSecondStageModelSnapshot(context, iteration_state.accepted_state) };
    const auto previous_objectives{ BuildObjectiveByKey(partition, iteration_state.objective_domain,
        context, previous_snapshot) };
    std::vector<char> changed_background_by_atom(context.atom_list.size(), 1);
    if (previous_background)
    {
        for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
            changed_background_by_atom.at(atom_index) =
                previous_background->response_by_atom.at(atom_index) != background->response_by_atom.at(atom_index);
    }
    for (const auto & [key, sample_refs] : partition.sample_id_list_by_key)
    {
        const bool changed{ std::ranges::any_of(sample_refs,
            [&](const auto & sample_ref) { return changed_background_by_atom.at(sample_ref.atom_index) != 0; }) };
        if (changed)
            iteration_state.cluster_objective_state[key] = ClusterObjectiveState{
                .best_objective = previous_objectives.at(key) };
    }
    RefreshBestAuditState(context, previous_snapshot, iteration_state);
    return false;
}

static IterationResult RunIteration(
    SecondStageContext & context,
    GraphTopology & graph_topology,
    const FitOptions & options,
    std::size_t attempt_number,
    IterationState & iteration_state,
    PerformanceCounters & performance_counters)
{
    // Prepare this attempt's frozen background, objectives, and active blocks.
    const bool background_partition_changed{ attempt_number > 1 && BeginFrozenBackgroundIteration(
        context, graph_topology, options, iteration_state, performance_counters) };
    const auto & previous_state{ iteration_state.accepted_state };
    const auto & selected_atom_index_list{ iteration_state.selected_atom_index_list };
    const auto & graph_partition{ iteration_state.graph_partition };
    const auto cluster_key_list{ BuildGraphClusterKeyList(graph_partition) };
    const auto & objective_domain{ iteration_state.objective_domain };

    const auto residual_baseline{
        BuildResidualBaseline(context, previous_state)
    };
    performance_counters.RecordGaussianCacheMisses();

    const auto previous_objective_by_key{
        BuildObjectiveByKey(graph_partition, objective_domain, residual_baseline)
    };
    ReconcileClusterObjectiveState(previous_objective_by_key, iteration_state.cluster_objective_state);
    iteration_state.trust_region_state.Reconcile(cluster_key_list);

    const auto quarantine_activity{
        iteration_state.quarantine_state.BeginIteration(
            iteration_state.accepted_iteration_count)
    };
    std::set<std::size_t> probation_atom_index_set;
    for (const auto & target : iteration_state.quarantine_state.probation_target_list)
    {
        probation_atom_index_set.insert(target.atom_index_list.begin(), target.atom_index_list.end());
    }
    std::vector<ClusterKey> probation_key_list;
    for (const auto & key : cluster_key_list)
    {
        if (std::ranges::any_of(
                key,
                [&](const auto atom_index)
                {
                    return probation_atom_index_set.contains(atom_index);
                }))
        {
            probation_key_list.emplace_back(key);
        }
    }
    iteration_state.trust_region_state.ResetToMinimum(probation_key_list);
    const auto joint_offset_ridge_multiplier_list{
        BuildSuspiciousJointOffsetRidgeMultiplierList(
            iteration_state.rollback_atom_mask,
            quarantine_activity,
            probation_atom_index_set)
    };
    // Build a constrained proposal while retaining unrestricted operator evidence.
    const auto iteration_phase_start{ std::chrono::steady_clock::now() };
    auto proposal_result{
        BuildIterationProposal(
            context,
            cluster_key_list,
            previous_state,
            options,
            joint_offset_ridge_multiplier_list,
            quarantine_activity,
            iteration_state.solver_workspace_by_key)
    };
    performance_counters.FinishIterationPhase(iteration_phase_start);
    performance_counters.RecordGaussianCacheHits();

    LogUnrestrictedOperatorAssessments(
        options.quiet_mode,
        proposal_result.assessment_by_atom,
        proposal_result.block_activity);

    const auto & proposal_state{
        proposal_result.proposal_state
    };
    const auto proposal_change_summary{
        SummarizeTransformedChanges(
            proposal_state, previous_state, selected_atom_index_list)
    };
    // Select cluster candidates, then reconcile their shared boundary samples.
    const CandidateSelectionInputs candidate_inputs{
        .context = context,
        .options = options,
        .residual_baseline = residual_baseline,
        .partition = graph_partition,
        .health_by_key = proposal_result.health_by_key,
        .previous_state = previous_state,
        .previous_polish_provenance = iteration_state.previous_polish_provenance,
        .proposal_state = proposal_state,
        .block_activity = proposal_result.block_activity,
        .ridge_multiplier_list = joint_offset_ridge_multiplier_list,
        .objective_domain = objective_domain,
        .previous_objective_by_key = previous_objective_by_key,
        .cluster_objective_state = iteration_state.cluster_objective_state,
        .best_audit_state = iteration_state.best_audit_state,
        .trust_region_state = iteration_state.trust_region_state,
        .solver_workspace_by_key = iteration_state.solver_workspace_by_key,
        .boundary_joint_correction_workspace_by_key = iteration_state.boundary_joint_correction_workspace_by_key,
        .performance_counters = performance_counters
    };
    auto selection{ SelectClusterCandidates(candidate_inputs) };
    const auto iteration_failure_atom_mask{
        BuildSuspiciousFailureAtomMask(
            selection.block_activity,
            proposal_result.assessment_by_atom)
    };
    const auto iteration_suspicious_atom_count{
        static_cast<std::size_t>(std::ranges::count_if(
            iteration_failure_atom_mask,
            [](char value) { return value != 0; }))
    };
    const auto has_suspicious_offset_fallback{ iteration_suspicious_atom_count > 0 };

    auto & assembled_state{ selection.assembled_state };
    auto & assembled_polish_provenance{ selection.assembled_polish_provenance };
    IterationResult result;
    // Quarantine transitions can change the assembled state and require a new audit.
    const auto has_quarantine_transition{
        iteration_state.quarantine_state.UpdateAfterIteration(
            selection.accepted_cluster_diagnostic_list,
            selection.rejected_cluster_diagnostic_list,
            selection.block_activity,
            proposal_result.assessment_by_atom,
            proposal_result.health_by_key,
            assembled_state,
            previous_state,
            iteration_state.previous_polish_provenance,
            assembled_polish_provenance,
            iteration_state.accepted_iteration_count + 1)
    };
    if (has_quarantine_transition)
    {
        selection.final_audit_objective.reset();
        ReauditFallbackSelection(candidate_inputs, selection);
    }
    // Publish audited history before the all-rejected exit, as on accepted attempts.
    iteration_state.cluster_objective_state = std::move(selection.cluster_objective_state);
    result.trust_region_update = iteration_state.trust_region_state.ApplyRadiusUpdates(
        selection.grow_trust_region_key_list, selection.shrink_trust_region_key_list,
        selection.rejected_key_list, selection.exhausted_key_list);
    const auto assembled_uses_polish{
        UsesPolish(assembled_polish_provenance)
    };
    result.accepted_cluster_diagnostic_list = std::move(selection.accepted_cluster_diagnostic_list);
    result.rejected_cluster_diagnostic_list = std::move(selection.rejected_cluster_diagnostic_list);
    result.boundary_reconciliation_diagnostic_list = std::move(selection.boundary_reconciliation_diagnostic_list);
    iteration_state.rollback_atom_mask = selection.block_activity.BuildCombinedFixedAtomMask();
    result.attempt_number = attempt_number;
    result.accepted_iteration_count = iteration_state.accepted_iteration_count;
    result.quarantine_atom_count = iteration_state.quarantine_state.AtomCount();
    result.active_atom_count = context.atom_list.size() - result.quarantine_atom_count;
    result.polish_progress = selection.polish_progress;
    result.suspicious_atom_count = iteration_suspicious_atom_count;
    result.proposal_maximum_transformed_change = std::ranges::max(proposal_change_summary.maximum_list);

    if (selection.accepted_key_list.empty())
    {
        result.stop_reason = attempt_number >= kMaximumIterations ?
            SecondStageStopReason::AllRejectedAtMaximumIterations :
            SecondStageStopReason::AllRejectedBacktrackingExhausted;
        return result;
    }

    // Certify accepted movement separately from the unrestricted operator residual.
    const auto active_population{
        BuildActiveCoordinatePopulation(
            selected_atom_index_list,
            selection.block_activity)
    };
    const auto transformed_change_summary{
        SummarizeTransformedChanges(
            assembled_state,
            previous_state,
            iteration_state.selected_atom_index_list)
    };
    auto certificate{
        SummarizeFixedPointOperator(
            proposal_result.fixed_point_operator,
            previous_state,
            iteration_state.selected_atom_index_list)
    };
    certificate.accepted_active_movement = SummarizeActiveDofChanges(
        assembled_state,
        previous_state,
        active_population);
    certificate.solver_qualified = AreActiveCoordinatesSolverQualified(
        iteration_state.selected_atom_index_list,
        cluster_key_list,
        selection.block_activity,
        proposal_result.local_refit_status_by_atom,
        proposal_result.health_by_key);
    // Advance accepted progress; a changed partition takes effect next attempt.
    iteration_state.accepted_iteration_count++;
    iteration_state.accepted_iterations_since_topology_rebuild++;
    result.objective_domain_changed = TryRebuildAdaptiveTopology(
        context,
        options,
        assembled_state,
        graph_topology,
        iteration_state,
        performance_counters) || background_partition_changed;
    if (result.objective_domain_changed)
    {
        iteration_state.quarantine_state.force_probation = true;
    }
    bool improved_best_audit{ false };
    bool improved_audit_baseline{ false };
    {
        auto candidate_audit_objective{ selection.final_audit_objective };
        if (!candidate_audit_objective.has_value())
        {
            const auto candidate_model_snapshot{
                BuildSecondStageModelSnapshot(context, assembled_state)
            };
            candidate_audit_objective = EvaluateAuditObjective(
                iteration_state.objective_domain,
                context, candidate_model_snapshot);
        }
        if (candidate_audit_objective.has_value())
        {
            const auto previous_audit_objective{ EvaluateAuditObjective(objective_domain, residual_baseline) };
            improved_audit_baseline = previous_audit_objective.has_value() && IsBetterAuditObjective(
                candidate_audit_objective->GetTotalObjective(), previous_audit_objective->GetTotalObjective(),
                kObjectiveStrictTolerance);
            improved_best_audit = TryUpdateBestAuditState(
                assembled_state,
                assembled_uses_polish,
                iteration_state.accepted_iteration_count,
                *candidate_audit_objective,
                iteration_state.best_audit_state);
        }
    }
    if (improved_best_audit)
    {
        performance_counters.RecordFullStateMaterialization();
    }
    const auto changed_rejected_trust_radius{
        std::ranges::any_of(
            result.rejected_cluster_diagnostic_list,
            [&](const auto & diagnostic)
            {
                return std::ranges::find(
                    result.trust_region_update.changed_key_list,
                    diagnostic.key) !=
                    result.trust_region_update.changed_key_list.end();
            })
    };
    const auto has_pending_quarantine_lifecycle{
        std::ranges::any_of(
            iteration_state.quarantine_state.state_by_target | std::views::values,
            [](const auto & state)
            {
                return state.lifecycle != QuarantineLifecycle::Exhausted;
            })
    };
    if (result.objective_domain_changed || improved_audit_baseline ||
        changed_rejected_trust_radius || has_pending_quarantine_lifecycle ||
        has_quarantine_transition)
    {
        iteration_state.audit_patience_count = 0;
    }
    else
    {
        iteration_state.audit_patience_count++;
    }

    result.accepted_iteration_count = iteration_state.accepted_iteration_count;
    result.accepted_maximum_transformed_change = std::ranges::max(transformed_change_summary.maximum_list);
    result.transformed_change_percentile = certificate.accepted_active_movement.percentile_list;
    certificate.objective_domain_changed = result.objective_domain_changed;
    certificate.quarantine_transition = has_quarantine_transition;
    certificate.suspicious_offset_fallback = has_suspicious_offset_fallback;
    certificate.rejected_cluster = !selection.rejected_key_list.empty();
    if (certificate.ProductionConverged())
    {
        result.stop_reason = SecondStageStopReason::Converged;
    }
    else if (iteration_state.audit_patience_count >= kAuditPatience)
    {
        result.stop_reason = SecondStageStopReason::AuditPatience;
    }

    LogConvergenceSafeguardAudit(options.quiet_mode, result, certificate);

    // Commit the accepted state even when this attempt reaches a stopping condition.
    iteration_state.accepted_state = std::move(assembled_state);
    iteration_state.previous_polish_provenance = std::move(assembled_polish_provenance);
    return result;
}

void ApplyFitState(
    ModelObject & model_object,
    const SecondStageContext & context,
    const FitState & iteration_state)
{
    const auto model_snapshot{
        BuildSecondStageModelSnapshot(context, iteration_state)
    };

    auto analysis{ model_object.EditAnalysis() };
    for (std::size_t i = 0; i < context.atom_list.size(); i++)
    {
        auto adjusted_sampling_entries{
            BuildSecondStageAdjustedSamples(context, i, model_snapshot)
        };
        analysis.ApplyAtomLocalSecondStageResult(
            *context.atom_list.at(i).atom,
            iteration_state.at(i),
            std::move(adjusted_sampling_entries));
    }
}

struct FinalPolishResidualSafetyResult
{
    FinalPolishResidualSafetyStatus status{ FinalPolishResidualSafetyStatus::NotEvaluated };
    std::optional<ConvergenceCertificate> base{};
    std::optional<ConvergenceCertificate> candidate{};
};

static std::optional<ConvergenceCertificate> EvaluateFinalPolishCertificate(
    const SecondStageContext & context,
    const FitOptions & options,
    IterationState & iteration_state,
    const SuspiciousBlockActivity & final_block_activity,
    const FitState & candidate_state)
{
    try
    {
        const auto cluster_key_list{
            BuildGraphClusterKeyList(iteration_state.graph_partition)
        };
        const auto joint_offset_ridge_multiplier_list{
            BuildSuspiciousJointOffsetRidgeMultiplierList(
                iteration_state.rollback_atom_mask,
                final_block_activity,
                {})
        };
        auto certificate_options{ options };
        certificate_options.quiet_mode = true;
        auto proposal_result{
            BuildIterationProposal(
                context,
                cluster_key_list,
                candidate_state,
                certificate_options,
                joint_offset_ridge_multiplier_list,
                final_block_activity,
                iteration_state.solver_workspace_by_key,
                "final-recertification")
        };
        auto certificate{ SummarizeFixedPointOperator(
            proposal_result.fixed_point_operator,
            candidate_state,
            iteration_state.selected_atom_index_list,
            "final-recertification") };
        certificate.solver_qualified = AreActiveCoordinatesSolverQualified(
            iteration_state.selected_atom_index_list,
            cluster_key_list,
            proposal_result.block_activity,
            proposal_result.local_refit_status_by_atom,
            proposal_result.health_by_key);
        return certificate;
    }
    catch (...)
    {
        return std::nullopt;
    }
}

static bool HasComparableFinalPolishOperatorEvidence(const ConvergenceCertificate & certificate)
{
    const auto & percentile_list{
        certificate.operator_nominal_residual.percentile_list
    };
    return certificate.solver_qualified && certificate.operator_complete &&
        std::ranges::all_of(
            percentile_list,
            [](double value) { return std::isfinite(value); });
}

static bool IsFinalPolishResidualNonWorsening(
    const ConvergenceCertificate & base_certificate,
    const ConvergenceCertificate & candidate_certificate)
{
    if (!HasComparableFinalPolishOperatorEvidence(base_certificate) ||
        !HasComparableFinalPolishOperatorEvidence(candidate_certificate))
    {
        return false;
    }
    const auto & base_percentile_list{
        base_certificate.operator_nominal_residual.percentile_list
    };
    const auto & candidate_percentile_list{
        candidate_certificate.operator_nominal_residual.percentile_list
    };
    for (std::size_t index = 0; index < base_percentile_list.size(); index++)
    {
        if (candidate_percentile_list.at(index) >
            std::max(base_percentile_list.at(index), kTransformedChangeTolerance))
        {
            return false;
        }
    }
    return true;
}

static FinalPolishResidualSafetyResult EvaluateFinalPolishResidualSafety(
    const SecondStageContext & context,
    const FitOptions & options,
    IterationState & iteration_state,
    const SuspiciousBlockActivity & final_block_activity,
    const FitState & base_state,
    const FitState & candidate_state,
    FinalPolishCertificationPolicy policy)
{
    FinalPolishResidualSafetyResult result;
    result.candidate = EvaluateFinalPolishCertificate(
        context,
        options,
        iteration_state,
        final_block_activity,
        candidate_state);
    if (!result.candidate.has_value())
    {
        result.status = FinalPolishResidualSafetyStatus::Error;
        return result;
    }
    if (result.candidate->StrictOperatorPassed())
    {
        result.status = FinalPolishResidualSafetyStatus::AbsolutePassed;
        return result;
    }
    if (policy == FinalPolishCertificationPolicy::RequireStrictFixedPoint)
    {
        result.status = FinalPolishResidualSafetyStatus::Failed;
        return result;
    }
    if (!HasComparableFinalPolishOperatorEvidence(*result.candidate))
    {
        result.status = FinalPolishResidualSafetyStatus::Failed;
        return result;
    }
    result.base = EvaluateFinalPolishCertificate(
        context,
        options,
        iteration_state,
        final_block_activity,
        base_state);
    if (!result.base.has_value())
    {
        result.status = FinalPolishResidualSafetyStatus::Error;
        return result;
    }
    result.status = IsFinalPolishResidualNonWorsening(
            *result.base,
            *result.candidate) ?
        FinalPolishResidualSafetyStatus::RelativePassed :
        FinalPolishResidualSafetyStatus::Failed;
    return result;
}

static const FitState & FinalizeSecondStageState(
    ModelObject & model_object,
    const SecondStageContext & context,
    const FitOptions & options,
    const GraphTopology & graph_topology,
    IterationState & iteration_state,
    bool use_best_audit_state,
    FinalPolishCertificationPolicy certification_policy,
    PerformanceCounters & performance_counters)
{
    const auto final_uses_best_audit{
        use_best_audit_state && iteration_state.best_audit_state.has_value()
    };
    auto & final_state{
        final_uses_best_audit ?
            iteration_state.best_audit_state->state :
            iteration_state.accepted_state
    };
    const auto final_block_activity{
        iteration_state.quarantine_state.BuildFinalActivity()
    };
    auto polish_result{
        RunFinalDependencyPolish(
            context,
            options,
            graph_topology,
            iteration_state.graph_partition,
            iteration_state.objective_domain,
            final_block_activity,
            iteration_state.trust_region_state,
            final_state,
            iteration_state.boundary_joint_correction_workspace_by_key,
            performance_counters)
    };
    FinalPolishResidualSafetyResult residual_safety;
    if (polish_result.accepted && polish_result.objective.has_value())
    {
        residual_safety = EvaluateFinalPolishResidualSafety(
            context,
            options,
            iteration_state,
            final_block_activity,
            final_state,
            polish_result.state,
            certification_policy);
    }
    const auto polish_applied{
        polish_result.accepted &&
        polish_result.objective.has_value() &&
        (residual_safety.status == FinalPolishResidualSafetyStatus::AbsolutePassed ||
            residual_safety.status == FinalPolishResidualSafetyStatus::RelativePassed)
    };
    LogFinalDependencyPolish(
        options.quiet_mode,
        polish_result,
        certification_policy,
        residual_safety.status,
        polish_applied,
        residual_safety.base.has_value() ? &*residual_safety.base : nullptr,
        residual_safety.candidate.has_value() ? &*residual_safety.candidate : nullptr);
    if (polish_applied)
    {
        if (final_uses_best_audit)
        {
            final_state = std::move(polish_result.state);
            iteration_state.best_audit_state->objective = *polish_result.objective;
            iteration_state.best_audit_state->uses_polish = true;
        }
        else
        {
            if (iteration_state.previous_polish_provenance.size() != context.atom_list.size())
            {
                iteration_state.previous_polish_provenance.resize(context.atom_list.size(), 0);
            }
            for (std::size_t atom_index = 0; atom_index < context.atom_list.size(); atom_index++)
            {
                if (IsTransformedChangeMaterial(
                        CalculateTransformedChange(
                            polish_result.state.at(atom_index).mdpde.GetModel(),
                            final_state.at(atom_index).mdpde.GetModel()),
                        kTransformedChangeTolerance))
                {
                    iteration_state.previous_polish_provenance.at(atom_index) = 1;
                }
            }
            final_state = std::move(polish_result.state);
            TryUpdateBestAuditState(
                final_state,
                true,
                iteration_state.accepted_iteration_count,
                *polish_result.objective,
                iteration_state.best_audit_state);
        }
    }
    ApplyFitState(model_object, context, final_state);
    return final_state;
}

} // namespace

bool RunSecondStageIterations(ModelObject & model_object, const FitOptions & options)
{
    if (options.enable_second_stage_dependency_polish &&
        options.second_stage_dependency_polish_max_iterations == 0)
    {
        throw std::invalid_argument(
            "Second-stage dependency polish maximum iterations must be positive when enabled.");
    }
    // Prepare seeds and sampling context before establishing the initial objective.
    SecondStageContext context;
    FitState initial_state;
    {
        auto initialization{
            BuildSecondStageInitialization(model_object, options)
        };
        LogSecondStageStart(options.quiet_mode);

        if (!initialization.has_value())
        {
            LogSecondStageInitializationFailure(options.quiet_mode);
            return false;
        }
        StoreSecondStageNeighborCounts(model_object, initialization->context,
            initialization->neighbor_count_list);
        LogSecondStageSeedSelections(
            initialization->selection_record_list,
            options.quiet_mode);
        context = std::move(initialization->context);
        initial_state = std::move(initialization->state);
    }
    auto graph_topology{
        BuildSecondStageGraphTopology(context, initial_state, options.quiet_mode)
    };
    LogGraphTopology(graph_topology, options.quiet_mode);
    auto iteration_state{
        BuildIterationState(context, graph_topology, std::move(initial_state), options)
    };
    PerformanceCounters performance_counters{
        options.quiet_mode,
        context,
        iteration_state.solver_workspace_by_key,
        iteration_state.boundary_joint_correction_workspace_by_key
    };
    if (iteration_state.best_audit_state.has_value())
    {
        performance_counters.RecordFullStateMaterialization();
    }
    LogObjectiveDomain(iteration_state.objective_domain, options.quiet_mode);
    const auto progress_column_widths{ BuildProgressColumnWidths(context.atom_list.size()) };
    LogProgressHeader(options.quiet_mode, progress_column_widths);

    const auto audit_comparison_objective_domain{ iteration_state.objective_domain };
    IterationResult terminal_result;
    if (context.atom_list.size() == 0)
    {
        terminal_result.attempt_number = 1;
        terminal_result.accepted_iteration_count =
            iteration_state.accepted_iteration_count;
        terminal_result.stop_reason = SecondStageStopReason::Quarantine;
    }
    else
    {
        for (std::size_t iter = 0; iter < kMaximumIterations; iter++)
        {
            terminal_result = RunIteration(
                context,
                graph_topology,
                options,
                iter + 1,
                iteration_state,
                performance_counters);
            LogAcceptedCandidateSearchDiagnostics(
                options.quiet_mode,
                terminal_result);
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
            LogTrustModelShadowDiagnostics(
                options.quiet_mode,
                terminal_result);
#endif
            LogRejectedClusterDiagnostics(
                options.quiet_mode,
                terminal_result.rejected_cluster_diagnostic_list);
            LogIterationProgress(
                options.quiet_mode,
                progress_column_widths,
                terminal_result);

            if (terminal_result.stop_reason ==
                    SecondStageStopReason::AllRejectedBacktrackingExhausted ||
                terminal_result.stop_reason ==
                    SecondStageStopReason::AllRejectedAtMaximumIterations)
            {
                LogAllRejectedResolution(
                    options.quiet_mode,
                    terminal_result);
            }
            if (terminal_result.stop_reason != SecondStageStopReason::None)
            {
                break;
            }
        }
    }
    if (terminal_result.stop_reason == SecondStageStopReason::None)
    {
        terminal_result.stop_reason = SecondStageStopReason::MaximumIterations;
    }

    // Finalize the chosen state and certify any dependency polish before persistence.
    const auto converged{
        terminal_result.stop_reason == SecondStageStopReason::Converged
    };
    const auto use_best_audit_state{
        !converged &&
        terminal_result.stop_reason != SecondStageStopReason::Quarantine &&
        iteration_state.best_audit_state.has_value()
    };
    const auto & finalized_state{
        FinalizeSecondStageState(
            model_object,
            context,
            options,
            graph_topology,
            iteration_state,
            use_best_audit_state,
            converged ?
                FinalPolishCertificationPolicy::RequireStrictFixedPoint :
                FinalPolishCertificationPolicy::RequireResidualNonRegression,
            performance_counters)
    };
    LogSecondStageAuditTerminal(
        options.quiet_mode,
        context,
        terminal_result.stop_reason,
        terminal_result.attempt_number,
        iteration_state.accepted_iteration_count,
        finalized_state,
        audit_comparison_objective_domain);

    if ((terminal_result.stop_reason == SecondStageStopReason::Quarantine || converged) &&
        iteration_state.quarantine_state.TargetCount() != 0)
    {
        LogQuarantineFallback(
            options.quiet_mode,
            iteration_state.accepted_iteration_count,
            iteration_state.quarantine_state.entered_target_count,
            iteration_state.quarantine_state.released_target_count,
            iteration_state.quarantine_state.failed_probation_count,
            iteration_state.quarantine_state.TargetCount(),
            finalized_state);
    }
    if (terminal_result.stop_reason == SecondStageStopReason::Quarantine)
    {
        LogNoSelectedAtoms(options.quiet_mode);
    }
    else if (converged)
    {
        if (iteration_state.quarantine_state.TargetCount() == 0)
        {
            LogConverged(
                options.quiet_mode,
                terminal_result,
                finalized_state);
        }
    }
    else if (terminal_result.stop_reason ==
        SecondStageStopReason::MaximumIterations)
    {
        LogMaximumIterations(
            options.quiet_mode,
            iteration_state.quarantine_state.entered_target_count,
            iteration_state.quarantine_state.released_target_count,
            iteration_state.quarantine_state.failed_probation_count,
            iteration_state.quarantine_state.TargetCount(),
            iteration_state.best_audit_state,
            iteration_state.accepted_state);
    }
    LogSecondStageSummary(
        options.quiet_mode,
        iteration_state.accepted_iteration_count,
        iteration_state.best_audit_state,
        iteration_state.previous_polish_provenance,
        terminal_result.stop_reason,
        use_best_audit_state);
    return true;
}

std::optional<SecondStageSeedSelection> SelectSecondStageSeed(
    const GaussianModel3DWithUncertainty & local_mdpde,
    const std::optional<GaussianModel3D> & global_median)
{
    if (IsValidSecondStageGaussianModel(local_mdpde.GetModel()))
    {
        return SecondStageSeedSelection{
            SecondStageSeedSource::LocalMdpde,
            local_mdpde
        };
    }
    if (global_median.has_value() &&
        IsValidSecondStageGaussianModel(*global_median))
    {
        return SecondStageSeedSelection{
            SecondStageSeedSource::GlobalMedian,
            GaussianModel3DWithUncertainty{
                *global_median,
                GaussianModel3DUncertainty{}
            }
        };
    }
    return std::nullopt;
}

AdaptiveTopologyRebuildDecision EvaluateAdaptiveTopologyRebuildTrigger(
    const FitState & accepted_state,
    const FittedGaussianSnapshot & topology_reference_state,
    const std::vector<std::size_t> & active_index_list,
    std::size_t accepted_iterations_since_rebuild)
{
    const auto drift_summary{
        SummarizeTransformedChanges(
            accepted_state,
            topology_reference_state,
            active_index_list)
    };
    const auto maximum_transformed_drift{ std::ranges::max(drift_summary.maximum_list) };
    if (maximum_transformed_drift >= kAdaptiveTopologyRebuildDriftThreshold)
    {
        return AdaptiveTopologyRebuildDecision{
            AdaptiveTopologyRebuildTrigger::Drift,
            maximum_transformed_drift
        };
    }
    if (accepted_iterations_since_rebuild >= kAdaptiveTopologyRebuildAcceptedIterationInterval)
    {
        return AdaptiveTopologyRebuildDecision{
            AdaptiveTopologyRebuildTrigger::Interval,
            maximum_transformed_drift
        };
    }
    return AdaptiveTopologyRebuildDecision{
        AdaptiveTopologyRebuildTrigger::None,
        maximum_transformed_drift
    };
}

bool ConvergenceCertificate::StrictOperatorPassed() const
{
    return solver_qualified && operator_complete &&
        IsTransformedPercentileConverged(operator_nominal_residual);
}

bool ConvergenceCertificate::ProductionConverged() const
{
    return StrictOperatorPassed() &&
        IsTransformedPercentileConverged(accepted_active_movement) &&
        !objective_domain_changed && !quarantine_transition &&
        !suspicious_offset_fallback && !rejected_cluster;
}

SuspiciousUpdateMask BuildSuspiciousFailureAtomMask(
    const SuspiciousBlockActivity & block_activity,
    std::span<const SuspiciousGaussianAssessment> assessment_by_atom)
{
    const auto atom_count{ assessment_by_atom.size() };
    if (block_activity.shape_fixed_atom_mask.size() != atom_count ||
        block_activity.offset_fixed_atom_mask.size() != atom_count ||
        block_activity.hard_failure_atom_mask.size() != atom_count)
    {
        throw std::invalid_argument(
            "Suspicious failure activity and assessment sizes are inconsistent.");
    }
    SuspiciousUpdateMask result(atom_count, 0);
    for (std::size_t atom_index = 0; atom_index < atom_count; atom_index++)
    {
        const auto has_fixed_endpoint{
            block_activity.shape_fixed_atom_mask.at(atom_index) != 0 ||
            block_activity.offset_fixed_atom_mask.at(atom_index) != 0
        };
        result.at(atom_index) =
            block_activity.hard_failure_atom_mask.at(atom_index) != 0 ||
            (has_fixed_endpoint &&
                assessment_by_atom[atom_index].reason !=
                    SuspiciousGaussianReason::None) ? 1 : 0;
    }
    return result;
}

ActiveCoordinatePopulation BuildActiveCoordinatePopulation(
    const std::vector<std::size_t> & atom_index_list,
    const SuspiciousBlockActivity & block_activity)
{
    ValidateBlockActivitySize(block_activity.shape_fixed_atom_mask.size(), block_activity);
    ActiveCoordinatePopulation result;
    for (const auto atom_index : atom_index_list)
    {
        if (block_activity.HasActiveShape(atom_index))
        {
            result.active_shape_atom_index_list.emplace_back(atom_index);
        }
        if (block_activity.HasActiveOffset(atom_index))
        {
            result.active_offset_atom_index_list.emplace_back(atom_index);
        }
    }
    return result;
}

TransformedChangeSummary SummarizeActiveDofChanges(
    const std::vector<TransformedChange> & change_list,
    const ActiveCoordinatePopulation & population)
{
    std::array<std::vector<double>, GaussianModel3D::TransformedCoordinateSize()> parameter_change_lists;
    for (const auto parameter_index : std::array<std::size_t, 2>{
        GaussianModel3D::LogPeakHeightCoordinateIndex(),
        GaussianModel3D::LogWidthCoordinateIndex() })
    {
        auto & values{ parameter_change_lists.at(parameter_index) };
        values.reserve(population.active_shape_atom_index_list.size());
        for (const auto atom_index : population.active_shape_atom_index_list)
        {
            if (atom_index >= change_list.size())
            {
                throw std::invalid_argument(
                    "Active-coordinate shape change input is inconsistent.");
            }
            values.emplace_back(change_list.at(atom_index).at(parameter_index));
        }
    }

    auto & offset_change_list{
        parameter_change_lists.at(GaussianModel3D::OffsetToPeakRatioCoordinateIndex())
    };
    offset_change_list.reserve(population.active_offset_atom_index_list.size());
    for (const auto atom_index : population.active_offset_atom_index_list)
    {
        if (atom_index >= change_list.size())
        {
            throw std::invalid_argument("Active-coordinate offset change input is inconsistent.");
        }
        const auto value{
            change_list.at(atom_index).at(GaussianModel3D::OffsetToPeakRatioCoordinateIndex())
        };
        offset_change_list.emplace_back(
            std::isfinite(value) ? std::abs(value) : std::numeric_limits<double>::infinity());
    }

    TransformedChangeSummary result;
    for (std::size_t parameter_index = 0; parameter_index < parameter_change_lists.size(); parameter_index++)
    {
        const auto & values{ parameter_change_lists.at(parameter_index) };
        result.population_size_list.at(parameter_index) = values.size();
        result.percentile_list.at(parameter_index) =
            array_helper::ComputePercentile(values, kConvergencePercentile);
        result.maximum_list.at(parameter_index) =
            values.empty() ? 0.0 : *std::ranges::max_element(values);
    }
    return result;
}

TransformedChangeSummary SummarizeActiveDofChanges(
    const FitState & current_state,
    const FitState & previous_state,
    const ActiveCoordinatePopulation & population)
{
    if (current_state.size() != previous_state.size())
    {
        throw std::invalid_argument(
            "Active-coordinate transformed state sizes are inconsistent.");
    }
    std::vector<TransformedChange> change_list;
    change_list.reserve(current_state.size());
    for (std::size_t atom_index = 0; atom_index < current_state.size(); atom_index++)
    {
        change_list.emplace_back(CalculateTransformedChange(
            GetFitModel(current_state, atom_index),
            GetFitModel(previous_state, atom_index)));
    }
    return SummarizeActiveDofChanges(change_list, population);
}

bool AreActiveCoordinatesSolverQualified(
    const std::vector<std::size_t> & atom_index_list,
    const std::vector<ClusterKey> & cluster_key_list,
    const SuspiciousBlockActivity & block_activity,
    std::span<const std::optional<RHBMEstimationStatus>> local_refit_status_by_atom,
    const ClusterHealthMap & health_by_key)
{
    const auto atom_count{ block_activity.shape_fixed_atom_mask.size() };
    ValidateBlockActivitySize(atom_count, block_activity);
    if (local_refit_status_by_atom.size() != atom_count)
    {
        throw std::invalid_argument(
            "Convergence audit qualification inputs are inconsistent.");
    }

    for (const auto atom_index : atom_index_list)
    {
        if (!block_activity.HasActiveShape(atom_index)) continue;
        const auto & status{ local_refit_status_by_atom[atom_index] };
        if (!status.has_value() || !IsLocalRefitStatusSolverQualified(*status))
        {
            return false;
        }
    }

    for (const auto atom_index : atom_index_list)
    {
        if (!block_activity.HasActiveOffset(atom_index)) continue;
        const auto owner{ std::ranges::find_if(cluster_key_list, [&](const auto & key)
        {
            return std::ranges::binary_search(key, atom_index);
        }) };
        if (owner == cluster_key_list.end()) return false;
        const auto health_iter{ health_by_key.find(*owner) };
        if (health_iter == health_by_key.end() ||
            health_iter->second.joint_offset_status != JointOffsetSolveStatus::Converged)
        {
            return false;
        }
    }

    return true;
}

} // namespace rhbm_gem::core::detail
