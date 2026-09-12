#include "core/detail/ComponentAssembly.hpp"
#include "core/detail/DependencyPolish.hpp"

#include "core/detail/CandidateSelection.hpp"
#include "core/detail/CandidateEvaluation.hpp"
#include "core/detail/Diagnosis.hpp"

#include <chrono>
#include <exception>
#include <utility>

#include <rhbm_gem/core/GaussianEstimator.hpp>

namespace rhbm_gem::core::detail {
namespace {
constexpr double kFinalDependencyPolishTrustRadius{ 1.0 };
}

FinalDependencyPolishResult RunFinalDependencyPolish(
    const SecondStageContext & context,
    const FitOptions & options,
    const GraphTopology & topology,
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & objective_domain,
    const SuspiciousBlockActivity & block_activity,
    const FitState & base_state,
    BoundaryJointCorrectionWorkspaceMap & workspace_by_key,
    PerformanceCounters & performance_counters)
{
    FinalDependencyPolishResult result{ .state = base_state };
    if (!options.enable_second_stage_dependency_polish) return result;

    const auto polish_start{ std::chrono::steady_clock::now() };
    const auto component_list{
        BuildUncutDependencyPolishComponents(
            topology,
            partition,
            objective_domain.owner_key_by_atom_index)
    };
    result.diagnostic.component_count = component_list.size();
    result.diagnostic.component_list.reserve(component_list.size());
    const auto base_baseline{ BuildResidualBaseline(context, base_state) };
    performance_counters.RecordGaussianCacheMisses();
    const auto base_objective{
        EvaluateAuditObjective(objective_domain, base_baseline)
    };
    if (base_objective.has_value())
    {
        result.objective = base_objective;
        result.diagnostic.objective_before = base_objective->GetTotalObjective();
        result.diagnostic.objective_after = base_objective->GetTotalObjective();
    }

    std::vector<std::optional<FitStatePatch>> accepted_patch_by_component(
        component_list.size());
    std::size_t accepted_component_count{ 0 };
    std::vector<double> ridge_multiplier_list(context.atom_list.size(), 1.0);

    for (std::size_t component_position = 0;
        component_position < component_list.size();
        component_position++)
    {
        const auto & component{ component_list.at(component_position) };
        std::vector<std::size_t> shape_active_index_list;
        std::vector<std::size_t> offset_active_index_list;
        for (const auto atom_index : component.atom_index_list)
        {
            if (block_activity.HasActiveShape(atom_index))
            {
                shape_active_index_list.emplace_back(atom_index);
            }
            if (block_activity.HasActiveOffset(atom_index))
            {
                offset_active_index_list.emplace_back(atom_index);
            }
        }
        const auto parameter_count{
            2 * shape_active_index_list.size() + offset_active_index_list.size()
        };
        result.diagnostic.atom_count += component.atom_index_list.size();
        result.diagnostic.parameter_count += parameter_count;
        auto & diagnostic{
            result.diagnostic.component_list.emplace_back(
                FinalDependencyPolishDiagnostic::Component{
                    .key_list = component.key_list,
                    .atom_count = component.atom_index_list.size(),
                    .parameter_count = parameter_count
                })
        };
        if (!base_objective.has_value()) continue;

        const auto component_start{ std::chrono::steady_clock::now() };
        if (shape_active_index_list.empty() && offset_active_index_list.empty())
        {
            diagnostic.elapsed_milliseconds =
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - component_start).count();
            continue;
        }
        diagnostic.objective_before = base_objective->GetTotalObjective();
        result.diagnostic.attempted_component_count++;

        try
        {
            FitStatePatch endpoint_patch{
                FitStatePatch::FromState(base_state, component.atom_index_list)
            };
            auto endpoint_objective{ *base_objective };
            const auto maximum_round_count{
                options.second_stage_dependency_polish_max_iterations
            };
            if (maximum_round_count != 0)
            {
                result.diagnostic.round_count++;
                diagnostic.round_count++;
                std::vector<BoundaryJointTrustRegion> trust_region_list;
                trust_region_list.reserve(component.key_list.size());
                for (const auto & key : component.key_list)
                {
                    trust_region_list.emplace_back(BoundaryJointTrustRegion{
                        key,
                        kFinalDependencyPolishTrustRadius
                    });
                }
                const BoundaryJointCorrectionWorkspaceKey workspace_key{
                    shape_active_index_list,
                    offset_active_index_list,
                    component.affected_sample_ref_list,
                };
                auto & solver{
                    workspace_by_key.try_emplace(workspace_key).first->second
                };
                for (std::size_t round = 0;
                    round < maximum_round_count;
                    round++)
                {
                    if (round != 0)
                    {
                        result.diagnostic.round_count++;
                        diagnostic.round_count++;
                    }
                    const FitStateView endpoint_state_view{
                        base_state,
                        endpoint_patch
                    };
                    const auto symbolic_analysis_count_before{
                        solver.GetSymbolicAnalysisCount()
                    };
                    const auto correction_result{
                        BuildBoundaryJointCorrection(
                            context,
                            endpoint_state_view,
                            shape_active_index_list,
                            offset_active_index_list,
                            component.affected_sample_ref_list,
                            ridge_multiplier_list,
                            trust_region_list,
                            solver,
                            "final-dependency-polish",
                            JointCorrectionTrustReference::Endpoint)
                    };
                    diagnostic.symbolic_analysis_count +=
                        solver.GetSymbolicAnalysisCount() -
                        symbolic_analysis_count_before;
                    if (correction_result.status !=
                            BoundaryJointCorrectionStatus::CandidateReady ||
                        !correction_result.patch.has_value())
                    {
                        break;
                    }
                    diagnostic.parameter_count = correction_result.parameter_count;
                    const CandidateEvaluationOverlay candidate_overlay{
                        context,
                        base_baseline,
                        base_state,
                        *correction_result.patch
                    };
                    const auto evaluation{ EvaluateCandidate(candidate_overlay,
                        FinalPolishCandidateReference{component, partition, objective_domain,
                            endpoint_state_view, *base_objective, endpoint_objective, performance_counters,
                            options.quiet_mode, diagnostic.objective_diagnostic_list, correction_result.damping, round + 1}) };
                    diagnostic.suspicious_candidate_atom_count += evaluation.suspicious_atom_count;
                    result.diagnostic.suspicious_candidate_atom_count += evaluation.suspicious_atom_count;
                    if (!evaluation.objective) break;
                    const auto & candidate_objective{ evaluation.objective };

                    endpoint_patch = *correction_result.patch;
                    endpoint_objective = *candidate_objective;
                    diagnostic.objective_after =
                        candidate_objective->GetTotalObjective();
                }
            }

            if (diagnostic.objective_after.has_value() &&
                IsBetterAuditObjective(
                    *diagnostic.objective_after,
                    *diagnostic.objective_before,
                    kObjectiveStrictTolerance))
            {
                diagnostic.accepted = true;
                accepted_patch_by_component.at(component_position) =
                    std::move(endpoint_patch);
                accepted_component_count++;
            }
        }
        catch (const std::exception &)
        {
            diagnostic.accepted = false;
        }
        diagnostic.elapsed_milliseconds =
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - component_start).count();
    }

    std::vector<const FitStatePatch *> selected_patches;
    selected_patches.reserve(accepted_patch_by_component.size());
    for (const auto & patch : accepted_patch_by_component)
    {
        selected_patches.emplace_back(patch.has_value() ? &*patch : nullptr);
    }
    auto assembled_state{ AssembleComponentState(base_state, selected_patches) };
    const auto evaluate_global_audit = [&](const FitState & state)
    {
        performance_counters.RecordFullStateMaterialization();
        const auto snapshot{ BuildSecondStageModelSnapshot(context, state) };
        return EvaluateAuditObjective(objective_domain, context, snapshot);
    };

    struct ComponentRemoval
    {
        std::size_t position;
        ObjectiveBreakdown objective;
        FitState state;
    };
    auto assembled_objective{ AuditAndSalvageComponents(
        [&]
        {
            return accepted_component_count == 0 ?
                base_objective : evaluate_global_audit(assembled_state);
        },
        [&](const auto & objective)
        {
            return !base_objective.has_value() || accepted_component_count == 0 ||
                (objective.has_value() && IsBetterAuditObjective(
                    objective->GetTotalObjective(),
                    base_objective->GetTotalObjective(),
                    kObjectiveStrictTolerance));
        },
        [&](const auto & assembled_objective) -> std::optional<ComponentRemoval>
        {
            std::optional<std::size_t> removal_position;
            std::optional<ObjectiveBreakdown> best_removal_objective;
            FitState best_removal_state;
            for (std::size_t candidate_position = 0;
                candidate_position < accepted_patch_by_component.size();
                candidate_position++)
            {
                if (!accepted_patch_by_component.at(candidate_position).has_value())
                {
                    continue;
                }
                auto candidate_state{
                    AssembleComponentState(base_state, selected_patches, candidate_position)
                };
                const auto candidate_objective{
                    evaluate_global_audit(candidate_state)
                };
                if (!candidate_objective.has_value() ||
                    (assembled_objective.has_value() &&
                        candidate_objective->GetTotalObjective() >=
                            assembled_objective->GetTotalObjective()))
                {
                    continue;
                }
                if (!best_removal_objective.has_value() ||
                    candidate_objective->GetTotalObjective() <
                        best_removal_objective->GetTotalObjective())
                {
                    removal_position = candidate_position;
                    best_removal_objective = candidate_objective;
                    best_removal_state = std::move(candidate_state);
                }
            }
            if (!removal_position.has_value()) return std::nullopt;
            return ComponentRemoval{
                *removal_position, *best_removal_objective, std::move(best_removal_state)
            };
        },
        [&](ComponentRemoval & removal) -> std::optional<ObjectiveBreakdown>
        {
            result.diagnostic.component_list.at(removal.position).accepted = false;
            selected_patches.at(removal.position) = nullptr;
            accepted_patch_by_component.at(removal.position).reset();
            accepted_component_count--;
            assembled_state = std::move(removal.state);
            return removal.objective;
        }) };

    if (base_objective.has_value() &&
        assembled_objective.has_value() &&
        accepted_component_count != 0 &&
        IsBetterAuditObjective(
            assembled_objective->GetTotalObjective(),
            base_objective->GetTotalObjective(),
            kObjectiveStrictTolerance))
    {
        result.state = std::move(assembled_state);
        result.objective = assembled_objective;
        result.accepted = true;
        result.diagnostic.objective_after =
            assembled_objective->GetTotalObjective();
    }
    else
    {
        for (auto & diagnostic : result.diagnostic.component_list)
        {
            diagnostic.accepted = false;
        }
    }
    result.diagnostic.accepted_component_count =
        result.accepted ? accepted_component_count : 0;
    result.diagnostic.elapsed_milliseconds =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - polish_start).count();
    performance_counters.RecordDependencyPolish(
        result.diagnostic.component_count,
        result.diagnostic.attempted_component_count,
        result.diagnostic.accepted_component_count,
        result.diagnostic.component_count -
            result.diagnostic.accepted_component_count,
        result.diagnostic.atom_count,
        result.diagnostic.parameter_count,
        result.diagnostic.round_count,
        result.diagnostic.elapsed_milliseconds);
    return result;
}

} // namespace rhbm_gem::core::detail
