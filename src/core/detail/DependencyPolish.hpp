#pragma once

#include "core/detail/ObjectiveEvaluation.hpp"
#include "core/detail/JointFitting.hpp"

namespace rhbm_gem::core {
struct FitOptions;
}

namespace rhbm_gem::core::detail {

struct FinalDependencyPolishDiagnostic
{
    std::size_t component_count{ 0 };
    std::size_t attempted_component_count{ 0 };
    std::size_t accepted_component_count{ 0 };
    std::size_t atom_count{ 0 };
    std::size_t parameter_count{ 0 };
    std::size_t round_count{ 0 };
    std::size_t suspicious_candidate_atom_count{ 0 };
    std::optional<double> objective_before{};
    std::optional<double> objective_after{};
    double elapsed_milliseconds{ 0.0 };
    struct Component
    {
        std::vector<ClusterKey> key_list{};
        std::size_t atom_count{ 0 };
        std::size_t parameter_count{ 0 };
        std::size_t round_count{ 0 };
        std::size_t suspicious_candidate_atom_count{ 0 };
        std::size_t symbolic_analysis_count{ 0 };
        std::optional<double> objective_before{};
        std::optional<double> objective_after{};
        double elapsed_milliseconds{ 0.0 };
        std::vector<JointCandidateObjectiveDiagnostic> objective_diagnostic_list{};
        bool accepted{ false };
    };
    std::vector<Component> component_list{};
};

struct FinalDependencyPolishResult
{
    FitState state{};
    std::optional<ObjectiveBreakdown> objective{};
    FinalDependencyPolishDiagnostic diagnostic{};
    bool accepted{ false };
};

FinalDependencyPolishResult RunFinalDependencyPolish(
    const SecondStageContext & context,
    const FitOptions & options,
    const GraphTopology & topology,
    const CouplingGraphPartition & partition,
    const ObjectiveDomain & objective_domain,
    const SuspiciousBlockActivity & block_activity,
    const FitState & base_state,
    BoundaryJointCorrectionWorkspaceMap & workspace_by_key,
    PerformanceCounters & performance_counters);

} // namespace rhbm_gem::core::detail
