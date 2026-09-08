#pragma once

#include "core/detail/JointFitting.hpp"
#include "core/detail/SuspiciousUpdate.hpp"

namespace rhbm_gem::core {
struct FitOptions;
}

namespace rhbm_gem::core::detail {

struct FixedPointOperatorEvidence
{
    FittedGaussianSnapshot state{};
    std::vector<char> shape_available_atom_mask{};
    std::vector<char> offset_available_atom_mask{};
};

struct IterationProposalResult
{
    FitState proposal_state{};
    FixedPointOperatorEvidence fixed_point_operator{};
    SuspiciousBlockActivity block_activity{};
    std::vector<SuspiciousGaussianAssessment> assessment_by_atom{};
    std::vector<std::optional<RHBMEstimationStatus>> local_refit_status_by_atom{};
    ClusterHealthMap health_by_key{};
};

IterationProposalResult BuildIterationProposal(
    const SecondStageContext & context,
    const std::vector<ClusterKey> & cluster_key_list,
    const FitState & previous_state,
    const FitOptions & options,
    const std::vector<double> & ridge_multiplier_list,
    const SuspiciousBlockActivity & quarantine_activity,
    ClusterSolverWorkspaceMap & solver_workspace_by_key,
    std::string_view diagnostic_phase = "outer-operator");

} // namespace rhbm_gem::core::detail
