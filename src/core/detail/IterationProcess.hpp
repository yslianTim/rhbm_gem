#pragma once

#include "core/detail/CandidateSelection.hpp"

#include <cstddef>
#include <map>
#include <optional>
#include <variant>
#include <vector>

namespace rhbm_gem::core::detail {

enum class SecondStageSeedSource
{
    GroupPosterior,
    GroupPrior,
    GroupMedian,
    GlobalMedian
};

struct SecondStageSeedCandidates
{
    std::optional<GaussianModel3DWithUncertainty> group_posterior{};
    std::optional<GaussianModel3DWithUncertainty> group_prior{};
    std::optional<GaussianModel3DWithUncertainty> group_median{};
    std::optional<GaussianModel3DWithUncertainty> global_median{};
};

struct SecondStageSeedSelection
{
    SecondStageSeedSource source{ SecondStageSeedSource::GlobalMedian };
    GaussianModel3DWithUncertainty model{};
};

std::optional<SecondStageSeedSelection> SelectSecondStageSeed(
    const SecondStageSeedCandidates & candidates);

constexpr std::size_t kPersistentTerminalFailureIterationLimit{ 5 };

using PersistentSuspiciousRollbackReason = std::vector<std::size_t>;
using PersistentTerminalFailureReason =
    std::variant<PersistentSuspiciousRollbackReason, JointOffsetSolveStatus>;

struct PersistentTerminalFailureState
{
    PersistentTerminalFailureReason reason{};
    std::size_t stable_iteration_count{ 0 };
};

using PersistentTerminalFailureStateMap =
    std::map<ClusterKey, PersistentTerminalFailureState>;
using TerminalPersistentFailureMap =
    std::map<ClusterKey, PersistentTerminalFailureReason>;

TerminalPersistentFailureMap UpdatePersistentTerminalFailureState(
    const std::vector<ClusterKey> & accepted_key_list,
    const SuspiciousUpdateMask & suspicious_atom_mask,
    const ClusterHealthMap & health_by_key,
    const FitState & assembled_state,
    const FitState & previous_state,
    PersistentTerminalFailureStateMap & state_by_key);

} // namespace rhbm_gem::core::detail
