#pragma once

#include "core/detail/JointFitting.hpp"
#include "core/detail/SuspiciousUpdate.hpp"

#include <span>
#include <variant>

namespace rhbm_gem::core::detail {

struct ClusterCandidateDiagnostic;
struct FixedPointOperatorEvidence;

constexpr std::size_t kPersistentQuarantineFailureIterationLimit{ 5 };

enum class QuarantineTargetKind
{
    ShapeAtom,
    OffsetAtom,
    HardFailureCluster
};

struct QuarantineTarget
{
    QuarantineTargetKind kind{ QuarantineTargetKind::ShapeAtom };
    std::vector<std::size_t> atom_index_list{};

    friend auto operator<=>(const QuarantineTarget &, const QuarantineTarget &) = default;
};

struct StabilizationTerminalFailure
{
    StabilizationTerminalReason category{ StabilizationTerminalReason::None };
    std::optional<SuspiciousGaussianReason> guard_reason{};

    friend auto operator<=>(
        const StabilizationTerminalFailure &,
        const StabilizationTerminalFailure &) = default;
};

using QuarantineFailureReason =
    std::variant<JointOffsetSolveStatus, StabilizationTerminalFailure>;

using QuarantineFailureReasonMap = std::map<QuarantineTarget, QuarantineFailureReason>;

enum class QuarantineLifecycle
{
    Active,
    Frozen
};

struct QuarantineFailureState
{
    QuarantineFailureReason reason{};
    std::size_t stable_iteration_count{ 0 };
    std::size_t last_domain_revision{ 0 };
    QuarantineLifecycle lifecycle{ QuarantineLifecycle::Active };
};

using QuarantineFailureStateMap = std::map<QuarantineTarget, QuarantineFailureState>;

struct QuarantineStateTransition
{
    std::vector<QuarantineTarget> entered_target_list{};
    std::vector<QuarantineTarget> released_target_list{};
    std::vector<QuarantineTarget> failed_retry_target_list{};
};

QuarantineStateTransition UpdateQuarantineFailureState(
    const QuarantineFailureReasonMap & failure_reason_by_target,
    const std::vector<QuarantineTarget> & retry_target_list,
    const std::vector<QuarantineTarget> & successful_retry_target_list,
    std::size_t domain_revision,
    QuarantineFailureStateMap & state_by_target);

struct QuarantineState
{
    QuarantineFailureStateMap state_by_target{};
    std::vector<QuarantineTarget> retry_target_list{};
    std::size_t entered_target_count{ 0 };
    std::size_t released_target_count{ 0 };
    std::size_t failed_retry_count{ 0 };
    QuarantineState() = default;

    explicit QuarantineState(std::size_t atom_count)
        : m_atom_count(atom_count)
    {
    }

    SuspiciousBlockActivity BeginIteration(std::size_t domain_revision);
    SuspiciousBlockActivity BuildFinalActivity() const;
    bool UpdateAfterIteration(
        std::span<const ClusterCandidateDiagnostic> accepted_diagnostic_list,
        std::span<const ClusterCandidateDiagnostic> rejected_diagnostic_list,
        const SuspiciousBlockActivity & block_activity,
        std::span<const SuspiciousGaussianAssessment> assessment_by_atom,
        const ClusterHealthMap & health_by_key,
        const FixedPointOperatorEvidence & operator_evidence,
        const FitState & assembled_state,
        const FitState & previous_state,
        std::size_t domain_revision);
    std::size_t AtomCount() const;
    std::size_t TargetCount() const
    {
        return static_cast<std::size_t>(std::ranges::count_if(
            state_by_target,
            [](const auto & entry)
            {
                return entry.second.lifecycle != QuarantineLifecycle::Active;
            }));
    }
private:
    std::size_t m_atom_count{ 0 };
};

} // namespace rhbm_gem::core::detail
