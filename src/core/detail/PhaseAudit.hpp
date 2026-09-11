#pragma once

#include "core/detail/IterationProposal.hpp"
#include "core/detail/ObjectiveEvaluation.hpp"

#include <atomic>
#include <mutex>
#include <string>

namespace rhbm_gem::core::detail {

struct CandidateSelectionInputs;

// Observation only: snapshots are owned here; production workspaces are never used.
class PhaseAudit
{
    struct BoundaryMember
    {
        ClusterKey key;
        std::vector<SampleRef> samples;
        std::optional<ObjectiveBreakdown> previous;
        ClusterObjectiveState history;
    };
    struct BoundaryGates
    {
        std::vector<BoundaryMember> members;
        ResidualBaseline residual_baseline;
        FitState previous_state;
        ObjectiveBreakdown improvement_reference;
    };
    struct Event
    {
        std::string id, parent_id, stage, disposition, reason;
        ClusterKey key;
        FittedGaussianSnapshot state;
        double factor{ 1.0 };
        bool probe{ false }, recertify{ true };
    };
    SecondStageContext m_context;
    ObjectiveDomain m_domain;
    FitState m_baseline;
    std::vector<ClusterKey> m_keys;
    std::size_t m_attempt, m_domain_id;
    std::vector<Event> m_events;
    std::map<std::string, BoundaryGates> m_boundary_gates;
    std::map<ClusterKey, std::vector<Event>> m_worker_events;
    void MergeWorkers();
    std::mutex m_mutex;
    std::atomic<std::size_t> m_capture_failures{ 0 };
    std::string Add(std::string stage, ClusterKey key, FittedGaussianSnapshot state,
        FittedGaussianSnapshot parent, double factor, std::string disposition,
        std::string reason, bool probe, bool recertify);

public:
    PhaseAudit(const SecondStageContext &, const ObjectiveDomain &, const FitState &,
        std::vector<ClusterKey>, std::size_t attempt, std::size_t domain_id);
    void Capture(std::string_view stage, const ClusterKey &, const FitStateView &,
        const FitStateView * parent = nullptr, double factor = 1.0,
        std::string_view disposition = "observed", std::string_view reason = "",
        bool probe = false, bool recertify = true) noexcept;
    void CaptureState(std::string_view stage, const FitState &, bool probe = false) noexcept;
    void CaptureCorrection(std::string_view stage, const ClusterKey &, const FitStateView &,
        const FitStateView & parent, double factor, std::string_view disposition,
        std::string_view reason, const CandidateSelectionInputs &, const std::vector<ClusterKey> &,
        const ObjectiveBreakdown & improvement_reference) noexcept;
    void CaptureIntermediate(std::string_view stage, const FittedGaussianSnapshot &) noexcept;
    void CaptureOperator(const FixedPointOperatorEvidence &) noexcept;
    void Missing(std::string_view stage, const ClusterKey &, std::string_view reason) noexcept;
    void CaptureSearchAssembly() noexcept;
    void Finish(const FitOptions &, const std::vector<double> & ridge,
        const SuspiciousBlockActivity &, const IterationProposalResult & production,
        const FitState & final_state) noexcept;
};

std::shared_ptr<PhaseAudit> BeginPhaseAudit(const SecondStageContext &, bool quiet,
    const ObjectiveDomain &, const FitState &, const std::vector<ClusterKey> &,
    std::size_t attempt, std::size_t domain_id) noexcept;
std::string_view PhaseAuditRejectionReason(const ObjectiveAttemptDiagnostic &);

} // namespace rhbm_gem::core::detail
