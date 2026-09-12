#pragma once

#include "core/detail/ObjectiveEvaluation.hpp"
#include "utils/hrl/EstimationAudit.hpp"

namespace rhbm_gem::core {
struct FitOptions;
}

namespace rhbm_gem::core::detail {

struct CandidateSelectionInputs;
struct CandidateSelection;
struct IterationResult;
struct IterationProposalResult;
struct TrustModelTrialRecord;
enum class BoundaryAcceptancePolicy;

// Observation entry points preserve the existing lifecycle and never choose candidates.
void BeginClusterHistoryObserver(SecondStageContext &, bool quiet) noexcept;
void ObserveHistoryPartition(const SecondStageContext &, const CouplingGraphPartition &,
    const ObjectiveDomain &, const FitState &) noexcept;
void ObserveHistoryBackground(const SecondStageContext &, const std::shared_ptr<const FrozenBackground> &,
    const CouplingGraphPartition &, const ObjectiveDomain &, const FitState &) noexcept;
void ObserveHistoryAttempt(SecondStageContext &, const ObjectiveByKey &, const FitState &,
    const CouplingGraphPartition &, const ObjectiveDomain &, std::size_t attempt, std::size_t accepted) noexcept;
void ObserveHistorySearch(const SecondStageContext &, const ClusterKey &) noexcept;
void ObserveLocalHistory(const CandidateEvaluationOverlay &, const ClusterKey &,
    const std::vector<SampleRef> &, const ObjectiveDomain &, std::string_view source,
    bool accepted, ObjectiveAttemptDiagnostic &) noexcept;
void ObserveBoundaryHistory(const SecondStageContext &, JointCandidateObjectiveDiagnostic *) noexcept;
void ObserveBoundaryMemberHistory(const CandidateEvaluationOverlay &, const ClusterKey &,
    const std::vector<SampleRef> &, const ObjectiveDomain &, bool accepted,
    const ObjectiveAttemptDiagnostic &, JointCandidateObjectiveDiagnostic *) noexcept;
void ObserveBoundaryHistoryAccepted(const SecondStageContext &, std::size_t observation) noexcept;
void ObserveHistoryRejected(const SecondStageContext &, const ClusterKey &) noexcept;
void ObserveHistoryPublication(const SecondStageContext &) noexcept;

void BeginPhaseObservation(SecondStageContext &, bool quiet, const ObjectiveDomain &, const FitState &,
    const std::vector<ClusterKey> &, std::size_t attempt, std::size_t domain_id) noexcept;

class ProductionObservationScope
{
    estimation_audit::Scope m_scope;
public:
    ProductionObservationScope(const SecondStageContext &, std::size_t attempt);
    ProductionObservationScope(const ProductionObservationScope &) = delete;
    ProductionObservationScope & operator=(const ProductionObservationScope &) = delete;
};

enum class BoundaryObservationStage { Endpoint, Correction, Backtracking };
std::string_view BoundaryDiagnosticName(BoundaryAcceptancePolicy, BoundaryObservationStage) noexcept;
std::string_view BoundaryPhaseName(BoundaryAcceptancePolicy, BoundaryObservationStage) noexcept;

class TrustModelTrialObserver
{
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
    const CandidateSelectionInputs & inputs;
    const ClusterKey & key;
    const std::vector<SampleRef> & samples;
    TrustModelTrialRecord & record;
    std::optional<std::size_t> final_trial{};
    std::size_t search_pass{ 0 };
public:
    TrustModelTrialObserver(const CandidateSelectionInputs &, const ClusterKey &, const std::vector<SampleRef> &);
    void BeginSearch() { ++search_pass; }
    void Generated();
    void Invalid();
    void Nonmaterial();
    void TrustSkipped();
    void GuardRejected();
    void Trial(const FitStatePatch &, const ObjectiveAttemptDiagnostic &, bool polish, double factor, bool accepted);
    void Finish(bool shrink_trust_region, std::optional<double>, const ObjectiveAttemptDiagnostic &);
#else
public:
    TrustModelTrialObserver(const CandidateSelectionInputs &, const ClusterKey &, const std::vector<SampleRef> &) {}
    void BeginSearch() {}
    void Generated() {}
    void Invalid() {}
    void Nonmaterial() {}
    void TrustSkipped() {}
    void GuardRejected() {}
    void Trial(const FitStatePatch &, const ObjectiveAttemptDiagnostic &, bool, double, bool) {}
    void Finish(bool, std::optional<double>, const ObjectiveAttemptDiagnostic &) {}
#endif
};

#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
void BeginTrustModelAudit(SecondStageContext &, const std::vector<ClusterKey> &);
void FinalizeTrustModelAudit(const SecondStageContext &, const CandidateSelection &);
void LogTrustModelAudit(const SecondStageContext &, bool, const IterationResult &);
#else
inline void BeginTrustModelAudit(SecondStageContext &, const std::vector<ClusterKey> &) {}
inline void FinalizeTrustModelAudit(const SecondStageContext &, const CandidateSelection &) {}
inline void LogTrustModelAudit(const SecondStageContext &, bool, const IterationResult &) {}
#endif

#ifdef RHBM_GEM_ENABLE_SECOND_STAGE_AUDIT_TRACE

void ObservePhaseMissing(const SecondStageContext &, std::string_view, const ClusterKey &, std::string_view) noexcept;
void ObservePhaseState(const SecondStageContext &, std::string_view, const FitState &, bool probe = false) noexcept;
void ObservePhaseCandidate(const SecondStageContext &, std::string_view, const ClusterKey &, const FitStateView &,
    const FitStateView * parent = nullptr, double factor = 1.0, std::string_view disposition = "observed",
    std::string_view reason = "", bool probe = false, bool recertify = true) noexcept;
void ObservePhaseCorrection(const SecondStageContext &, std::string_view, const ClusterKey &, const FitStateView &,
    const FitStateView &, double, std::string_view, std::string_view, const CandidateSelectionInputs &,
    const std::vector<ClusterKey> &, const ObjectiveBreakdown &) noexcept;
void ObservePhaseLocalPolish(const SecondStageContext &, const ClusterKey &, const FitStateView &,
    const FitStateView &, double, bool, const ObjectiveAttemptDiagnostic &) noexcept;
void ObservePhaseSearchAssembly(const SecondStageContext &, const FitState &) noexcept;
void ObservePhaseProposal(const SecondStageContext &, const IterationProposalResult &) noexcept;
void ObservePhaseFinish(const SecondStageContext &, const FitOptions &, const std::vector<double> &,
    const SuspiciousBlockActivity &, const IterationProposalResult &, const FitState &) noexcept;
void ObservePhaseIntermediate(const SecondStageContext &, const FittedGaussianSnapshot &) noexcept;
#else

inline void ObservePhaseMissing(const SecondStageContext &, std::string_view, const ClusterKey &, std::string_view) noexcept {}
inline void ObservePhaseState(const SecondStageContext &, std::string_view, const FitState &, bool = false) noexcept {}
inline void ObservePhaseCandidate(const SecondStageContext &, std::string_view, const ClusterKey &, const FitStateView &,
    const FitStateView * = nullptr, double = 1.0, std::string_view = "observed",
    std::string_view = "", bool = false, bool = true) noexcept {}
inline void ObservePhaseCorrection(const SecondStageContext &, std::string_view, const ClusterKey &, const FitStateView &,
    const FitStateView &, double, std::string_view, std::string_view, const CandidateSelectionInputs &,
    const std::vector<ClusterKey> &, const ObjectiveBreakdown &) noexcept {}
inline void ObservePhaseLocalPolish(const SecondStageContext &, const ClusterKey &, const FitStateView &,
    const FitStateView &, double, bool, const ObjectiveAttemptDiagnostic &) noexcept {}
inline void ObservePhaseSearchAssembly(const SecondStageContext &, const FitState &) noexcept {}
inline void ObservePhaseProposal(const SecondStageContext &, const IterationProposalResult &) noexcept {}
inline void ObservePhaseFinish(const SecondStageContext &, const FitOptions &, const std::vector<double> &,
    const SuspiciousBlockActivity &, const IterationProposalResult &, const FitState &) noexcept {}
inline void ObservePhaseIntermediate(const SecondStageContext &, const FittedGaussianSnapshot &) noexcept {}
#endif

} // namespace rhbm_gem::core::detail
