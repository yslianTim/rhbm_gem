#pragma once

#include "core/detail/CandidateSelection.hpp"

namespace rhbm_gem::core::detail {
struct IterationResult;
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
enum class TrustModelPredictionStatus
{
    Available,
    NonmaterialStep,
    ObjectiveUnavailable,
    ModelUnavailable,
    ResidualUnavailable,
    Nonfinite,
    NonpositivePrediction,
    NonmaterialPrediction
};

enum class TrustModelCandidateSource
{
    Base,
    Polish
};

enum class TrustModelTrialDisposition
{
    Accepted,
    ObjectiveRejected
};

struct TrustModelCandidateFunnel
{
    std::size_t generated_count{ 0 };
    std::size_t invalid_count{ 0 };
    std::size_t trust_skipped_count{ 0 };
    std::size_t guard_rejected_count{ 0 };
    std::size_t nonmaterial_count{ 0 };
    std::size_t objective_evaluated_count{ 0 };
    std::size_t polish_objective_evaluated_count{ 0 };
};

struct TrustModelShadowDiagnostic
{
    TrustModelPredictionStatus status{ TrustModelPredictionStatus::ModelUnavailable };
    TrustModelCandidateSource candidate_source{ TrustModelCandidateSource::Base };
    TrustModelTrialDisposition trial_disposition{
        TrustModelTrialDisposition::ObjectiveRejected };
    std::size_t search_pass{ 0 };
    std::size_t trial_number{ 0 };
    double factor{ 0.0 };
    double step_norm{ 0.0 };
    std::optional<double> actual_reduction{};
    std::optional<double> polish_reduction{};
    std::optional<double> predicted_residual_reduction{};
    std::optional<double> predicted_penalty_reduction{};
    std::optional<double> predicted_reduction{};
    std::optional<double> rho{};
    double boundary_utilization{ 0.0 };
    TrustRegionRadiusAction current_action{ TrustRegionRadiusAction::Keep };
    std::optional<TrustRegionRadiusAction> shadow_action{};
    bool objective_backtracked{ false };
    bool rejected_by_previous{ false };
    bool rejected_by_best{ false };
    bool rejected_by_strict_polish{ false };
    bool final_local_candidate{ false };
    bool readiness_eligible{ false };
    double elapsed_milliseconds{ 0.0 };
};
#endif
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
TrustRegionRadiusAction DetermineTrustModelShadowAction(
    const TrustModelShadowDiagnostic & diagnostic);

TrustModelShadowDiagnostic EvaluateTrustModelShadow(
    const SecondStageContext & context,
    const ResidualBaseline & residual_baseline,
    const FitState & previous_state,
    const FitStatePatch & candidate_patch,
    const ClusterKey & key,
    const std::vector<SampleRef> & objective_sample_ref_list,
    const ObjectiveDomain & objective_domain,
    const std::optional<ObjectiveBreakdown> & previous_objective,
    const std::optional<ObjectiveBreakdown> & candidate_objective,
    double trust_region_radius,
    TrustRegionRadiusAction current_action,
    TrustModelCandidateSource candidate_source,
    bool objective_backtracked);
#endif

class TrustModelAudit
{
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
    struct Record
    {
        std::vector<TrustModelShadowDiagnostic> trials{};
        TrustModelCandidateFunnel funnel{};
        bool boundary_touched{ false };
    };
    std::map<ClusterKey, Record> records;
    friend class TrustModelTrialObserver;
public:
    explicit TrustModelAudit(const std::vector<ClusterKey> & keys);
    void Finalize(const CandidateSelection & selection);
    void Log(bool quiet_mode, const IterationResult & result) const;
#endif
};

class TrustModelTrialObserver
{
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
    const CandidateSelectionInputs & inputs;
    const ClusterKey & key;
    const std::vector<SampleRef> & samples;
    TrustModelAudit::Record & record;
    std::optional<std::size_t> final_trial{};
    std::size_t search_pass{ 0 };
public:
    TrustModelTrialObserver(const CandidateSelectionInputs &, const ClusterKey &, const std::vector<SampleRef> &);
    void BeginSearch() { ++search_pass; }
    void Generated() { ++record.funnel.generated_count; }
    void Invalid() { ++record.funnel.invalid_count; }
    void Nonmaterial() { ++record.funnel.nonmaterial_count; }
    void TrustSkipped() { ++record.funnel.trust_skipped_count; }
    void GuardRejected() { ++record.funnel.guard_rejected_count; }
    void Trial(const FitStatePatch &, const ObjectiveAttemptDiagnostic &, bool polish, double factor, bool accepted);
    void Finish(TrustRegionRadiusAction, std::optional<double>, const ObjectiveAttemptDiagnostic &);
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
    void Finish(TrustRegionRadiusAction, std::optional<double>, const ObjectiveAttemptDiagnostic &) {}
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
} // namespace rhbm_gem::core::detail
