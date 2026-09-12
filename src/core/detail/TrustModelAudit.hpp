#pragma once

#include "core/detail/CandidateSelection.hpp"
#include "core/detail/SecondStageObservation.hpp"

namespace rhbm_gem::core::detail {
struct IterationResult;
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
enum class TrustRegionRadiusAction
{
    Keep,
    Grow, // Diagnostic shadow action only; production never grows radii.
    Shrink
};

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

#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
struct TrustModelTrialRecord
{
    std::vector<TrustModelShadowDiagnostic> trials{};
    TrustModelCandidateFunnel funnel{};
    bool boundary_touched{ false };
};
#endif

class TrustModelAudit
{
#ifdef RHBM_GEM_ENABLE_TRUST_MODEL_EXPERIMENT
    std::map<ClusterKey, TrustModelTrialRecord> records;
    friend class TrustModelTrialObserver;
public:
    explicit TrustModelAudit(const std::vector<ClusterKey> & keys);
    void Finalize(const CandidateSelection & selection);
    void Log(bool quiet_mode, const IterationResult & result) const;
#endif
};

} // namespace rhbm_gem::core::detail
