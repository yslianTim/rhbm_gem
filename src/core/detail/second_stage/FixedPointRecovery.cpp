#include "core/detail/second_stage/FixedPointRecovery.hpp"
#include <cmath>

namespace rhbm_gem::core::detail {

bool IsRecoveryProgressAcceptable(double current, double trial, double best,
    double objective, double factor)
{
    return std::isfinite(current) && current > 0.0 && std::isfinite(trial) && trial >= 0.0 &&
        std::isfinite(best) && std::isfinite(objective) && factor > 0.0 && factor <= 1.0 &&
        trial <= (1.0 - 1.0e-3 * factor) * current &&
        !IsObjectiveDeteriorated(objective, best, kObjectiveProgressTolerance);
}

FixedPointRecoveryResult RunFixedPointRecovery(
    const SecondStageContext & context, const std::vector<ClusterKey> & keys,
    const FitState & current, const FitOptions & options, const std::vector<double> & ridge,
    const ObjectiveDomain & domain, const BestAuditState & best,
    const TrustRegionStateSet & radii, const SuspiciousBlockActivity & activity)
{
    FixedPointRecoveryResult result;
    auto & diagnostic{ result.diagnostics };
    diagnostic.attempted = true;
    diagnostic.reason = "current-operator-unqualified";
    ++diagnostic.operator_evaluations;
    result.current_operator = EvaluateNominalOperator(context, keys, current, options, ridge);
    diagnostic.current_residual = QualifiedNominalResidualMeanSquare(result.current_operator, current);
    if (!diagnostic.current_residual) return result;
    diagnostic.reason = "best-objective-unavailable";
    if (!best) return result;
    const auto reference{ EvaluateAuditObjective(domain, context, BuildSecondStageModelSnapshot(context, best->state)) };
    if (!reference) return result;
    diagnostic.best_objective = reference->GetTotalObjective();
    diagnostic.reason = "inactive-coordinates";
    for (std::size_t atom = 0; atom < current.size(); ++atom)
        if (!activity.HasActiveShape(atom) || !activity.HasActiveOffset(atom)) return result;

    const auto current_assessment{ AssessNominalOperator(result.current_operator, current) };
    const auto current_objective{ EvaluateAuditObjective(domain, context, BuildSecondStageModelSnapshot(context, current)) };
    if (current_assessment.certificate.StrictOperatorPassed() && current_objective &&
        !IsObjectiveDeteriorated(current_objective->GetTotalObjective(),
            *diagnostic.best_objective, kObjectiveProgressTolerance))
    {
        diagnostic.accepted = true;
        diagnostic.reason = "already-fixed";
        result.state = current;
        result.assessment = current_assessment;
        result.accepted_operator = result.current_operator;
        return result;
    }
    const auto baseline{ BuildResidualBaseline(context, current) };
    std::vector<GaussianModel3D> previous;
    for (const auto & entry : current) previous.push_back(entry.mdpde.GetModel());
    diagnostic.reason = "search-exhausted";
    for (std::size_t trial = 0; trial < 8; ++trial)
    {
        const double factor{ std::ldexp(1.0, -static_cast<int>(trial)) };
        auto & detail{ diagnostic.trials.emplace_back(RecoveryTrial{ factor, "invalid-candidate" }) };
        const auto models{ BuildDampedModelList(previous, result.current_operator.state, factor) };
        if (!models) continue;
        auto candidate{ current };
        for (std::size_t atom = 0; atom < current.size(); ++atom)
            candidate[atom].mdpde = GaussianModel3DWithUncertainty{
                (*models)[atom], current[atom].mdpde.GetStandardDeviationModel() };
        bool valid{ true };
        for (const auto & key : keys)
        {
            std::vector<GaussianModel3D> before, after;
            for (const auto atom : key) { before.push_back(previous[atom]); after.push_back((*models)[atom]); }
            const auto step_norm{ CalculateModelTrustRegionStepNorm(before, after) };
            if (!step_norm || !IsTrustRegionStepWithinRadius(*step_norm, radii.GetRadius(key)))
            { detail.reason = "trust-region"; valid = false; break; }
            const auto patch{ FitStatePatch::FromState(candidate, key) };
            const FitStateView view{ candidate, patch };
            if (EvaluateClusterCandidateGuard(context, baseline.model_snapshot, key, view, activity))
            { detail.reason = "candidate-guard"; valid = false; break; }
        }
        if (!valid) continue;
        const auto objective{ EvaluateAuditObjective(domain, context, BuildSecondStageModelSnapshot(context, candidate)) };
        detail.reason = "objective-unavailable";
        if (!objective) continue;
        detail.objective = objective->GetTotalObjective();
        detail.reason = "best-objective-bound";
        if (IsObjectiveDeteriorated(*detail.objective, *diagnostic.best_objective, kObjectiveProgressTolerance)) continue;
        ++diagnostic.operator_evaluations;
        auto endpoint{ EvaluateNominalOperator(context, keys, candidate, options, ridge) };
        detail.residual = QualifiedNominalResidualMeanSquare(endpoint, candidate);
        detail.reason = "trial-operator-unqualified";
        if (!detail.residual) continue;
        detail.reason = "insufficient-residual-progress";
        if (!IsRecoveryProgressAcceptable(*diagnostic.current_residual, *detail.residual,
                *diagnostic.best_objective, *detail.objective, factor)) continue;
        detail.reason = diagnostic.reason = "accepted";
        diagnostic.accepted = true;
        result.assessment = AssessNominalOperator(endpoint, candidate);
        result.accepted_operator = std::move(endpoint);
        result.state = std::move(candidate);
        return result;
    }
    return result;
}

} // namespace rhbm_gem::core::detail
