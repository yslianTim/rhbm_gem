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
    const auto current_operator{ EvaluateNominalOperator(context, keys, current, options, ridge) };
    const auto current_residual{ QualifiedNominalResidualMeanSquare(current_operator, current) };
    if (!current_residual) return result;
    if (!best) return result;
    const auto reference{ EvaluateAuditObjective(domain, context, BuildSecondStageModelSnapshot(context, best->state)) };
    if (!reference) return result;
    const auto best_objective{ reference->GetTotalObjective() };
    for (std::size_t atom = 0; atom < current.size(); ++atom)
        if (!activity.HasActiveShape(atom) || !activity.HasActiveOffset(atom)) return result;

    const auto current_assessment{ AssessNominalOperator(current_operator, current) };
    const auto current_objective{ EvaluateAuditObjective(domain, context, BuildSecondStageModelSnapshot(context, current)) };
    if (current_assessment.certificate.StrictOperatorPassed() && current_objective &&
        !IsObjectiveDeteriorated(current_objective->GetTotalObjective(),
            best_objective, kObjectiveProgressTolerance))
    {
        result.state = current;
        result.assessment = current_assessment;
        return result;
    }
    const auto baseline{ BuildResidualBaseline(context, current) };
    std::vector<GaussianModel3D> previous;
    for (const auto & entry : current) previous.push_back(entry.mdpde.GetModel());
    for (std::size_t trial = 0; trial < 8; ++trial)
    {
        const double factor{ std::ldexp(1.0, -static_cast<int>(trial)) };
        const auto models{ BuildDampedModelList(previous, current_operator.state, factor) };
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
            { valid = false; break; }
            const auto patch{ FitStatePatch::FromState(candidate, key) };
            const FitStateView view{ candidate, patch };
            if (EvaluateClusterCandidateGuard(context, baseline.model_snapshot, key, view, activity))
            { valid = false; break; }
        }
        if (!valid) continue;
        const auto objective{ EvaluateAuditObjective(domain, context, BuildSecondStageModelSnapshot(context, candidate)) };
        if (!objective) continue;
        const auto trial_objective{ objective->GetTotalObjective() };
        if (IsObjectiveDeteriorated(trial_objective, best_objective, kObjectiveProgressTolerance)) continue;
        auto endpoint{ EvaluateNominalOperator(context, keys, candidate, options, ridge) };
        const auto trial_residual{ QualifiedNominalResidualMeanSquare(endpoint, candidate) };
        if (!trial_residual) continue;
        if (!IsRecoveryProgressAcceptable(*current_residual, *trial_residual,
                best_objective, trial_objective, factor)) continue;
        result.assessment = AssessNominalOperator(endpoint, candidate);
        result.state = std::move(candidate);
        return result;
    }
    return result;
}

} // namespace rhbm_gem::core::detail
