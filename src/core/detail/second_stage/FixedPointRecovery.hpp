#pragma once

#include "core/detail/second_stage/ConvergenceCertificate.hpp"
#include "core/detail/second_stage/IterationProposal.hpp"
#include "core/detail/second_stage/CandidateTransaction.hpp"

#include <vector>

namespace rhbm_gem::core::detail {

struct FixedPointRecoveryResult
{
    std::optional<FitState> state{};
    std::optional<ConvergenceAssessment> assessment{};
};

bool IsRecoveryProgressAcceptable(double current_residual, double trial_residual,
    double best_objective, double trial_objective, double factor);

FixedPointRecoveryResult RunFixedPointRecovery(
    const SecondStageContext &, const std::vector<ClusterKey> &, const FitState &,
    const FitOptions &, const std::vector<double> & ridge, const ObjectiveDomain &,
    const BestAuditState &, const TrustRegionStateSet &, const SuspiciousBlockActivity &);

} // namespace rhbm_gem::core::detail
