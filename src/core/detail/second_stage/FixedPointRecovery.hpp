#pragma once

#include "core/detail/second_stage/ConvergenceCertificate.hpp"
#include "core/detail/second_stage/IterationProposal.hpp"
#include "core/detail/second_stage/CandidateTransaction.hpp"

namespace rhbm_gem::core::detail {

struct RecoveryTrial
{
    double factor{ 0.0 };
    std::string_view reason{ "not-evaluated" };
    std::optional<double> residual{}, objective{};
};

struct RecoveryDiagnostics
{
    bool attempted{ false }, accepted{ false };
    std::size_t operator_evaluations{ 0 };
    std::string_view reason{ "not-attempted" };
    std::optional<double> current_residual{}, best_objective{};
    std::vector<RecoveryTrial> trials{};
};

struct FixedPointRecoveryResult
{
    RecoveryDiagnostics diagnostics{};
    std::optional<FitState> state{};
    std::optional<ConvergenceAssessment> assessment{};
    FixedPointOperatorEvidence current_operator{}, accepted_operator{};
};

bool IsRecoveryProgressAcceptable(double current_residual, double trial_residual,
    double best_objective, double trial_objective, double factor);

FixedPointRecoveryResult RunFixedPointRecovery(
    const SecondStageContext &, const std::vector<ClusterKey> &, const FitState &,
    const FitOptions &, const std::vector<double> & ridge, const ObjectiveDomain &,
    const BestAuditState &, const TrustRegionStateSet &, const SuspiciousBlockActivity &);

} // namespace rhbm_gem::core::detail
