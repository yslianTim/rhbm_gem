#pragma once

#include "core/detail/SecondStageFitting.hpp"

#include <limits>

namespace rhbm_gem::core { struct FitOptions; }

namespace rhbm_gem::core::detail {

enum class SuspiciousGaussianReason
{
    None,
    InvalidModel,
    NonFiniteResponse,
    OffsetMagnitude,
    CenterSignFlip,
    RadialRebound,
    WidthGrowth,
    AmplitudeOffsetCompensation
};

using SuspiciousUpdateMask = std::vector<char>;

struct SuspiciousBlockActivity
{
    SuspiciousUpdateMask shape_fixed_atom_mask{};
    SuspiciousUpdateMask offset_fixed_atom_mask{};
    SuspiciousUpdateMask hard_failure_atom_mask{};

    SuspiciousUpdateMask BuildCombinedFixedAtomMask() const;
    bool HasActiveShape(std::size_t atom_index) const;
    bool HasActiveOffset(std::size_t atom_index) const;
};

struct SuspiciousProfileAnalysis
{
    bool all_responses_finite{ true };
    double distance_range{ 0.0 };
    double max_abs_response{ 0.0 };
    double robust_residual_scale{ 0.0 };
    std::vector<double> radius_response_median_list{};
};

struct SuspiciousUpdateBaseline
{
    GaussianModel3D previous_model{};
    SuspiciousProfileAnalysis previous_analysis{};
};

enum class SuspiciousUpdateMode
{
    OffsetOnly,
    PostRefit
};

struct SuspiciousGaussianAssessment
{
    SuspiciousGaussianReason reason{ SuspiciousGaussianReason::None };
    double normalized_margin{ -std::numeric_limits<double>::infinity() };

    bool IsSuspicious() const { return reason != SuspiciousGaussianReason::None; }
};

SuspiciousGaussianAssessment AssessSuspiciousGaussianUpdate(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & candidate_model,
    const FitOptions & options,
    const SuspiciousUpdateBaseline & previous_baseline,
    SuspiciousUpdateMode mode);

SuspiciousUpdateBaseline BuildPreviousSuspiciousProfileBaseline(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & previous_model,
    const FitOptions & options);

enum class StabilizationTerminalReason
{
    None,
    GuardInfeasible,
    ObjectiveExhausted,
    InvalidCandidate
};

struct StabilizationTerminalDiagnostic
{
    StabilizationTerminalReason reason{ StabilizationTerminalReason::None };
    std::optional<std::size_t> guard_atom_index{};
    std::optional<SuspiciousUpdateMode> guard_mode{};
    std::optional<SuspiciousGaussianReason> guard_reason{};
};

std::optional<StabilizationTerminalDiagnostic> EvaluateClusterCandidateGuard(
    const SecondStageContext & context,
    const FitOptions & options,
    const SecondStageModelSnapshot & previous_snapshot,
    const ClusterKey & key,
    const FitStateView & candidate_state,
    const SuspiciousBlockActivity & block_activity);

std::size_t CountSuspiciousPolishAtoms(
    const SecondStageContext & context,
    const FitOptions & options,
    const std::vector<std::size_t> & atom_index_list,
    const FitStateView & endpoint_state,
    const FitStateView & candidate_state);

} // namespace rhbm_gem::core::detail
