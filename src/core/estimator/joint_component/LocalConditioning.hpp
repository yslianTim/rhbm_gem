#pragma once
#include "Numerics.hpp"

namespace rhbm_gem::core::joint_component {
// Provisional roundoff-amplification warning; this is diagnostic evidence only.
inline constexpr double LocalConditioningRoundoffAmplificationWarning=1e-8;
struct LocalConditioningAssessment
{
    JointBlockConditioning evidence;
    MemoryEstimate memory_estimate;
    double seconds{};
};
LocalConditioningAssessment AssessBcdLocalConditioning(
    const Domain &,VectorRef,const Vector &,const Vector &,const EvaluationContext &,
    std::string,std::size_t=0);

struct AcCompensationAssessment
{
    std::string reason,raw_reason,conditional_reason,cancellation_reason,weak_reason;
    bool amplitude_active{},width_upper_bound_active{},conditional_on_active_face{};
    std::optional<double> raw_correlation,conditional_correlation,cancellation_ratio;
    std::optional<double> gaussian_norm,charge_norm,residualized_gaussian_norm,residualized_charge_norm;
    std::optional<Eigen::Index> free_rank,nuisance_rank;
    bool free_rank_threshold_sensitive{},nuisance_rank_threshold_sensitive{};
    Vector singular_values;
    Eigen::Index weak_dimension{};
    Vector weak_atom_shares;
    std::optional<double> weak_cancellation_alignment;
    // Available for an isolated weakest direction. Entries follow beta order
    // in physical units, with zeros at active amplitudes.
    Vector physical_weak_direction;
    MemoryEstimate memory_estimate;
};

// Explicit research diagnostic. Supply the accepted block's Fitting Domain,
// endpoint parameters and parent scale. Other blocks are held fixed.
AcCompensationAssessment AssessAcCompensation(
    const Domain &,VectorRef,const Vector &,const Vector &,const EvaluationContext &,
    Eigen::Index target,std::size_t resident_bytes=0);
}
