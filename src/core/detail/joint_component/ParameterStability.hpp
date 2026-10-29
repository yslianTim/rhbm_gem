#pragma once
#include "Numerics.hpp"
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <functional>
#endif

namespace rhbm_gem::core::joint_component {
inline constexpr std::size_t MaximumLocalPerturbationDirectionsPerBlock=3;
inline constexpr std::size_t LocalPerturbationSeedsPerDirection=2;
inline constexpr std::size_t MaximumLocalPerturbationProfileEvaluations=12;
inline constexpr std::size_t MaximumPerturbedBlocksPerFit=8;
inline constexpr std::size_t MaximumRepresentativeLocalStabilityProbesPerComponent=3;
inline constexpr double LocalWidthPerturbation=1e-3;
inline constexpr double LocalParameterChangeWarning=1e-8;
inline constexpr double LocalPredictionChangeWarning=1e-8;
inline constexpr double MaximumLocalPerturbationCpuSeconds=2.0;
struct LocalParameterStabilityClassification
{
    std::string status,prediction_status;
};
// Pure classification of completed perturbation measurements; does not perform solver or replay work.
LocalParameterStabilityClassification ClassifyLocalParameterStability(
    double,double,double,bool,bool,double=LocalParameterChangeWarning,
    double=LocalPredictionChangeWarning);
struct LocalParameterStabilityAssessment
{
    JointBlockParameterStability evidence;
    MemoryEstimate memory_estimate;
    double seconds{},cpu_seconds{};
};
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct LocalStabilityTraceForTesting
{
    std::string stage,block_id,direction,sign,stop_reason,qualification_status;
    std::vector<std::string> core_atom_ids;
    std::vector<double> direction_vector;
    Vector baseline_eta,baseline_beta,seed_eta,seed_beta,requested_displacement,effective_displacement,
        endpoint_eta,endpoint_beta;
    Certificate baseline_certificate,certificate;
    std::optional<IndependentReplayEvidence> baseline_replay,replay;
    std::optional<double> baseline_width_kkt,width_kkt,parameter_change,width_parameter_change,
        scaled_ac_parameter_change,prediction_change,objective_change,objective_difference;
    std::vector<double> objective_trajectory;
    int evaluations{},accepted_updates{};
    bool feasible{},width_bound_feasible{true},budget_exhausted{};
    double seconds{},cpu_seconds{};
};
struct LocalParameterStabilityResearchHooks
{
    JointConvergenceCriteria eligibility_criteria{PracticalJointConvergenceCriteria};
    std::function<void(const LocalStabilityTraceForTesting &)> trace_observer;
};
#endif
LocalParameterStabilityAssessment AssessLocalParameterStability(
    const Domain &,VectorRef,const Vector &,const Vector &,const EvaluationContext &,
    std::string,std::size_t=0
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    ,const LocalParameterStabilityResearchHooks * = nullptr
#endif
    );
}
