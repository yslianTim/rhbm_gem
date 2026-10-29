#pragma once
#include "Numerics.hpp"
#include "SplitDiagnostics.hpp"

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
    std::optional<LocalSplitDiagnostics> split_diagnostics;
    std::size_t perturbation_replay_checks{};
    double seconds{},cpu_seconds{};
};
LocalParameterStabilityAssessment AssessLocalParameterStability(
    const Domain &,VectorRef,const Vector &,const Vector &,const EvaluationContext &,
    std::string,std::size_t=0,SplitDiagnosticsMode=SplitDiagnosticsMode::Disabled);
}
