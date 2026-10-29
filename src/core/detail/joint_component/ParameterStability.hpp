#pragma once
#include "Numerics.hpp"

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
LocalParameterStabilityClassification ClassifyLocalParameterStability(
    double,double,double,bool,bool,double=LocalParameterChangeWarning,
    double=LocalPredictionChangeWarning);
struct LocalParameterStabilityAssessment
{
    JointBlockParameterStability evidence;
    MemoryEstimate memory_estimate;
    double seconds{},cpu_seconds{};
};
LocalParameterStabilityAssessment AssessLocalParameterStability(
    const Domain &,VectorRef,const Vector &,const Vector &,const EvaluationContext &,
    std::string,std::size_t=0);
}
