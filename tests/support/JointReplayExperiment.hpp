#pragma once
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"

namespace second_stage_test {
namespace joint=rhbm_gem::core::joint_component;

enum class ReplayPredictionStrategy {DeltaAddDifference,DeltaSubtractThenAdd,AffectedRowExact};
enum class ReplayObjectiveStrategy {FullScan,Incremental};

struct ReplayContributor
{
    std::size_t atom_position{};
    double square{};
};
struct ReplayIncidence
{
    joint::Indices atoms;
    std::vector<std::vector<ReplayContributor>> rows;
    std::size_t StorageBytes() const;
};
struct ReplayCandidateResult
{
    joint::Vector prediction,residual;
    double objective{};
    std::size_t scratch_bytes{};
};

ReplayIncidence PrepareReplayIncidence(const joint::Domain &,const rhbm_gem::JointParameterLayout &,
    const joint::Indices & core_atoms,const joint::Indices & neighbor_atoms,const joint::Indices & affected_rows);
bool VerifyReplayIncidenceCoverage(const joint::Domain &,const rhbm_gem::JointParameterLayout &,
    const joint::Indices & affected_rows,const ReplayIncidence &);
ReplayCandidateResult EvaluateReplayCandidate(const joint::Domain &,const rhbm_gem::JointParameterLayout &,
    joint::VectorRef observations,const joint::BlockCoordinateState &,joint::VectorRef candidate_eta,
    joint::VectorRef candidate_beta,double scale,const joint::Indices & core_atoms,
    const joint::Indices & affected_rows,ReplayPredictionStrategy,ReplayObjectiveStrategy,
    const ReplayIncidence * =nullptr);

}
