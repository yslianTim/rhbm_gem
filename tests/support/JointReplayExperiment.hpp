#pragma once
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"

namespace second_stage_test {
namespace joint=rhbm_gem::core::joint_component;

enum class ReplayPredictionStrategy {DeltaAddDifference,DeltaSubtractThenAdd,AffectedRowExact};
enum class ReplayObjectiveStrategy
{
    FullScan,Incremental,IncrementalScaled,IncrementalKahan,IncrementalNeumaier,
    IncrementalKahanScaled,IncrementalNeumaierScaled,CachedExternally
};

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

struct ObjectiveCacheCandidate
{
    double objective{};
    std::size_t tiles_recomputed{};
    std::vector<std::pair<std::size_t,double>> updates;
};

class TiledObjectiveCache
{
public:
    TiledObjectiveCache(const std::vector<std::size_t> & informative_rows,joint::VectorRef residuals,
        std::size_t tile_size,bool hierarchical);
    ObjectiveCacheCandidate Evaluate(joint::VectorRef accepted_residuals,
        const joint::Indices & affected_rows,joint::VectorRef candidate_residuals,double scale) const;
    void Commit(const ObjectiveCacheCandidate &);
    std::size_t StorageBytes() const;
    std::size_t TileCount() const {return tile_count_;}

private:
    const std::vector<std::size_t> * informative_rows_{};
    std::size_t tile_size_{},tile_count_{},leaf_capacity_{};
    bool hierarchical_{};
    std::vector<double> values_;
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
