#pragma once
#include "StructuralPartition.hpp"
#include "Numerics.hpp"

namespace rhbm_gem::core::joint_component {
enum class FixedBBlockOrder {Forward,Reverse};
struct FixedBBlockPolicy
{
    std::size_t core_atoms{128},maximum_sweeps{200};
    FixedBBlockOrder order{FixedBBlockOrder::Forward};
    bool capture_diagnostics{};
};
struct FixedBBlockRecord
{
    std::size_t sweep{},block{};
    Indices core_atoms,affected_rows,neighbor_atoms;
    double objective_before{},objective_after{},objective_reduction{},max_scaled_ac_change{},factor_seconds{};
    double local_objective_before{},local_objective_after{},local_objective_reduction{};
    double local_kkt{},local_feasibility{},global_kkt_before{},global_kkt_after{};
    double raw_ac_change{},solver_beta_norm_difference{},global_replay_delta{},local_global_delta_error{},objective_replay_enclosure{};
    int local_iterations{},active_A{},free_columns{};
    int linear_solves{},rank{},factorizations{};
    std::size_t active_atoms{};
    bool diagnostics_captured{},solver_beta_equal_old{},accepted{};
    bool local_certificate_available{},local_certificate_feasible{},local_certificate_kkt_passed{};
    double local_certificate_projected_kkt{},local_certificate_rss{},local_certificate_objective{};
    std::string status,reason;
};
struct FixedBBlockSweep
{
    double objective_before{},objective_after{},relative_objective_change{},max_scaled_ac_change{};
    std::size_t accepted_blocks{},unchanged_blocks{},failed_blocks{};
    double global_a_feasibility{},global_ac_kkt{},cache_replay_error{},objective_replay_error{},wall_seconds{};
};
struct BlockCoordinateState
{
    Vector eta,beta,prediction,residual;
    double objective{};
};
struct FixedBBlockResult
{
    BlockCoordinateState state;
    std::vector<FixedBBlockRecord> blocks;
    std::vector<FixedBBlockSweep> sweeps;
    bool success{};
    std::string reason;
    std::size_t sweeps_to_global_kkt{};
};

FixedBBlockResult SearchFixedBBlocks(const JointProblemInput &,const JointParameterLayout &,VectorRef observations,
    VectorRef eta,const EvaluationContext &,const FixedBBlockPolicy & = {});
double BlockObjectiveReplayEnclosure(double reference);
bool WithinBlockObjectiveReplay(double error,double reference);
}
