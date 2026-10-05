#pragma once
#include "FixedBBlockCoordinate.hpp"
#include <functional>

namespace rhbm_gem::core::joint_component {
enum class FixedNeighborBlockOrder {Forward,Reverse};
struct FixedNeighborBlockRecord
{
    std::size_t sweep{},block{},affected_rows{};
    Indices core_atoms;
    double objective_before{},objective_after{},local_objective_before{},local_objective_after{};
    double global_replay_delta{},local_global_delta_error{},objective_replay_enclosure{};
    double local_final_gradient_inf_norm{},local_final_ac_kkt{};
    double search_seconds{};
    int profile_evaluations{},accepted_updates{};
    bool accepted{};
    std::string status,reason,local_search_stop_reason;
};
struct FixedNeighborBlockSweep
{
    double objective_before{},objective_after{},global_a_feasibility{},global_ac_kkt{};
    double global_width_gradient_inf_norm{},cache_replay_error{},objective_replay_error{},wall_seconds{};
    std::size_t block_solves{},profile_evaluations{},maximum_block_rows{},maximum_block_columns{};
};
struct FixedNeighborPolicy
{
    std::size_t core_atoms{128},maximum_sweeps{30};
    FixedNeighborBlockOrder order{FixedNeighborBlockOrder::Forward};
    bool stop_after_stationarity{true};
    std::function<void(const FixedNeighborBlockSweep &)> sweep_observer;
    std::function<void(std::size_t,const BlockCoordinateState &,const FixedNeighborBlockSweep &,
        const std::vector<FixedNeighborBlockRecord> &)> state_observer;
};
struct FixedNeighborResult
{
    BlockCoordinateState state;
    std::vector<FixedNeighborBlockRecord> blocks;
    std::vector<FixedNeighborBlockSweep> sweeps;
    bool search_converged{},endpoint_certified{};
    std::string reason;
    std::size_t sweeps_to_stationarity{};
    Assessment assessment;
    TrustEvidence endpoint_trust;
    JointFitResult fit;
};

FixedNeighborResult SearchFixedNeighbor(const JointProblem &,VectorRef initial_eta,
    const FixedNeighborPolicy & = {});
}
