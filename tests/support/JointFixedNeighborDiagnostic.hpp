#pragma once
#include "core/detail/joint_component/GaussSeidelBcd.hpp"

namespace rhbm_gem::core::joint_component {
struct BcdDiagnosticResult
{
    BlockCoordinateState state;
    std::vector<BcdBlockRecord> blocks;
    std::vector<BcdSweepRecord> sweeps;
    bool search_converged{},endpoint_certified{};
    std::string reason;
    std::size_t first_order_stationarity_sweep{},confirmed_stationarity_sweep{};
    std::size_t sweep_count{},total_block_solves{},total_profile_evaluations{},total_accepted_local_updates{};
    std::size_t block_preparations{},domain_preparations{},mapping_preparations{};
    std::size_t prepared_mapping_entries{},prepared_mapping_bytes{};
    MemoryEstimate memory_estimate,global_state_memory_estimate;
    BcdWork work;
    double search_seconds{};
    std::optional<BcdSweepRecord> final_sweep;
    std::vector<JointBlockConditioning> local_conditioning;
    bool certification_attempted{};
    Assessment assessment;
    TrustEvidence endpoint_trust;
    JointFitResult fit;
};

BcdDiagnosticResult RunBcdDiagnostic(const JointProblem &,VectorRef initial_eta,
    const BcdSearchPolicy & = {});
}
