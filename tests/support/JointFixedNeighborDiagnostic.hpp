#pragma once
#include "core/detail/joint_component/FixedNeighborBlockCoordinate.hpp"

namespace rhbm_gem::core::joint_component {
struct FixedNeighborDiagnosticResult
{
    BlockCoordinateState state;
    std::vector<FixedNeighborBlockRecord> blocks;
    std::vector<FixedNeighborBlockSweep> sweeps;
    bool search_converged{},endpoint_certified{};
    std::string reason;
    std::size_t first_order_stationarity_sweep{},confirmed_stationarity_sweep{};
    std::size_t block_preparations{},domain_preparations{},mapping_preparations{};
    FixedNeighborWork work;
    Assessment assessment;
    TrustEvidence endpoint_trust;
    JointFitResult fit;
};

FixedNeighborDiagnosticResult RunFixedNeighborDiagnostic(const JointProblem &,VectorRef initial_eta,
    const FixedNeighborPolicy & = {},Vector initial_beta = {});
}
