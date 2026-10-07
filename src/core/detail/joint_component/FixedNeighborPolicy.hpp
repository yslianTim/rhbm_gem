#pragma once

#include <cstddef>

namespace rhbm_gem::core::joint_component {

enum class FixedNeighborBlockOrder {Forward,Reverse};
enum class FixedNeighborLocalWork {Full,OneAcceptedUpdate,TwoAcceptedUpdates};
enum class FixedNeighborLocalSearch {LegacyCompact,OperatorPcg};
enum class FixedNeighborLocalPreconditioner {Identity,Diagonal,Schwarz};

// The production policy contains only the qualified route controls. Numerical
// tolerances and diagnostic switches remain implementation details.
struct FixedNeighborSearchPolicy
{
    std::size_t core_atoms{64};
    std::size_t maximum_sweeps{30};
    FixedNeighborBlockOrder order{FixedNeighborBlockOrder::Forward};
    FixedNeighborLocalWork local_work{FixedNeighborLocalWork::OneAcceptedUpdate};
};

} // namespace rhbm_gem::core::joint_component
