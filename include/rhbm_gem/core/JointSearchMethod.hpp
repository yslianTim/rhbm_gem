#pragma once

namespace rhbm_gem::core {

enum class JointSearchMethod : int
{
    LegacyCompact = 0,
    OperatorPcg = 1,
    FixedNeighbor = 2
};

} // namespace rhbm_gem::core
