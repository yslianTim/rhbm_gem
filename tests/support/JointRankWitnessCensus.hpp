#pragma once
#include "core/detail/joint_component/FreeDesignRank.hpp"

namespace second_stage_test {
namespace n=rhbm_gem::core::joint_component;
using LocalRankWitnessCensus=n::FreeDesignLocalWitness;
inline LocalRankWitnessCensus DiagnoseLocalRankWitnesses(const n::Sparse & design,double threshold_upper)
{
    return n::DiagnoseFreeDesignLocalWitnesses(design,threshold_upper);
}
}
