#pragma once
#include "CompactSvd.hpp"
#include "SparseFactor.hpp"

namespace rhbm_gem::core::joint_component {
enum class FreeDesignRankBackend {Dense,SpqrBounds};
enum class FreeDesignRankStatus {Unavailable,FullRank,Deficient};
struct RankBudget
{
    double seconds{120};
    std::size_t entries{100000000},workspace_bytes{256*1024*1024};
};
struct FreeDesignRankResult
{
    FreeDesignRankStatus status{FreeDesignRankStatus::Unavailable};
    std::string reason{"rank-unavailable"};
    Eigen::Index rank_lower{},rank_upper{};
    double minimum_lower{},maximum_lower{},maximum_upper{},threshold_lower{},threshold_upper{};
    double reconstruction_error{unavailable},orthogonal_minimum{unavailable},witness_upper{unavailable},seconds{};
    std::size_t entries{},workspace_bytes{};
};
// No dense fallback. A null factor permits structural checks only.
FreeDesignRankResult EvaluateFreeDesignRank(const Sparse &,const FreeDesignFactor *,const RankRequest &,const RankBudget & = {});
}
