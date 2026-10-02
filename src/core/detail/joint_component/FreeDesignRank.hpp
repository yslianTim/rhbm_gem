#pragma once
#include "CompactSvd.hpp"
#include "SparseFactor.hpp"
#include <optional>
#include <string_view>

namespace rhbm_gem::core::joint_component {
enum class FreeDesignRankStatus {Unavailable,FullRank,Deficient};
enum class FreeDesignRankWorkStage
{
    None,StructuralScan,DuplicateCheck,FactorInspection,WeakDirection,InverseBound,OrthogonalBound,Reconstruction
};
std::string_view FreeDesignRankWorkStageName(FreeDesignRankWorkStage);
struct FreeDesignRankResult
{
    FreeDesignRankStatus status{FreeDesignRankStatus::Unavailable};
    std::string reason{"rank-unavailable"};
    FreeDesignRankWorkStage work_stage{FreeDesignRankWorkStage::None};
    Eigen::Index rank_lower{},rank_upper{};
    double minimum_lower{},maximum_lower{},maximum_upper{},threshold_lower{},threshold_upper{};
    double reconstruction_error{unavailable},orthogonal_minimum{unavailable},witness_upper{unavailable},seconds{};
    std::optional<std::size_t> estimated_total_entries,estimated_remaining_entries,estimated_reconstruction_entries;
    std::size_t design_nonzeros{},factor_r_nonzeros{},reflector_nonzeros{},reflector_count{};
    std::size_t entries{},workspace_bytes{};
};
// No dense fallback. A null factor permits structural checks only.
FreeDesignRankResult EvaluateFreeDesignRank(const Sparse &,const FreeDesignFactor *,const RankRequest &,const RankBudget & = {});
}
