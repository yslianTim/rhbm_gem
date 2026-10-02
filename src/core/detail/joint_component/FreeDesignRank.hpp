#pragma once
#include "CompactSvd.hpp"
#include "SparseFactor.hpp"
#include <optional>
#include <string_view>

namespace rhbm_gem::core::joint_component {
enum class FreeDesignRankStatus {Unavailable,FullRank,Deficient};
enum class FreeDesignRankWorkStage
{
    None,StructuralScan,DuplicateCheck,LocalWitness,FactorInspection,WeakDirection,InverseBound,OrthogonalBound,Reconstruction
};
enum class FreeDesignRankCertificate {None,Structural,WeakDirection,LocalSupport,SpqrReconstruction,DenseOracle};
std::string_view FreeDesignRankWorkStageName(FreeDesignRankWorkStage);
std::string_view FreeDesignRankCertificateName(FreeDesignRankCertificate);
std::optional<double> CertifiedMagnitudeLowerBound(double);
std::optional<double> CertifiedSmallestSingularLowerBound2x2(double,double,double,double);
struct FreeDesignLocalWitness
{
    std::size_t groups{},covered_columns{},total_columns{},exclusive_rows{},max_group_size{};
    double coverage_fraction{},threshold_upper{};
    std::optional<double> minimum_lower;
    bool exclusive_rows_disjoint{true},would_certify{};
    std::string_view reason{"partial-column-coverage"};
};
FreeDesignLocalWitness DiagnoseFreeDesignLocalWitnesses(const Sparse &,double threshold_upper);
struct FreeDesignRankResult
{
    FreeDesignRankStatus status{FreeDesignRankStatus::Unavailable};
    FreeDesignRankCertificate certificate{FreeDesignRankCertificate::None};
    std::string reason{"rank-unavailable"};
    FreeDesignRankWorkStage work_stage{FreeDesignRankWorkStage::None};
    Eigen::Index rank_lower{},rank_upper{};
    double minimum_lower{},maximum_lower{},maximum_upper{},threshold_lower{},threshold_upper{};
    double reconstruction_error{unavailable},orthogonal_minimum{unavailable},witness_upper{unavailable},seconds{};
    double local_witness_seconds{};
    FreeDesignLocalWitness local_witness;
    std::optional<std::size_t> estimated_total_entries,estimated_remaining_entries,estimated_reconstruction_entries;
    std::size_t design_nonzeros{},factor_r_nonzeros{},reflector_nonzeros{},reflector_count{};
    std::size_t entries{},workspace_bytes{};
};
// No dense fallback. A null factor permits structural and local-support certificates only.
FreeDesignRankResult EvaluateFreeDesignRank(const Sparse &,const FreeDesignFactor *,const RankRequest &,const RankBudget & = {});
}
