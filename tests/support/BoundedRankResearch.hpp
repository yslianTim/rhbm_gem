#pragma once
#include "core/detail/joint_component/Numerics.hpp"

namespace second_stage_test::bounded_rank {
namespace runtime=rhbm_gem::core::joint_component;
using Vector=Eigen::VectorXd;

struct StreamingRankEstimate
{
    bool available{},gershgorin_suggests_full_rank{},certified_full_rank{};
    Eigen::Index rows{},columns{},nonzero_entries{};
    double lambda_min_lower{runtime::unavailable},lambda_max_upper{runtime::unavailable},bound_gap{runtime::unavailable};
    double relative_threshold{},threshold_upper{runtime::unavailable};
    std::size_t estimated_peak_bytes{};
    std::optional<std::size_t> current_rss_before_bytes,current_rss_after_bytes,process_peak_rss_bytes;
    double seconds{};
    std::string reason;
};

enum class SparseQROrdering {Natural,COLAMD};
const char * SparseQROrderingName(SparseQROrdering);
struct SparseQRRankResult
{
    bool preflight_passed{},valid{};
    std::string status{"Unavailable"},reason,ordering;
    Eigen::Index rows{},columns{},input_nonzeros{},rank{},r_nonzeros{};
    double threshold{},fill_in_ratio{},symbolic_seconds{},numeric_seconds{};
    std::size_t estimated_peak_bytes{},estimated_temporary_bytes{},factor_memory_bytes{};
    std::optional<std::size_t> current_rss_before_bytes,current_rss_after_bytes,process_peak_rss_bytes;
    std::optional<std::size_t> observed_rss_delta_bytes;
};

StreamingRankEstimate ScalarGershgorin(const runtime::Sparse &,const runtime::RankPolicy &,
    std::size_t resident_bytes=0,std::size_t budget=runtime::JointMemoryPlanningBudgetBytes);
runtime::Sparse SyntheticDesign(const std::string & topology,int atoms,double width=.55);
Vector NormalizeColumns(runtime::Sparse &);
// The matrix is already column-normalized and sigma_max must refer to it.
SparseQRRankResult SparseQRRank(const runtime::Sparse &,const runtime::RankPolicy &,double sigma_max,
    SparseQROrdering,std::size_t resident_bytes=0,
    std::size_t budget=runtime::JointMemoryPlanningBudgetBytes);
}
