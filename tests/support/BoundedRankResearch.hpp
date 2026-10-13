#pragma once
#include "core/detail/joint_component/Numerics.hpp"

namespace second_stage_test::bounded_rank {
namespace runtime=rhbm_gem::core::joint_component;

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

StreamingRankEstimate ScalarGershgorin(const runtime::Sparse &,const runtime::RankPolicy &,
    std::size_t resident_bytes=0,std::size_t budget=runtime::JointMemoryPlanningBudgetBytes);
}
