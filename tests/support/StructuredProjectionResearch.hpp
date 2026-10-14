#pragma once
#include "ProjectedWidthResearch.hpp"

namespace second_stage_test::projected_width {
struct StructuredProjectionResult
{
    Vector projection;
    Eigen::Index rows{},columns{},input_nonzeros{},gram_nonzeros{},factor_nonzeros{},
        exact_blocks{},articulation_columns{},separator_columns{},maximum_interior_columns{},
        interface_dimension{},maximum_factor_size{};
    std::size_t estimated_peak_bytes{};
    std::optional<std::size_t> current_rss_before_bytes,current_rss_after_bytes,process_peak_rss_bytes;
    double projection_residual{std::numeric_limits<double>::quiet_NaN()},
        idempotence_error{std::numeric_limits<double>::quiet_NaN()},
        diagonal_condition_estimate{std::numeric_limits<double>::quiet_NaN()},seconds{};
    bool preflight_passed{},valid{};
    std::string status{"Unavailable"},reason;
    std::string rank_status{"Unavailable"};
};

struct StructuredProjectionPreflight
{
    std::size_t estimated_peak_bytes{},budget_bytes{};
    std::string status{"Unavailable"},reason;
};
StructuredProjectionPreflight PreflightStructuredProjection(std::size_t rows,std::size_t columns,
    std::size_t nonzeros,std::size_t resident_bytes,
    std::size_t budget=runtime::JointMemoryPlanningBudgetBytes);

StructuredProjectionResult ProjectWithSparseNormalEquations(const runtime::Sparse & design,
    const Eigen::Ref<const Vector> & observation_direction,std::size_t resident_bytes=0,
    std::size_t budget=runtime::JointMemoryPlanningBudgetBytes);
}
