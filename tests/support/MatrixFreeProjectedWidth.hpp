#pragma once
#include "ProjectedWidthResearch.hpp"

namespace second_stage_test::projected_width {
struct OperatorResult
{
    runtime::Sparse normalized_design,derivative;
    Matrix basis;
    Vector column_norms;
    Eigen::Index design_rank{},width_columns{};
    double design_threshold{};
    std::size_t estimated_peak_bytes{},basis_bytes{};
    mutable std::size_t forward_calls{},adjoint_calls{},projection_calls{};
    bool available{};
    std::string status{"Unavailable"},reason;

    Vector Project(const Eigen::Ref<const Vector> &) const;
    Vector Apply(const Eigen::Ref<const Vector> &) const;
    Vector ApplyAdjoint(const Eigen::Ref<const Vector> &) const;
    double OrthogonalityResidual(const Eigen::Ref<const Vector> &) const;
};

OperatorResult BuildExactOracleOperator(const runtime::Sparse & design,const runtime::Sparse & derivative,
    const runtime::RankPolicy &,std::size_t resident_bytes=0,
    std::size_t budget=runtime::JointMemoryPlanningBudgetBytes);

struct BoundedProjectionPlan
{
    std::size_t estimated_action_bytes{},estimated_peak_bytes{},budget_bytes{};
    std::string resource_status{"Unavailable"},numerical_status{"Unavailable"},reason;
};
BoundedProjectionPlan PreflightBoundedProjection(std::size_t rows,std::size_t free_columns,
    std::size_t width_columns,std::size_t design_nonzeros,std::size_t derivative_nonzeros,
    std::size_t resident_bytes,std::size_t budget=runtime::JointMemoryPlanningBudgetBytes);
}
