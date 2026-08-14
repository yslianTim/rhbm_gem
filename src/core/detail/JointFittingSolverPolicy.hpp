#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace rhbm_gem::core::detail {

constexpr double kRobustLossCutoffMultiplier{ 1.345 };
constexpr double kJointOffsetRidgeRatio{ 1.0e-3 };
constexpr double kCollinearJointOffsetRidgeMultiplier{ 10.0 };
constexpr double kJointOffsetConditioningPivotRatioThreshold{ 1.0e-8 };

struct JointOffsetConditioning
{
    bool guard_required{ false };
    double pivot_ratio{ 0.0 };
};

inline JointOffsetConditioning EvaluateJointOffsetConditioning(
    const Eigen::SparseMatrix<double> & design_matrix)
{
    if (design_matrix.cols() == 0)
    {
        return JointOffsetConditioning{ true, 0.0 };
    }

    Eigen::SparseMatrix<double> normalized_design{ design_matrix };
    for (Eigen::Index column = 0;
        column < normalized_design.outerSize();
        column++)
    {
        double square_sum{ 0.0 };
        for (Eigen::SparseMatrix<double>::InnerIterator entry(
                normalized_design,
                column);
            entry;
            ++entry)
        {
            square_sum += entry.value() * entry.value();
        }
        if (!std::isfinite(square_sum) ||
            square_sum <= std::numeric_limits<double>::epsilon())
        {
            return JointOffsetConditioning{ true, 0.0 };
        }
        const auto scale{ std::sqrt(square_sum) };
        for (Eigen::SparseMatrix<double>::InnerIterator entry(
                normalized_design,
                column);
            entry;
            ++entry)
        {
            entry.valueRef() /= scale;
        }
    }

    Eigen::SparseMatrix<double> normalized_gram{
        normalized_design.transpose() * normalized_design
    };
    normalized_gram.makeCompressed();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(normalized_gram);
    if (solver.info() != Eigen::Success)
    {
        return JointOffsetConditioning{ true, 0.0 };
    }

    const auto diagonal{ solver.vectorD().eval() };
    if (diagonal.size() == 0 || !diagonal.allFinite())
    {
        return JointOffsetConditioning{ true, 0.0 };
    }
    if (diagonal.minCoeff() <= 0.0)
    {
        return JointOffsetConditioning{ true, 0.0 };
    }
    const auto maximum_pivot{ diagonal.maxCoeff() };
    const auto minimum_pivot{ diagonal.minCoeff() };
    if (!std::isfinite(maximum_pivot) ||
        maximum_pivot <= std::numeric_limits<double>::epsilon())
    {
        return JointOffsetConditioning{ true, 0.0 };
    }

    const auto pivot_ratio{ minimum_pivot / maximum_pivot };
    return JointOffsetConditioning{
        !std::isfinite(pivot_ratio) ||
            pivot_ratio <= kJointOffsetConditioningPivotRatioThreshold,
        std::isfinite(pivot_ratio) ? pivot_ratio : 0.0
    };
}

inline double CalculateJointOffsetRidgeDiagonal(
    double column_square_sum,
    double multiplier)
{
    if (!std::isfinite(multiplier) || multiplier <= 0.0)
    {
        throw std::invalid_argument(
            "Local fitting ridge multiplier must be positive and finite.");
    }
    const auto base_ridge{
        column_square_sum > std::numeric_limits<double>::epsilon() ?
            kJointOffsetRidgeRatio * column_square_sum : 1.0
    };
    return multiplier * base_ridge;
}

} // namespace rhbm_gem::core::detail
