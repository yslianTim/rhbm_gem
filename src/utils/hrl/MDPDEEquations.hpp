#pragma once

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

namespace rhbm_gem::mdpde_detail {

struct MDPDEEquationEvidence
{
    Eigen::VectorXd raw, scaled, weights, singular_values;
    double denominator{}, condition{};
    int rank{}, floor_count{};
    bool valid{};
    std::string reason;
};

MDPDEEquationEvidence EvaluateMDPDEEquations(const RHBMMemberDataset &,
    double alpha, const Eigen::VectorXd & beta, double variance, double floor);

} // namespace rhbm_gem::mdpde_detail
