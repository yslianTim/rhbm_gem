#pragma once

#include "utils/hrl/MDPDEEndpointRefinement.hpp"
#include <string>

namespace second_stage_test {

using rhbm_gem::mdpde_detail::MDPDEEquationEvidence;
using rhbm_gem::mdpde_detail::EvaluateMDPDEEquations;
Eigen::VectorXd MDPDETestBeta(const rhbm_gem::RHBMMemberDataset &,
    const Eigen::VectorXd & weights, const std::string & backend);
double MDPDETestVariance(const rhbm_gem::RHBMMemberDataset &, double alpha,
    const Eigen::VectorXd & weights, const Eigen::VectorXd & beta);

} // namespace second_stage_test
