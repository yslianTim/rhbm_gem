#pragma once

#include "utils/hrl/MDPDEEndpointRefinement.hpp"
#include <boost/json.hpp>
#include <string>

namespace second_stage_test {

using rhbm_gem::mdpde_detail::MDPDEEquationEvidence;
using rhbm_gem::mdpde_detail::EvaluateMDPDEEquations;
Eigen::VectorXd MDPDETestBeta(const rhbm_gem::RHBMMemberDataset &,
    const Eigen::VectorXd & weights, const std::string & backend);
double MDPDETestVariance(const rhbm_gem::RHBMMemberDataset &, double alpha,
    const Eigen::VectorXd & weights, const Eigen::VectorXd & beta);

struct ShapeFixture
{
    rhbm_gem::RHBMMemberDataset dataset;
    double alpha{};
    rhbm_gem::RHBMExecutionOptions options;
    rhbm_gem::RHBMBetaEstimateResult expected;
};
ShapeFixture ReadShapeFixture(const std::string & path);
boost::json::object EquationJSON(const MDPDEEquationEvidence &);

} // namespace second_stage_test
