#pragma once

#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <boost/json.hpp>
#include <string>

namespace second_stage_test {

struct MDPDEEquationEvidence
{
    Eigen::VectorXd raw, scaled, weights, singular_values;
    double denominator{}, condition{};
    int rank{}, floor_count{};
    bool valid{};
    std::string reason;
};

// Implemented beside the production primitives, only with BUILD_TESTING.
MDPDEEquationEvidence EvaluateMDPDEEquations(const rhbm_gem::RHBMMemberDataset &,
    double alpha, const Eigen::VectorXd & beta, double variance, double floor);
Eigen::VectorXd MDPDETestBeta(const rhbm_gem::RHBMMemberDataset &,
    const Eigen::VectorXd & weights, const std::string & backend);
double MDPDETestVariance(const rhbm_gem::RHBMMemberDataset &, double alpha,
    const Eigen::VectorXd & weights, const Eigen::VectorXd & beta);
void RecordMDPDEIteration(const Eigen::VectorXd &, double, double, double);

struct ShapeFixture
{
    rhbm_gem::RHBMMemberDataset dataset;
    double alpha{};
    rhbm_gem::RHBMExecutionOptions options;
    rhbm_gem::RHBMBetaEstimateResult expected;
};
ShapeFixture ReadShapeFixture(const std::string & path);
boost::json::object CompareMDPDE(const ShapeFixture &, bool perturb = true);
boost::json::object EquationJSON(const MDPDEEquationEvidence &);

} // namespace second_stage_test
