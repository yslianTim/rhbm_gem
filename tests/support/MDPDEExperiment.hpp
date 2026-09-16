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
rhbm_gem::RHBMDiagonalMatrix MDPDETestCovariance(double variance, const Eigen::VectorXd & weights);
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

struct EndpointRefinementResult
{
    rhbm_gem::RHBMBetaEstimateResult result;
    bool accepted{};
    boost::json::object evidence;
};

// The reference always continues the original endpoint; it is never a fallback.
EndpointRefinementResult RefineMDPDEEndpoint(const ShapeFixture &, int equation_budget,
    const MDPDEEquationEvidence * initial = nullptr);
boost::json::object CompareMDPDEBranches(const Eigen::VectorXd & candidate_beta, double candidate_variance,
    const MDPDEEquationEvidence & candidate, const Eigen::VectorXd & reference_beta, double reference_variance,
    const MDPDEEquationEvidence & reference, double floor);

} // namespace second_stage_test
