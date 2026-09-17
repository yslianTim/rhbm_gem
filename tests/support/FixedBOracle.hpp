#pragma once
#include "support/MatchedJointAC.hpp"

namespace second_stage_test::matched::fixed_b {
// Mean-only LS certificate: exact residuals require no positive variance.
boost::json::object Certificate(const Eigen::SparseMatrix<double> &, const Eigen::VectorXd &,
    const Eigen::VectorXd &);
boost::json::object Fit(const Eigen::SparseMatrix<double> &, const Eigen::VectorXd &,
    const boost::json::object & spectrum);
void Run(const std::string & manifest, const std::string & map,
    const std::string & checkpoint, const std::string & output);
} // namespace second_stage_test::matched::fixed_b
