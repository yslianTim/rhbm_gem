#pragma once
#include "support/UniqueStencilGrid.hpp"
#include <rhbm_gem/data/object/MapObject.hpp>
#include <memory>

namespace second_stage_test::matched::fixed_b {
struct Data
{
    unique_grid::Grid grid;
    std::vector<Atom> atoms;
    std::unique_ptr<rhbm_gem::MapObject> generation;
    Eigen::VectorXd y64, y32, checkpoint_b;
    boost::json::array identities;
};
// Truth is used only inside preparation for forward validation and scoring output.
Data Prepare(const std::string & manifest, const std::string & map,
    const std::string & checkpoint, const std::string & output,
    const std::string & experiment = "fixed-b-oracle");
// Mean-only LS certificate: exact residuals require no positive variance.
boost::json::object Certificate(const Eigen::SparseMatrix<double> &, const Eigen::VectorXd &,
    const Eigen::VectorXd &);
boost::json::object Fit(const Eigen::SparseMatrix<double> &, const Eigen::VectorXd &,
    const boost::json::object & spectrum);
void Run(const std::string & manifest, const std::string & map,
    const std::string & checkpoint, const std::string & output);
} // namespace second_stage_test::matched::fixed_b
