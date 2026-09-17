#pragma once
#include "support/MatchedJointAC.hpp"

namespace second_stage_test::matched::unique_grid {
inline constexpr std::array<double,4> alphas{0,.1,.5,1};
inline constexpr int iteration_budget{100}, refinement_budget{100};
struct Voxel
{
    std::size_t index{}, multiplicity{};
    Position position{};
    double observed{}, nearest_distance{};
    bool in_stencil{};
};
struct Grid
{
    std::vector<Voxel> voxels;
    std::vector<std::array<std::size_t,64>> sample_rows;
};
Grid BuildGrid(const std::vector<Stencil> &, const rhbm_gem::MapObject &);
Eigen::MatrixXd BuildDesign(const Grid &, const std::vector<Atom> &, double cutoff);
Eigen::VectorXd Project(const Grid &, const std::vector<Stencil> &, const Eigen::VectorXd &);
double Direct(const Position &, const std::vector<Atom> &, double cutoff, double * absolute = nullptr);
joint_ac::Blocks GlobalBlock(Eigen::Index rows, double alpha);
boost::json::object Fit(const Eigen::MatrixXd &, const Eigen::VectorXd &, const Eigen::VectorXd &, double alpha,
    const Eigen::SparseMatrix<double> * sparse_design = nullptr);
boost::json::object Fit(const Eigen::SparseMatrix<double> &, const Eigen::VectorXd &, const Eigen::VectorXd &, double alpha);
void RunUnion(const std::string &, const std::string &, const std::string &, const std::string &);
void Run(const std::string & manifest, const std::string & map,
    const std::string & state_index, const std::string & output);
} // namespace second_stage_test::matched::unique_grid
