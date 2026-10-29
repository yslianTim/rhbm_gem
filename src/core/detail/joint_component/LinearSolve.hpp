#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <vector>

namespace rhbm_gem::core::joint_component {
struct LinearSolveWorkspace
{
    using Sparse=Eigen::SparseMatrix<double>;
    using RowMajorSparse=Eigen::SparseMatrix<double,Eigen::RowMajor>;

    Sparse normalized,selected,reduced;
    RowMajorSparse rows;
    Eigen::SparseQR<Sparse,Eigen::COLAMDOrdering<int>> qr;
    Eigen::MatrixXd local;
    Eigen::VectorXd scales,rhs,b,candidate,gradient,solution,transformed,reduced_rhs;
    std::vector<Eigen::Index> columns,previous_columns;
    std::vector<bool> free;
    std::vector<Eigen::Triplet<double>> entries;
    std::vector<double> response;
    std::vector<int> local_column,tile_columns;
};
}
