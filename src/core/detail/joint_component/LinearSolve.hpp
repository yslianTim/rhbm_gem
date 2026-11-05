#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <vector>

namespace rhbm_gem::core::joint_component {
struct ConstrainedAcWorkspace
{
    using Sparse=Eigen::SparseMatrix<double>;
    using RowMajorSparse=Eigen::SparseMatrix<double,Eigen::RowMajor>;

    Sparse normalized,free_design,reduced;
    RowMajorSparse rows;
    Eigen::SparseQR<Sparse,Eigen::COLAMDOrdering<int>> qr;
    Eigen::MatrixXd local;
    Eigen::VectorXd scales,rhs,b,candidate,gradient,solution,transformed,reduced_rhs;
    std::vector<Eigen::Index> free_column_indices,previous_free_column_indices;
    // Every C is free; an A is free exactly when its constraint is outside the working set.
    std::vector<bool> is_free_variable;
    std::vector<Eigen::Triplet<double>> entries;
    std::vector<double> response;
    std::vector<int> local_column,tile_columns;
};
}
