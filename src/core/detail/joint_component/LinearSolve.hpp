#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <cstdint>
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
#ifdef RHBM_GEM_TEST_INSTRUMENTATION
    bool has_previous_pattern{},has_previous_reduced_pattern{};
    Eigen::Index previous_pattern_rows{},previous_pattern_columns{};
    std::size_t previous_pattern_nonzeros{};
    std::uint64_t previous_face_hash1{},previous_face_hash2{},previous_pattern_hash1{},previous_pattern_hash2{};
    Eigen::Index previous_reduced_rows{},previous_reduced_columns{};
    std::size_t previous_reduced_nonzeros{};
    std::uint64_t previous_reduced_pattern_hash1{},previous_reduced_pattern_hash2{},
        previous_reduced_numeric_hash1{},previous_reduced_numeric_hash2{};
#endif
};
}
