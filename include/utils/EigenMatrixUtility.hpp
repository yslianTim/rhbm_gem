#pragma once

#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif

#include <algorithm>
#include <iostream>
#include <Eigen/Dense>

namespace EigenMatrixUtility
{

template<typename Derived>
typename Derived::Scalar GetMedian(Eigen::DenseBase<Derived>& d)
{
    auto array{ d.reshaped() };
    auto array_size{ array.size() };
    if (array_size == 0)
    {
        std::cout <<"[Warning] 'EigenMatrixUtility::GetMedian' The input data size is zero, return 0..."<< std::endl;
        return 0.0;
    }
    std::sort(array.begin(), array.end());
    return (array_size % 2 == 0) ? array.segment((array_size - 2)/2, 2).mean() : array(array_size/2);
}

template<typename Derived>
typename Derived::Scalar GetMedian(const Eigen::DenseBase<Derived>& d)
{
    typename Derived::PlainObject m{ d.replicate(1,1) };
    return GetMedian(m);
}

template<typename Derived>
typename Derived::Scalar GetStandardDeviation(Eigen::DenseBase<Derived>& d)
{
    auto array{ d.reshaped() };
    auto array_size{ array.size() };
    if (array_size == 0)
    {
        std::cout <<"[Warning] 'EigenMatrixUtility::GetStandardDeviation' The input data size is zero, return 0..."<< std::endl;
        return 0.0;
    }
    if (array_size <= 1) return 0.0f;
    return (array - array.mean()).matrix().norm()/std::sqrt(array_size - 1);
}

template<typename Derived>
typename Derived::Scalar GetStandardDeviation(const Eigen::DenseBase<Derived>& d)
{
    typename Derived::PlainObject m{ d.replicate(1,1) };
    return GetStandardDeviation(m);
}

template<typename Derived>
typename Derived::PlainObject GetInverseMatrix(const Eigen::MatrixBase<Derived>& matrix)
{
    Eigen::FullPivLU<typename Derived::PlainObject> lu_decomposition(matrix);
    if (!lu_decomposition.isInvertible())
    {
        return matrix.completeOrthogonalDecomposition().pseudoInverse();
    }
    return lu_decomposition.inverse();
}

template<typename Scalar, int SizeAtCompileTime, int MaxSizeAtCompileTime = SizeAtCompileTime>
Eigen::DiagonalMatrix<Scalar, SizeAtCompileTime, MaxSizeAtCompileTime>
GetInverseDiagonalMatrix(const Eigen::DiagonalMatrix<Scalar, SizeAtCompileTime, MaxSizeAtCompileTime>& matrix)
{
    Eigen::DiagonalMatrix<Scalar, SizeAtCompileTime, MaxSizeAtCompileTime> inverseDiagMatrix;
    inverseDiagMatrix.diagonal() = matrix.diagonal().unaryExpr([](const Scalar& val) {
        if (val == 0.0) {
            return 0.0;
        }
        return 1.0 / val;
    });
    return inverseDiagMatrix;
}

}
