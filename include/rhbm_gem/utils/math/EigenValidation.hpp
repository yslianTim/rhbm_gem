#pragma once

#include <stdexcept>
#include <string>

#include <Eigen/Dense>

namespace rhbm_gem::EigenValidation {
namespace detail {

template <typename NameType>
inline std::string GetName(const NameType & name)
{
    return std::string(name);
}

} // namespace detail

template <typename Derived, typename NameType>
inline void RequireFinite(
    const Eigen::DenseBase<Derived> & dense,
    const NameType & name)
{
    if (!dense.allFinite())
    {
        throw std::invalid_argument(
            detail::GetName(name) + " must contain only finite values.");
    }
}

template <typename Derived, typename NameType>
inline void RequireNonEmpty(
    const Eigen::DenseBase<Derived> & dense,
    const NameType & name)
{
    if (dense.size() <= 0)
    {
        throw std::invalid_argument(detail::GetName(name) + " must not be empty.");
    }
}

template <typename Derived, typename NameType>
inline void RequireVectorSize(
    const Eigen::DenseBase<Derived> & vector_like,
    Eigen::Index expected_size,
    const NameType & name)
{
    if (vector_like.size() != expected_size)
    {
        throw std::invalid_argument(
            detail::GetName(name) + " must have size " +
            std::to_string(static_cast<long long>(expected_size)) + ".");
    }
}

template <typename LeftDerived, typename RightDerived, typename NameType>
inline void RequireSameSize(
    const Eigen::DenseBase<LeftDerived> & lhs,
    const Eigen::DenseBase<RightDerived> & rhs,
    const NameType & name)
{
    if (lhs.size() != rhs.size())
    {
        throw std::invalid_argument(detail::GetName(name) + " must have matching sizes.");
    }
}

template <typename Derived, typename NameType>
inline void RequireShape(
    const Eigen::DenseBase<Derived> & matrix,
    Eigen::Index expected_rows,
    Eigen::Index expected_cols,
    const NameType & name)
{
    if (matrix.rows() != expected_rows || matrix.cols() != expected_cols)
    {
        throw std::invalid_argument(
            detail::GetName(name) + " must have shape " +
            std::to_string(static_cast<long long>(expected_rows)) + "x" +
            std::to_string(static_cast<long long>(expected_cols)) + ".");
    }
}

template <typename Derived, typename NameType>
inline void RequireSquare(
    const Eigen::DenseBase<Derived> & matrix,
    const NameType & name)
{
    if (matrix.rows() != matrix.cols())
    {
        throw std::invalid_argument(detail::GetName(name) + " must be square.");
    }
}

template <typename Derived, typename NameType>
inline Eigen::Vector3d RequireVector3d(
    const Eigen::DenseBase<Derived> & vector_like,
    const NameType & name)
{
    RequireVectorSize(vector_like, 3, name);
    RequireFinite(vector_like, name);
    return Eigen::Vector3d{
        static_cast<double>(vector_like(0)),
        static_cast<double>(vector_like(1)),
        static_cast<double>(vector_like(2))
    };
}

} // namespace rhbm_gem::EigenValidation
