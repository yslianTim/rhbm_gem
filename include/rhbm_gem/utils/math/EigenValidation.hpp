#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

#include <Eigen/Dense>

namespace rhbm_gem::eigen_validation {
namespace detail {

template <typename NameType>
inline std::string GetName(const NameType & name)
{
    return std::string(name);
}

inline std::string ApplyContext(
    std::string message,
    std::string_view context)
{
    if (context.empty())
    {
        return message;
    }
    return std::string(context) + " Details: " + message;
}

} // namespace detail

template <typename Derived, typename NameType>
inline void RequireFinite(
    const Eigen::DenseBase<Derived> & dense,
    const NameType & name,
    std::string_view context = {})
{
    if (!dense.allFinite())
    {
        throw std::invalid_argument(detail::ApplyContext(
            detail::GetName(name) + " must contain only finite values.",
            context));
    }
}

template <typename Derived, typename NameType>
inline void RequireNonEmpty(
    const Eigen::DenseBase<Derived> & dense,
    const NameType & name,
    std::string_view context = {})
{
    if (dense.size() <= 0)
    {
        throw std::invalid_argument(detail::ApplyContext(
            detail::GetName(name) + " must not be empty.",
            context));
    }
}

template <typename Derived, typename NameType>
inline void RequireVectorSize(
    const Eigen::DenseBase<Derived> & vector_like,
    Eigen::Index expected_size,
    const NameType & name,
    std::string_view context = {})
{
    if (vector_like.size() != expected_size)
    {
        throw std::invalid_argument(detail::ApplyContext(
            detail::GetName(name) + " must have size " +
            std::to_string(static_cast<long long>(expected_size)) + ".",
            context));
    }
}

template <typename Derived, typename NameType>
inline void RequireRows(
    const Eigen::DenseBase<Derived> & matrix,
    Eigen::Index expected_rows,
    const NameType & name,
    std::string_view context = {})
{
    if (matrix.rows() != expected_rows)
    {
        throw std::invalid_argument(detail::ApplyContext(
            detail::GetName(name) + " must have " +
            std::to_string(static_cast<long long>(expected_rows)) + " rows.",
            context));
    }
}

template <typename Derived, typename NameType>
inline void RequireCols(
    const Eigen::DenseBase<Derived> & matrix,
    Eigen::Index expected_cols,
    const NameType & name,
    std::string_view context = {})
{
    if (matrix.cols() != expected_cols)
    {
        throw std::invalid_argument(detail::ApplyContext(
            detail::GetName(name) + " must have " +
            std::to_string(static_cast<long long>(expected_cols)) + " columns.",
            context));
    }
}

template <typename LeftDerived, typename RightDerived, typename NameType>
inline void RequireSameSize(
    const Eigen::DenseBase<LeftDerived> & lhs,
    const Eigen::DenseBase<RightDerived> & rhs,
    const NameType & name,
    std::string_view context = {})
{
    if (lhs.size() != rhs.size())
    {
        throw std::invalid_argument(detail::ApplyContext(
            detail::GetName(name) + " must have matching sizes.",
            context));
    }
}

template <typename Derived, typename NameType>
inline void RequireShape(
    const Eigen::DenseBase<Derived> & matrix,
    Eigen::Index expected_rows,
    Eigen::Index expected_cols,
    const NameType & name,
    std::string_view context = {})
{
    if (matrix.rows() != expected_rows || matrix.cols() != expected_cols)
    {
        throw std::invalid_argument(detail::ApplyContext(
            detail::GetName(name) + " must have shape " +
            std::to_string(static_cast<long long>(expected_rows)) + "x" +
            std::to_string(static_cast<long long>(expected_cols)) + ".",
            context));
    }
}

template <typename Derived, typename NameType>
inline void RequireSquare(
    const Eigen::DenseBase<Derived> & matrix,
    const NameType & name,
    std::string_view context = {})
{
    if (matrix.rows() != matrix.cols())
    {
        throw std::invalid_argument(detail::ApplyContext(
            detail::GetName(name) + " must be square.",
            context));
    }
}

template <typename Derived, typename NameType>
inline Eigen::Vector3d RequireVector3d(
    const Eigen::DenseBase<Derived> & vector_like,
    const NameType & name,
    std::string_view context = {})
{
    RequireVectorSize(vector_like, 3, name, context);
    RequireFinite(vector_like, name, context);
    return Eigen::Vector3d{
        static_cast<double>(vector_like(0)),
        static_cast<double>(vector_like(1)),
        static_cast<double>(vector_like(2))
    };
}

} // namespace rhbm_gem::eigen_validation
