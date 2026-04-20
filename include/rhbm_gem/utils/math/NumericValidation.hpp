#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace rhbm_gem::NumericValidation {
namespace detail {

template <typename Type>
inline bool IsFinite(Type value)
{
    if constexpr (std::is_floating_point_v<Type>)
    {
        return std::isfinite(value);
    }
    return true;
}

template <typename Type>
inline std::string ToString(Type value)
{
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

template <typename NameType>
inline std::string GetName(const NameType & name)
{
    return std::string(name);
}

} // namespace detail

template <typename Type, typename NameType>
inline Type RequirePositive(Type value, const NameType & name)
{
    if (!detail::IsFinite(value) || value <= static_cast<Type>(0))
    {
        throw std::invalid_argument(detail::GetName(name) + " must be positive.");
    }
    return value;
}

template <typename Type, typename NameType>
inline Type RequireNonNegative(Type value, const NameType & name)
{
    if (!detail::IsFinite(value) || value < static_cast<Type>(0))
    {
        throw std::invalid_argument(detail::GetName(name) + " must be non-negative.");
    }
    return value;
}

template <typename Type, typename NameType>
inline Type RequireFinite(Type value, const NameType & name)
{
    if (!detail::IsFinite(value))
    {
        throw std::invalid_argument(detail::GetName(name) + " must be finite.");
    }
    return value;
}

template <typename Type, typename NameType>
inline Type RequireFinitePositive(Type value, const NameType & name)
{
    if (!detail::IsFinite(value) || value <= static_cast<Type>(0))
    {
        throw std::invalid_argument(detail::GetName(name) + " must be finite and positive.");
    }
    return value;
}

template <typename Type, typename NameType>
inline Type RequireFiniteNonNegative(Type value, const NameType & name)
{
    if (!detail::IsFinite(value) || value < static_cast<Type>(0))
    {
        throw std::invalid_argument(detail::GetName(name) + " must be finite and non-negative.");
    }
    return value;
}

template <typename Type, typename NameType>
inline std::pair<Type, Type> RequireFiniteNonNegativeRange(
    Type min_value,
    Type max_value,
    const NameType & name)
{
    if (!detail::IsFinite(min_value) || !detail::IsFinite(max_value) ||
        min_value < static_cast<Type>(0) || max_value < static_cast<Type>(0) ||
        max_value < min_value)
    {
        throw std::invalid_argument(
            detail::GetName(name) + " must be finite and satisfy 0 <= min <= max.");
    }
    return { min_value, max_value };
}

template <typename Type, typename LowerType, typename UpperType, typename NameType>
inline Type RequireFiniteInclusiveRange(
    Type value,
    LowerType lower,
    UpperType upper,
    const NameType & name)
{
    using CommonType = std::common_type_t<Type, LowerType, UpperType>;
    const auto lower_bound{ static_cast<CommonType>(lower) };
    const auto upper_bound{ static_cast<CommonType>(upper) };
    const auto converted_value{ static_cast<CommonType>(value) };

    if (!detail::IsFinite(value) || !detail::IsFinite(lower) || !detail::IsFinite(upper) ||
        converted_value < lower_bound || converted_value > upper_bound)
    {
        throw std::invalid_argument(
            detail::GetName(name) + " must be finite and within [" +
            detail::ToString(lower) + ", " + detail::ToString(upper) + "].");
    }
    return value;
}

} // namespace rhbm_gem::NumericValidation
