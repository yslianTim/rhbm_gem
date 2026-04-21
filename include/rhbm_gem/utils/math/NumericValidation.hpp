#pragma once

#include <array>
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

template <typename Type, typename Predicate>
inline void RequireArrayValues(
    const Type & values,
    Predicate predicate,
    const std::string & message)
{
    for (const auto & value : values)
    {
        if (!predicate(value))
        {
            throw std::invalid_argument(message);
        }
    }
}

} // namespace detail

template <typename Type>
inline bool IsPositive(Type value)
{
    return value > static_cast<Type>(0);
}

template <typename Type>
inline bool IsNonNegative(Type value)
{
    return value >= static_cast<Type>(0);
}

template <typename Type>
inline bool IsFinite(Type value)
{
    return detail::IsFinite(value);
}

template <typename Type>
inline bool IsFinitePositive(Type value)
{
    return IsFinite(value) && IsPositive(value);
}

template <typename Type>
inline bool IsFiniteNonNegative(Type value)
{
    return IsFinite(value) && IsNonNegative(value);
}

template <typename Type, typename LowerType>
inline bool IsAtLeast(Type value, LowerType lower)
{
    using CommonType = std::common_type_t<Type, LowerType>;
    return static_cast<CommonType>(value) >= static_cast<CommonType>(lower);
}

template <typename Type, typename UpperType>
inline bool IsAtMost(Type value, UpperType upper)
{
    using CommonType = std::common_type_t<Type, UpperType>;
    return static_cast<CommonType>(value) <= static_cast<CommonType>(upper);
}

template <typename Type, typename LowerType, typename UpperType>
inline bool IsFiniteInclusiveRange(
    Type value,
    LowerType lower,
    UpperType upper)
{
    using CommonType = std::common_type_t<Type, LowerType, UpperType>;
    const auto lower_bound{ static_cast<CommonType>(lower) };
    const auto upper_bound{ static_cast<CommonType>(upper) };
    const auto converted_value{ static_cast<CommonType>(value) };

    return IsFinite(value) && IsFinite(lower) && IsFinite(upper) &&
           converted_value >= lower_bound && converted_value <= upper_bound;
}

template <typename Type, typename LowerType, typename UpperType>
inline bool IsFiniteExclusiveInclusiveRange(
    Type value,
    LowerType lower,
    UpperType upper)
{
    using CommonType = std::common_type_t<Type, LowerType, UpperType>;
    const auto lower_bound{ static_cast<CommonType>(lower) };
    const auto upper_bound{ static_cast<CommonType>(upper) };
    const auto converted_value{ static_cast<CommonType>(value) };

    return IsFinite(value) && IsFinite(lower) && IsFinite(upper) &&
           converted_value > lower_bound && converted_value <= upper_bound;
}

template <typename Type, typename NameType>
inline Type RequirePositive(Type value, const NameType & name)
{
    if (!IsPositive(value))
    {
        throw std::invalid_argument(detail::GetName(name) + " must be positive.");
    }
    return value;
}

template <typename Type, typename NameType>
inline Type RequireNonNegative(Type value, const NameType & name)
{
    if (!IsNonNegative(value))
    {
        throw std::invalid_argument(detail::GetName(name) + " must be non-negative.");
    }
    return value;
}

template <typename Type, typename NameType>
inline Type RequireFinite(Type value, const NameType & name)
{
    if (!IsFinite(value))
    {
        throw std::invalid_argument(detail::GetName(name) + " must be finite.");
    }
    return value;
}

template <typename Type, typename NameType>
inline Type RequireFinitePositive(Type value, const NameType & name)
{
    if (!IsFinitePositive(value))
    {
        throw std::invalid_argument(detail::GetName(name) + " must be finite and positive.");
    }
    return value;
}

template <typename Type, typename NameType>
inline Type RequireFiniteNonNegative(Type value, const NameType & name)
{
    if (!IsFiniteNonNegative(value))
    {
        throw std::invalid_argument(detail::GetName(name) + " must be finite and non-negative.");
    }
    return value;
}

template <typename Type, typename LowerType, typename NameType>
inline Type RequireAtLeast(Type value, LowerType lower, const NameType & name)
{
    if (!IsAtLeast(value, lower))
    {
        throw std::invalid_argument(
            detail::GetName(name) + " must be at least " + detail::ToString(lower) + ".");
    }
    return value;
}

template <typename Type, typename UpperType, typename NameType>
inline Type RequireAtMost(Type value, UpperType upper, const NameType & name)
{
    if (!IsAtMost(value, upper))
    {
        throw std::invalid_argument(
            detail::GetName(name) + " must be at most " + detail::ToString(upper) + ".");
    }
    return value;
}

template <typename Type, typename NameType>
inline std::pair<Type, Type> RequireFiniteNonNegativeRange(
    Type min_value,
    Type max_value,
    const NameType & name)
{
    if (!IsFiniteNonNegative(min_value) || !IsFiniteNonNegative(max_value) ||
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
    if (!IsFiniteInclusiveRange(value, lower, upper))
    {
        throw std::invalid_argument(
            detail::GetName(name) + " must be finite and within [" +
            detail::ToString(lower) + ", " + detail::ToString(upper) + "].");
    }
    return value;
}

template <typename Type, typename LowerType, typename UpperType, typename NameType>
inline Type RequireFiniteExclusiveInclusiveRange(
    Type value,
    LowerType lower,
    UpperType upper,
    const NameType & name)
{
    if (!IsFiniteExclusiveInclusiveRange(value, lower, upper))
    {
        throw std::invalid_argument(
            detail::GetName(name) + " must be finite and within (" +
            detail::ToString(lower) + ", " + detail::ToString(upper) + "].");
    }
    return value;
}

template <typename Type, std::size_t Size, typename NameType>
inline const std::array<Type, Size> & RequireAllPositive(
    const std::array<Type, Size> & values,
    const NameType & name)
{
    detail::RequireArrayValues(
        values,
        [](const auto value) { return IsPositive(value); },
        detail::GetName(name) + " must contain only positive values.");
    return values;
}

template <typename Type, std::size_t Size, typename NameType>
inline const std::array<Type, Size> & RequireAllFinitePositive(
    const std::array<Type, Size> & values,
    const NameType & name)
{
    detail::RequireArrayValues(
        values,
        [](const auto value) { return IsFinitePositive(value); },
        detail::GetName(name) + " must contain only finite positive values.");
    return values;
}

} // namespace rhbm_gem::NumericValidation
