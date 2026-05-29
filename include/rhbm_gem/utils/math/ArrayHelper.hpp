#pragma once

#include <cstddef>
#include <cmath>
#include <limits>
#include <vector>
#include <tuple>
#include <array>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <type_traits>

namespace rhbm_gem {
namespace array_helper {

template <typename Type>
Type ComputeMin(const Type * data, size_t size)
{
    if (size <= 0) return std::numeric_limits<Type>::quiet_NaN();
    Type min_value{ std::numeric_limits<Type>::max() };
    for (size_t i = 0; i < size; i++)
    {
        min_value = std::min(min_value, data[i]);
    }
    return min_value;
}

template <typename Type>
Type ComputeMax(const Type * data, size_t size)
{
    if (size <= 0) return std::numeric_limits<Type>::quiet_NaN();
    Type max_value{ std::numeric_limits<Type>::lowest() };
    for (size_t i = 0; i < size; i++)
    {
        max_value = std::max(max_value, data[i]);
    }
    return max_value;
}

template <typename Type>
Type ComputePercentile(const std::vector<Type> & data, double percentile)
{
    if (data.empty() || percentile <= 0.0) return static_cast<Type>(0.0);

    if (percentile >= 1.0)
    {
        return data.empty() ? static_cast<Type>(0.0) : *std::max_element(data.begin(), data.end());
    }
    std::vector<Type> buffer(data);
    auto n{ buffer.size() };
    auto pos{ percentile * static_cast<double>(n - 1) };
    auto idx{ static_cast<size_t>(std::floor(pos)) };
    double frac{ pos - static_cast<double>(idx) };

    std::nth_element(
        buffer.begin(),
        buffer.begin() + static_cast<typename std::vector<Type>::difference_type>(idx),
        buffer.end()
    );
    Type lower{ buffer[idx] };
    if (idx + 1 < n)
    {
        std::nth_element(
            buffer.begin(),
            buffer.begin() + static_cast<typename std::vector<Type>::difference_type>(idx + 1),
            buffer.end()
        );
        Type upper{ buffer[idx + 1] };
        return static_cast<Type>(lower * (1 - frac) + upper * frac);
    }
    else
    {
        return lower;
    }
}

template <typename Type>
Type ComputeMean(const Type * data, size_t size)
{
    if (size <= 0) return static_cast<Type>(0.0);

    const auto sum{ std::accumulate(data, data + size, static_cast<Type>(0.0)) };
    return sum / static_cast<Type>(size);
}

template <typename Type>
Type ComputeStandardDeviation(const Type * data, size_t size, Type mean)
{
    if (size <= 1) return static_cast<Type>(0.0);

    const auto sum_sq_diff{
        std::accumulate(
            data,
            data + size,
            static_cast<Type>(0.0),
            [mean](Type sum, Type value)
            {
                const auto diff{ static_cast<Type>(value - mean) };
                return static_cast<Type>(sum + diff * diff);
            })
    };
    return std::sqrt(sum_sq_diff / static_cast<Type>(size - 1));
}

template <typename Type>
Type ComputeMedian(const std::vector<Type> & data)
{
    if (data.empty()) return static_cast<Type>(0.0);

    std::vector<Type> buffer(data);
    size_t n{ buffer.size() };
    size_t mid{ n / 2 };

    std::nth_element(
        buffer.begin(),
        buffer.begin() + static_cast<typename std::vector<Type>::difference_type>(mid),
        buffer.end()
    );
    Type upper_mid{ buffer[mid] };
    if (n % 2 == 1)
    {
        return upper_mid;
    }
    else
    {
        std::nth_element(
            buffer.begin(),
            buffer.begin() + static_cast<typename std::vector<Type>::difference_type>(mid - 1),
            buffer.begin() + static_cast<typename std::vector<Type>::difference_type>(mid)
        );
        Type lower_mid{ buffer[mid - 1] };
        return static_cast<Type>((lower_mid + upper_mid) / 2.0);
    }
}

inline std::size_t ComputeProportionValueCount(std::size_t data_size, double ratio)
{
    if (data_size == 0 || !std::isfinite(ratio) || ratio <= 0.0) return 0;

    const auto clamped_ratio{ std::min(ratio, 1.0) };
    return std::max<std::size_t>(
        1,
        static_cast<std::size_t>(std::ceil(static_cast<double>(data_size) * clamped_ratio)));
}

template <typename Type>
std::vector<Type> ComputeSmallestProportionValues(const std::vector<Type> & data, double ratio)
{
    const auto value_count{ ComputeProportionValueCount(data.size(), ratio) };
    if (value_count == 0) return {};

    std::vector<Type> buffer(data);
    std::sort(buffer.begin(), buffer.end());
    buffer.resize(value_count);
    return buffer;
}

template <typename Type>
std::vector<Type> ComputeLargestProportionValues(const std::vector<Type> & data, double ratio)
{
    const auto value_count{ ComputeProportionValueCount(data.size(), ratio) };
    if (value_count == 0) return {};

    std::vector<Type> buffer(data);
    std::sort(
        buffer.begin(),
        buffer.end(),
        [](const Type & lhs, const Type & rhs)
        {
            return lhs > rhs;
        });
    buffer.resize(value_count);
    return buffer;
}

template <typename Type>
std::tuple<Type, Type> ComputePercentileRangeTuple(
    const std::vector<Type> & data, double percentile_min=0.01, double percentile_max=0.99)
{
    if (data.empty()) return std::make_tuple(static_cast<Type>(0.0), static_cast<Type>(0.0));
    auto min_value{ ComputePercentile(data, percentile_min) };
    auto max_value{ ComputePercentile(data, percentile_max) };
    return std::make_tuple(min_value, max_value);
}

template <typename Type>
std::tuple<Type, Type> ComputeScalingPercentileRangeTuple(
    const std::vector<Type> & data, Type scaling,
    double percentile_min=0.01, double percentile_max=0.99)
{
    auto range_tuple{ ComputePercentileRangeTuple(data, percentile_min, percentile_max) };
    auto range{ std::get<1>(range_tuple) - std::get<0>(range_tuple) };
    auto min_value{ std::get<0>(range_tuple) - scaling * range };
    auto max_value{ std::get<1>(range_tuple) + scaling * range };
    return std::make_tuple(min_value, max_value);
}

template <typename Type>
std::tuple<Type, Type> ComputeRangeTuple(const std::vector<Type> & data)
{
    if (data.empty()) return std::make_tuple(static_cast<Type>(0.0), static_cast<Type>(0.0));
    auto min_value{ ComputeMin(data.data(), data.size()) };
    auto max_value{ ComputeMax(data.data(), data.size()) };
    return std::make_tuple(min_value, max_value);
}

template <typename Type>
std::tuple<Type, Type> ComputeScalingRangeTuple(
    const std::vector<Type> & data,
    Type scaling,
    Type min_range = static_cast<Type>(0.1))
{
    if (data.empty()) return std::make_tuple(static_cast<Type>(0.0), static_cast<Type>(0.0));
    auto range_tuple{ ComputeRangeTuple(data) };
    auto range{ std::get<1>(range_tuple) - std::get<0>(range_tuple) };
    auto min_value{ std::get<0>(range_tuple) - scaling * range };
    auto max_value{ std::get<1>(range_tuple) + scaling * range };
    auto range_mean{ (min_value + max_value) / static_cast<Type>(2.0) };
    auto half_range{ (max_value - min_value) / static_cast<Type>(2.0) };
    if (half_range < min_range)
    {
        min_value = range_mean - min_range;
        max_value = range_mean + min_range;
    }
    return std::make_tuple(min_value, max_value);
}

template <typename Type, std::size_t N>
Type ComputeDotProduct(const std::array<Type, N> & v1, const std::array<Type, N> & v2)
{
    return std::inner_product(v1.begin(), v1.end(), v2.begin(), static_cast<Type>(0.0));
}

template <typename Type, std::size_t N>
Type ComputeNorm(const std::array<Type, N> & vec)
{
    return static_cast<Type>(std::sqrt(ComputeDotProduct(vec, vec)));
}

template <typename Type, std::size_t N>
Type ComputeNorm(const std::array<Type, N> & v1, const std::array<Type, N> & v2)
{
    Type norm_square{ static_cast<Type>(0.0) };
    for (std::size_t i = 0; i < N; i++)
    {
        const auto diff{ static_cast<Type>(v1[i] - v2[i]) };
        norm_square += diff * diff;
    }
    return static_cast<Type>(std::sqrt(norm_square));
}

template <typename Type, std::size_t N>
std::array<Type, N> ComputeVector(
    const std::array<Type, N> & p1,
    const std::array<Type, N> & p2,
    bool normalize)
{
    static_assert(
        std::is_floating_point<Type>::value,
        "array_helper::ComputeVector requires a floating point Type.");

    std::array<Type, N> vector{};
    for (std::size_t i = 0; i < N; i++)
    {
        vector[i] = static_cast<Type>(p2[i] - p1[i]);
    }
    if (!normalize) return vector;

    const auto norm{ ComputeNorm(vector) };
    if (norm == static_cast<Type>(0.0)) return vector;

    for (auto & value : vector)
    {
        value = static_cast<Type>(value / norm);
    }
    return vector;
}

template <typename Type, std::size_t N>
int ComputeRank(const std::array<Type, N> & values, std::size_t index_to_rank)
{
    if (index_to_rank >= N)
    {
        throw std::out_of_range("GetRank index_to_rank out of range");
    }

    std::array<Type, N> sorted_values{ values };
    std::sort(sorted_values.begin(), sorted_values.end());

    auto it{ std::find(sorted_values.begin(), sorted_values.end(), values[index_to_rank]) };
    return static_cast<int>(std::distance(sorted_values.begin(), it)) + 1;
}

} // namespace array_helper
} // namespace rhbm_gem
