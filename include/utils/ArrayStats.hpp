#pragma once

#include <cmath>
#include <limits>
#include <vector>
#include <tuple>
#include <array>
#include <algorithm>

#ifdef USE_OPENMP
#include <omp.h>
#endif

template <typename Type>
class ArrayStats
{
public:
    static Type ComputeMin(
        const Type * data, size_t size, [[maybe_unused]] int thread_size = 1)
    {
        if (size <= 0) return std::numeric_limits<Type>::quiet_NaN();
        Type min_value{ std::numeric_limits<Type>::max() };
        #ifdef USE_OPENMP
        #pragma omp parallel for reduction(min:min_value) num_threads(thread_size)
        #endif
        for (size_t i = 0; i < size; i++)
        {
            min_value = std::min(min_value, data[i]);
        }
        return min_value;
    }

    static Type ComputeMax(
        const Type * data, size_t size, [[maybe_unused]] int thread_size = 1)
    {
        if (size <= 0) return std::numeric_limits<Type>::quiet_NaN();
        Type max_value{ std::numeric_limits<Type>::lowest() };
        #ifdef USE_OPENMP
        #pragma omp parallel for reduction(max:max_value) num_threads(thread_size)
        #endif
        for (size_t i = 0; i < size; i++)
        {
            max_value = std::max(max_value, data[i]);
        }
        return max_value;
    }

    static Type ComputePercentile(const std::vector<Type> & data, double percentile)
    {
        if (data.empty() || percentile <= 0.0)
        {
            return static_cast<Type>(0.0);
        }
        if (percentile >= 1.0)
        {
            return data.empty() ? static_cast<Type>(0.0) :
                   *std::max_element(data.begin(), data.end());
        }
        std::vector<Type> sorted_data(data);
        std::sort(sorted_data.begin(), sorted_data.end());
        auto n{ sorted_data.size() };
        auto pos{ percentile * static_cast<double>(n - 1) };
        auto idx{ static_cast<size_t>(std::floor(pos)) };
        double frac{ pos - static_cast<double>(idx) };
        if (idx + 1 < n)
        {
            return static_cast<Type>(sorted_data[idx] * (1 - frac) + sorted_data[idx + 1] * frac);
        }
        else
        {
            return sorted_data[idx];
        }
    }

    static Type ComputeMean(
        const Type * data, size_t size, [[maybe_unused]] int thread_size = 1)
    {
        if (size <= 0)
        {
            return static_cast<Type>(0.0);
        }
        auto sum{ static_cast<Type>(0.0) };
        #ifdef USE_OPENMP
        #pragma omp parallel for reduction(+:sum) num_threads(thread_size)
        #endif
        for (size_t i = 0; i < size; i++)
        {
            sum += data[i];
        }
        return sum / static_cast<Type>(size);
    }

    static Type ComputeStandardDeviation(
        const Type * data, size_t size, Type mean, [[maybe_unused]] int thread_size = 1)
    {
        if (size <= 1)
        {
            return static_cast<Type>(0.0);
        }
        auto sum_sq_diff{ static_cast<Type>(0.0) };
        #ifdef USE_OPENMP
        #pragma omp parallel for reduction(+:sum_sq_diff) num_threads(thread_size)
        #endif
        for (size_t i = 0; i < size; i++)
        {
            auto diff{ static_cast<Type>(data[i]) - static_cast<Type>(mean) };
            sum_sq_diff += diff * diff;
        }
        return std::sqrt(sum_sq_diff / static_cast<Type>(size - 1));
    }

    static Type ComputeMedian(const std::vector<Type> & data)
    {
        if (data.empty())
        {
            return static_cast<Type>(0.0);
        }
        
        std::vector<Type> sorted_data(data);
        std::sort(sorted_data.begin(), sorted_data.end());
        
        size_t n{ sorted_data.size() };
        if (n % 2 == 1)
        {
            return sorted_data[n / 2];
        }
        else
        {
            return static_cast<Type>((sorted_data[n / 2 - 1] + sorted_data[n / 2]) / 2.0);
        }
    }

    static std::tuple<Type, Type> ComputePercentileRangeTuple(
        const std::vector<Type> & data, double percentile_min=0.01, double percentile_max=0.99)
    {
        if (data.empty())
        {
            return std::make_tuple(static_cast<Type>(0.0), static_cast<Type>(0.0));
        }
        auto min_value{ ComputePercentile(data, percentile_min) };
        auto max_value{ ComputePercentile(data, percentile_max) };
        return std::make_tuple(min_value, max_value);
    }

    static std::tuple<Type, Type> ComputeScalingPercentileRangeTuple(
        const std::vector<Type> & data, Type scaling,
        double percentile_min=0.01, double percentile_max=0.99)
    {
        auto range_tuple{ ComputePercentileRangeTuple(data, percentile_min, percentile_max) };
        auto range{ std::get<1>(range_tuple) - std::get<0>(range_tuple) };
        auto min_value{ std::get<0>(range_tuple) - scaling * range };
        auto max_value{ std::get<1>(range_tuple) + scaling * range };
        return std::make_tuple(min_value, max_value);
    }

    static std::tuple<Type, Type> ComputeRangeTuple(
        const std::vector<Type> & data, int thread_size = 1)
    {
        if (data.empty())
        {
            return std::make_tuple(static_cast<Type>(0.0), static_cast<Type>(0.0));
        }
        auto min_value{ ComputeMin(data.data(), data.size(), thread_size) };
        auto max_value{ ComputeMax(data.data(), data.size(), thread_size) };
        return std::make_tuple(min_value, max_value);
    }

    static std::tuple<Type, Type> ComputeScalingRangeTuple(
        const std::vector<Type> & data, Type scaling, int thread_size = 1)
    {
        auto range_tuple{ ComputeRangeTuple(data, thread_size) };
        auto range{ std::get<1>(range_tuple) - std::get<0>(range_tuple) };
        auto min_value{ std::get<0>(range_tuple) - scaling * range };
        auto max_value{ std::get<1>(range_tuple) + scaling * range };
        return std::make_tuple(min_value, max_value);
    }

    static Type ComputeNorm(const std::array<Type, 3> & v1, const std::array<Type, 3> & v2)
    {
        auto diff_x{ static_cast<Type>(v1[0] - v2[0]) };
        auto diff_y{ static_cast<Type>(v1[1] - v2[1]) };
        auto diff_z{ static_cast<Type>(v1[2] - v2[2]) };
        return static_cast<Type>(std::sqrt(std::pow(diff_x, 2) + std::pow(diff_y, 2) + std::pow(diff_z, 2)));
    }
};