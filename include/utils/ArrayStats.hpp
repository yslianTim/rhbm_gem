#pragma once

#include <cmath>
#include <limits>
#include <vector>
#include <tuple>
#include <algorithm>

#ifdef USE_OPENMP
#include <omp.h>
#endif

template <typename Type>
class ArrayStats
{
public:
    static Type ComputeMin(const Type * data, int size, int thread_size = 1)
    {
        if (size <= 0) return std::numeric_limits<Type>::quiet_NaN();
        Type min_value{ std::numeric_limits<Type>::max() };
        #ifdef USE_OPENMP
        #pragma omp parallel for reduction(min:min_value) num_threads(thread_size)
        #endif
        for (int i = 0; i < size; i++)
        {
            min_value = std::min(min_value, data[i]);
        }
        return min_value;
    }

    static Type ComputeMax(const Type * data, int size, int thread_size = 1)
    {
        if (size <= 0) return std::numeric_limits<Type>::quiet_NaN();
        Type max_value{ std::numeric_limits<Type>::lowest() };
        #ifdef USE_OPENMP
        #pragma omp parallel for reduction(max:max_value) num_threads(thread_size)
        #endif
        for (int i = 0; i < size; i++)
        {
            max_value = std::max(max_value, data[i]);
        }
        return max_value;
    }

    static Type ComputeMean(const Type * data, int size, int thread_size = 1)
    {
        if (size <= 0) return Type(0);
        double sum{ 0.0 };
        #ifdef USE_OPENMP
        #pragma omp parallel for reduction(+:sum) num_threads(thread_size)
        #endif
        for (int i = 0; i < size; i++)
        {
            sum += data[i];
        }
        return static_cast<Type>(sum / size);
    }

    static Type ComputeStandardDeviation(const Type * data, int size, Type mean, int thread_size = 1)
    {
        if (size <= 1) return Type(0);
        double sum_sq_diff{ 0.0 };
        #ifdef USE_OPENMP
        #pragma omp parallel for reduction(+:sum_sq_diff) num_threads(thread_size)
        #endif
        for (int i = 0; i < size; i++)
        {
            double diff{ static_cast<double>(data[i]) - static_cast<double>(mean) };
            sum_sq_diff += diff * diff;
        }
        return static_cast<Type>(std::sqrt(sum_sq_diff / (size - 1)));
    }

    static Type ComputeMedian(const std::vector<Type>& data)
    {
        if (data.empty())
        {
            return Type(0);
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

    static std::tuple<Type, Type> ComputeRangeTuple(const std::vector<Type> & data, int thread_size = 1)
    {
        if (data.empty()) return std::make_tuple(0.0, 0.0);
        auto min_value{ ComputeMin(data.data(), static_cast<int>(data.size()), thread_size) };
        auto max_value{ ComputeMax(data.data(), static_cast<int>(data.size()), thread_size) };
        return std::make_tuple(min_value, max_value);
    }

    static std::tuple<Type, Type> ComputeScalingRangeTuple(
        const std::vector<Type> & data, double scaling, int thread_size = 1)
    {
        auto range_tuple{ ComputeRangeTuple(data, thread_size) };
        auto range{ std::get<1>(range_tuple) - std::get<0>(range_tuple) };
        auto min_value{ std::get<0>(range_tuple) - scaling * range };
        auto max_value{ std::get<1>(range_tuple) + scaling * range };
        return std::make_tuple(min_value, max_value);
    }
};