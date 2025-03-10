#pragma once

#include <cmath>
#include <limits>
#include <algorithm>
#include <omp.h>

template <typename Type>
class ArrayStats
{
public:
    static Type ComputeMin(const Type * data, int size, int thread_size = 1)
    {
        if (size <= 0) return std::numeric_limits<Type>::quiet_NaN();
        Type min_val{ std::numeric_limits<Type>::max() };
        #pragma omp parallel for reduction(min:min_val) num_threads(thread_size)
        for (int i = 0; i < size; i++)
        {
            min_val = std::min(min_val, data[i]);
        }
        return min_val;
    }

    static Type ComputeMax(const Type * data, int size, int thread_size = 1)
    {
        if (size <= 0) return std::numeric_limits<Type>::quiet_NaN();
        Type max_val{ std::numeric_limits<Type>::lowest() };
        #pragma omp parallel for reduction(max:max_val) num_threads(thread_size)
        for (int i = 0; i < size; i++)
        {
            max_val = std::max(max_val, data[i]);
        }
        return max_val;
    }

    static Type ComputeMean(const Type * data, int size, int thread_size = 1)
    {
        if (size <= 0) return Type(0);
        double sum{ 0.0 };
        #pragma omp parallel for reduction(+:sum) num_threads(thread_size)
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
        #pragma omp parallel for reduction(+:sum_sq_diff) num_threads(thread_size)
        for (int i = 0; i < size; i++)
        {
            double diff{ static_cast<double>(data[i]) - static_cast<double>(mean) };
            sum_sq_diff += diff * diff;
        }
        return static_cast<Type>(std::sqrt(sum_sq_diff / (size - 1)));
    }
};