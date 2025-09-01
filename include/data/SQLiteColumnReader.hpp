#pragma once

#include <cstddef>
#include <sqlite3.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <tuple>
#include <cstring>

template<typename T>
struct SQLiteColumnReader
{
};

// int 32-bits specialization
template<>
struct SQLiteColumnReader<int>
{
    static int Get(sqlite3_stmt * stmt, int index)
    {
        return sqlite3_column_int(stmt, index);
    }
};

template<>
struct SQLiteColumnReader<uint16_t>
{
    static uint16_t Get(sqlite3_stmt * stmt, int index)
    {
        return static_cast<uint16_t>(sqlite3_column_int(stmt, index));
    }
};

// int 64-bits specialization
template<>
struct SQLiteColumnReader<int64_t>
{
    static int64_t Get(sqlite3_stmt * stmt, int index)
    {
        return sqlite3_column_int64(stmt, index);
    }
};

template<>
struct SQLiteColumnReader<uint64_t>
{
    static uint64_t Get(sqlite3_stmt * stmt, int index)
    {
        return static_cast<uint64_t>(sqlite3_column_int64(stmt, index));
    }
};

// double specialization
template<>
struct SQLiteColumnReader<double>
{
    static double Get(sqlite3_stmt * stmt, int index)
    {
        return sqlite3_column_double(stmt, index);
    }
};

// std::string specialization
template<>
struct SQLiteColumnReader<std::string>
{
    static std::string Get(sqlite3_stmt * stmt, int index)
    {
        const unsigned char * text{ sqlite3_column_text(stmt, index) };
        if (!text)
        {
            throw std::runtime_error("NULL column string");
        }
        return std::string(reinterpret_cast<const char*>(text));
    }
};

// std::vector<float> specialization
template<>
struct SQLiteColumnReader<std::vector<float>>
{
    static std::vector<float> Get(sqlite3_stmt * stmt, int index)
    {
        const void * blob_data = sqlite3_column_blob(stmt, index);
        int blob_size = sqlite3_column_bytes(stmt, index);
        if (!blob_data || blob_size <= 0)
        {
            return {};
        }
        auto count{ blob_size / static_cast<int>(sizeof(float)) };
        std::vector<float> result;
        result.resize(static_cast<size_t>(count));
        std::memcpy(result.data(), blob_data, static_cast<size_t>(blob_size));
        return result;
    }
};

// std::vector<double> specialization
template<>
struct SQLiteColumnReader<std::vector<double>>
{
    static std::vector<double> Get(sqlite3_stmt * stmt, int index)
    {
        const void * blob_data = sqlite3_column_blob(stmt, index);
        int blob_size{ sqlite3_column_bytes(stmt, index) };
        if (!blob_data || blob_size <= 0)
        {
            return {};
        }
        int count{ blob_size / static_cast<int>(sizeof(double)) };
        std::vector<double> result;
        result.resize(static_cast<size_t>(count));
        std::memcpy(result.data(), blob_data, static_cast<size_t>(blob_size));
        return result;
    }
};

// std::vector<std::tuple<float, float>> specialization
template <>
struct SQLiteColumnReader<std::vector<std::tuple<float, float>>>
{
    static std::vector<std::tuple<float, float>> Get(sqlite3_stmt* stmt, int index)
    {
        const void * blob_data = sqlite3_column_blob(stmt, index);
        int blob_size{ sqlite3_column_bytes(stmt, index) };
        if (!blob_data || blob_size <= 0)
        {
            return {};
        }
        int count{ blob_size / (2 * static_cast<int>(sizeof(float))) };
        std::vector<std::tuple<float, float>> result;
        result.reserve(static_cast<size_t>(count));
        const float * blob_floats{ reinterpret_cast<const float *>(blob_data) };
        for (int i = 0; i < count; ++i)
        {
            float first = blob_floats[2 * i];
            float second = blob_floats[2 * i + 1];
            result.emplace_back(std::make_tuple(first, second));
        }
        return result;
    }
};

// std::vector<std::tuple<double, double>> specialization
template <>
struct SQLiteColumnReader<std::vector<std::tuple<double, double>>>
{
    static std::vector<std::tuple<double, double>> Get(sqlite3_stmt* stmt, int index)
    {
        const void * blob_data = sqlite3_column_blob(stmt, index);
        int blob_size{ sqlite3_column_bytes(stmt, index) };
        if (!blob_data || blob_size <= 0)
        {
            return {};
        }
        int count{ blob_size / (2 * static_cast<int>(sizeof(double))) };
        std::vector<std::tuple<double, double>> result;
        result.reserve(static_cast<size_t>(count));
        const double * blob_doubles{ reinterpret_cast<const double *>(blob_data) };
        for (int i = 0; i < count; ++i)
        {
            double first = blob_doubles[2 * i];
            double second = blob_doubles[2 * i + 1];
            result.emplace_back(std::make_tuple(first, second));
        }
        return result;
    }
};
