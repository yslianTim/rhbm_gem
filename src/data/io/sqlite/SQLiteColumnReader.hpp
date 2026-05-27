#pragma once

#include <cstddef>
#include <sqlite3.h>
#include <string>
#include <cstring>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {

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
struct SQLiteColumnReader<uint8_t>
{
    static uint8_t Get(sqlite3_stmt * stmt, int index)
    {
        return static_cast<uint8_t>(sqlite3_column_int(stmt, index));
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

template<>
struct SQLiteColumnReader<uint32_t>
{
    static uint32_t Get(sqlite3_stmt * stmt, int index)
    {
        return static_cast<uint32_t>(sqlite3_column_int(stmt, index));
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
        if (sqlite3_column_type(stmt, index) == SQLITE_NULL)
        {
            return {};
        }

        const unsigned char * text{ sqlite3_column_text(stmt, index) };
        // According to sqlite3 documentation sqlite3_column_text returns a pointer to a buffer
        // that is valid until the next call to sqlite3_step() for the same statement. It may
        // return nullptr only if the column type is SQLITE_NULL, which we check above.
        // However, perform an extra safety check here just in case.
        if (!text)
        {
            return {};
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

// LocalPotentialSampleList specialization
template <>
struct SQLiteColumnReader<LocalPotentialSampleList>
{
    static LocalPotentialSampleList Get(sqlite3_stmt* stmt, int index)
    {
        const void * blob_data = sqlite3_column_blob(stmt, index);
        int blob_size{ sqlite3_column_bytes(stmt, index) };
        if (!blob_data || blob_size <= 0)
        {
            return {};
        }
        if (blob_size % (3 * static_cast<int>(sizeof(float))) != 0)
        {
            throw std::runtime_error("Local potential sample blob requires sampling size for legacy decoding.");
        }
        int count{ blob_size / (3 * static_cast<int>(sizeof(float))) };
        LocalPotentialSampleList result;
        result.reserve(static_cast<size_t>(count));
        const float * blob_floats{ reinterpret_cast<const float *>(blob_data) };
        for (int i = 0; i < count; ++i)
        {
            const int offset{ 3 * i };
            result.emplace_back(LocalPotentialSample{
                blob_floats[offset + 1],
                SamplingPoint{
                    blob_floats[offset],
                    { 0.0f, 0.0f, 0.0f },
                    blob_floats[offset + 2] != 0.0f
                }
            });
        }
        return result;
    }
};

inline LocalPotentialSampleList ReadLocalPotentialSampleListColumn(
    sqlite3_stmt * stmt,
    int index,
    int sampling_size)
{
    const void * blob_data = sqlite3_column_blob(stmt, index);
    const int blob_size{ sqlite3_column_bytes(stmt, index) };
    if (!blob_data || blob_size <= 0)
    {
        return {};
    }
    if (blob_size % static_cast<int>(sizeof(float)) != 0)
    {
        throw std::runtime_error("Invalid local potential sample blob size.");
    }

    const int float_count{ blob_size / static_cast<int>(sizeof(float)) };
    const int legacy_float_count{ sampling_size * 2 };
    const int current_float_count{ sampling_size * 3 };
    int stride{ 0 };
    if (float_count == current_float_count)
    {
        stride = 3;
    }
    else if (float_count == legacy_float_count)
    {
        stride = 2;
    }
    else
    {
        throw std::runtime_error("Local potential sample blob size is inconsistent with sampling size.");
    }

    LocalPotentialSampleList result;
    result.reserve(static_cast<size_t>(sampling_size));
    const float * blob_floats{ reinterpret_cast<const float *>(blob_data) };
    for (int i = 0; i < sampling_size; ++i)
    {
        const int offset{ stride * i };
        const auto is_selected{
            stride == 2 ? true : blob_floats[offset + 2] != 0.0f
        };
        result.emplace_back(LocalPotentialSample{
            blob_floats[offset + 1],
            SamplingPoint{
                blob_floats[offset],
                { 0.0f, 0.0f, 0.0f },
                is_selected
            }
        });
    }
    return result;
}

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

} // namespace rhbm_gem
