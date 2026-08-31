#pragma once

#include <array>
#include <cstddef>
#include <cstring>
#include <limits>
#include <sqlite3.h>
#include <string>
#include <stdexcept>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {

static_assert(sizeof(double) == 8, "SQLite sampling BLOB requires 64-bit double.");
static_assert(
    std::numeric_limits<double>::is_iec559,
    "SQLite sampling BLOB requires IEEE-754 double.");

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
        if (blob_size % (3 * static_cast<int>(sizeof(double))) != 0)
        {
            throw std::runtime_error("Invalid local potential sample blob size.");
        }
        int count{ blob_size / (3 * static_cast<int>(sizeof(double))) };
        LocalPotentialSampleList result;
        result.reserve(static_cast<size_t>(count));
        const auto * blob_bytes{ static_cast<const std::byte *>(blob_data) };
        for (int i = 0; i < count; ++i)
        {
            std::array<double, 3> values{};
            std::memcpy(
                values.data(),
                blob_bytes + static_cast<size_t>(i) * 3 * sizeof(double),
                3 * sizeof(double));
            result.emplace_back(LocalPotentialSample{
                values[1],
                SamplingPoint{
                    values[0],
                    { 0.0, 0.0, 0.0 },
                    values[2] != 0.0
                }
            });
        }
        return result;
    }
};

} // namespace rhbm_gem
