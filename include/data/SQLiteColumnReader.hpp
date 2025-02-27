#pragma once

#include <sqlite3.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstring>

template<typename T>
struct SQLiteColumnReader
{
};

// int specialization
template<>
struct SQLiteColumnReader<int>
{
    static int Get(sqlite3_stmt * stmt, int index)
    {
        return sqlite3_column_int(stmt, index);
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
        auto text{ sqlite3_column_text(stmt, index) };
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
        auto blob_data{ sqlite3_column_blob(stmt, index) };
        auto blob_size{ sqlite3_column_bytes(stmt, index) };
        if (!blob_data || blob_size <= 0)
        {
            return {};
        }
        auto count{ blob_size / static_cast<int>(sizeof(float)) };
        std::vector<float> result(count);
        std::memcpy(result.data(), blob_data, blob_size);
        return result;
    }
};

// std::vector<double> specialization
template<>
struct SQLiteColumnReader<std::vector<double>>
{
    static std::vector<double> Get(sqlite3_stmt * stmt, int index)
    {
        auto blob_data{ sqlite3_column_blob(stmt, index) };
        auto blob_size{ sqlite3_column_bytes(stmt, index) };
        if (!blob_data || blob_size <= 0)
        {
            return {};
        }
        auto count{ blob_size / static_cast<int>(sizeof(double)) };
        std::vector<double> result(count);
        std::memcpy(result.data(), blob_data, blob_size);
        return result;
    }
};