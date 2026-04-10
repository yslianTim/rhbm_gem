#pragma once

#include <sqlite3.h>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace rhbm_gem {

template<typename T>
struct SQLiteBinder
{
};

// int 32-bits specialization
template<>
struct SQLiteBinder<int>
{
    static int Bind(sqlite3_stmt * stmt, int index, int value)
    {
        return sqlite3_bind_int(stmt, index, value);
    }
};

template<>
struct SQLiteBinder<uint8_t>
{
    static int Bind(sqlite3_stmt * stmt, int index, uint8_t value)
    {
        return sqlite3_bind_int(stmt, index, static_cast<int>(value));
    }
};

template<>
struct SQLiteBinder<uint16_t>
{
    static int Bind(sqlite3_stmt * stmt, int index, uint16_t value)
    {
        return sqlite3_bind_int(stmt, index, static_cast<int>(value));
    }
};

template<>
struct SQLiteBinder<uint32_t>
{
    static int Bind(sqlite3_stmt * stmt, int index, uint32_t value)
    {
        return sqlite3_bind_int(stmt, index, static_cast<int>(value));
    }
};

// int 64-bits specialization
template<>
struct SQLiteBinder<int64_t>
{
    static int Bind(sqlite3_stmt * stmt, int index, int64_t value)
    {
        return sqlite3_bind_int64(stmt, index, value);
    }
};

// uint64_t specialization
template<>
struct SQLiteBinder<uint64_t>
{
    static int Bind(sqlite3_stmt * stmt, int index, uint64_t value)
    {
        return sqlite3_bind_int64(stmt, index,
                                  static_cast<sqlite3_int64>(value));
    }
};

// double specialization
template<>
struct SQLiteBinder<double>
{
    static int Bind(sqlite3_stmt * stmt, int index, double value)
    {
        return sqlite3_bind_double(stmt, index, value);
    }
};

// std::string specialization
template<>
struct SQLiteBinder<std::string>
{
    static int Bind(sqlite3_stmt * stmt, int index, const std::string & value)
    {
        // Use SQLITE_TRANSIENT so SQLite makes its own copy of the string
        return sqlite3_bind_text(stmt, index, value.c_str(), -1, SQLITE_TRANSIENT);
    }
};

// std::string_view specialization (since C++17)
template<>
struct SQLiteBinder<std::string_view>
{
    static int Bind(sqlite3_stmt * stmt, int index, std::string_view value)
    {
        return sqlite3_bind_text(stmt, index, value.data(),
                                 static_cast<int>(value.size()), SQLITE_TRANSIENT);
    }
};

// std::vector<float> specialization
template<>
struct SQLiteBinder<std::vector<float>>
{
    static int Bind(sqlite3_stmt * stmt, int index, const std::vector<float> & value)
    {
        return sqlite3_bind_blob(
            stmt,
            index,
            reinterpret_cast<const void*>(value.data()),
            static_cast<int>(value.size() * sizeof(float)),
            SQLITE_TRANSIENT
        );
    }
};

// std::vector<double> specialization
template<>
struct SQLiteBinder<std::vector<double>>
{
    static int Bind(sqlite3_stmt * stmt, int index, const std::vector<double> & value)
    {
        return sqlite3_bind_blob(
            stmt,
            index,
            reinterpret_cast<const void*>(value.data()),
            static_cast<int>(value.size() * sizeof(double)),
            SQLITE_TRANSIENT
        );
    }
};

// LocalPotentialSampleList specialization
template<>
struct SQLiteBinder<LocalPotentialSampleList>
{
    static int Bind(sqlite3_stmt * stmt, int index, const LocalPotentialSampleList & value)
    {
        if(value.empty())
        {
            return sqlite3_bind_blob(stmt, index, nullptr, 0, SQLITE_STATIC);
        }
        
        std::vector<float> contiguous;
        contiguous.reserve(value.size() * 2);
        for(const auto & sample : value)
        {
            contiguous.push_back(sample.distance);
            contiguous.push_back(sample.response);
        }
        
        return sqlite3_bind_blob(
            stmt,
            index,
            reinterpret_cast<const void*>(contiguous.data()),
            static_cast<int>(contiguous.size() * sizeof(float)),
            SQLITE_TRANSIENT
        );
    }
};

// std::vector<std::tuple<double, double>> specialization
template<>
struct SQLiteBinder<std::vector<std::tuple<double, double>>>
{
    static int Bind(sqlite3_stmt * stmt, int index, const std::vector<std::tuple<double, double>> & value)
    {
        if(value.empty())
        {
            return sqlite3_bind_blob(stmt, index, nullptr, 0, SQLITE_STATIC);
        }
        
        std::vector<double> contiguous;
        contiguous.reserve(value.size() * 2);
        for(const auto & tup : value)
        {
            contiguous.push_back(std::get<0>(tup));
            contiguous.push_back(std::get<1>(tup));
        }
        
        return sqlite3_bind_blob(
            stmt,
            index,
            reinterpret_cast<const void*>(contiguous.data()),
            static_cast<int>(contiguous.size() * sizeof(double)),
            SQLITE_TRANSIENT
        );
    }
};

} // namespace rhbm_gem
