#pragma once

#include <sqlite3.h>
#include <string>
#include <string_view>
#include <vector>
#include <tuple>

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

// int 64-bits specialization
template<>
struct SQLiteBinder<int64_t>
{
    static int Bind(sqlite3_stmt * stmt, int index, int64_t value)
    {
        return sqlite3_bind_int64(stmt, index, value);
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
        // SQLITE_TRANSIENT or SQLITE_STATIC
        return sqlite3_bind_text(stmt, index, value.c_str(), -1, SQLITE_STATIC);
    }
};

// std::string_view specialization (since C++17)
template<>
struct SQLiteBinder<std::string_view>
{
    static int Bind(sqlite3_stmt * stmt, int index, std::string_view value)
    {
        return sqlite3_bind_text(stmt, index, value.data(),
                                 static_cast<int>(value.size()), SQLITE_STATIC);
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
            SQLITE_STATIC
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
            SQLITE_STATIC
        );
    }
};

// std::vector<std::tuple<float, float>> specialization
template<>
struct SQLiteBinder<std::vector<std::tuple<float, float>>>
{
    static int Bind(sqlite3_stmt * stmt, int index, const std::vector<std::tuple<float, float>> & value)
    {
        if(value.empty())
        {
            return sqlite3_bind_blob(stmt, index, nullptr, 0, SQLITE_STATIC);
        }
        
        std::vector<float> contiguous;
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
            static_cast<int>(contiguous.size() * sizeof(float)),
            SQLITE_TRANSIENT
        );
    }
};

// std::vector<std::tuple<float, float, float>> specialization
template<>
struct SQLiteBinder<std::vector<std::tuple<float, float, float>>>
{
    static int Bind(sqlite3_stmt * stmt, int index, const std::vector<std::tuple<float, float, float>> & value)
    {
        if(value.empty())
        {
            return sqlite3_bind_blob(stmt, index, nullptr, 0, SQLITE_STATIC);
        }
        
        std::vector<float> contiguous;
        contiguous.reserve(value.size() * 3);
        for(const auto & tup : value)
        {
            contiguous.push_back(std::get<0>(tup));
            contiguous.push_back(std::get<1>(tup));
            contiguous.push_back(std::get<2>(tup));
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
