#pragma once

#include <iostream>
#include <tuple>

template <typename T>
struct SingleHash
{
    std::size_t operator()(const T & value) const
    {
        return std::hash<T>()(value);
    }
};

template <typename... Args>
struct TupleHash
{
    std::size_t operator()(const std::tuple<Args...> & tuple) const
    {
        return GetHashTuple(tuple, std::index_sequence_for<Args...>());
    }

private:
    template <std::size_t... Is>
    std::size_t GetHashTuple(const std::tuple<Args...> & tuple, std::index_sequence<Is...>) const
    {
        return (GetHashElement(std::get<Is>(tuple)) ^ ...);
    }
    // Recursively compute the hash value of each tuple element
    template <typename T>
    std::size_t GetHashElement(const T & element) const
    {
        return SingleHash<T>()(element);
    }
};