#pragma once

#include <iostream>
#include <tuple>
#include <functional>

// 輔助函式：類似 boost::hash_combine 的實作
inline void hash_combine(std::size_t & seed, std::size_t hash_value)
{
    seed ^= hash_value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// 輔助函式：利用折疊表達式計算 tuple 各元素 hash 的組合值
template <typename Tuple, std::size_t... Is>
std::size_t hash_tuple_impl(const Tuple & t, std::index_sequence<Is...>)
{
    std::size_t seed = 0;
    // 依序將 tuple 中每個元素的 hash 組合到 seed 中
    ((hash_combine(seed, std::hash<std::tuple_element_t<Is, Tuple>>{}(std::get<Is>(t)))), ...);
    return seed;
}

// 通用的 Tuple hash functor：適用於任何 std::tuple<Ts...>
template <typename... Ts>
struct TupleHash
{
    std::size_t operator()(const std::tuple<Ts...>& t) const
    {
        return hash_tuple_impl(t, std::index_sequence_for<Ts...>{});
    }
};

// 或者定義一個泛用的模板特化，讓使用者只需要傳入整個 tuple 型態
template <typename T>
struct GenericTupleHash;

template <typename... Ts>
struct GenericTupleHash<std::tuple<Ts...>>
{
    std::size_t operator()(const std::tuple<Ts...>& t) const
    {
        return hash_tuple_impl(t, std::index_sequence_for<Ts...>{});
    }
};