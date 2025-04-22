#pragma once

#include <unordered_map>
#include <optional>

template<typename KeyPacker, typename ValueType>
class CompositeKeyMap
{
public:
    template<typename... Args>
    void Set(Args... args, ValueType const & value)
    {
        m_data_map[ KeyPacker::Pack(args...) ] = value;
    }

    template<typename... Args>
    std::optional<std::reference_wrapper<ValueType>> Get(Args... args)
    {
        uint64_t key{ KeyPacker::Pack(args...) };
        auto it{ m_data_map.find(key) };
        if (it == m_data_map.end()) return std::nullopt;
        return it->second;
    }

    template<typename... Args>
    bool Remove(Args... args)
    {
        return m_data_map.erase( KeyPacker::Pack(args...) ) > 0;
    }

private:
    std::unordered_map<uint64_t, ValueType> m_data_map;

};