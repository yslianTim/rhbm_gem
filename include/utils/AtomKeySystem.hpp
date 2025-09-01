#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <mutex>

using AtomKey = uint16_t;

class AtomKeySystem
{
    static constexpr AtomKey kDynamicBase{ 1000 };
    std::mutex m_mutex;
    AtomKey m_next_dynamic_key{ kDynamicBase };
    std::unordered_map<std::string, AtomKey> m_id_to_key_map;
    std::unordered_map<AtomKey, std::string> m_key_to_id_map;
    
public:
    static AtomKeySystem & Instance(void);
    void RegisterAtom(const std::string & id);
    AtomKey GetAtomKey(const std::string & id);
    std::string GetComponentId(AtomKey key);

private:
    AtomKeySystem(void);
    ~AtomKeySystem(void);

};
