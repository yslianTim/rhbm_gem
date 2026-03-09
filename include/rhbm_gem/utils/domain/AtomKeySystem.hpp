#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <mutex>

#include <rhbm_gem/utils/GlobalEnumClass.hpp>

using AtomKey = uint16_t;

class AtomKeySystem
{
    static const AtomKey k_dynamic_base;
    static const AtomKey k_max_key;
    std::mutex m_mutex;
    AtomKey m_next_dynamic_key;
    std::unordered_map<std::string, AtomKey> m_id_to_key_map;
    std::unordered_map<AtomKey, std::string> m_key_to_id_map;
    
public:
    AtomKeySystem();
    ~AtomKeySystem();
    AtomKeySystem(const AtomKeySystem & other);

    void RegisterAtom(const std::string & atom_id);
    void RegisterAtom(const std::string & atom_id, AtomKey atom_key);
    AtomKey GetAtomKey(const std::string & atom_id);
    std::string GetAtomId(AtomKey atom_key);
    bool IsBuildInAtom(const std::string & atom_id) const;
    bool IsBuildInAtom(AtomKey atom_key) const;

};
