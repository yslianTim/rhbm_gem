#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <mutex>

#include "GlobalEnumClass.hpp"

using AtomKey = uint32_t;

class AtomKeySystem
{
    static const AtomKey kDynamicBase;
    std::mutex m_mutex;
    AtomKey m_next_dynamic_key{ kDynamicBase };
    std::unordered_map<std::string, AtomKey> m_id_to_key_map;
    std::unordered_map<AtomKey, std::string> m_key_to_id_map;
    
public:
    static AtomKeySystem & Instance(void);
    void RegisterAtom(const std::string & atom_id);
    AtomKey GetAtomKey(const std::string & atom_id);
    AtomKey GetAtomKey(Element element, Remoteness remoteness, Branch branch);
    std::string GetAtomId(AtomKey atom_key);
    bool IsBuildInAtom(const std::string & atom_id) const;
    bool IsBuildInAtom(AtomKey atom_key) const;
    void ParseAtomId(
        const std::string & atom_id,
        const std::string & element_type,
        bool is_standard_monomer,
        Element & element,
        Remoteness & remoteness,
        Branch & branch
    );

private:
    AtomKeySystem(void);
    ~AtomKeySystem(void);

};
