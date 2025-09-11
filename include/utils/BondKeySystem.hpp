#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <mutex>
#include <utility>

#include "GlobalEnumClass.hpp"

using BondKey = uint16_t;

class BondKeySystem
{
    static const BondKey k_dynamic_base;
    static const BondKey k_max_key;
    std::mutex m_mutex;
    BondKey m_next_dynamic_key;
    std::unordered_map<std::string, BondKey> m_id_to_key_map;
    std::unordered_map<BondKey, std::string> m_key_to_id_map;
    
public:
    BondKeySystem(void);
    ~BondKeySystem();
    BondKeySystem(const BondKeySystem & other);

    void RegisterBond(const std::string & atom_id_1, const std::string & atom_id_2);
    void RegisterBond(const std::string & atom_id_1, const std::string & atom_id_2, BondKey bond_key);
    BondKey GetBondKey(const std::string & atom_id_1, const std::string & atom_id_2);
    std::string GetBondId(BondKey bond_key);
    bool IsBuildInBond(const std::string & atom_id_1, const std::string & atom_id_2) const;
    bool IsBuildInBond(BondKey bond_key) const;

private:
    std::string BuildBondIdFromAtomIdPair(const std::string & atom_id_1, const std::string & atom_id_2) const;
    std::pair<std::string, std::string> BuildAtomIdPairFromBondId(const std::string & bond_id) const;

};
