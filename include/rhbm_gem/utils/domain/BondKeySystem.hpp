#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <mutex>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

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
    BondKeySystem();
    ~BondKeySystem();
    BondKeySystem(const BondKeySystem & other);

    void RegisterBond(const std::string & bond_id);
    void RegisterBond(const std::string & atom_id_1, const std::string & atom_id_2);
    void RegisterBond(const std::string & bond_id, BondKey bond_key);
    void RegisterBond(const std::string & atom_id_1, const std::string & atom_id_2, BondKey bond_key);
    BondKey GetBondKey(const std::string & bond_id);
    BondKey GetBondKey(const std::string & atom_id_1, const std::string & atom_id_2);
    std::string GetBondId(BondKey bond_key);
    bool IsRegistedBond(const std::string & atom_id_1, const std::string & atom_id_2) const;

private:
    std::string BuildBondIdFromAtomIdPair(const std::string & atom_id_1, const std::string & atom_id_2) const;

};
