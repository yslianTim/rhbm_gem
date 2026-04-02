#pragma once

#include <memory>

#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/data/object/BondObject.hpp>

namespace rhbm_gem {

class BondObjectAccess
{
public:
    static void SetSelectedFlag(BondObject & bond_object, bool value) { bond_object.m_is_selected = value; }
    static void SetSpecialBondFlag(BondObject & bond_object, bool value) { bond_object.m_is_special_bond = value; }
    static void SetBondKey(BondObject & bond_object, BondKey value) { bond_object.m_bond_key = value; }
    static void SetBondType(BondObject & bond_object, BondType value) { bond_object.m_bond_type = value; }
    static void SetBondOrder(BondObject & bond_object, BondOrder value) { bond_object.m_bond_order = value; }
    static void SetLocalPotentialEntry(
        BondObject & bond_object,
        std::unique_ptr<LocalPotentialEntry> entry)
    {
        bond_object.SetLocalPotentialEntry(std::move(entry));
    }
    static LocalPotentialEntry * LocalPotentialEntryPtr(BondObject & bond_object)
    {
        return bond_object.m_local_potential_entry.get();
    }
    static const LocalPotentialEntry * LocalPotentialEntryPtr(const BondObject & bond_object)
    {
        return bond_object.m_local_potential_entry.get();
    }
};

} // namespace rhbm_gem
