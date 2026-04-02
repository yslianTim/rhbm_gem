#pragma once

#include <memory>
#include <stdexcept>
#include <string>

#include "data/detail/AtomObjectAccess.hpp"
#include "data/detail/BondObjectAccess.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>

namespace rhbm_gem {

inline LocalPotentialEntry * GetLocalPotentialEntry(AtomObject & atom_object)
{
    return AtomObjectAccess::LocalPotentialEntryPtr(atom_object);
}

inline const LocalPotentialEntry * GetLocalPotentialEntry(const AtomObject & atom_object)
{
    return AtomObjectAccess::LocalPotentialEntryPtr(atom_object);
}

inline LocalPotentialEntry * GetLocalPotentialEntry(BondObject & bond_object)
{
    return BondObjectAccess::LocalPotentialEntryPtr(bond_object);
}

inline const LocalPotentialEntry * GetLocalPotentialEntry(const BondObject & bond_object)
{
    return BondObjectAccess::LocalPotentialEntryPtr(bond_object);
}

inline void SetLocalPotentialEntry(AtomObject & atom_object, std::unique_ptr<LocalPotentialEntry> entry)
{
    AtomObjectAccess::SetLocalPotentialEntry(atom_object, std::move(entry));
}

inline void SetLocalPotentialEntry(BondObject & bond_object, std::unique_ptr<LocalPotentialEntry> entry)
{
    BondObjectAccess::SetLocalPotentialEntry(bond_object, std::move(entry));
}

inline const LocalPotentialEntry & RequireLocalPotentialEntry(
    const LocalPotentialEntry * entry,
    const char * context)
{
    if (entry == nullptr)
    {
        throw std::runtime_error(std::string(context) + " is not available.");
    }
    return *entry;
}

inline const LocalPotentialEntry & RequireLocalPotentialEntry(const AtomObject & atom_object)
{
    return RequireLocalPotentialEntry(GetLocalPotentialEntry(atom_object), "Atom local entry");
}

inline const LocalPotentialEntry & RequireLocalPotentialEntry(const BondObject & bond_object)
{
    return RequireLocalPotentialEntry(GetLocalPotentialEntry(bond_object), "Bond local entry");
}

} // namespace rhbm_gem
