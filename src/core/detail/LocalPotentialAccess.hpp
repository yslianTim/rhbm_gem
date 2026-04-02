#pragma once

#include <stdexcept>
#include <string>

#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>

namespace rhbm_gem {

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
    return RequireLocalPotentialEntry(atom_object.GetLocalPotentialEntry(), "Atom local entry");
}

inline const LocalPotentialEntry & RequireLocalPotentialEntry(const BondObject & bond_object)
{
    return RequireLocalPotentialEntry(bond_object.GetLocalPotentialEntry(), "Bond local entry");
}

} // namespace rhbm_gem
