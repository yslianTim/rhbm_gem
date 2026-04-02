#pragma once

#include <stdexcept>
#include <string>

#include "data/detail/ModelAnalysisAccess.hpp"
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

inline const LocalPotentialEntry * FindLocalPotentialEntry(const AtomObject & atom_object)
{
    return ModelAnalysisAccess::FindLocalEntryConst(atom_object);
}

inline const LocalPotentialEntry * FindLocalPotentialEntry(const BondObject & bond_object)
{
    return ModelAnalysisAccess::FindLocalEntryConst(bond_object);
}

inline const LocalPotentialEntry & RequireLocalPotentialEntry(const AtomObject & atom_object)
{
    return RequireLocalPotentialEntry(FindLocalPotentialEntry(atom_object), "Atom local entry");
}

inline const LocalPotentialEntry & RequireLocalPotentialEntry(const BondObject & bond_object)
{
    return RequireLocalPotentialEntry(FindLocalPotentialEntry(bond_object), "Bond local entry");
}

} // namespace rhbm_gem
