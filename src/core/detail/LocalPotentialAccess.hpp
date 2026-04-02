#pragma once

#include <memory>
#include <stdexcept>
#include <string>

#include "data/detail/ModelAnalysisAccess.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

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

inline LocalPotentialEntry & EnsureLocalPotentialEntry(
    ModelObject & model_object,
    const AtomObject & atom_object)
{
    return ModelAnalysisAccess::Mutable(model_object).Atoms().EnsureLocalEntry(atom_object);
}

inline LocalPotentialEntry & EnsureLocalPotentialEntry(
    ModelObject & model_object,
    const BondObject & bond_object)
{
    return ModelAnalysisAccess::Mutable(model_object).Bonds().EnsureLocalEntry(bond_object);
}

inline void SetLocalPotentialEntry(
    ModelObject & model_object,
    const AtomObject & atom_object,
    std::unique_ptr<LocalPotentialEntry> entry)
{
    ModelAnalysisAccess::Mutable(model_object).Atoms().SetLocalEntry(atom_object, std::move(entry));
}

inline void SetLocalPotentialEntry(
    ModelObject & model_object,
    const BondObject & bond_object,
    std::unique_ptr<LocalPotentialEntry> entry)
{
    ModelAnalysisAccess::Mutable(model_object).Bonds().SetLocalEntry(bond_object, std::move(entry));
}

inline LocalPotentialEntry * FindLocalPotentialEntry(
    ModelObject & model_object,
    const AtomObject & atom_object)
{
    return ModelAnalysisAccess::Mutable(model_object).Atoms().FindLocalEntry(atom_object);
}

inline const LocalPotentialEntry * FindLocalPotentialEntry(
    const ModelObject & model_object,
    const AtomObject & atom_object)
{
    return ModelAnalysisAccess::Read(model_object).Atoms().FindLocalEntry(atom_object);
}

inline LocalPotentialEntry * FindLocalPotentialEntry(
    ModelObject & model_object,
    const BondObject & bond_object)
{
    return ModelAnalysisAccess::Mutable(model_object).Bonds().FindLocalEntry(bond_object);
}

inline const LocalPotentialEntry * FindLocalPotentialEntry(
    const ModelObject & model_object,
    const BondObject & bond_object)
{
    return ModelAnalysisAccess::Read(model_object).Bonds().FindLocalEntry(bond_object);
}

inline const LocalPotentialEntry * FindLocalPotentialEntry(const AtomObject & atom_object)
{
    auto * owner{ ModelAnalysisAccess::OwnerOf(atom_object) };
    return owner == nullptr ? nullptr
                            : ModelAnalysisAccess::Read(*owner).Atoms().FindLocalEntry(atom_object);
}

inline const LocalPotentialEntry * FindLocalPotentialEntry(const BondObject & bond_object)
{
    auto * owner{ ModelAnalysisAccess::OwnerOf(bond_object) };
    return owner == nullptr ? nullptr
                            : ModelAnalysisAccess::Read(*owner).Bonds().FindLocalEntry(bond_object);
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
