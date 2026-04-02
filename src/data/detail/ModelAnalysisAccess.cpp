#include "data/detail/ModelAnalysisAccess.hpp"

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/ModelSpatialCache.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>

#include <stdexcept>
#include <string>

namespace rhbm_gem {

ModelAnalysisData & ModelAnalysisAccess::Mutable(ModelObject & model_object)
{
    return *model_object.m_analysis_data;
}

const ModelAnalysisData & ModelAnalysisAccess::Read(const ModelObject & model_object)
{
    return *model_object.m_analysis_data;
}

void ModelAnalysisAccess::ClearFitStates(ModelObject & model_object)
{
    Mutable(model_object).ClearFitStates();
}

ModelObject * ModelAnalysisAccess::OwnerOf(const AtomObject & atom_object)
{
    return atom_object.m_owner_model;
}

ModelObject * ModelAnalysisAccess::OwnerOf(const BondObject & bond_object)
{
    return bond_object.m_owner_model;
}

LocalPotentialEntry & ModelAnalysisAccess::EnsureLocalEntry(
    ModelObject & model_object,
    const AtomObject & atom_object)
{
    return Mutable(model_object).Atoms().EnsureLocalEntry(atom_object);
}

LocalPotentialEntry & ModelAnalysisAccess::EnsureLocalEntry(
    ModelObject & model_object,
    const BondObject & bond_object)
{
    return Mutable(model_object).Bonds().EnsureLocalEntry(bond_object);
}

void ModelAnalysisAccess::SetLocalEntry(
    ModelObject & model_object,
    const AtomObject & atom_object,
    std::unique_ptr<LocalPotentialEntry> entry)
{
    Mutable(model_object).Atoms().SetLocalEntry(atom_object, std::move(entry));
}

void ModelAnalysisAccess::SetLocalEntry(
    ModelObject & model_object,
    const BondObject & bond_object,
    std::unique_ptr<LocalPotentialEntry> entry)
{
    Mutable(model_object).Bonds().SetLocalEntry(bond_object, std::move(entry));
}

LocalPotentialEntry * ModelAnalysisAccess::FindLocalEntry(
    ModelObject & model_object,
    const AtomObject & atom_object)
{
    return Mutable(model_object).Atoms().FindLocalEntry(atom_object);
}

const LocalPotentialEntry * ModelAnalysisAccess::FindLocalEntry(
    const ModelObject & model_object,
    const AtomObject & atom_object)
{
    return Read(model_object).Atoms().FindLocalEntry(atom_object);
}

LocalPotentialEntry * ModelAnalysisAccess::FindLocalEntry(
    ModelObject & model_object,
    const BondObject & bond_object)
{
    return Mutable(model_object).Bonds().FindLocalEntry(bond_object);
}

const LocalPotentialEntry * ModelAnalysisAccess::FindLocalEntry(
    const ModelObject & model_object,
    const BondObject & bond_object)
{
    return Read(model_object).Bonds().FindLocalEntry(bond_object);
}

const LocalPotentialEntry * ModelAnalysisAccess::FindLocalEntry(const AtomObject & atom_object)
{
    auto * owner{ OwnerOf(atom_object) };
    return owner == nullptr ? nullptr : Read(*owner).Atoms().FindLocalEntry(atom_object);
}

const LocalPotentialEntry * ModelAnalysisAccess::FindLocalEntry(const BondObject & bond_object)
{
    auto * owner{ OwnerOf(bond_object) };
    return owner == nullptr ? nullptr : Read(*owner).Bonds().FindLocalEntry(bond_object);
}

const LocalPotentialEntry & ModelAnalysisAccess::RequireLocalEntry(
    const LocalPotentialEntry * entry,
    const char * context)
{
    if (entry == nullptr)
    {
        throw std::runtime_error(std::string(context) + " is not available.");
    }
    return *entry;
}

const LocalPotentialEntry & ModelAnalysisAccess::RequireLocalEntry(const AtomObject & atom_object)
{
    return RequireLocalEntry(FindLocalEntry(atom_object), "Atom local entry");
}

const LocalPotentialEntry & ModelAnalysisAccess::RequireLocalEntry(const BondObject & bond_object)
{
    return RequireLocalEntry(FindLocalEntry(bond_object), "Bond local entry");
}

GroupPotentialEntry & ModelAnalysisAccess::EnsureAtomGroupEntry(
    ModelObject & model_object,
    const std::string & class_key)
{
    return Mutable(model_object).Atoms().EnsureGroupEntry(class_key);
}

GroupPotentialEntry & ModelAnalysisAccess::EnsureBondGroupEntry(
    ModelObject & model_object,
    const std::string & class_key)
{
    return Mutable(model_object).Bonds().EnsureGroupEntry(class_key);
}

GroupPotentialEntry * ModelAnalysisAccess::FindAtomGroupEntry(
    ModelObject & model_object,
    const std::string & class_key)
{
    return Mutable(model_object).Atoms().FindGroupEntry(class_key);
}

const GroupPotentialEntry * ModelAnalysisAccess::FindAtomGroupEntry(
    const ModelObject & model_object,
    const std::string & class_key)
{
    return Read(model_object).Atoms().FindGroupEntry(class_key);
}

GroupPotentialEntry * ModelAnalysisAccess::FindBondGroupEntry(
    ModelObject & model_object,
    const std::string & class_key)
{
    return Mutable(model_object).Bonds().FindGroupEntry(class_key);
}

const GroupPotentialEntry * ModelAnalysisAccess::FindBondGroupEntry(
    const ModelObject & model_object,
    const std::string & class_key)
{
    return Read(model_object).Bonds().FindGroupEntry(class_key);
}

std::vector<AtomObject *> ModelAnalysisAccess::FindAtomsInRange(
    ModelObject & model_object,
    const AtomObject & center_atom,
    double range)
{
    model_object.EnsureKDTreeRoot();
    if (model_object.m_spatial_cache == nullptr ||
        model_object.m_spatial_cache->kd_tree_root == nullptr)
    {
        return {};
    }

    return KDTreeAlgorithm<AtomObject>::RangeSearch(
        model_object.m_spatial_cache->kd_tree_root.get(),
        const_cast<AtomObject *>(&center_atom),
        range);
}

} // namespace rhbm_gem
