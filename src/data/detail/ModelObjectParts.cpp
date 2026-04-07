#include "data/detail/ModelObjectParts.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem {

ModelObjectParts::ModelObjectParts() = default;
ModelObjectParts::~ModelObjectParts() = default;
ModelObjectParts::ModelObjectParts(ModelObjectParts && other) noexcept = default;
ModelObjectParts & ModelObjectParts::operator=(ModelObjectParts && other) noexcept = default;

ModelObject AssembleModelObject(ModelObjectParts parts)
{
    if (parts.component_key_system == nullptr)
    {
        parts.component_key_system = std::make_unique<ComponentKeySystem>();
    }
    if (parts.atom_key_system == nullptr)
    {
        parts.atom_key_system = std::make_unique<AtomKeySystem>();
    }
    if (parts.bond_key_system == nullptr)
    {
        parts.bond_key_system = std::make_unique<BondKeySystem>();
    }

    ModelObject model;
    model.m_parts = std::make_unique<ModelObjectParts>(std::move(parts));
    model.AttachOwnedObjects();
    model.InvalidateDerivedState();
    model.RebuildObjectIndex();
    model.RebuildSelection();
    return model;
}

} // namespace rhbm_gem
