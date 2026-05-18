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
    model.m_atom_list = std::move(parts.atom_list);
    model.m_bond_list = std::move(parts.bond_list);
    model.m_chain_id_list_map = std::move(parts.chain_id_list_map);
    model.m_chemical_component_entry_map = std::move(parts.chemical_component_entry_map);
    model.m_component_key_system = std::move(parts.component_key_system);
    model.m_atom_key_system = std::move(parts.atom_key_system);
    model.m_bond_key_system = std::move(parts.bond_key_system);
    model.AttachOwnedObjects();
    model.InvalidateDerivedState();
    model.RebuildObjectIndex();
    model.RebuildSelection();
    return model;
}

} // namespace rhbm_gem
