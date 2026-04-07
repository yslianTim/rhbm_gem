#include "data/detail/ModelObjectAssembly.hpp"

#include "data/detail/ModelAnalysisData.hpp"

namespace rhbm_gem {

ModelObject AssembleModelObject(ModelObjectAssembly assembly)
{
    ModelObject model;
    model.m_atom_list = std::move(assembly.atom_list);
    model.m_bond_list = std::move(assembly.bond_list);
    model.m_chain_id_list_map = std::move(assembly.chain_id_list_map);
    model.m_chemical_component_entry_map = std::move(assembly.chemical_component_entry_map);
    model.m_component_key_system =
        assembly.component_key_system != nullptr ?
            std::move(assembly.component_key_system) :
            std::make_unique<ComponentKeySystem>();
    model.m_atom_key_system =
        assembly.atom_key_system != nullptr ?
            std::move(assembly.atom_key_system) :
            std::make_unique<AtomKeySystem>();
    model.m_bond_key_system =
        assembly.bond_key_system != nullptr ?
            std::move(assembly.bond_key_system) :
            std::make_unique<BondKeySystem>();

    model.AttachOwnedObjects();
    model.SyncDerivedState();
    return model;
}

} // namespace rhbm_gem
