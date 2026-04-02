#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/BondKeySystem.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class ChemicalComponentEntry;

class ModelObjectBuilder
{
    ModelObject m_model_object;

public:
    ModelObjectBuilder() = default;

    explicit ModelObjectBuilder(std::vector<std::unique_ptr<AtomObject>> atom_list)
    {
        SetAtomList(std::move(atom_list));
    }

    void AddAtom(std::unique_ptr<AtomObject> atom)
    {
        m_model_object.m_atom_list.emplace_back(std::move(atom));
    }

    void AddBond(std::unique_ptr<BondObject> bond)
    {
        m_model_object.m_bond_list.emplace_back(std::move(bond));
    }

    void SetAtomList(std::vector<std::unique_ptr<AtomObject>> atom_list)
    {
        m_model_object.m_atom_list = std::move(atom_list);
    }

    void SetBondList(std::vector<std::unique_ptr<BondObject>> bond_list)
    {
        m_model_object.m_bond_list = std::move(bond_list);
    }

    void SetChainIDListMap(
        const std::unordered_map<std::string, std::vector<std::string>> & value)
    {
        m_model_object.m_chain_id_list_map = value;
    }

    void AddChemicalComponentEntry(
        ComponentKey component_key,
        std::unique_ptr<ChemicalComponentEntry> entry)
    {
        m_model_object.m_chemical_component_entry_map[component_key] = std::move(entry);
    }

    void SetChemicalComponentEntryMap(
        std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> entry_map)
    {
        m_model_object.m_chemical_component_entry_map = std::move(entry_map);
    }

    ChemicalComponentEntry * FindChemicalComponentEntry(ComponentKey component_key)
    {
        const auto iter{ m_model_object.m_chemical_component_entry_map.find(component_key) };
        return iter == m_model_object.m_chemical_component_entry_map.end() ?
            nullptr : iter->second.get();
    }

    ComponentKeySystem & MutableComponentKeySystem()
    {
        return *m_model_object.m_component_key_system;
    }

    AtomKeySystem & MutableAtomKeySystem()
    {
        return *m_model_object.m_atom_key_system;
    }

    BondKeySystem & MutableBondKeySystem()
    {
        return *m_model_object.m_bond_key_system;
    }

    void SetComponentKeySystem(std::unique_ptr<ComponentKeySystem> component_key_system)
    {
        m_model_object.m_component_key_system = std::move(component_key_system);
    }

    void SetAtomKeySystem(std::unique_ptr<AtomKeySystem> atom_key_system)
    {
        m_model_object.m_atom_key_system = std::move(atom_key_system);
    }

    void SetBondKeySystem(std::unique_ptr<BondKeySystem> bond_key_system)
    {
        m_model_object.m_bond_key_system = std::move(bond_key_system);
    }

    ModelObject Build()
    {
        ModelObject built_model{ std::move(m_model_object) };
        built_model.AttachOwnedObjects();
        built_model.SyncDerivedState();
        return built_model;
    }
};

} // namespace rhbm_gem
