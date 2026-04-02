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
class ModelBuilder
{
public:
    static void AddAtom(ModelObject & model_object, std::unique_ptr<AtomObject> atom)
    {
        model_object.AddAtom(std::move(atom));
    }

    static void AddBond(ModelObject & model_object, std::unique_ptr<BondObject> bond)
    {
        model_object.AddBond(std::move(bond));
    }

    static void SetAtomList(
        ModelObject & model_object,
        std::vector<std::unique_ptr<AtomObject>> atom_list)
    {
        model_object.SetAtomList(std::move(atom_list));
    }

    static void SetBondList(
        ModelObject & model_object,
        std::vector<std::unique_ptr<BondObject>> bond_list)
    {
        model_object.SetBondList(std::move(bond_list));
    }

    static void SetChainIDListMap(
        ModelObject & model_object,
        const std::unordered_map<std::string, std::vector<std::string>> & value)
    {
        model_object.SetChainIDListMap(value);
    }

    static void AddChemicalComponentEntry(
        ModelObject & model_object,
        ComponentKey component_key,
        std::unique_ptr<ChemicalComponentEntry> entry)
    {
        model_object.AddChemicalComponentEntry(component_key, std::move(entry));
    }

    static void SetChemicalComponentEntryMap(
        ModelObject & model_object,
        std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> entry_map)
    {
        model_object.SetChemicalComponentEntryMap(std::move(entry_map));
    }

    static std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> &
    ChemicalComponentEntryMapRef(ModelObject & model_object)
    {
        return model_object.m_chemical_component_entry_map;
    }

    static void SetComponentKeySystem(
        ModelObject & model_object,
        std::unique_ptr<ComponentKeySystem> component_key_system)
    {
        model_object.SetComponentKeySystem(std::move(component_key_system));
    }

    static void SetAtomKeySystem(
        ModelObject & model_object,
        std::unique_ptr<AtomKeySystem> atom_key_system)
    {
        model_object.SetAtomKeySystem(std::move(atom_key_system));
    }

    static void SetBondKeySystem(
        ModelObject & model_object,
        std::unique_ptr<BondKeySystem> bond_key_system)
    {
        model_object.SetBondKeySystem(std::move(bond_key_system));
    }

    static ComponentKeySystem & ComponentKeySystemRef(ModelObject & model_object)
    {
        return *model_object.m_component_key_system;
    }

    static AtomKeySystem & AtomKeySystemRef(ModelObject & model_object)
    {
        return *model_object.m_atom_key_system;
    }

    static BondKeySystem & BondKeySystemRef(ModelObject & model_object)
    {
        return *model_object.m_bond_key_system;
    }

    static void FinalizeLoad(ModelObject & model_object)
    {
        model_object.FinalizeLoad();
    }
};

} // namespace rhbm_gem
