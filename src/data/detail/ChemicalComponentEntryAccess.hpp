#pragma once

#include <map>
#include <memory>
#include <string>

#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>

namespace rhbm_gem {

class ChemicalComponentEntryAccess
{
public:
    static void SetComponentId(ChemicalComponentEntry & entry, const std::string & comp_id)
    {
        entry.m_component_id = comp_id;
    }
    static void SetComponentName(ChemicalComponentEntry & entry, const std::string & name)
    {
        entry.m_component_name = name;
    }
    static void SetComponentType(ChemicalComponentEntry & entry, const std::string & type)
    {
        entry.m_component_type = type;
    }
    static void SetComponentFormula(ChemicalComponentEntry & entry, const std::string & formula)
    {
        entry.m_component_formula = formula;
    }
    static void SetComponentMolecularWeight(ChemicalComponentEntry & entry, float weight)
    {
        entry.m_component_molecular_weight = weight;
    }
    static void SetStandardMonomerFlag(ChemicalComponentEntry & entry, bool flag)
    {
        entry.m_standard_monomer_flag = flag;
    }
    static void SetComponentAtomEntryMap(
        ChemicalComponentEntry & entry,
        std::map<AtomKey, ComponentAtomEntry> atom_entry_map)
    {
        entry.m_component_atom_entry_map = std::move(atom_entry_map);
    }
    static void SetComponentBondEntryMap(
        ChemicalComponentEntry & entry,
        std::map<BondKey, ComponentBondEntry> bond_entry_map)
    {
        entry.m_component_bond_entry_map = std::move(bond_entry_map);
    }
    static void AddComponentAtomEntry(
        ChemicalComponentEntry & entry,
        AtomKey atom_key,
        const ComponentAtomEntry & atom_info)
    {
        entry.AddComponentAtomEntry(atom_key, atom_info);
    }
    static void AddComponentBondEntry(
        ChemicalComponentEntry & entry,
        BondKey bond_key,
        const ComponentBondEntry & bond_info)
    {
        entry.AddComponentBondEntry(bond_key, bond_info);
    }
    static const ComponentAtomEntry * FindComponentAtomEntry(
        const ChemicalComponentEntry & entry,
        AtomKey atom_key)
    {
        return entry.FindComponentAtomEntry(atom_key);
    }
    static const ComponentBondEntry * FindComponentBondEntry(
        const ChemicalComponentEntry & entry,
        BondKey bond_key)
    {
        return entry.FindComponentBondEntry(bond_key);
    }
};

} // namespace rhbm_gem
