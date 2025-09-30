#pragma once

#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include <map>

#include "GlobalEnumClass.hpp"
#include "AtomKeySystem.hpp"
#include "BondKeySystem.hpp"

struct ComponentAtomEntry
{
    std::string atom_id;
    Element element_type;
    bool aromatic_atom_flag;
    char chiral_config; // 'N', 'R', 'S'
};

struct ComponentBondEntry
{
    std::pair<std::string, std::string> atom_id_pair;
    BondType bond_type;
    BondOrder bond_order;
    bool aromatic_atom_flag;
    char chiral_config; // 'N', 'R', 'S'
};

class ChemicalComponentEntry
{
    std::string m_component_id;
    std::string m_component_name;
    std::string m_component_type;
    std::string m_component_formula;
    float m_component_molecular_weight;
    bool m_standard_monomer_flag;
    std::map<AtomKey, ComponentAtomEntry> m_component_atom_entry_map;
    std::map<BondKey, ComponentBondEntry> m_component_bond_entry_map;

public:
    ChemicalComponentEntry(void);
    ~ChemicalComponentEntry();

    void SetComponentId(const std::string & comp_id) { m_component_id = comp_id; }
    void SetComponentName(const std::string & name) { m_component_name = name; }
    void SetComponentType(const std::string & type) { m_component_type = type; }
    void SetComponentFormula(const std::string & formula) { m_component_formula = formula; }
    void SetComponentMolecularWeight(float weight) { m_component_molecular_weight = weight; }
    void SetStandardMonomerFlag(bool flag) { m_standard_monomer_flag = flag; }
    void SetComponentAtomEntryMap(std::map<AtomKey, ComponentAtomEntry> & atom_entry_map);
    void SetComponentBondEntryMap(std::map<BondKey, ComponentBondEntry> & bond_entry_map);
    void AddComponentAtomEntry(AtomKey atom_key, const ComponentAtomEntry & atom_info);
    void AddComponentBondEntry(BondKey bond_key, const ComponentBondEntry & bond_info);

    bool HasComponentAtomEntry(AtomKey atom_key) const;
    bool HasComponentBondEntry(BondKey bond_key) const;
    std::string GetComponentId(void) const { return m_component_id; }
    std::string GetComponentName(void) const { return m_component_name; }
    std::string GetComponentType(void) const { return m_component_type; }
    std::string GetComponentFormula(void) const { return m_component_formula; }
    float GetComponentMolecularWeight(void) const { return m_component_molecular_weight; }
    bool IsStandardMonomer(void) const { return m_standard_monomer_flag; }
    const std::map<AtomKey, ComponentAtomEntry> &
    GetComponentAtomEntryMap(void) const { return m_component_atom_entry_map; }
    const std::map<BondKey, ComponentBondEntry> &
    GetComponentBondEntryMap(void) const { return m_component_bond_entry_map; }
    const ComponentAtomEntry * GetComponentAtomEntryPtr(AtomKey atom_key) const;
    const ComponentBondEntry * GetComponentBondEntryPtr(BondKey bond_key) const;

};
