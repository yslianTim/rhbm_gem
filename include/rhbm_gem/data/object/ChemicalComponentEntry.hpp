#pragma once

#include <string>
#include <vector>
#include <utility>
#include <unordered_map>
#include <map>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/BondKeySystem.hpp>

namespace rhbm_gem {

struct ComponentAtomEntry
{
    std::string atom_id;
    Element element_type;
    bool aromatic_atom_flag;
    StereoChemistry stereo_config;
};

struct ComponentBondEntry
{
    std::string bond_id;
    BondType bond_type;
    BondOrder bond_order;
    bool aromatic_atom_flag;
    StereoChemistry stereo_config;
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
    ChemicalComponentEntry();
    ~ChemicalComponentEntry();

    void SetComponentId(const std::string & comp_id) { m_component_id = comp_id; }
    void SetComponentName(const std::string & name) { m_component_name = name; }
    void SetComponentType(const std::string & type) { m_component_type = type; }
    void SetComponentFormula(const std::string & formula) { m_component_formula = formula; }
    void SetComponentMolecularWeight(float weight) { m_component_molecular_weight = weight; }
    void SetStandardMonomerFlag(bool flag) { m_standard_monomer_flag = flag; }

    std::string GetComponentId() const { return m_component_id; }
    std::string GetComponentName() const { return m_component_name; }
    std::string GetComponentType() const { return m_component_type; }
    std::string GetComponentFormula() const { return m_component_formula; }
    float GetComponentMolecularWeight() const { return m_component_molecular_weight; }
    bool IsStandardMonomer() const { return m_standard_monomer_flag; }
    const std::map<AtomKey, ComponentAtomEntry> & AtomEntries() const { return m_component_atom_entry_map; }
    const std::map<BondKey, ComponentBondEntry> & BondEntries() const { return m_component_bond_entry_map; }
    void AddComponentAtomEntry(AtomKey atom_key, const ComponentAtomEntry & atom_info);
    void AddComponentBondEntry(BondKey bond_key, const ComponentBondEntry & bond_info);
    bool HasComponentBondEntry(BondKey bond_key) const;
    const ComponentAtomEntry * FindComponentAtomEntry(AtomKey atom_key) const;
    const ComponentBondEntry * FindComponentBondEntry(BondKey bond_key) const;

};

} // namespace rhbm_gem
