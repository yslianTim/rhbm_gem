#pragma once

#include <array>
#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <unordered_map>

#include "GlobalEnumClass.hpp"
#include "ComponentKeySystem.hpp"
#include "AtomKeySystem.hpp"

class AtomObject;
class BondObject;
class ChemicalComponentEntry;
struct ComponentAtomEntry;
struct ComponentBondEntry;

class AtomicModelDataBlock
{
    std::string m_map_id, m_model_id;
    std::string m_resolution, m_resolution_method;
    std::unique_ptr<ComponentKeySystem> m_component_key_system;
    std::unique_ptr<AtomKeySystem> m_atom_key_system;

    std::unordered_map<int, std::vector<std::unique_ptr<AtomObject>>> m_atom_object_list_map;
    std::unordered_map<int, std::vector<std::unique_ptr<BondObject>>> m_bond_object_list_map;
    
    std::unordered_map<std::string, Entity> m_entity_type_map; // key : entity_id
    std::unordered_map<std::string, int> m_molecules_size_map; // key : entity_id
    std::unordered_map<std::string, std::vector<std::string>> m_chain_id_list_map; // key : entity_id
    std::unordered_map<Entity, std::vector<std::string>> m_entity_id_list_map; // key: entity type
    
    std::unordered_map<std::string, int> m_struct_sheet_strand_map;
    std::unordered_map<std::string, std::array<std::string, 5>> m_struct_helix_range_map;
    std::unordered_map<std::string, std::array<std::string, 4>> m_struct_sheet_range_map;
    std::vector<Element> m_element_type_list;

    std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> m_chemical_component_entry_map;

public:
    AtomicModelDataBlock(void);
    ~AtomicModelDataBlock();

    void AddAtomObject(int model_number, std::unique_ptr<AtomObject> atom_object);
    void AddBondObject(int model_number, std::unique_ptr<BondObject> bond_object);
    void AddEntityTypeInEntityMap(const std::string & entity_id, Entity entity);
    void AddChainIDInEntityMap(const std::string & entity_id, const std::string & chain_id);
    void AddMoleculesSizeInEntityMap(const std::string & entity_id, int molecules_size);
    void AddSheetStrands(const std::string & sheet_id, int strands_size);
    void AddSheetRange(const std::string & composite_sheet_id, const std::array<std::string, 4> & range);
    void AddHelixRange(const std::string & helix_id, const std::array<std::string, 5> & range);
    void AddElementType(const Element & element);
    void AddChemicalComponentEntry(ComponentKey comp_key, std::unique_ptr<ChemicalComponentEntry> entry);
    void AddComponentAtomEntry(
        ComponentKey comp_key,
        AtomKey atom_key,
        const ComponentAtomEntry & atom_entry
    );
    void AddComponentBondEntry(
        ComponentKey comp_key,
        const std::pair<AtomKey, AtomKey> & atom_key_pair,
        const ComponentBondEntry & bond_entry
    );

    void SetPdbID(const std::string & label) { m_model_id = label; }
    void SetEmdID(const std::string & label) { m_map_id = label; }
    void SetResolution(const std::string & value) { m_resolution = value; }
    void SetResolutionMethod(const std::string & value) { m_resolution_method = value; }
    void SetStructureInfo(AtomObject * atom_object);

    const std::string & GetPdbID(void) const;
    const std::string & GetEmdID(void) const;
    double GetResolution(void) const;
    const std::string & GetResolutionMethod(void) const;
    const std::vector<Element> & GetElementTypeList(void) const;
    const std::unordered_map<std::string, Entity> & GetEntityTypeMap(void) const;
    const std::unordered_map<std::string, int> & GetMoleculesSizeMap(void) const;
    const std::unordered_map<Entity, std::vector<std::string>> & GetEntityIDListMap(void) const;
    const std::unordered_map<std::string, std::vector<std::string>> & GetChainIDListMap(void) const;
    ChemicalComponentEntry * GetChemicalComponentEntryPtr(ComponentKey key);
    std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> & GetChemicalComponentEntryMap(void);
    ComponentKeySystem * GetComponentKeySystemPtr(void) { return m_component_key_system.get(); }
    AtomKeySystem * GetAtomKeySystemPtr(void) { return m_atom_key_system.get(); }
    std::vector<std::unique_ptr<AtomObject>> MoveAtomObjectList(int model_number=1);
    std::vector<std::unique_ptr<BondObject>> MoveBondObjectList(int model_number=1);
    std::unique_ptr<ComponentKeySystem> MoveComponentKeySystem(void) { return std::move(m_component_key_system); }
    std::unique_ptr<AtomKeySystem> MoveAtomKeySystem(void) { return std::move(m_atom_key_system); }

};
