#pragma once

#include <array>
#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <tuple>
#include <map>
#include <unordered_map>

#include <rhbm_gem/utils/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/ComponentKeySystem.hpp>
#include <rhbm_gem/utils/AtomKeySystem.hpp>
#include <rhbm_gem/utils/BondKeySystem.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class ChemicalComponentEntry;
struct ComponentAtomEntry;
struct ComponentBondEntry;
struct SheetRange
{
    std::string chain_id_beg;
    int seq_id_beg;
    std::string chain_id_end;
    int seq_id_end;
};

struct HelixRange
{
    std::string chain_id_beg;
    int seq_id_beg;
    std::string chain_id_end;
    int seq_id_end;
    std::string conf_type;
};

class AtomicModelDataBlock
{
    std::string m_map_id, m_model_id;
    std::string m_resolution, m_resolution_method;
    std::unique_ptr<ComponentKeySystem> m_component_key_system;
    std::unique_ptr<AtomKeySystem> m_atom_key_system;
    std::unique_ptr<BondKeySystem> m_bond_key_system;

    std::unordered_map<int, std::vector<std::unique_ptr<AtomObject>>> m_atom_object_list_map; // key: model number
    std::unordered_map<int, std::map<std::tuple<std::string, std::string, std::string, std::string>, AtomObject *>> m_atom_object_in_tuple_map; // key1: model number ; key2: (chain_id, comp_id, seq_id, atom_id)
    std::vector<std::unique_ptr<BondObject>> m_bond_object_list;
    
    std::unordered_map<std::string, Entity> m_entity_type_map; // key: entity_id
    std::unordered_map<std::string, int> m_molecules_size_map; // key: entity_id
    std::unordered_map<std::string, std::vector<std::string>> m_chain_id_list_map; // key: entity_id
    std::unordered_map<Entity, std::vector<std::string>> m_entity_id_list_map; // key: entity type
    
    std::unordered_map<std::string, int> m_struct_sheet_strand_map;
    std::unordered_map<std::string, HelixRange> m_struct_helix_range_map;
    std::unordered_map<std::string, SheetRange> m_struct_sheet_range_map;
    std::vector<Element> m_element_type_list;

    std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> m_chemical_component_entry_map;

public:
    AtomicModelDataBlock();
    ~AtomicModelDataBlock();

    void AddAtomObject(
        int model_number,
        std::unique_ptr<AtomObject> atom_object,
        const std::string & raw_sequence_id_token = "");
    void AddBondObject(std::unique_ptr<BondObject> bond_object);
    void AddEntityTypeInEntityMap(const std::string & entity_id, Entity entity);
    void AddChainIDInEntityMap(const std::string & entity_id, const std::string & chain_id);
    void AddMoleculesSizeInEntityMap(const std::string & entity_id, int molecules_size);
    void AddSheetStrands(const std::string & sheet_id, int strands_size);
    void AddSheetRange(const std::string & composite_sheet_id, const SheetRange & range);
    void AddHelixRange(const std::string & helix_id, const HelixRange & range);
    void AddElementType(const Element & element);
    void AddChemicalComponentEntry(ComponentKey comp_key, std::unique_ptr<ChemicalComponentEntry> entry);
    void AddComponentAtomEntry(
        ComponentKey comp_key,
        AtomKey atom_key,
        const ComponentAtomEntry & atom_entry
    );
    void AddComponentBondEntry(
        ComponentKey comp_key,
        BondKey bond_key,
        const ComponentBondEntry & bond_entry
    );

    void SetPdbID(const std::string & label) { m_model_id = label; }
    void SetEmdID(const std::string & label) { m_map_id = label; }
    void SetResolution(const std::string & value) { m_resolution = value; }
    void SetResolutionMethod(const std::string & value) { m_resolution_method = value; }
    void SetStructureInfo(AtomObject * atom_object);

    const std::string & GetPdbID() const;
    const std::string & GetEmdID() const;
    double GetResolution() const;
    const std::string & GetResolutionMethod() const;
    const std::vector<Element> & GetElementTypeList() const;
    const std::unordered_map<std::string, Entity> & GetEntityTypeMap() const;
    const std::unordered_map<std::string, int> & GetMoleculesSizeMap() const;
    const std::unordered_map<Entity, std::vector<std::string>> & GetEntityIDListMap() const;
    const std::unordered_map<std::string, std::vector<std::string>> & GetChainIDListMap() const;
    const std::unordered_map<int, std::vector<std::unique_ptr<AtomObject>>> & GetAtomObjectMap() const;
    const std::vector<std::unique_ptr<BondObject>> & GetBondObjectList() const;
    std::vector<int> GetModelNumberList() const;
    bool HasModelNumber(int model_number) const;
    ChemicalComponentEntry * GetChemicalComponentEntryPtr(ComponentKey key) const;
    const ComponentBondEntry * GetComponentBondEntryPtr(ComponentKey comp_key, BondKey bond_key) const;
    std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> & GetChemicalComponentEntryMap();
    bool HasChemicalComponentEntry(ComponentKey comp_key) const;
    bool HasComponentBondEntry(ComponentKey comp_key, BondKey bond_key) const;
    AtomObject * GetAtomObjectPtrInTuple(
        int model_number,
        const std::string & chain_id,
        const std::string & comp_id,
        const std::string & seq_id,
        const std::string & atom_id
    ) const;
    AtomObject * GetAtomObjectPtrInAnyModel(
        const std::string & chain_id,
        const std::string & comp_id,
        const std::string & seq_id,
        const std::string & atom_id,
        int * model_number = nullptr
    ) const;
    ComponentKeySystem * GetComponentKeySystemPtr() { return m_component_key_system.get(); }
    AtomKeySystem * GetAtomKeySystemPtr() { return m_atom_key_system.get(); }
    BondKeySystem * GetBondKeySystemPtr() { return m_bond_key_system.get(); }
    std::vector<std::unique_ptr<AtomObject>> MoveAtomObjectList(int model_number=1);
    std::vector<std::unique_ptr<BondObject>> MoveBondObjectList();
    std::unique_ptr<ComponentKeySystem> MoveComponentKeySystem() { return std::move(m_component_key_system); }
    std::unique_ptr<AtomKeySystem> MoveAtomKeySystem() { return std::move(m_atom_key_system); }
    std::unique_ptr<BondKeySystem> MoveBondKeySystem() { return std::move(m_bond_key_system); }

};

} // namespace rhbm_gem
