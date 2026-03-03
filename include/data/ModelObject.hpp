#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "AtomKeySystem.hpp"
#include "BondKeySystem.hpp"
#include "ComponentKeySystem.hpp"
#include "DataObjectBase.hpp"

template <typename T> struct KDNode;

namespace rhbm_gem {

class AtomObject;
class BondObject;
class ChemicalComponentEntry;
class GroupPotentialEntry;

class ModelObject : public DataObjectBase
{
    std::vector<std::unique_ptr<AtomObject>> m_atom_list;
    std::vector<std::unique_ptr<BondObject>> m_bond_list;
    std::vector<AtomObject *> m_selected_atom_list;
    std::vector<BondObject *> m_selected_bond_list;
    std::string m_key_tag, m_pdb_id, m_emd_id;
    std::string m_resolution_method;
    double m_resolution;
    std::map<int, AtomObject*> m_serial_id_atom_map;
    std::unordered_map<std::string, std::vector<std::string>> m_chain_id_list_map;
    std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> m_chemical_component_entry_map;
    std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> m_atom_group_potential_entry_map;
    std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> m_bond_group_potential_entry_map;
    std::unique_ptr<::KDNode<AtomObject>> m_kd_tree_root;
    std::unique_ptr<std::array<float, 3>> m_center_of_mass_position;
    std::unique_ptr<std::tuple<double, double>> m_model_position_range[3];
    std::unique_ptr<ComponentKeySystem> m_component_key_system;
    std::unique_ptr<AtomKeySystem> m_atom_key_system;
    std::unique_ptr<BondKeySystem> m_bond_key_system;

public:
    ModelObject(void);
    explicit ModelObject(std::vector<std::unique_ptr<AtomObject>> atom_object_list);
    ~ModelObject();
    ModelObject(const ModelObject & other);
    std::unique_ptr<DataObjectBase> Clone(void) const override;
    void Display(void) const override;
    void Update(void) override;
    void Accept(DataObjectVisitorBase * visitor) override;
    void SetKeyTag(const std::string & label) override { m_key_tag = label; }
    std::string GetKeyTag(void) const override { return m_key_tag; }

    void AddAtom(std::unique_ptr<AtomObject> atom);
    void AddBond(std::unique_ptr<BondObject> bond);
    void SetAtomList(std::vector<std::unique_ptr<AtomObject>> atom_list);
    void SetBondList(std::vector<std::unique_ptr<BondObject>> bond_list);
    void SetPdbID(const std::string & label) { m_pdb_id = label; }
    void SetEmdID(const std::string & label) { m_emd_id = label; }
    void SetResolution(double value) { m_resolution = value; }
    void SetResolutionMethod(const std::string & value) { m_resolution_method = value; }
    void SetChainIDListMap(
        const std::unordered_map<std::string, std::vector<std::string>> & value) { m_chain_id_list_map = value; }
    void AddAtomGroupPotentialEntry(
        const std::string & class_key,
        std::unique_ptr<GroupPotentialEntry> & entry);
    void AddBondGroupPotentialEntry(
        const std::string & class_key,
        std::unique_ptr<GroupPotentialEntry> & entry);
    void AddChemicalComponentEntry(
        ComponentKey component_key,
        std::unique_ptr<ChemicalComponentEntry> entry);
    void SetChemicalComponentEntryMap(
        std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> & entry_map);
    void SetComponentKeySystem(std::unique_ptr<ComponentKeySystem> component_key_system);
    void SetAtomKeySystem(std::unique_ptr<AtomKeySystem> atom_key_system);
    void SetBondKeySystem(std::unique_ptr<BondKeySystem> bond_key_system);
    void BuildKDTreeRoot(void);
    void FilterAtomFromSymmetry(bool is_asymmetry);
    void FilterBondFromSymmetry(bool is_asymmetry);

    bool HasStandardRNAComponent(void) const;
    bool HasStandardDNAComponent(void) const;
    size_t GetNumberOfAtom(void) const { return m_atom_list.size(); }
    size_t GetNumberOfBond(void) const { return m_bond_list.size(); }
    size_t GetNumberOfSelectedAtom(void) const { return m_selected_atom_list.size(); }
    size_t GetNumberOfSelectedBond(void) const { return m_selected_bond_list.size(); }
    const std::vector<std::unique_ptr<AtomObject>> & GetAtomList(void) const { return m_atom_list; }
    const std::vector<std::unique_ptr<BondObject>> & GetBondList(void) const { return m_bond_list; }
    const std::vector<AtomObject *> & GetSelectedAtomList(void) const { return m_selected_atom_list; }
    const std::vector<BondObject *> & GetSelectedBondList(void) const { return m_selected_bond_list; }
    std::string GetPdbID(void) const { return m_pdb_id; }
    std::string GetEmdID(void) const { return m_emd_id; }
    double GetResolution(void) const { return m_resolution; }
    std::string GetResolutionMethod(void) const { return m_resolution_method; }
    std::array<float, 3> GetCenterOfMassPosition(void);
    std::tuple<double, double> GetModelPositionRange(int axis);
    double GetModelPosition(int axis, double normalized_pos);
    double GetModelLength(int axis);
    AtomObject * GetAtomPtr(int serial_id) const { return m_serial_id_atom_map.at(serial_id); }
    GroupPotentialEntry * GetAtomGroupPotentialEntry(const std::string & class_key) const;
    GroupPotentialEntry * GetBondGroupPotentialEntry(const std::string & class_key) const;
    ChemicalComponentEntry * GetChemicalComponentEntry(ComponentKey key) const;
    const std::map<int, AtomObject *> & GetSerialIDAtomMap(void) const;
    const std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> &
    GetAtomGroupPotentialEntryMap(void) const;
    const std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> &
    GetBondGroupPotentialEntryMap(void) const;
    const std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> &
    GetChemicalComponentEntryMap(void) const;
    ::KDNode<AtomObject> * GetKDTreeRoot(void) const { return m_kd_tree_root.get(); }
    ComponentKeySystem * GetComponentKeySystemPtr(void) { return m_component_key_system.get(); }
    AtomKeySystem * GetAtomKeySystemPtr(void) { return m_atom_key_system.get(); }
    BondKeySystem * GetBondKeySystemPtr(void) { return m_bond_key_system.get(); }
    std::vector<ComponentKey> GetComponentKeyList(void) const;

private:
    void BuildSelectedAtomList(void);
    void BuildSelectedBondList(void);
};

} // namespace rhbm_gem
