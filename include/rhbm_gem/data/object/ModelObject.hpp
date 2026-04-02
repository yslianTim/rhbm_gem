#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/BondKeySystem.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>

template <typename T> struct KDNode;

namespace rhbm_gem {

class AtomObject;
class BondObject;
class ChemicalComponentEntry;
class ModelObjectAccess;
class ModelSelectionView;
class ModelAnalysisState;

class ModelObject
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
    std::unique_ptr<::KDNode<AtomObject>> m_kd_tree_root;
    std::unique_ptr<std::array<float, 3>> m_center_of_mass_position;
    std::unique_ptr<std::tuple<double, double>> m_model_position_range[3];
    std::unique_ptr<ComponentKeySystem> m_component_key_system;
    std::unique_ptr<AtomKeySystem> m_atom_key_system;
    std::unique_ptr<BondKeySystem> m_bond_key_system;
    std::unique_ptr<ModelAnalysisState> m_analysis_state;

public:
    ModelObject();
    explicit ModelObject(std::vector<std::unique_ptr<AtomObject>> atom_object_list);
    ~ModelObject();
    ModelObject(const ModelObject & other);
    void SetKeyTag(const std::string & label) { m_key_tag = label; }
    std::string GetKeyTag() const { return m_key_tag; }
    void SetPdbID(const std::string & label) { m_pdb_id = label; }
    void SetEmdID(const std::string & label) { m_emd_id = label; }
    void SetResolution(double value) { m_resolution = value; }
    void SetResolutionMethod(const std::string & value) { m_resolution_method = value; }

    bool HasStandardRNAComponent() const;
    bool HasStandardDNAComponent() const;
    size_t GetNumberOfAtom() const { return m_atom_list.size(); }
    size_t GetNumberOfBond() const { return m_bond_list.size(); }
    const std::vector<std::unique_ptr<AtomObject>> & GetAtomList() const { return m_atom_list; }
    const std::vector<std::unique_ptr<BondObject>> & GetBondList() const { return m_bond_list; }
    std::string GetPdbID() const { return m_pdb_id; }
    std::string GetEmdID() const { return m_emd_id; }
    double GetResolution() const { return m_resolution; }
    std::string GetResolutionMethod() const { return m_resolution_method; }
    const std::unordered_map<std::string, std::vector<std::string>> &
    GetChainIDListMap() const { return m_chain_id_list_map; }
    std::array<float, 3> GetCenterOfMassPosition();
    std::tuple<double, double> GetModelPositionRange(int axis);
    double GetModelPosition(int axis, double normalized_pos);
    double GetModelLength(int axis);
    const std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> &
    GetChemicalComponentEntryMap() const;
    std::vector<ComponentKey> GetComponentKeyList() const;

private:
    friend class ModelObjectAccess;
    friend class ModelSelectionView;
    AtomObject * GetAtomPtr(int serial_id) const { return m_serial_id_atom_map.at(serial_id); }

    void AddAtom(std::unique_ptr<AtomObject> atom);
    void AddBond(std::unique_ptr<BondObject> bond);
    void SetAtomList(std::vector<std::unique_ptr<AtomObject>> atom_list);
    void SetBondList(std::vector<std::unique_ptr<BondObject>> bond_list);
    void SetChainIDListMap(
        const std::unordered_map<std::string, std::vector<std::string>> & value) { m_chain_id_list_map = value; }
    void AddChemicalComponentEntry(
        ComponentKey component_key,
        std::unique_ptr<ChemicalComponentEntry> entry);
    void SetChemicalComponentEntryMap(
        std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> entry_map);
    void SetComponentKeySystem(std::unique_ptr<ComponentKeySystem> component_key_system);
    void SetAtomKeySystem(std::unique_ptr<AtomKeySystem> atom_key_system);
    void SetBondKeySystem(std::unique_ptr<BondKeySystem> bond_key_system);
    void BuildKDTreeRoot();
    void RebuildSelectionIndex();
    void FinalizeLoad();
    void ApplySymmetrySelection(bool is_asymmetry);
    void RebuildSelectionState();
    void InvalidateDerivedCaches();
    void SyncDerivedState();
    void BuildSelectedAtomList();
    void BuildSelectedBondList();
    void FilterSelectionFromSymmetry(bool is_asymmetry);
};

} // namespace rhbm_gem
