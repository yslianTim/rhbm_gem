#pragma once

#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/BondKeySystem.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>

namespace rhbm_gem {

class ModelObject;
class AtomObject;
class BondObject;
class ChemicalComponentEntry;
class ModelObjectBuilder;
class ModelAnalysisData;
class ModelObjectStorage;
class ModelSelectionAccess;
struct ModelSpatialCache;
ModelAnalysisData & MutableAnalysisData(ModelObject & model_object);
const ModelAnalysisData & ReadAnalysisData(const ModelObject & model_object);
void ClearAnalysisFitStates(ModelObject & model_object);
std::vector<AtomObject *> FindAtomsInRange(
    ModelObject & model_object,
    const AtomObject & center_atom,
    double range);

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
    std::unique_ptr<ModelSpatialCache> m_spatial_cache;
    std::unique_ptr<std::array<float, 3>> m_center_of_mass_position;
    std::unique_ptr<std::tuple<double, double>> m_model_position_range[3];
    std::unique_ptr<ComponentKeySystem> m_component_key_system;
    std::unique_ptr<AtomKeySystem> m_atom_key_system;
    std::unique_ptr<BondKeySystem> m_bond_key_system;
    std::unique_ptr<ModelAnalysisData> m_analysis_data;

public:
    ModelObject();
    explicit ModelObject(std::vector<std::unique_ptr<AtomObject>> atom_object_list);
    ~ModelObject();
    ModelObject(const ModelObject & other);
    ModelObject(ModelObject && other) noexcept;
    ModelObject & operator=(ModelObject && other) noexcept;
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
    const std::vector<AtomObject *> & GetSelectedAtoms() const { return m_selected_atom_list; }
    const std::vector<BondObject *> & GetSelectedBonds() const { return m_selected_bond_list; }
    size_t GetSelectedAtomCount() const { return m_selected_atom_list.size(); }
    size_t GetSelectedBondCount() const { return m_selected_bond_list.size(); }
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
    AtomObject * FindAtomPtr(int serial_id) const;
    std::string FindComponentID(ComponentKey component_key) const;
    std::string FindAtomID(AtomKey atom_key) const;
    std::string FindBondID(BondKey bond_key) const;
    bool HasChemicalComponentEntry(ComponentKey component_key) const;
    const ChemicalComponentEntry * FindChemicalComponentEntry(ComponentKey component_key) const;
    std::vector<ComponentKey> GetComponentKeyList() const;
    void SelectAllAtoms(bool selected = true);
    void SelectAllBonds(bool selected = true);
    void SelectAtoms(const std::function<bool(const AtomObject &)> & predicate);
    void SelectBonds(const std::function<bool(const BondObject &)> & predicate);
    void SetAtomSelected(int serial_id, bool selected);
    void SetBondSelected(int atom_serial_id_1, int atom_serial_id_2, bool selected);
    void RebuildSelection();
    void ApplySymmetrySelection(bool is_asymmetry);

private:
    friend class ModelObjectBuilder;
    friend class ModelObjectStorage;
    friend class ModelSelectionAccess;
    friend ModelAnalysisData & MutableAnalysisData(ModelObject & model_object);
    friend const ModelAnalysisData & ReadAnalysisData(const ModelObject & model_object);
    friend void ClearAnalysisFitStates(ModelObject & model_object);
    friend std::vector<AtomObject *> FindAtomsInRange(
        ModelObject & model_object,
        const AtomObject & center_atom,
        double range);

    void RebuildSelectionState();
    void RebuildObjectIndex();
    void AttachOwnedObjects();
    void InvalidateDerivedCaches();
    void SyncDerivedState();
    void EnsureKDTreeRoot();
    void SetAtomSelectionBulk(const std::unordered_set<int> & selected_serial_ids);
    void SetBondSelectionBulk(const std::set<std::pair<int, int>> & selected_serial_pairs);
    void BuildSelectedAtomList();
    void BuildSelectedBondList();
    void FilterSelectionFromSymmetry(bool is_asymmetry);
};

} // namespace rhbm_gem
