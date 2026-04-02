#pragma once

#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "data/detail/AtomObjectAccess.hpp"
#include "data/detail/BondObjectAccess.hpp"
#include "data/detail/ModelAnalysisState.hpp"
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/BondKeySystem.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>

namespace rhbm_gem {

class ModelAnalysisState;
template <typename T> struct KDNode;

class ModelObjectAccess
{
public:
    static ModelAnalysisState & AnalysisState(ModelObject & model_object)
    {
        return *model_object.m_analysis_state;
    }

    static const ModelAnalysisState & AnalysisState(const ModelObject & model_object)
    {
        return *model_object.m_analysis_state;
    }

    static ::KDNode<AtomObject> * KDTreeRoot(const ModelObject & model_object)
    {
        return model_object.m_kd_tree_root.get();
    }

    static void BuildKDTreeRoot(ModelObject & model_object)
    {
        model_object.BuildKDTreeRoot();
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

    static void SetAtomList(
        ModelObject & model_object,
        std::vector<std::unique_ptr<AtomObject>> atom_list)
    {
        model_object.m_atom_list = std::move(atom_list);
        model_object.m_analysis_state->Clear();
        model_object.SyncDerivedState();
    }

    static void SetBondList(
        ModelObject & model_object,
        std::vector<std::unique_ptr<BondObject>> bond_list)
    {
        model_object.m_bond_list = std::move(bond_list);
        model_object.m_analysis_state->Clear();
        model_object.SyncDerivedState();
    }

    static void SetChainIDListMap(
        ModelObject & model_object,
        const std::unordered_map<std::string, std::vector<std::string>> & value)
    {
        model_object.m_chain_id_list_map = value;
    }

    static void AddAtom(
        ModelObject & model_object,
        std::unique_ptr<AtomObject> atom)
    {
        model_object.m_atom_list.emplace_back(std::move(atom));
        model_object.m_analysis_state->Clear();
        model_object.SyncDerivedState();
    }

    static void AddBond(
        ModelObject & model_object,
        std::unique_ptr<BondObject> bond)
    {
        model_object.m_bond_list.emplace_back(std::move(bond));
        model_object.m_analysis_state->Clear();
        model_object.SyncDerivedState();
    }

    static void AddChemicalComponentEntry(
        ModelObject & model_object,
        ComponentKey component_key,
        std::unique_ptr<ChemicalComponentEntry> entry)
    {
        model_object.m_chemical_component_entry_map[component_key] = std::move(entry);
    }

    static void SetChemicalComponentEntryMap(
        ModelObject & model_object,
        std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> entry_map)
    {
        model_object.m_chemical_component_entry_map = std::move(entry_map);
    }

    static void SetComponentKeySystem(
        ModelObject & model_object,
        std::unique_ptr<ComponentKeySystem> component_key_system)
    {
        model_object.m_component_key_system = std::move(component_key_system);
    }

    static void SetAtomKeySystem(
        ModelObject & model_object,
        std::unique_ptr<AtomKeySystem> atom_key_system)
    {
        model_object.m_atom_key_system = std::move(atom_key_system);
    }

    static void SetBondKeySystem(
        ModelObject & model_object,
        std::unique_ptr<BondKeySystem> bond_key_system)
    {
        model_object.m_bond_key_system = std::move(bond_key_system);
    }

    static void FinalizeLoad(ModelObject & model_object)
    {
        model_object.SyncDerivedState();
    }

    static void ApplySymmetrySelection(ModelObject & model_object, bool is_asymmetry)
    {
        model_object.FilterSelectionFromSymmetry(is_asymmetry);
        model_object.SyncDerivedState();
    }

    static void RebuildSelectionIndex(ModelObject & model_object)
    {
        model_object.SyncDerivedState();
    }

    static AtomObject * FindAtomPtr(const ModelObject & model_object, int serial_id)
    {
        return model_object.GetAtomPtr(serial_id);
    }

    static std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> &
    ChemicalComponentEntryMap(ModelObject & model_object)
    {
        return model_object.m_chemical_component_entry_map;
    }

    static void ClearAnalysisFitStates(ModelObject & model_object);
};

} // namespace rhbm_gem
