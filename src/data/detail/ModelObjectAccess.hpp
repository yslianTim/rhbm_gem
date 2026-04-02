#pragma once

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
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

    static void AddAtom(
        ModelObject & model_object,
        std::unique_ptr<AtomObject> atom)
    {
        model_object.AddAtom(std::move(atom));
    }

    static void AddBond(
        ModelObject & model_object,
        std::unique_ptr<BondObject> bond)
    {
        model_object.AddBond(std::move(bond));
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

    static void FinalizeLoad(ModelObject & model_object)
    {
        model_object.FinalizeLoad();
    }

    static void ApplySymmetrySelection(ModelObject & model_object, bool is_asymmetry)
    {
        model_object.ApplySymmetrySelection(is_asymmetry);
    }

    static void RebuildSelectionIndex(ModelObject & model_object)
    {
        model_object.RebuildSelectionIndex();
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

    static const std::vector<AtomObject *> & SelectedAtoms(const ModelObject & model_object)
    {
        return model_object.m_selected_atom_list;
    }

    static const std::vector<BondObject *> & SelectedBonds(const ModelObject & model_object)
    {
        return model_object.m_selected_bond_list;
    }

    static size_t SelectedAtomCount(const ModelObject & model_object)
    {
        return model_object.m_selected_atom_list.size();
    }

    static size_t SelectedBondCount(const ModelObject & model_object)
    {
        return model_object.m_selected_bond_list.size();
    }

    static void ClearAnalysisFitStates(ModelObject & model_object);
};

} // namespace rhbm_gem
