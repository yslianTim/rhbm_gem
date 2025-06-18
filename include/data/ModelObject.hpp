#pragma once

#include <memory>
#include <string>
#include <vector>
#include <array>
#include <tuple>
#include <unordered_map>
#include "DataObjectBase.hpp"
#include "KDTreeAlgorithm.hpp"

class AtomObject;
class GroupPotentialEntry;
class GroupChargeEntry;
class ModelObject : public DataObjectBase
{
    std::vector<std::unique_ptr<AtomObject>> m_atom_list;
    std::vector<AtomObject *> m_selected_atom_list;
    std::string m_key_tag, m_pdb_id, m_emd_id;
    std::string m_resolution_method;
    double m_resolution;
    std::unordered_map<std::string, std::vector<std::string>> m_chain_id_list_map; // key : entity_id
    std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> m_group_potential_entry_map;
    std::unordered_map<std::string, std::unique_ptr<GroupChargeEntry>> m_group_charge_entry_map;
    std::unique_ptr<KDNode<AtomObject>> m_kd_tree_root;
    std::unique_ptr<std::array<float, 3>> m_center_of_mass_position;
    std::unique_ptr<std::tuple<double, double>> m_model_position_range[3];

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

    void AddAtom(std::unique_ptr<AtomObject> component);
    void SetPdbID(const std::string & label) { m_pdb_id = label; }
    void SetEmdID(const std::string & label) { m_emd_id = label; }
    void SetResolution(double value) { m_resolution = value; }
    void SetResolutionMethod(const std::string & value) { m_resolution_method = value; }
    void SetChainIDListMap(const std::unordered_map<std::string, std::vector<std::string>> & value) { m_chain_id_list_map = value; }
    void AddGroupPotentialEntry(const std::string & class_key, std::unique_ptr<GroupPotentialEntry> & entry);
    void AddGroupChargeEntry(const std::string & class_key, std::unique_ptr<GroupChargeEntry> & entry);
    void BuildKDTreeRoot(void);
    void FilterAtomFromSymmetry(bool is_asymmetry);

    size_t GetNumberOfAtom(void) const { return m_atom_list.size(); }
    size_t GetNumberOfSelectedAtom(void) const { return m_selected_atom_list.size(); }
    const std::vector<std::unique_ptr<AtomObject>> & GetComponentsList(void) const { return m_atom_list; }
    const std::vector<AtomObject *> & GetSelectedAtomList(void) const { return m_selected_atom_list; }
    std::string GetPdbID(void) const { return m_pdb_id; }
    std::string GetEmdID(void) const { return m_emd_id; }
    double GetResolution(void) const { return m_resolution; }
    std::string GetResolutionMethod(void) const { return m_resolution_method; }
    std::array<float, 3> GetCenterOfMassPosition(void);
    std::tuple<double, double> GetModelPositionRange(int axis);
    double GetModelPosition(int axis, double normalized_pos);
    double GetModelLength(int axis);
    GroupPotentialEntry * GetGroupPotentialEntry(const std::string & class_key) const;
    GroupChargeEntry * GetGroupChargeEntry(const std::string & class_key) const;
    const std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> & GetGroupPotentialEntryMap(void) const;
    const std::unordered_map<std::string, std::unique_ptr<GroupChargeEntry>> & GetGroupChargeEntryMap(void) const;
    KDNode<AtomObject> * GetKDTreeRoot(void) const { return m_kd_tree_root.get(); }

private:
    void BuildSelectedAtomList(void);
    
};
