#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <tuple>
#include <unordered_map>
#include "DataObjectBase.hpp"
#include "KDTreeAlgorithm.hpp"

class AtomObject;
class GroupPotentialEntry;
class ModelObject : public DataObjectBase
{
    std::vector<std::unique_ptr<AtomObject>> m_atom_list;
    std::string m_key_tag, m_pdb_id, m_emd_id;
    std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> m_group_potential_entry_map;
    std::unique_ptr<KDNode<AtomObject>> m_kd_tree_root;

public:
    ModelObject(void);
    explicit ModelObject(std::vector<std::unique_ptr<AtomObject>> atom_object_list);
    ~ModelObject();
    ModelObject(const ModelObject & other);
    std::unique_ptr<DataObjectBase> Clone(void) const override;
    void Display(void) const override;
    void Accept(DataObjectVisitorBase * visitor) override;
    void SetKeyTag(const std::string & label) override { m_key_tag = label; }
    std::string GetKeyTag(void) const override { return m_key_tag; }

    void AddAtom(std::unique_ptr<AtomObject> component);
    void SetPdbID(const std::string & label) { m_pdb_id = label; }
    void SetEmdID(const std::string & label) { m_emd_id = label; }
    void AddGroupPotentialEntry(const std::string & class_key, std::unique_ptr<GroupPotentialEntry> & entry);
    void BuildKDTreeRoot(void);

    int GetNumberOfAtom(void) const { return static_cast<int>(m_atom_list.size()); }
    size_t GetNumberOfSelectedAtom(void) const;
    const std::vector<std::unique_ptr<AtomObject>> & GetComponentsList(void) const { return m_atom_list; }
    std::string GetPdbID(void) const { return m_pdb_id; }
    std::string GetEmdID(void) const { return m_emd_id; }
    std::tuple<double, double> GetModelPositionRange(int axis, double margin=0.0) const;
    GroupPotentialEntry * GetGroupPotentialEntry(const std::string & class_key) const;
    const std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> & GetGroupPotentialEntryMap(void) const;
    KDNode<AtomObject> * GetKDTreeRoot(void) const { return m_kd_tree_root.get(); }

private:
    
};