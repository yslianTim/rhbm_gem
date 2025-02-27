#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "DataObjectBase.hpp"

class AtomObject;

class ModelObject : public DataObjectBase
{
    std::vector<std::unique_ptr<AtomObject>> m_atom_list;
    std::string m_key_tag, m_pdb_id, m_emd_id;

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

    int GetNumberOfAtom(void) const { return m_atom_list.size(); }
    int GetNumberOfSelectedAtom(void) const;
    const std::vector<std::unique_ptr<AtomObject>> & GetComponentsList(void) const { return m_atom_list; }
    std::string GetPdbID(void) const { return m_pdb_id; }
    std::string GetEmdID(void) const { return m_emd_id; }

private:
    
};