#pragma once

#include <array>
#include <string>

#include "DataObjectBase.hpp"
#include "GlobalEnumClass.hpp"
#include "ComponentKeySystem.hpp"
#include "AtomKeySystem.hpp"

class AtomObject;
class AtomicPotentialEntry;

class BondObject : public DataObjectBase
{
    std::string m_key_tag;
    bool m_is_selected;
    AtomObject * m_atom_object_1;
    AtomObject * m_atom_object_2;
    std::array<float, 3> m_position, m_bond_vector;
    std::unique_ptr<AtomicPotentialEntry> m_atomic_potential_entry;

public:
    BondObject(void);
    explicit BondObject(AtomObject * atom_object_1, AtomObject * atom_object_2);
    ~BondObject();
    BondObject(const BondObject & other);
    std::unique_ptr<DataObjectBase> Clone(void) const override;
    void Display(void) const override;
    void Update(void) override;
    void Accept(DataObjectVisitorBase * visitor) override;
    void SetKeyTag(const std::string & label) override { m_key_tag = label; }
    std::string GetKeyTag(void) const override { return m_key_tag; }

    std::unique_ptr<BondObject> BondObjectClone(void) const;
    void AddAtomicPotentialEntry(std::unique_ptr<AtomicPotentialEntry> entry);

    std::string GetInfo(void) const;
    std::array<float, 3> GetPosition(void) const { return m_position; }
    std::array<float, 3> GetBondVector(void) const { return m_bond_vector; }
    AtomicPotentialEntry * GetAtomicPotentialEntry(void) const { return m_atomic_potential_entry.get(); }

};
