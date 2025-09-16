#pragma once

#include <array>
#include <string>

#include "DataObjectBase.hpp"
#include "GlobalEnumClass.hpp"
#include "ComponentKeySystem.hpp"
#include "AtomKeySystem.hpp"
#include "BondKeySystem.hpp"

class AtomObject;
class AtomicPotentialEntry;

class BondObject : public DataObjectBase
{
    std::string m_key_tag;
    bool m_is_selected;
    bool m_is_special_bond;
    BondKey m_bond_key;
    int m_atom_serial_id_1;
    int m_atom_serial_id_2;
    AtomObject * m_atom_object_1;
    AtomObject * m_atom_object_2;
    std::array<float, 3> m_position, m_bond_vector;
    std::unique_ptr<AtomicPotentialEntry> m_local_potential_entry;

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
    void SetSelectedFlag(bool value) { m_is_selected = value; }
    void SetSpecialBondFlag(bool value) { m_is_special_bond = value; }
    void SetBondKey(BondKey value) { m_bond_key = value; }
    void AddLocalPotentialEntry(std::unique_ptr<AtomicPotentialEntry> entry);

    std::string GetInfo(void) const;
    int GetAtomSerialID1(void) const { return m_atom_serial_id_1; }
    int GetAtomSerialID2(void) const { return m_atom_serial_id_2; }
    ComponentKey GetComponentKey(void) const;
    AtomKey GetAtomKey1(void) const;
    AtomKey GetAtomKey2(void) const;
    BondKey GetBondKey(void) const { return m_bond_key; }
    bool GetSelectedFlag(void) const { return m_is_selected; }
    bool GetSpecialBondFlag(void) const { return m_is_special_bond; }
    std::array<float, 3> GetPosition(void) const { return m_position; }
    std::array<float, 3> GetBondVector(void) const { return m_bond_vector; }
    AtomObject * GetAtomObject1(void) const { return m_atom_object_1; }
    AtomObject * GetAtomObject2(void) const { return m_atom_object_2; }
    AtomicPotentialEntry * GetLocalPotentialEntry(void) const { return m_local_potential_entry.get(); }

};
