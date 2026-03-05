#pragma once

#include <array>
#include <string>

#include "DataObjectBase.hpp"
#include "GlobalEnumClass.hpp"
#include "ComponentKeySystem.hpp"
#include "AtomKeySystem.hpp"
#include "BondKeySystem.hpp"

namespace rhbm_gem {

class AtomObject;
class LocalPotentialEntry;

class BondObject : public DataObjectBase
{
    std::string m_key_tag;
    bool m_is_selected;
    bool m_is_special_bond;
    BondKey m_bond_key;
    BondType m_bond_type;
    BondOrder m_bond_order;
    int m_atom_serial_id_1;
    int m_atom_serial_id_2;
    AtomObject * m_atom_object_1;
    AtomObject * m_atom_object_2;
    std::array<float, 3> m_position, m_bond_vector, m_unit_vector;
    std::unique_ptr<LocalPotentialEntry> m_local_potential_entry;

public:
    BondObject();
    explicit BondObject(AtomObject * atom_object_1, AtomObject * atom_object_2);
    ~BondObject();
    BondObject(const BondObject & other);
    std::unique_ptr<DataObjectBase> Clone() const override;
    void Display() const override;
    void Update() override;
    void Accept(DataObjectVisitor & visitor) override;
    void Accept(ConstDataObjectVisitor & visitor) const override;
    void Accept(DataObjectVisitor & visitor, ModelVisitMode model_mode) override;
    void Accept(ConstDataObjectVisitor & visitor, ModelVisitMode model_mode) const override;
    void SetKeyTag(const std::string & label) override { m_key_tag = label; }
    std::string GetKeyTag() const override { return m_key_tag; }

    std::unique_ptr<BondObject> BondObjectClone() const;
    void SetSelectedFlag(bool value) { m_is_selected = value; }
    void SetSpecialBondFlag(bool value) { m_is_special_bond = value; }
    void SetBondKey(BondKey value) { m_bond_key = value; }
    void SetBondType(BondType value) { m_bond_type = value; }
    void SetBondOrder(BondOrder value) { m_bond_order = value; }
    void AddLocalPotentialEntry(std::unique_ptr<LocalPotentialEntry> entry);

    std::string GetInfo() const;
    int GetAtomSerialID1() const { return m_atom_serial_id_1; }
    int GetAtomSerialID2() const { return m_atom_serial_id_2; }
    float GetBondLength() const;
    ComponentKey GetComponentKey() const;
    AtomKey GetAtomKey1() const;
    AtomKey GetAtomKey2() const;
    BondKey GetBondKey() const { return m_bond_key; }
    BondType GetBondType() const { return m_bond_type; }
    BondOrder GetBondOrder() const { return m_bond_order; }
    bool GetSelectedFlag() const { return m_is_selected; }
    bool GetSpecialBondFlag() const { return m_is_special_bond; }
    std::array<float, 3> GetPosition() const { return m_position; }
    std::array<float, 3> GetBondVector() const { return m_bond_vector; }
    std::array<float, 3> GetUnitVector() const { return m_unit_vector; }
    AtomObject * GetAtomObject1() const { return m_atom_object_1; }
    AtomObject * GetAtomObject2() const { return m_atom_object_2; }
    LocalPotentialEntry * GetLocalPotentialEntry() const { return m_local_potential_entry.get(); }

};

} // namespace rhbm_gem
