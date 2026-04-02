#pragma once

#include <array>
#include <string>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/BondKeySystem.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class ModelObject;
class ModelObjectBuilder;
ModelObject * OwnerModelOf(const BondObject & bond_object);

class BondObject
{
    bool m_is_selected;
    bool m_is_special_bond;
    BondKey m_bond_key;
    BondType m_bond_type;
    BondOrder m_bond_order;
    int m_atom_serial_id_1;
    int m_atom_serial_id_2;
    AtomObject * m_atom_object_1;
    AtomObject * m_atom_object_2;
    ModelObject * m_owner_model;
    std::array<float, 3> m_position, m_bond_vector, m_unit_vector;

public:
    BondObject();
    explicit BondObject(AtomObject * atom_object_1, AtomObject * atom_object_2);
    ~BondObject();
    BondObject(const BondObject & other);

    void SetSpecialBondFlag(bool value) { m_is_special_bond = value; }
    void SetBondKey(BondKey value) { m_bond_key = value; }
    void SetBondType(BondType value) { m_bond_type = value; }
    void SetBondOrder(BondOrder value) { m_bond_order = value; }
    std::string GetInfo() const;
    int GetAtomSerialID1() const { return m_atom_serial_id_1; }
    int GetAtomSerialID2() const { return m_atom_serial_id_2; }
    BondKey GetBondKey() const { return m_bond_key; }
    BondType GetBondType() const { return m_bond_type; }
    BondOrder GetBondOrder() const { return m_bond_order; }
    bool GetSpecialBondFlag() const { return m_is_special_bond; }
    std::array<float, 3> GetPosition() const { return m_position; }
    std::array<float, 3> GetBondVector() const { return m_bond_vector; }
    AtomObject * GetAtomObject1() const { return m_atom_object_1; }
    AtomObject * GetAtomObject2() const { return m_atom_object_2; }

private:
    friend class ModelObject;
    friend class ModelObjectBuilder;
    friend ModelObject * OwnerModelOf(const BondObject & bond_object);

    void SetSelectedFlag(bool value) { m_is_selected = value; }
    void SetOwnerModel(ModelObject * value) { m_owner_model = value; }
};

} // namespace rhbm_gem
