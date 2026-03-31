#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>

#include <stdexcept>

namespace rhbm_gem {

BondObject::BondObject() :
    m_key_tag{ "" },
    m_is_selected{ false }, m_bond_key{ 0 },
    m_bond_type{ BondType::UNK }, m_bond_order{ BondOrder::UNK },
    m_atom_serial_id_1{ 0 }, m_atom_serial_id_2{ 0 },
    m_atom_object_1{ nullptr }, m_atom_object_2{ nullptr },
    m_position{ 0.0, 0.0, 0.0 }, m_bond_vector{ 0.0, 0.0, 0.0 }, m_unit_vector{ 0.0, 0.0, 0.0 },
    m_local_potential_entry{ nullptr }
{

}

BondObject::BondObject(AtomObject * atom_object_1, AtomObject * atom_object_2) :
    m_key_tag{ "" },
    m_is_selected{ false }, m_bond_key{ 0 },
    m_bond_type{ BondType::UNK }, m_bond_order{ BondOrder::UNK },
    m_atom_serial_id_1{ atom_object_1->GetSerialID() },
    m_atom_serial_id_2{ atom_object_2->GetSerialID() },
    m_atom_object_1{ atom_object_1 }, m_atom_object_2{ atom_object_2 },
    m_position{ 0.0, 0.0, 0.0 }, m_bond_vector{ 0.0, 0.0, 0.0 }, m_unit_vector{ 0.0, 0.0, 0.0 },
    m_local_potential_entry{ nullptr }
{
    auto position_1{ m_atom_object_1->GetPosition() };
    auto position_2{ m_atom_object_2->GetPosition() };
    m_position = {
        (position_1[0] + position_2[0]) / 2.0f,
        (position_1[1] + position_2[1]) / 2.0f,
        (position_1[2] + position_2[2]) / 2.0f
    };
    m_bond_vector = {
        position_2[0] - position_1[0],
        position_2[1] - position_1[1],
        position_2[2] - position_1[2]
    };
    auto norm{ ArrayStats<float>::ComputeNorm(m_bond_vector) };
    if (norm > 0.0f)
    {
        m_unit_vector = {
            m_bond_vector[0] / norm,
            m_bond_vector[1] / norm,
            m_bond_vector[2] / norm
        };
    }
}

BondObject::~BondObject()
{

}

BondObject::BondObject(const BondObject & other) :
    m_key_tag{ other.m_key_tag },
    m_is_selected{ other.m_is_selected }, m_bond_key{ other.m_bond_key },
    m_bond_type{ other.m_bond_type }, m_bond_order{ other.m_bond_order },
    m_atom_serial_id_1{ other.m_atom_serial_id_1 },
    m_atom_serial_id_2{ other.m_atom_serial_id_2 },
    m_atom_object_1{ other.m_atom_object_1 }, m_atom_object_2{ other.m_atom_object_2 },
    m_position{ other.m_position }, m_bond_vector{ other.m_bond_vector },
    m_unit_vector{ other.m_unit_vector },
    m_local_potential_entry{ nullptr }
{

}

std::unique_ptr<DataObjectBase> BondObject::Clone() const
{
    return std::make_unique<BondObject>(*this);
}

std::unique_ptr<BondObject> BondObject::BondObjectClone() const
{
    auto base_clone{ this->Clone() };
    auto derived_ptr{ dynamic_cast<BondObject*>(base_clone.release()) };
    if (!derived_ptr)
    {
        throw std::runtime_error("Clone() did not return an BondObject instance!");
    }
    return std::unique_ptr<BondObject>(derived_ptr);
}

std::string BondObject::GetInfo() const
{
    return "[Component ID1] " + m_atom_object_1->GetComponentID() + " " +
           "[Component ID2] " + m_atom_object_2->GetComponentID() + " " +
           "[Atom ID1] " + m_atom_object_1->GetAtomID() + " " +
           "[Atom ID2] " + m_atom_object_2->GetAtomID() + " " +
           "Position = (" +
           std::to_string(m_position.at(0)) + ", " +
           std::to_string(m_position.at(1)) + ", " +
           std::to_string(m_position.at(2)) + ") " +
           "Axis Vector = (" +
           std::to_string(m_bond_vector.at(0)) + ", " +
           std::to_string(m_bond_vector.at(1)) + ", " +
           std::to_string(m_bond_vector.at(2)) + ")";
}

ComponentKey BondObject::GetComponentKey() const
{
    return m_atom_object_1->GetComponentKey();
}

AtomKey BondObject::GetAtomKey1() const
{
    return m_atom_object_1->GetAtomKey();
}

AtomKey BondObject::GetAtomKey2() const
{
    return m_atom_object_2->GetAtomKey();
}

void BondObject::AddLocalPotentialEntry(std::unique_ptr<LocalPotentialEntry> entry)
{
    m_local_potential_entry = std::move(entry);
}

float BondObject::GetBondLength() const
{
    if (m_atom_object_1 == nullptr || m_atom_object_2 == nullptr) return 0.0f;
    return ArrayStats<float>::ComputeNorm(
        m_atom_object_1->GetPositionRef(), m_atom_object_2->GetPositionRef());
}

} // namespace rhbm_gem
