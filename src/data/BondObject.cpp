#include "BondObject.hpp"
#include "AtomObject.hpp"
#include "AtomicInfoHelper.hpp"
#include "DataObjectVisitorBase.hpp"
#include "AtomicPotentialEntry.hpp"
#include "GlobalEnumClass.hpp"
#include "AtomClassifier.hpp"
#include "Logger.hpp"

#include <stdexcept>

BondObject::BondObject(void) :
    m_key_tag{ "" },
    m_is_selected{ false }, m_bond_key{ 0 },
    m_atom_serial_id_1{ 0 }, m_atom_serial_id_2{ 0 },
    m_atom_object_1{ nullptr }, m_atom_object_2{ nullptr },
    m_position{ 0.0, 0.0, 0.0 }, m_bond_vector{ 0.0, 0.0, 0.0 },
    m_local_potential_entry{ nullptr }
{

}

BondObject::BondObject(AtomObject * atom_object_1, AtomObject * atom_object_2) :
    m_key_tag{ "" },
    m_is_selected{ false }, m_bond_key{ 0 },
    m_atom_serial_id_1{ atom_object_1->GetSerialID() },
    m_atom_serial_id_2{ atom_object_2->GetSerialID() },
    m_atom_object_1{ atom_object_1 }, m_atom_object_2{ atom_object_2 },
    m_position{ 0.0, 0.0, 0.0 }, m_bond_vector{ 0.0, 0.0, 0.0 },
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
}

BondObject::~BondObject()
{

}

BondObject::BondObject(const BondObject & other) :
    m_key_tag{ other.m_key_tag },
    m_is_selected{ other.m_is_selected }, m_bond_key{ other.m_bond_key },
    m_atom_serial_id_1{ other.m_atom_serial_id_1 },
    m_atom_serial_id_2{ other.m_atom_serial_id_2 },
    m_atom_object_1{ other.m_atom_object_1 }, m_atom_object_2{ other.m_atom_object_2 },
    m_position{ other.m_position }, m_bond_vector{ other.m_bond_vector }
{

}

std::unique_ptr<DataObjectBase> BondObject::Clone() const
{
    return std::make_unique<BondObject>(*this);
}

void BondObject::Display(void) const
{
    Logger::Log(LogLevel::Info, "BondObject Display: " + GetInfo());
}

void BondObject::Update(void)
{
    Logger::Log(LogLevel::Info, "BondObject Update: " + GetInfo());
}

void BondObject::Accept(DataObjectVisitorBase * visitor)
{
    visitor->VisitBondObject(this);
}

std::unique_ptr<BondObject> BondObject::BondObjectClone(void) const
{
    auto base_clone{ this->Clone() };
    auto derived_ptr{ dynamic_cast<BondObject*>(base_clone.release()) };
    if (!derived_ptr)
    {
        throw std::runtime_error("Clone() did not return an BondObject instance!");
    }
    return std::unique_ptr<BondObject>(derived_ptr);
}

std::string BondObject::GetInfo(void) const
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

ComponentKey BondObject::GetComponentKey(void) const
{
    return m_atom_object_1->GetComponentKey();
}

AtomKey BondObject::GetAtomKey1(void) const
{
    return m_atom_object_1->GetAtomKey();
}

AtomKey BondObject::GetAtomKey2(void) const
{
    return m_atom_object_2->GetAtomKey();
}

void BondObject::AddLocalPotentialEntry(std::unique_ptr<AtomicPotentialEntry> entry)
{
    m_local_potential_entry = std::move(entry);
}

