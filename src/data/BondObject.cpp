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
    m_is_selected{ false },
    m_atom_object_1{ nullptr }, m_atom_object_2{ nullptr },
    m_position{ 0.0, 0.0, 0.0 }, m_bond_vector{ 0.0, 0.0, 0.0 },
    m_atomic_potential_entry{ nullptr }
{

}

BondObject::BondObject(AtomObject * atom_object_1, AtomObject * atom_object_2) :
    m_key_tag{ "" },
    m_is_selected{ false },
    m_atom_object_1{ atom_object_1 }, m_atom_object_2{ atom_object_2 },
    m_position{ 0.0, 0.0, 0.0 }, m_bond_vector{ 0.0, 0.0, 0.0 },
    m_atomic_potential_entry{ nullptr }
{

}

BondObject::~BondObject()
{

}

BondObject::BondObject(const BondObject & other) :
    m_key_tag{ other.m_key_tag },
    m_is_selected{ other.m_is_selected },
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

void BondObject::AddAtomicPotentialEntry(std::unique_ptr<AtomicPotentialEntry> entry)
{
    m_atomic_potential_entry = std::move(entry);
}

