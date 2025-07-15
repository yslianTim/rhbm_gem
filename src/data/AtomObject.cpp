#include "AtomObject.hpp"
#include "AtomicInfoHelper.hpp"
#include "DataObjectVisitorBase.hpp"
#include "AtomicPotentialEntry.hpp"
#include "GlobalEnumClass.hpp"
#include "AtomClassifier.hpp"
#include "Logger.hpp"

#include <stdexcept>

AtomObject::AtomObject(void) :
    m_key_tag{ "" },
    m_is_selected{ false }, m_is_special_atom{ false },
    m_serial_id{ 0 }, m_residue_id{ 0 },
    m_chain_id{ "" }, m_indicator{ "" },
    m_occupancy{ 0.0 }, m_temperature{ 0.0 },
    m_residue{ Residue::UNK }, m_element{ Element::UNK },
    m_remoteness{ Remoteness::UNK }, m_branch{ Branch::UNK },
    m_structure{ Structure::UNK },
    m_position{ 0.0, 0.0, 0.0 },
    m_atomic_potential_entry{ nullptr }
{

}

AtomObject::~AtomObject()
{

}

AtomObject::AtomObject(const AtomObject & other) :
    m_key_tag{ other.m_key_tag },
    m_is_selected{ other.m_is_selected }, m_is_special_atom{ other.m_is_special_atom },
    m_serial_id{ other.m_serial_id }, m_residue_id{ other.m_residue_id },
    m_chain_id{ other.m_chain_id }, m_indicator{ other.m_indicator },
    m_occupancy{ other.m_occupancy }, m_temperature{ other.m_temperature },
    m_residue{ other.m_residue }, m_element{ other.m_element },
    m_remoteness{ other.m_remoteness }, m_branch{ other.m_branch },
    m_structure{ other.m_structure },
    m_position{ other.m_position }
{

}

std::unique_ptr<DataObjectBase> AtomObject::Clone() const
{
    return std::make_unique<AtomObject>(*this);
}

void AtomObject::Display(void) const
{
    Logger::Log(LogLevel::Info, "AtomObject Display: " + GetInfo());
}

void AtomObject::Update(void)
{
    Logger::Log(LogLevel::Info, "AtomObject Update: " + GetInfo());
}

void AtomObject::Accept(DataObjectVisitorBase * visitor)
{
    visitor->VisitAtomObject(this);
}

std::unique_ptr<AtomObject> AtomObject::AtomObjectClone(void) const
{
    auto base_clone{ this->Clone() };
    auto derived_ptr{ dynamic_cast<AtomObject*>(base_clone.release()) };
    if (!derived_ptr)
    {
        throw std::runtime_error("Clone() did not return an AtomObject instance!");
    }
    return std::unique_ptr<AtomObject>(derived_ptr);
}

void AtomObject::SetResidue(Residue value) { m_residue = value; }
void AtomObject::SetElement(Element value) { m_element = value; }
void AtomObject::SetRemoteness(Remoteness value) { m_remoteness = value; }
void AtomObject::SetBranch(Branch value) { m_branch = value; }
void AtomObject::SetStructure(Structure value) { m_structure = value; }

void AtomObject::SetResidue(const std::string & name)
{
    m_residue = AtomicInfoHelper::GetResidueFromString(name);
}

void AtomObject::SetElement(const std::string & name)
{
    m_element = AtomicInfoHelper::GetElementFromString(name);
}

void AtomObject::SetRemoteness(const std::string & name)
{
    m_remoteness = AtomicInfoHelper::GetRemotenessFromString(name);
}

void AtomObject::SetBranch(const std::string & name)
{
    m_branch = AtomicInfoHelper::GetBranchFromString(name);
}

void AtomObject::SetPosition(float x, float y, float z)
{
    m_position.at(0) = x;
    m_position.at(1) = y;
    m_position.at(2) = z;
}

std::string AtomObject::GetInfo(void) const
{
    return "";
}

void AtomObject::AddAtomicPotentialEntry(std::unique_ptr<AtomicPotentialEntry> entry)
{
    m_atomic_potential_entry = std::move(entry);
}

void AtomObject::AddAlternatePosition(
    const std::string & indicator, const std::array<float, 3> & value)
{
    m_alternate_position_map[indicator] = value;
}

void AtomObject::AddAlternateOccupancy(
    const std::string & indicator, float value)
{
    m_alternate_occupancy_map[indicator] = value;
}

void AtomObject::AddAlternateTemperature(
    const std::string & indicator, float value)
{
    m_alternate_temperature_map[indicator] = value;
}

Element AtomObject::GetElement(void) const { return m_element; }
Residue AtomObject::GetResidue(void) const { return m_residue; }
Remoteness AtomObject::GetRemoteness(void) const { return m_remoteness; }
Branch AtomObject::GetBranch(void) const { return m_branch; }
Structure AtomObject::GetStructure(void) const { return m_structure; }

bool AtomObject::IsUnknownAtom(void) const
{
    auto element{ GetElement() };
    auto residue{ GetResidue() };
    auto remoteness{ GetRemoteness() };
    auto branch{ GetBranch() };
    if (element == Element::UNK ||
        residue == Residue::UNK ||
        remoteness == Remoteness::UNK ||
        branch == Branch::UNK)
    {
        return true;
    }
    return false;
}
