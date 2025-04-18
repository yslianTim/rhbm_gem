#include "AtomObject.hpp"
#include "AtomicInfoHelper.hpp"
#include "DataObjectVisitorBase.hpp"
#include "AtomicPotentialEntry.hpp"

AtomObject::AtomObject(void) :
    m_key_tag{ "" },
    m_is_selected{ true }, m_is_special_atom{ false },
    m_serial_id{ 0 }, m_residue_id{ 0 },
    m_chain_id{ "" }, m_indicator{ "" },
    m_occupancy{ 0.0 }, m_temperature{ 0.0 },
    m_residue_type{ 0 }, m_element_type{ 0 }, m_remoteness_type{ 0 }, m_branch_type{ 0 },
    m_status{ 0 },
    m_position{ 0.0, 0.0, 0.0 }, m_atomic_potential_entry{ nullptr }
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
    m_residue_type{ other.m_residue_type }, m_element_type{ other.m_element_type },
    m_remoteness_type{ other.m_remoteness_type }, m_branch_type{ other.m_branch_type },
    m_status{ other.m_status },
    m_position{ other.m_position }
{

}

std::unique_ptr<DataObjectBase> AtomObject::Clone() const
{
    return std::make_unique<AtomObject>(*this);
}

void AtomObject::Display(void) const
{
    std::cout << "This is AtomObject." << std::endl;
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

void AtomObject::SetResidue(const std::string & name)
{
    m_residue_type = AtomicInfoHelper::AtomInfoMapping<AtomicInfoHelper::ResidueTag>(name);
}

void AtomObject::SetElement(const std::string & name)
{
    m_element_type = AtomicInfoHelper::AtomInfoMapping<AtomicInfoHelper::ElementTag>(name);
}

void AtomObject::SetRemoteness(const std::string & name)
{
    m_remoteness_type = AtomicInfoHelper::AtomInfoMapping<AtomicInfoHelper::RemotenessTag>(name);
}

void AtomObject::SetBranch(const std::string & name)
{
    m_branch_type = AtomicInfoHelper::AtomInfoMapping<AtomicInfoHelper::BranchTag>(name);
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

Element AtomObject::GetElementType(void) const
{
    return static_cast<Element>(m_element_type);
}

Residue AtomObject::GetResidueType(void) const
{
    return static_cast<Residue>(m_residue_type);
}

Remoteness AtomObject::GetRemotenessType(void) const
{
    return static_cast<Remoteness>(m_remoteness_type);
}

Branch AtomObject::GetBranchType(void) const
{
    return static_cast<Branch>(m_branch_type);
}

bool AtomObject::IsUnknownAtom(void) const
{
    auto element{ GetElementType() };
    auto residue{ GetResidueType() };
    auto remoteness{ GetRemotenessType() };
    auto branch{ GetBranchType() };
    if (element == Element::UNK ||
        residue == Residue::UNK ||
        remoteness == Remoteness::UNK ||
        branch == Branch::UNK)
    {
        return true;
    }
    return false;
}