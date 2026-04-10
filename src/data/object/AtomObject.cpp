#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>

#include <stdexcept>

namespace rhbm_gem {

AtomObject::AtomObject() :
    m_is_selected{ false }, m_is_special_atom{ false },
    m_serial_id{ 0 }, m_residue_id{ 0 },
    m_component_id{ "" }, m_atom_id{ "" }, m_chain_id{ "" }, m_indicator{ "" },
    m_occupancy{ 0.0 }, m_temperature{ 0.0 },
    m_component_key{ 0 }, m_atom_key{ 0 },
    m_residue{ Residue::UNK }, m_element{ Element::UNK }, m_spot{ Spot::UNK },
    m_structure{ Structure::UNK }, m_owner_model{ nullptr },
    m_position{ 0.0, 0.0, 0.0 }
{

}

AtomObject::~AtomObject()
{

}

AtomObject::AtomObject(const AtomObject & other) :
    m_is_selected{ other.m_is_selected }, m_is_special_atom{ other.m_is_special_atom },
    m_serial_id{ other.m_serial_id }, m_residue_id{ other.m_residue_id },
    m_component_id{ other.m_component_id }, m_atom_id{ other.m_atom_id },
    m_chain_id{ other.m_chain_id }, m_indicator{ other.m_indicator },
    m_occupancy{ other.m_occupancy }, m_temperature{ other.m_temperature },
    m_component_key{ other.m_component_key }, m_atom_key{ other.m_atom_key },
    m_residue{ other.m_residue }, m_element{ other.m_element }, m_spot{ other.m_spot },
    m_structure{ other.m_structure }, m_owner_model{ nullptr },
    m_position{ other.m_position },
    m_alternate_position_map{ other.m_alternate_position_map },
    m_alternate_occupancy_map{ other.m_alternate_occupancy_map },
    m_alternate_temperature_map{ other.m_alternate_temperature_map }
{
}

void AtomObject::SetComponentID(const std::string & value)
{
    m_component_id = value;
    m_residue = ChemicalDataHelper::GetResidueFromString(m_component_id, false);
}

void AtomObject::SetAtomID(const std::string & value)
{
    m_atom_id = value;
    m_spot = ChemicalDataHelper::GetSpotFromString(m_atom_id, false);
}

void AtomObject::SetPosition(float x, float y, float z)
{
    m_position.at(0) = x;
    m_position.at(1) = y;
    m_position.at(2) = z;
    if (m_owner_model != nullptr)
    {
        m_owner_model->InvalidateDerivedState();
    }
}

void AtomObject::SetPosition(const std::array<float, 3> & value)
{
    m_position = value;
    if (m_owner_model != nullptr)
    {
        m_owner_model->InvalidateDerivedState();
    }
}

std::string AtomObject::GetInfo() const
{
    return "[Serial ID] " + std::to_string(m_serial_id) + " " +
           "[Residue ID] " + std::to_string(m_residue_id) + " " +
           "[Atom ID] " + m_atom_id + " " +
           "[Chain ID] " + m_chain_id + " " +
           ChemicalDataHelper::GetLabel(m_residue) + " " +
           ChemicalDataHelper::GetLabel(m_element) + " " +
           "[Component Key] " + std::to_string(m_component_key) + " " +
           "[Atom Key] " + std::to_string(m_atom_key) + " " +
           m_indicator + " " +
           "Position = (" +
           std::to_string(m_position.at(0)) + ", " +
           std::to_string(m_position.at(1)) + ", " +
           std::to_string(m_position.at(2)) + ")";
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

const std::unordered_map<std::string, std::array<float, 3>> &
AtomObject::GetAlternatePositions() const
{
    return m_alternate_position_map;
}

const std::unordered_map<std::string, float> &
AtomObject::GetAlternateOccupancies() const
{
    return m_alternate_occupancy_map;
}

const std::unordered_map<std::string, float> &
AtomObject::GetAlternateTemperatures() const
{
    return m_alternate_temperature_map;
}

std::vector<AtomObject *> AtomObject::FindNeighborAtoms(double radius, bool include_self) const
{
    if (m_owner_model == nullptr)
    {
        throw std::runtime_error("AtomObject::FindNeighborAtoms requires an attached owner model.");
    }
    return m_owner_model->FindNeighborAtoms(*this, radius, include_self);
}

bool AtomObject::IsUnknownAtom() const
{
    if (m_element == Element::UNK || m_residue == Residue::UNK || m_spot == Spot::UNK)
    {
        return true;
    }
    return false;
}

bool AtomObject::IsMainChainAtom() const
{
    size_t dummy_id;
    return AtomClassifier::IsMainChainMember(m_spot, dummy_id);
}

} // namespace rhbm_gem
