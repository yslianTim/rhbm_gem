#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <stdexcept>

namespace rhbm_gem {

AtomObject::AtomObject() :
    m_key_tag{ "" },
    m_is_selected{ false }, m_is_special_atom{ false },
    m_serial_id{ 0 }, m_residue_id{ 0 },
    m_component_id{ "" }, m_atom_id{ "" }, m_chain_id{ "" }, m_indicator{ "" },
    m_occupancy{ 0.0 }, m_temperature{ 0.0 },
    m_component_key{ 0 }, m_atom_key{ 0 },
    m_residue{ Residue::UNK }, m_element{ Element::UNK }, m_spot{ Spot::UNK },
    m_structure{ Structure::UNK },
    m_position{ 0.0, 0.0, 0.0 },
    m_local_potential_entry{ nullptr }
{

}

AtomObject::~AtomObject()
{

}

AtomObject::AtomObject(const AtomObject & other) :
    m_key_tag{ other.m_key_tag },
    m_is_selected{ other.m_is_selected }, m_is_special_atom{ other.m_is_special_atom },
    m_serial_id{ other.m_serial_id }, m_residue_id{ other.m_residue_id },
    m_component_id{ other.m_component_id }, m_atom_id{ other.m_atom_id },
    m_chain_id{ other.m_chain_id }, m_indicator{ other.m_indicator },
    m_occupancy{ other.m_occupancy }, m_temperature{ other.m_temperature },
    m_component_key{ other.m_component_key }, m_atom_key{ other.m_atom_key },
    m_residue{ other.m_residue }, m_element{ other.m_element }, m_spot{ other.m_spot },
    m_structure{ other.m_structure },
    m_position{ other.m_position }
{

}

std::unique_ptr<DataObjectBase> AtomObject::Clone() const
{
    return std::make_unique<AtomObject>(*this);
}

void AtomObject::Display() const
{
    Logger::Log(LogLevel::Info, "AtomObject Display: " + GetInfo());
}

void AtomObject::Update()
{
    Logger::Log(LogLevel::Info, "AtomObject Update: " + GetInfo());
}

std::unique_ptr<AtomObject> AtomObject::AtomObjectClone() const
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
void AtomObject::SetSpot(Spot value) { m_spot = value; }
void AtomObject::SetStructure(Structure value) { m_structure = value; }

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

void AtomObject::SetResidue(const std::string & name)
{
    m_residue = ChemicalDataHelper::GetResidueFromString(name, false);
}

void AtomObject::SetElement(const std::string & name)
{
    m_element = ChemicalDataHelper::GetElementFromString(name);
}

void AtomObject::SetSpot(const std::string & name)
{
    m_spot = ChemicalDataHelper::GetSpotFromString(name, false);
}

void AtomObject::SetPosition(float x, float y, float z)
{
    m_position.at(0) = x;
    m_position.at(1) = y;
    m_position.at(2) = z;
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

void AtomObject::AddLocalPotentialEntry(std::unique_ptr<LocalPotentialEntry> entry)
{
    m_local_potential_entry = std::move(entry);
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

Element AtomObject::GetElement() const { return m_element; }
Residue AtomObject::GetResidue() const { return m_residue; }
Spot AtomObject::GetSpot() const { return m_spot; }
Structure AtomObject::GetStructure() const { return m_structure; }

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
