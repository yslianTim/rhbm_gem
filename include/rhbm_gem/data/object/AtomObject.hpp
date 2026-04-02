#pragma once

#include <array>
#include <string>
#include <unordered_map>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>

namespace rhbm_gem {

class LocalPotentialEntry;

class AtomObject
{
    bool m_is_selected;
    bool m_is_special_atom;
    int m_serial_id, m_residue_id;
    std::string m_component_id, m_atom_id, m_chain_id, m_indicator;
    float m_occupancy, m_temperature;
    ComponentKey m_component_key;
    AtomKey m_atom_key;
    Residue m_residue;
    Element m_element;
    Spot m_spot;
    Structure m_structure;
    std::array<float, 3> m_position;
    std::unordered_map<std::string, std::array<float, 3>> m_alternate_position_map;
    std::unordered_map<std::string, float> m_alternate_occupancy_map;
    std::unordered_map<std::string, float> m_alternate_temperature_map;
    std::unique_ptr<LocalPotentialEntry> m_local_potential_entry;

public:
    AtomObject();
    ~AtomObject();
    AtomObject(const AtomObject & other);

    void SetSelectedFlag(bool value) { m_is_selected = value; }
    void SetSpecialAtomFlag(bool value) { m_is_special_atom = value; }
    void SetSerialID(int value) { m_serial_id = value; }
    void SetSequenceID(int value) { m_residue_id = value; }
    void SetComponentID(const std::string & value);
    void SetComponentKey(ComponentKey value) { m_component_key = value; }
    void SetAtomID(const std::string & value);
    void SetAtomKey(AtomKey value) { m_atom_key = value; }
    void SetChainID(const std::string & value) { m_chain_id = value; }
    void SetIndicator(const std::string & value) { m_indicator = value; }
    void SetOccupancy(float value) { m_occupancy = value; }
    void SetTemperature(float value) { m_temperature = value; }
    void SetResidue(Residue value);
    void SetElement(Element value);
    void SetSpot(Spot value);
    void SetStructure(Structure value);
    void SetResidue(const std::string & name);
    void SetElement(const std::string & name);
    void SetSpot(const std::string & name);
    void SetPosition(float x, float y, float z);
    void SetPosition(const std::array<float, 3> & value) { m_position = value; }
    void SetLocalPotentialEntry(std::unique_ptr<LocalPotentialEntry> entry);
    void AddAlternatePosition(const std::string & indicator, const std::array<float, 3> & value);
    void AddAlternateOccupancy(const std::string & indicator, float value);
    void AddAlternateTemperature(const std::string & indicator, float value);
    std::string GetInfo() const;
    ComponentKey GetComponentKey() const { return m_component_key; }
    AtomKey GetAtomKey() const { return m_atom_key; }
    Element GetElement() const;
    Residue GetResidue() const;
    Spot GetSpot() const;
    Structure GetStructure() const;
    bool IsUnknownAtom() const;
    bool IsMainChainAtom() const;
    bool GetSelectedFlag() const { return m_is_selected; }
    bool GetSpecialAtomFlag() const { return m_is_special_atom; }
    int GetSerialID() const { return m_serial_id; }
    int GetSequenceID() const { return m_residue_id; }
    std::string GetComponentID() const { return m_component_id; }
    std::string GetAtomID() const { return m_atom_id; }
    std::string GetChainID() const { return m_chain_id; }
    std::string GetIndicator() const { return m_indicator; }
    float GetOccupancy() const { return m_occupancy; }
    float GetTemperature() const { return m_temperature; }
    std::array<float, 3> GetPosition() const { return m_position; }
    const std::array<float, 3> & GetPositionRef() const { return m_position; }
    const std::unordered_map<std::string, std::array<float, 3>> & GetAlternatePositions() const;
    const std::unordered_map<std::string, float> & GetAlternateOccupancies() const;
    const std::unordered_map<std::string, float> & GetAlternateTemperatures() const;
    LocalPotentialEntry * GetLocalPotentialEntry() const { return m_local_potential_entry.get(); }

};

} // namespace rhbm_gem
