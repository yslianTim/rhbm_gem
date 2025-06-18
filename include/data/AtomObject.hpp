#pragma once

#include <array>
#include <string>
#include <unordered_map>
#include "DataObjectBase.hpp"

class AtomicPotentialEntry;
class AtomicChargeEntry;
enum class Residue : uint16_t;
enum class Element : uint16_t;
enum class Remoteness : uint8_t;
enum class Branch : uint8_t;
enum class Structure : uint8_t;

class AtomObject : public DataObjectBase
{
    std::string m_key_tag;
    bool m_is_selected;
    bool m_is_special_atom;
    int m_serial_id, m_residue_id;
    std::string m_chain_id, m_indicator;
    float m_occupancy, m_temperature;
    Residue m_residue;
    Element m_element;
    Remoteness m_remoteness;
    Branch m_branch;
    Structure m_structure;
    std::array<float, 3> m_position;
    std::unordered_map<std::string, std::array<float, 3>> m_alternate_position_map;
    std::unordered_map<std::string, float> m_alternate_occupancy_map;
    std::unordered_map<std::string, float> m_alternate_temperature_map;
    std::unique_ptr<AtomicPotentialEntry> m_atomic_potential_entry;
    std::unique_ptr<AtomicChargeEntry> m_atomic_charge_entry;

public:
    AtomObject(void);
    ~AtomObject();
    AtomObject(const AtomObject & other);
    std::unique_ptr<DataObjectBase> Clone(void) const override;
    void Display(void) const override;
    void Update(void) override;
    void Accept(DataObjectVisitorBase * visitor) override;
    void SetKeyTag(const std::string & label) override { m_key_tag = label; }
    std::string GetKeyTag(void) const override { return m_key_tag; }

    std::unique_ptr<AtomObject> AtomObjectClone(void) const;
    void SetSelectedFlag(bool value) { m_is_selected = value; }
    void SetSpecialAtomFlag(bool value) { m_is_special_atom = value; }
    void SetSerialID(int value) { m_serial_id = value; }
    void SetResidueID(int value) { m_residue_id = value; }
    void SetChainID(const std::string & value) { m_chain_id = value; }
    void SetIndicator(const std::string & value) { m_indicator = value; }
    void SetOccupancy(float value) { m_occupancy = value; }
    void SetTemperature(float value) { m_temperature = value; }
    void SetResidue(Residue value);
    void SetElement(Element value);
    void SetRemoteness(Remoteness value);
    void SetBranch(Branch value);
    void SetStructure(Structure value);
    void SetResidue(const std::string & name);
    void SetElement(const std::string & name);
    void SetRemoteness(const std::string & name);
    void SetBranch(const std::string & name);
    void SetPosition(float x, float y, float z);
    void SetPosition(const std::array<float, 3> & value) { m_position = value; }
    void AddAtomicPotentialEntry(std::unique_ptr<AtomicPotentialEntry> entry);
    void AddAtomicChargeEntry(std::unique_ptr<AtomicChargeEntry> entry);
    void AddAlternatePosition(const std::string & indicator, const std::array<float, 3> & value);
    void AddAlternateOccupancy(const std::string & indicator, float value);
    void AddAlternateTemperature(const std::string & indicator, float value);

    std::string GetInfo(void) const;
    Element GetElement(void) const;
    Residue GetResidue(void) const;
    Remoteness GetRemoteness(void) const;
    Branch GetBranch(void) const;
    Structure GetStructure(void) const;
    bool IsUnknownAtom(void) const;
    bool GetSelectedFlag(void) const { return m_is_selected; }
    bool GetSpecialAtomFlag(void) const { return m_is_special_atom; }
    int GetSerialID(void) const { return m_serial_id; }
    int GetResidueID(void) const { return m_residue_id; }
    std::string GetChainID(void) const { return m_chain_id; }
    std::string GetIndicator(void) const { return m_indicator; }
    float GetOccupancy(void) const { return m_occupancy; }
    float GetTemperature(void) const { return m_temperature; }
    std::array<float, 3> GetPosition(void) const { return m_position; }
    const std::array<float, 3> & GetPositionRef(void) const { return m_position; }
    AtomicPotentialEntry * GetAtomicPotentialEntry(void) const { return m_atomic_potential_entry.get(); }
    AtomicChargeEntry * GetAtomicChargeEntry(void) const { return m_atomic_charge_entry.get(); }

private:
    
};
