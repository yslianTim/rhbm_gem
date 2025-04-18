#pragma once

#include <array>
#include "DataObjectBase.hpp"

class AtomicPotentialEntry;
enum class Element : int;
enum class Residue : int;
enum class Remoteness : int;
enum class Branch : int;

class AtomObject : public DataObjectBase
{
    std::string m_key_tag;
    bool m_is_selected;
    bool m_is_special_atom;
    int m_serial_id, m_residue_id;
    std::string m_chain_id, m_indicator;
    float m_occupancy, m_temperature;
    int m_residue_type;
    int m_element_type;
    int m_remoteness_type;
    int m_branch_type;
    int m_status;
    std::array<float, 3> m_position;
    std::unique_ptr<AtomicPotentialEntry> m_atomic_potential_entry;

public:
    AtomObject(void);
    ~AtomObject();
    AtomObject(const AtomObject & other);
    std::unique_ptr<DataObjectBase> Clone(void) const override;
    void Display(void) const override;
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
    void SetResidue(int value) { m_residue_type = value; }
    void SetElement(int value) { m_element_type = value; }
    void SetRemoteness(int value) { m_remoteness_type = value; }
    void SetBranch(int value) { m_branch_type = value; }
    void SetResidue(const std::string & name);
    void SetElement(const std::string & name);
    void SetRemoteness(const std::string & name);
    void SetBranch(const std::string & name);
    void SetStatus(int value) { m_status = value; }
    void SetPosition(float x, float y, float z);
    void SetPosition(const std::array<float, 3> & value) { m_position = value; }
    void AddAtomicPotentialEntry(std::unique_ptr<AtomicPotentialEntry> entry);

    std::string GetInfo(void) const;
    Element GetElementType(void) const;
    Residue GetResidueType(void) const;
    Remoteness GetRemotenessType(void) const;
    Branch GetBranchType(void) const;
    bool IsUnknownAtom(void) const;
    bool GetSelectedFlag(void) const { return m_is_selected; }
    bool GetSpecialAtomFlag(void) const { return m_is_special_atom; }
    int GetSerialID(void) const { return m_serial_id; }
    int GetResidueID(void) const { return m_residue_id; }
    std::string GetChainID(void) const { return m_chain_id; }
    std::string GetIndicator(void) const { return m_indicator; }
    float GetOccupancy(void) const { return m_occupancy; }
    float GetTemperature(void) const { return m_temperature; }
    int GetResidue(void) const { return m_residue_type; }
    int GetElement(void) const { return m_element_type; }
    int GetRemoteness(void) const { return m_remoteness_type; }
    int GetBranch(void) const { return m_branch_type; }
    int GetStatus(void) const { return m_status; }
    std::array<float, 3> GetPosition(void) const { return m_position; }
    AtomicPotentialEntry * GetAtomicPotentialEntry(void) const { return m_atomic_potential_entry.get(); }

private:
    
};