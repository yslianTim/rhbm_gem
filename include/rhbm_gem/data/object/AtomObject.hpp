#pragma once

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>

namespace rhbm_gem {

class AtomObject;
class AtomLocalPotentialView;
class ModelAnalysisData;
class ModelObject;

class AtomObject
{
    bool m_is_selected;
    bool m_is_special_atom;
    int m_serial_id, m_residue_id;
    std::string m_component_id, m_atom_id, m_chain_id, m_indicator;
    double m_occupancy, m_temperature;
    double m_standard_qscore;
    ComponentKey m_component_key;
    AtomKey m_atom_key;
    Residue m_residue;
    Element m_element;
    Spot m_spot;
    Structure m_structure;
    ModelObject * m_owner_model;
    std::array<double, 3> m_position;
    std::unordered_map<std::string, std::array<double, 3>> m_alternate_position_map;
    std::unordered_map<std::string, double> m_alternate_occupancy_map;
    std::unordered_map<std::string, double> m_alternate_temperature_map;

public:
    AtomObject();
    ~AtomObject();
    AtomObject(const AtomObject & other);

    void SetSpecialAtomFlag(bool value) { m_is_special_atom = value; }
    void SetSerialID(int value) { m_serial_id = value; }
    void SetSequenceID(int value) { m_residue_id = value; }
    void SetComponentID(const std::string & value);
    void SetComponentKey(ComponentKey value) { m_component_key = value; }
    void SetAtomID(const std::string & value);
    void SetAtomKey(AtomKey value) { m_atom_key = value; }
    void SetChainID(const std::string & value) { m_chain_id = value; }
    void SetIndicator(const std::string & value) { m_indicator = value; }
    void SetOccupancy(double value) { m_occupancy = value; }
    void SetTemperature(double value) { m_temperature = value; }
    void SetStandardQScore(double value) { m_standard_qscore = value; }
    void SetResidue(Residue value) { m_residue = value; }
    void SetElement(Element value) { m_element = value; }
    void SetSpot(Spot value) { m_spot = value; }
    void SetStructure(Structure value) { m_structure = value; }
    void SetPosition(double x, double y, double z);
    void SetPosition(const std::array<double, 3> & value);
    void AddAlternatePosition(const std::string & indicator, const std::array<double, 3> & value);
    void AddAlternateOccupancy(const std::string & indicator, double value);
    void AddAlternateTemperature(const std::string & indicator, double value);
    std::string GetInfo() const;
    ComponentKey GetComponentKey() const { return m_component_key; }
    AtomKey GetAtomKey() const { return m_atom_key; }
    Element GetElement() const { return m_element; }
    Residue GetResidue() const { return m_residue; }
    Spot GetSpot() const { return m_spot; }
    Structure GetStructure() const { return m_structure; }
    bool IsUnknownAtom() const;
    bool IsMainChainAtom() const;
    bool GetSpecialAtomFlag() const { return m_is_special_atom; }
    int GetSerialID() const { return m_serial_id; }
    int GetSequenceID() const { return m_residue_id; }
    std::string GetComponentID() const { return m_component_id; }
    std::string GetAtomID() const { return m_atom_id; }
    std::string GetChainID() const { return m_chain_id; }
    std::string GetIndicator() const { return m_indicator; }
    double GetOccupancy() const { return m_occupancy; }
    double GetTemperature() const { return m_temperature; }
    double GetStandardQScore() const { return m_standard_qscore; }
    std::array<double, 3> GetPosition() const { return m_position; }
    const std::array<double, 3> & GetPositionRef() const { return m_position; }
    const std::unordered_map<std::string, std::array<double, 3>> & GetAlternatePositions() const;
    const std::unordered_map<std::string, double> & GetAlternateOccupancies() const;
    const std::unordered_map<std::string, double> & GetAlternateTemperatures() const;
    std::vector<AtomObject *> FindNeighborAtoms(double radius = 2.5, bool include_self = false) const;

private:
    friend class AtomLocalPotentialView;
    friend class ModelObject;

    void SetSelectedFlag(bool value) { m_is_selected = value; }
    void SetOwnerModel(ModelObject * value) { m_owner_model = value; }
};

} // namespace rhbm_gem
