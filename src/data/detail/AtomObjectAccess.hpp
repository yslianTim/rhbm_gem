#pragma once

#include <array>
#include <memory>
#include <string>

#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>

namespace rhbm_gem {

class AtomObjectAccess
{
public:
    static void SetSelectedFlag(AtomObject & atom_object, bool value) { atom_object.m_is_selected = value; }
    static void SetSpecialAtomFlag(AtomObject & atom_object, bool value) { atom_object.m_is_special_atom = value; }
    static void SetSerialID(AtomObject & atom_object, int value) { atom_object.m_serial_id = value; }
    static void SetSequenceID(AtomObject & atom_object, int value) { atom_object.m_residue_id = value; }
    static void SetComponentID(AtomObject & atom_object, const std::string & value)
    {
        atom_object.SetComponentID(value);
    }
    static void SetComponentKey(AtomObject & atom_object, ComponentKey value) { atom_object.m_component_key = value; }
    static void SetAtomID(AtomObject & atom_object, const std::string & value) { atom_object.SetAtomID(value); }
    static void SetAtomKey(AtomObject & atom_object, AtomKey value) { atom_object.m_atom_key = value; }
    static void SetChainID(AtomObject & atom_object, const std::string & value) { atom_object.m_chain_id = value; }
    static void SetIndicator(AtomObject & atom_object, const std::string & value) { atom_object.m_indicator = value; }
    static void SetOccupancy(AtomObject & atom_object, float value) { atom_object.m_occupancy = value; }
    static void SetTemperature(AtomObject & atom_object, float value) { atom_object.m_temperature = value; }
    static void SetResidue(AtomObject & atom_object, Residue value) { atom_object.SetResidue(value); }
    static void SetElement(AtomObject & atom_object, Element value) { atom_object.SetElement(value); }
    static void SetSpot(AtomObject & atom_object, Spot value) { atom_object.SetSpot(value); }
    static void SetStructure(AtomObject & atom_object, Structure value) { atom_object.SetStructure(value); }
    static void SetResidue(AtomObject & atom_object, const std::string & name) { atom_object.SetResidue(name); }
    static void SetElement(AtomObject & atom_object, const std::string & name) { atom_object.SetElement(name); }
    static void SetSpot(AtomObject & atom_object, const std::string & name) { atom_object.SetSpot(name); }
    static void SetPosition(AtomObject & atom_object, float x, float y, float z) { atom_object.SetPosition(x, y, z); }
    static void SetPosition(AtomObject & atom_object, const std::array<float, 3> & value)
    {
        atom_object.m_position = value;
    }
    static void AddAlternatePosition(
        AtomObject & atom_object,
        const std::string & indicator,
        const std::array<float, 3> & value)
    {
        atom_object.AddAlternatePosition(indicator, value);
    }
    static void AddAlternateOccupancy(AtomObject & atom_object, const std::string & indicator, float value)
    {
        atom_object.AddAlternateOccupancy(indicator, value);
    }
    static void AddAlternateTemperature(AtomObject & atom_object, const std::string & indicator, float value)
    {
        atom_object.AddAlternateTemperature(indicator, value);
    }
    static void SetLocalPotentialEntry(
        AtomObject & atom_object,
        std::unique_ptr<LocalPotentialEntry> entry)
    {
        atom_object.SetLocalPotentialEntry(std::move(entry));
    }
    static LocalPotentialEntry * LocalPotentialEntryPtr(AtomObject & atom_object)
    {
        return atom_object.m_local_potential_entry.get();
    }
    static const LocalPotentialEntry * LocalPotentialEntryPtr(const AtomObject & atom_object)
    {
        return atom_object.m_local_potential_entry.get();
    }
};

} // namespace rhbm_gem
