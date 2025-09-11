#include "ChemicalComponentEntry.hpp"
#include "AtomicInfoHelper.hpp"
#include "Logger.hpp"

ChemicalComponentEntry::ChemicalComponentEntry(void)
{
    Logger::Log(LogLevel::Debug, "ChemicalComponentEntry::ChemicalComponentEntry() called");
}

ChemicalComponentEntry::~ChemicalComponentEntry()
{
    Logger::Log(LogLevel::Debug, "ChemicalComponentEntry::~ChemicalComponentEntry() called");
}

void ChemicalComponentEntry::SetComponentAtomEntryMap(
    std::map<AtomKey, ComponentAtomEntry> & atom_entry_map)
{
    m_component_atom_entry_map = std::move(atom_entry_map);
}

void ChemicalComponentEntry::SetComponentBondEntryMap(
    std::map<BondKey, ComponentBondEntry> & bond_entry_map)
{
    m_component_bond_entry_map = std::move(bond_entry_map);
}

void ChemicalComponentEntry::AddComponentAtomEntry(
    AtomKey atom_key, const ComponentAtomEntry & atom_info)
{
    m_component_atom_entry_map[atom_key] = atom_info;
}

void ChemicalComponentEntry::AddComponentBondEntry(
    BondKey bond_key, const ComponentBondEntry & bond_info)
{
    m_component_bond_entry_map[bond_key] = bond_info;
}
