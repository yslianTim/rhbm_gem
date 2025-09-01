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

void ChemicalComponentEntry::AddComponentAtomEntry(
    AtomKey atom_key, const ComponentAtomEntry & atom_info)
{
    m_component_atom_entry_map[atom_key] = atom_info;
}

void ChemicalComponentEntry::AddComponentBondEntry(
    const std::pair<AtomKey, AtomKey> & atom_key_pair,
    const ComponentBondEntry & bond_info)
{
    m_component_bond_entry_map[atom_key_pair] = bond_info;
}
