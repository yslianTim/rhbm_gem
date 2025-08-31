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
    const std::string & atom_id, const ComponentAtomEntry & atom_info)
{
    m_component_atom_entry_map[atom_id] = atom_info;
}

void ChemicalComponentEntry::AddComponentBondEntry(
    const std::pair<std::string, std::string> & atom_id_pair,
    const ComponentBondEntry & bond_info)
{
    m_component_bond_entry_map[atom_id_pair] = bond_info;
}
