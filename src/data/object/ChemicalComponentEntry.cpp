#include <rhbm_gem/data/ChemicalComponentEntry.hpp>
#include <rhbm_gem/utils/Logger.hpp>

namespace rhbm_gem {

ChemicalComponentEntry::ChemicalComponentEntry()
{
}

ChemicalComponentEntry::~ChemicalComponentEntry()
{
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

bool ChemicalComponentEntry::HasComponentAtomEntry(AtomKey atom_key) const
{
    return (m_component_atom_entry_map.find(atom_key) != m_component_atom_entry_map.end());
}

bool ChemicalComponentEntry::HasComponentBondEntry(BondKey bond_key) const
{
    return (m_component_bond_entry_map.find(bond_key) != m_component_bond_entry_map.end());
}

const ComponentAtomEntry * ChemicalComponentEntry::GetComponentAtomEntryPtr(AtomKey atom_key) const
{
    if (m_component_atom_entry_map.find(atom_key) == m_component_atom_entry_map.end())
    {
        Logger::Log(LogLevel::Warning,
            "Component atom entry (atom_key: "
            + std::to_string(atom_key)
            + ") not found in chemical component map.");
        return nullptr;
    }
    return &m_component_atom_entry_map.at(atom_key);
}

const ComponentBondEntry * ChemicalComponentEntry::GetComponentBondEntryPtr(BondKey bond_key) const
{
    if (m_component_bond_entry_map.find(bond_key) == m_component_bond_entry_map.end())
    {
        Logger::Log(LogLevel::Warning,
            "Component bond entry (bond_key: "
            + std::to_string(bond_key)
            + ") not found in chemical component map.");
        return nullptr;
    }
    return &m_component_bond_entry_map.at(bond_key);
}   

} // namespace rhbm_gem
