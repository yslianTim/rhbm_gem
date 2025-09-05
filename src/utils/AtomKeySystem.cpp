#include "AtomKeySystem.hpp"
#include "GlobalEnumClass.hpp"
#include "AtomicInfoHelper.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"

#include <algorithm>

const AtomKey AtomKeySystem::k_dynamic_base{ 10000 };

AtomKeySystem::AtomKeySystem(void) :
    m_next_dynamic_key{ k_dynamic_base }
{
    Logger::Log(LogLevel::Debug, "AtomKeySystem::AtomKeySystem() called");
    const auto & build_in_spot_map{ AtomicInfoHelper::GetSpotMap() };
    for (const auto & [atom_id, spot] : build_in_spot_map)
    {
        auto atom_key{ static_cast<AtomKey>(spot) };
        m_id_to_key_map.emplace(atom_id, atom_key);
        m_key_to_id_map.emplace(atom_key, std::string{atom_id});
    }
}

AtomKeySystem::~AtomKeySystem(void)
{
    Logger::Log(LogLevel::Debug, "AtomKeySystem::~AtomKeySystem() called");
}

AtomKeySystem & AtomKeySystem::Instance(void)
{
    Logger::Log(LogLevel::Debug, "AtomKeySystem::Instance() called");
    static AtomKeySystem instance;
    return instance;
}

void AtomKeySystem::RegisterAtom(const std::string & atom_id)
{
    Logger::Log(LogLevel::Debug, "AtomKeySystem::RegisterAtom() called for " + atom_id);
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(atom_id) != m_id_to_key_map.end()) return;
    AtomKey new_atom_key{ m_next_dynamic_key++ };
    m_id_to_key_map[atom_id] = new_atom_key;
    m_key_to_id_map[new_atom_key] = atom_id;
}

void AtomKeySystem::RegisterAtom(const std::string & atom_id, AtomKey atom_key)
{
    Logger::Log(LogLevel::Debug, "AtomKeySystem::RegisterAtom() called for " + atom_id);
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(atom_id) != m_id_to_key_map.end()) return;
    m_id_to_key_map[atom_id] = atom_key;
    m_key_to_id_map[atom_key] = atom_id;
}

AtomKey AtomKeySystem::GetAtomKey(const std::string & atom_id)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(atom_id) == m_id_to_key_map.end())
    {
        Logger::Log(LogLevel::Warning, 
            "AtomKeySystem::GetAtomKey() - Unknown atom id: " + atom_id);
        return static_cast<AtomKey>(0);
    }
    return m_id_to_key_map.at(atom_id);
}

std::string AtomKeySystem::GetAtomId(AtomKey atom_key)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_key_to_id_map.find(atom_key) == m_key_to_id_map.end())
    {
        Logger::Log(LogLevel::Warning, 
            "AtomKeySystem::GetAtomId() - Unknown atom key: "+ std::to_string(atom_key));
        return "UNK";
    }
    return m_key_to_id_map.at(atom_key);
}

bool AtomKeySystem::IsBuildInAtom(const std::string & atom_id) const
{
    if (m_id_to_key_map.find(atom_id) == m_id_to_key_map.end()) return false;
    return m_id_to_key_map.at(atom_id) < k_dynamic_base;
}

bool AtomKeySystem::IsBuildInAtom(AtomKey atom_key) const
{
    return atom_key < k_dynamic_base;
}
