#include "AtomKeySystem.hpp"
#include "GlobalEnumClass.hpp"
#include "AtomicInfoHelper.hpp"
#include "Logger.hpp"

AtomKeySystem::AtomKeySystem(void)
{
    Logger::Log(LogLevel::Debug, "AtomKeySystem::AtomKeySystem() called");
    //const auto & build_in_residue_map{ AtomicInfoHelper::GetResidueMap() };
    //for (const auto & [id, residue] : build_in_residue_map)
    //{
    //    auto component_key{ static_cast<AtomKey>(residue) };
    //    m_id_to_key_map.emplace(id, component_key);
    //    m_key_to_id_map.emplace(component_key, std::string{id});
    //};
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

void AtomKeySystem::RegisterAtom(const std::string & id)
{
    Logger::Log(LogLevel::Debug, "AtomKeySystem::RegisterAtom() called for " + id);
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(id) != m_id_to_key_map.end()) return;
    AtomKey new_key{ m_next_dynamic_key++ };
    m_id_to_key_map[id] = new_key;
    m_key_to_id_map[new_key] = id;
}

AtomKey AtomKeySystem::GetAtomKey(const std::string & id)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(id) == m_id_to_key_map.end())
    {
        Logger::Log(LogLevel::Warning, 
            "AtomKeySystem::GetAtomKey() - Unknown component id: " + id);
        return static_cast<AtomKey>(0);
    }
    return m_id_to_key_map.at(id);
}

std::string AtomKeySystem::GetComponentId(AtomKey key)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_key_to_id_map.find(key) == m_key_to_id_map.end())
    {
        Logger::Log(LogLevel::Warning, 
            "AtomKeySystem::GetComponentId() - Unknown component key: "+ std::to_string(key));
        return "UNK";
    }
    return m_key_to_id_map.at(key);
}
