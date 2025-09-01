#include "ComponentKeySystem.hpp"
#include "GlobalEnumClass.hpp"
#include "AtomicInfoHelper.hpp"
#include "Logger.hpp"

ComponentKeySystem::ComponentKeySystem(void)
{
    Logger::Log(LogLevel::Debug, "ComponentKeySystem::ComponentKeySystem() called");
    const auto & build_in_residue_map{ AtomicInfoHelper::GetResidueMap() };
    for (const auto & [id, residue] : build_in_residue_map)
    {
        auto component_key{ static_cast<ComponentKey>(residue) };
        m_id_to_key_map.emplace(id, component_key);
        m_key_to_id_map.emplace(component_key, std::string{id});
    };
}

ComponentKeySystem::~ComponentKeySystem(void)
{
    Logger::Log(LogLevel::Debug, "ComponentKeySystem::~ComponentKeySystem() called");
}

ComponentKeySystem & ComponentKeySystem::Instance(void)
{
    Logger::Log(LogLevel::Debug, "ComponentKeySystem::Instance() called");
    static ComponentKeySystem instance;
    return instance;
}

void ComponentKeySystem::RegisterComponent(const std::string & id)
{
    Logger::Log(LogLevel::Debug, "ComponentKeySystem::RegisterComponent() called for " + id);
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(id) != m_id_to_key_map.end()) return;
    ComponentKey new_key{ m_next_dynamic_key++ };
    m_id_to_key_map[id] = new_key;
    m_key_to_id_map[new_key] = id;
}

ComponentKey ComponentKeySystem::GetComponentKey(const std::string & id)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(id) == m_id_to_key_map.end())
    {
        Logger::Log(LogLevel::Warning, 
            "ComponentKeySystem::GetComponentKey() - Unknown component id: " + id);
        return static_cast<ComponentKey>(Residue::UNK);
    }
    return m_id_to_key_map.at(id);
}

std::string ComponentKeySystem::GetComponentId(ComponentKey key)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_key_to_id_map.find(key) == m_key_to_id_map.end())
    {
        Logger::Log(LogLevel::Warning, 
            "ComponentKeySystem::GetComponentId() - Unknown component key: "+ std::to_string(key));
        return "UNK";
    }
    return m_key_to_id_map.at(key);
}
