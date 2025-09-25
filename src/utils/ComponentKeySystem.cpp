#include "ComponentKeySystem.hpp"
#include "GlobalEnumClass.hpp"
#include "ChemicalDataHelper.hpp"
#include "Logger.hpp"

#include <limits>

const ComponentKey ComponentKeySystem::k_dynamic_base{ 30 };
const ComponentKey ComponentKeySystem::k_max_key{ std::numeric_limits<ComponentKey>::max() };

ComponentKeySystem::ComponentKeySystem(void) :
    m_next_dynamic_key{ k_dynamic_base }
{
    Logger::Log(LogLevel::Debug, "ComponentKeySystem::ComponentKeySystem() called");
    const auto & build_in_residue_map{ ChemicalDataHelper::GetResidueMap() };
    for (const auto & [component_id, residue] : build_in_residue_map)
    {
        auto component_key{ static_cast<ComponentKey>(residue) };
        m_id_to_key_map.emplace(component_id, component_key);
        m_key_to_id_map.emplace(component_key, std::string{component_id});
    };
}

ComponentKeySystem::ComponentKeySystem(const ComponentKeySystem & other) :
    m_next_dynamic_key{ other.m_next_dynamic_key },
    m_id_to_key_map{ other.m_id_to_key_map },
    m_key_to_id_map{ other.m_key_to_id_map }
{
    Logger::Log(LogLevel::Debug, "ComponentKeySystem copy ctor called");
}

ComponentKeySystem::~ComponentKeySystem(void)
{
    Logger::Log(LogLevel::Debug, "ComponentKeySystem::~ComponentKeySystem() called");
}

void ComponentKeySystem::RegisterComponent(const std::string & component_id)
{
    Logger::Log(LogLevel::Debug, "ComponentKeySystem::RegisterComponent() called for " + component_id);
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(component_id) != m_id_to_key_map.end()) return;
    if (m_next_dynamic_key >= k_max_key)
    {
        Logger::Log(LogLevel::Warning,
            "ComponentKeySystem::RegisterComponent() - Component key overflow for " + component_id
            + ", avoiding register new component key duw to maximum key reached."
        );
        return;
    }
    ComponentKey new_component_key{ m_next_dynamic_key++ };
    m_id_to_key_map[component_id] = new_component_key;
    m_key_to_id_map[new_component_key] = component_id;
}

void ComponentKeySystem::RegisterComponent(
    const std::string & component_id, ComponentKey component_key)
{
    Logger::Log(LogLevel::Debug, "ComponentKeySystem::RegisterComponent() called for " + component_id);
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(component_id) != m_id_to_key_map.end()) return;
    m_id_to_key_map[component_id] = component_key;
    m_key_to_id_map[component_key] = component_id;
}

ComponentKey ComponentKeySystem::GetComponentKey(const std::string & component_id)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_id_to_key_map.find(component_id) == m_id_to_key_map.end())
    {
        Logger::Log(LogLevel::Warning, 
            "ComponentKeySystem::GetComponentKey() - Unknown component id: " + component_id);
        return static_cast<ComponentKey>(Residue::UNK);
    }
    return m_id_to_key_map.at(component_id);
}

std::string ComponentKeySystem::GetComponentId(ComponentKey component_key)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_key_to_id_map.find(component_key) == m_key_to_id_map.end())
    {
        Logger::Log(LogLevel::Warning, 
            "ComponentKeySystem::GetComponentId() - Unknown component key: "
            + std::to_string(component_key));
        return "UNK";
    }
    return m_key_to_id_map.at(component_key);
}

bool ComponentKeySystem::IsBuildInComponent(const std::string & component_id) const
{
    if (m_id_to_key_map.find(component_id) == m_id_to_key_map.end()) return false;
    return m_id_to_key_map.at(component_id) < k_dynamic_base;
}

bool ComponentKeySystem::IsBuildInComponent(ComponentKey component_key) const
{
    return component_key < k_dynamic_base;
}
