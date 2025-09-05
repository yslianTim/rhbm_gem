#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <mutex>

using ComponentKey = uint8_t;

class ComponentKeySystem
{
    static const ComponentKey k_dynamic_base;
    std::mutex m_mutex;
    ComponentKey m_next_dynamic_key;
    std::unordered_map<std::string, ComponentKey> m_id_to_key_map;
    std::unordered_map<ComponentKey, std::string> m_key_to_id_map;
    
public:
    ComponentKeySystem(void);
    ~ComponentKeySystem();
    ComponentKeySystem(const ComponentKeySystem & other);
    ComponentKeySystem & operator=(const ComponentKeySystem & other);

    void RegisterComponent(const std::string & component_id);
    void RegisterComponent(const std::string & component_id, ComponentKey component_key);
    ComponentKey GetComponentKey(const std::string & component_id);
    std::string GetComponentId(ComponentKey component_key);
    bool IsBuildInComponent(const std::string & component_id) const;
    bool IsBuildInComponent(ComponentKey component_key) const;

};
