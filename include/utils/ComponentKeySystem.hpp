#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <mutex>

using ComponentKey = uint16_t;

class ComponentKeySystem
{
    static constexpr ComponentKey kDynamicBase{ 1000 };
    std::mutex m_mutex;
    ComponentKey m_next_dynamic_key{ kDynamicBase };
    std::unordered_map<std::string, ComponentKey> m_id_to_key_map;
    std::unordered_map<ComponentKey, std::string> m_key_to_id_map;
    
public:
    static ComponentKeySystem & Instance(void);
    void RegisterComponent(const std::string & id);
    ComponentKey GetComponentKey(const std::string & id);
    std::string GetComponentId(ComponentKey key);

private:
    ComponentKeySystem(void);
    ~ComponentKeySystem(void);

};
