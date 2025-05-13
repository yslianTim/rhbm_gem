#pragma once

#include <string>
#include <vector>
#include <unordered_map>

enum class Element : uint16_t;

class AtomicModelDataBlock
{
    std::string m_map_id, m_model_id;
    std::string m_resolution, m_resolution_method;

    std::vector<Element> m_element_type_list;

public:
    AtomicModelDataBlock(void);
    ~AtomicModelDataBlock();

};