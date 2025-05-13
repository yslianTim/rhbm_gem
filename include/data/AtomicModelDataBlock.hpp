#pragma once

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>

enum class Element : uint16_t;

class AtomObject;

class AtomicModelDataBlock
{
    std::string m_map_id, m_model_id;
    std::string m_resolution, m_resolution_method;

    std::vector<std::unique_ptr<AtomObject>> m_atom_object_list;
    std::vector<Element> m_element_type_list;

public:
    AtomicModelDataBlock(void);
    ~AtomicModelDataBlock();

    void AddAtomObject(std::unique_ptr<AtomObject> atom_object);
    void AddElementType(const Element & element);
    void SetPdbID(const std::string & label) { m_model_id = label; }
    void SetEmdID(const std::string & label) { m_map_id = label; }
    void SetResolution(const std::string & value) { m_resolution = value; }
    void SetResolutionMethod(const std::string & value) { m_resolution_method = value; }

    std::vector<std::unique_ptr<AtomObject>> GetAtomObjectList(void);
    std::string GetPdbID(void) const;
    std::string GetEmdID(void) const;
    double GetResolution(void) const;
    std::string GetResolutionMethod(void) const;

};