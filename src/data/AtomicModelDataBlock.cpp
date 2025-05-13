#include "AtomicModelDataBlock.hpp"
#include "AtomObject.hpp"
#include "GlobalEnumClass.hpp"

AtomicModelDataBlock::AtomicModelDataBlock(void)
{

}

AtomicModelDataBlock::~AtomicModelDataBlock()
{

}

void AtomicModelDataBlock::AddAtomObject(
    int model_number, std::unique_ptr<AtomObject> atom_object)
{
    m_atom_object_list_map[model_number].emplace_back(std::move(atom_object));
}

void AtomicModelDataBlock::AddElementType(const Element & element)
{
    m_element_type_list.emplace_back(element);
}

std::vector<std::unique_ptr<AtomObject>> AtomicModelDataBlock::GetAtomObjectList(int model_number)
{
    return std::move(m_atom_object_list_map.at(model_number));
}

std::string AtomicModelDataBlock::GetPdbID(void) const
{
    return m_model_id;
}

std::string AtomicModelDataBlock::GetEmdID(void) const
{
    return m_map_id;
}

double AtomicModelDataBlock::GetResolution(void) const
{
    return std::stod(m_resolution);
}

std::string AtomicModelDataBlock::GetResolutionMethod(void) const
{
    return m_resolution_method;
}