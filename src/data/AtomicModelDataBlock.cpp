#include "AtomicModelDataBlock.hpp"
#include "AtomObject.hpp"
#include "GlobalEnumClass.hpp"
#include "AtomicInfoHelper.hpp"

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

void AtomicModelDataBlock::AddChainEntityType(const std::string & entity_id, Entity entity)
{
    m_chain_entity_type_map[entity_id] = entity;
}

void AtomicModelDataBlock::AddSheetStrands(
    const std::string & sheet_id, int strands_size)
{
    m_struct_sheet_strand_map[sheet_id] = strands_size;
}

void AtomicModelDataBlock::AddSheetRange(
    const std::string & composite_sheet_id, const std::array<std::string, 4> & range)
{
    m_struct_sheet_range_map[composite_sheet_id] = range;
}

void AtomicModelDataBlock::AddHelixRange(
    const std::string & helix_id, const std::array<std::string, 5> & range)
{
    m_struct_helix_range_map[helix_id] = range;
}

void AtomicModelDataBlock::AddElementType(const Element & element)
{
    m_element_type_list.emplace_back(element);
}

void AtomicModelDataBlock::SetStructureInfo(AtomObject * atom_object)
{
    auto chain_id{ atom_object->GetChainID() };
    auto residue_id{ atom_object->GetResidueID() };

    for (auto & [helix_id, range] : m_struct_helix_range_map)
    {
        if (chain_id == range.at(0) || chain_id == range.at(2))
        {
            auto beg{ std::stoi(range.at(1)) };
            auto end{ std::stoi(range.at(3)) };
            if (residue_id >= beg && residue_id <= end)
            {
                atom_object->SetStructure(AtomicInfoHelper::GetStructureFromString(range.at(4)));
                return;
            }
        }
    }

    for (auto & [composite_sheet_id, range] : m_struct_sheet_range_map)
    {
        if (chain_id == range.at(0) || chain_id == range.at(2))
        {
            auto beg{ std::stoi(range.at(1)) };
            auto end{ std::stoi(range.at(3)) };
            if (residue_id >= beg && residue_id <= end)
            {
                atom_object->SetStructure(Structure::STRN);
                return;
            }
        }
    }

    atom_object->SetStructure(Structure::FREE);
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