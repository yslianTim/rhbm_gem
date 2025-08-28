#include "AtomicModelDataBlock.hpp"
#include "AtomObject.hpp"
#include "GlobalEnumClass.hpp"
#include "AtomicInfoHelper.hpp"
#include "Logger.hpp"

AtomicModelDataBlock::AtomicModelDataBlock(void)
{
    Logger::Log(LogLevel::Debug, "AtomicModelDataBlock::AtomicModelDataBlock() called");
}

AtomicModelDataBlock::~AtomicModelDataBlock()
{
    Logger::Log(LogLevel::Debug, "AtomicModelDataBlock::~AtomicModelDataBlock() called");
}

void AtomicModelDataBlock::AddAtomObject(
    int model_number, std::unique_ptr<AtomObject> atom_object)
{
    m_atom_object_list_map[model_number].emplace_back(std::move(atom_object));
}

void AtomicModelDataBlock::AddEntityTypeInEntityMap(
    const std::string & entity_id, Entity entity)
{
    m_entity_type_map[entity_id] = entity;
    m_entity_id_list_map[entity].emplace_back(entity_id);
}

void AtomicModelDataBlock::AddChainIDInEntityMap(
    const std::string & entity_id, const std::string & chain_id)
{
    m_chain_id_list_map[entity_id].emplace_back(chain_id);
}

void AtomicModelDataBlock::AddMoleculesSizeInEntityMap(
    const std::string & entity_id, int molecules_size)
{
    m_molecules_size_map[entity_id] = molecules_size;
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

void AtomicModelDataBlock::AddChemicalComponentFormula(
    const std::string & comp_id, const std::string & formula)
{
    m_chemical_component_formula_map[comp_id] = formula;
}

void AtomicModelDataBlock::AddChemicalComponentName(
    const std::string & comp_id, const std::string & name)
{
    m_chemical_component_name_map[comp_id] = name;
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
                atom_object->SetStructure(Structure::SHEET);
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

const std::string & AtomicModelDataBlock::GetPdbID(void) const
{
    return m_model_id;
}

const std::string & AtomicModelDataBlock::GetEmdID(void) const
{
    return m_map_id;
}

double AtomicModelDataBlock::GetResolution(void) const
{
    double resolution_value;
    try
    {
        resolution_value = std::stod(m_resolution);
    }
    catch (const std::exception & e)
    {
        Logger::Log(LogLevel::Error, "AtomicModelDataBlock::GetResolution()" + std::string(e.what()));
        Logger::Log(LogLevel::Error, "Invalid resolution value: " + m_resolution
                    + ". Setting resolution to -1.0");
        resolution_value = -1.0;
    }
    return resolution_value;
}

const std::string & AtomicModelDataBlock::GetResolutionMethod(void) const
{
    return m_resolution_method;
}

const std::vector<Element> & AtomicModelDataBlock::GetElementTypeList(void) const
{
    return m_element_type_list;
}

const std::unordered_map<std::string, Entity> &
AtomicModelDataBlock::GetEntityTypeMap(void) const
{
    return m_entity_type_map;
}

const std::unordered_map<std::string, int> &
AtomicModelDataBlock::GetMoleculesSizeMap(void) const
{
    return m_molecules_size_map;
}

const std::unordered_map<Entity, std::vector<std::string>> &
AtomicModelDataBlock::GetEntityIDListMap(void) const
{
    return m_entity_id_list_map;
}

const std::unordered_map<std::string, std::vector<std::string>> &
AtomicModelDataBlock::GetChainIDListMap(void) const
{
    return m_chain_id_list_map;
}
