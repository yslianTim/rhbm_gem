#include "AtomicModelDataBlock.hpp"
#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "GlobalEnumClass.hpp"
#include "ChemicalDataHelper.hpp"
#include "ChemicalComponentEntry.hpp"
#include "Logger.hpp"

AtomicModelDataBlock::AtomicModelDataBlock(void) :
    m_component_key_system{ std::make_unique<ComponentKeySystem>() },
    m_atom_key_system{ std::make_unique<AtomKeySystem>() },
    m_bond_key_system{ std::make_unique<BondKeySystem>() }
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
    auto sequence_id_str{ std::to_string(atom_object->GetSequenceID()) };
    sequence_id_str = (sequence_id_str == "-1") ? "." : sequence_id_str;
    auto atom_tuple_key{ std::make_tuple(
        atom_object->GetChainID(),
        atom_object->GetComponentID(),
        sequence_id_str,
        atom_object->GetAtomID())
    };
    m_atom_object_in_tuple_map[model_number][atom_tuple_key] = atom_object.get();
    m_atom_object_list_map[model_number].emplace_back(std::move(atom_object));
}

void AtomicModelDataBlock::AddBondObject(
    std::unique_ptr<BondObject> bond_object)
{
    m_bond_object_list.emplace_back(std::move(bond_object));
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

void AtomicModelDataBlock::AddChemicalComponentEntry(
    ComponentKey comp_key, std::unique_ptr<ChemicalComponentEntry> entry)
{
    m_chemical_component_entry_map[comp_key] = std::move(entry);
}

void AtomicModelDataBlock::AddComponentAtomEntry(
    ComponentKey comp_key,
    AtomKey atom_key,
    const ComponentAtomEntry & atom_entry)
{
    Logger::Log(LogLevel::Debug, "AtomicModelDataBlock::AddComponentAtomEntry() called");
    if (HasChemicalComponentEntry(comp_key) == false)
    {
        Logger::Log(LogLevel::Warning,
            "Chemical component key " + std::to_string(comp_key) +
            " not found in chemical component map.");
        return;
    }
    m_chemical_component_entry_map.at(comp_key)->AddComponentAtomEntry(atom_key, atom_entry);
}

void AtomicModelDataBlock::AddComponentBondEntry(
    ComponentKey comp_key,
    BondKey bond_key,
    const ComponentBondEntry & bond_entry)
{
    Logger::Log(LogLevel::Debug, "AtomicModelDataBlock::AddComponentBondEntry() called");
    if (HasChemicalComponentEntry(comp_key) == false)
    {
        Logger::Log(LogLevel::Warning,
            "Chemical component key " + std::to_string(comp_key) +
            " not found in chemical component map.");
        return;
    }
    m_chemical_component_entry_map.at(comp_key)->AddComponentBondEntry(bond_key, bond_entry);
}

void AtomicModelDataBlock::SetStructureInfo(AtomObject * atom_object)
{
    auto chain_id{ atom_object->GetChainID() };
    auto sequence_id{ atom_object->GetSequenceID() };

    for (auto & [helix_id, range] : m_struct_helix_range_map)
    {
        if (chain_id == range.at(0) || chain_id == range.at(2))
        {
            auto beg{ std::stoi(range.at(1)) };
            auto end{ std::stoi(range.at(3)) };
            if (sequence_id >= beg && sequence_id <= end)
            {
                atom_object->SetStructure(ChemicalDataHelper::GetStructureFromString(range.at(4)));
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
            if (sequence_id >= beg && sequence_id <= end)
            {
                atom_object->SetStructure(Structure::SHEET);
                return;
            }
        }
    }

    atom_object->SetStructure(Structure::FREE);
}

std::vector<std::unique_ptr<AtomObject>> AtomicModelDataBlock::MoveAtomObjectList(int model_number)
{
    return std::move(m_atom_object_list_map.at(model_number));
}

std::vector<std::unique_ptr<BondObject>> AtomicModelDataBlock::MoveBondObjectList(void)
{
    if (m_bond_object_list.empty())
    {
        return {};
    }
    return std::move(m_bond_object_list);
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
        Logger::Log(LogLevel::Error,
            "AtomicModelDataBlock::GetResolution()" + std::string(e.what()));
        Logger::Log(LogLevel::Error,
            "Invalid resolution value: " + m_resolution
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

const std::unordered_map<int, std::vector<std::unique_ptr<AtomObject>>> &
AtomicModelDataBlock::GetAtomObjectMap(void) const
{
    return m_atom_object_list_map;
}

const std::vector<std::unique_ptr<BondObject>> & AtomicModelDataBlock::GetBondObjectList(void) const
{
    return m_bond_object_list;
}

AtomObject * AtomicModelDataBlock::GetAtomObjectPtrInTuple(
    int model_number,
    const std::string & chain_id,
    const std::string & comp_id,
    const std::string & seq_id,
    const std::string & atom_id) const
{
    if (m_atom_object_in_tuple_map.find(model_number) == m_atom_object_in_tuple_map.end())
    {
        Logger::Log(LogLevel::Warning,
            "Model number " + std::to_string(model_number) +
            " not found in atom object map.");
        return nullptr;
    }
    auto atom_tuple_key{ std::make_tuple(chain_id, comp_id, seq_id, atom_id) };
    if (m_atom_object_in_tuple_map.at(model_number).find(atom_tuple_key) ==
        m_atom_object_in_tuple_map.at(model_number).end())
    {
        Logger::Log(LogLevel::Warning,
            "Atom object ("
            + chain_id + ", "
            + comp_id + ", "
            + seq_id + ", "
            + atom_id + ") not found in atom object map.");
        return nullptr;
    }
    return m_atom_object_in_tuple_map.at(model_number).at(atom_tuple_key);
}

ChemicalComponentEntry * AtomicModelDataBlock::GetChemicalComponentEntryPtr(ComponentKey key) const
{
    if (HasChemicalComponentEntry(key) == false)
    {
        Logger::Log(LogLevel::Warning,
            "Chemical component key " + std::to_string(key) +
            " not found in chemical component map.");
        return nullptr;
    }
    return m_chemical_component_entry_map.at(key).get();
}

const ComponentBondEntry * AtomicModelDataBlock::GetComponentBondEntryPtr(
    ComponentKey comp_key, BondKey bond_key) const
{
    auto comp_entry_ptr{ GetChemicalComponentEntryPtr(comp_key) };
    if (comp_entry_ptr == nullptr || comp_entry_ptr->HasComponentBondEntry(bond_key) == false)
    {
        Logger::Log(LogLevel::Warning,
            "Component bond entry (comp_key: "+ std::to_string(comp_key)
            + ", bond_key: "+ std::to_string(bond_key)
            + ") not found in chemical component map.");
        return nullptr;
    }
    return comp_entry_ptr->GetComponentBondEntryPtr(bond_key);
}

std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> &
AtomicModelDataBlock::GetChemicalComponentEntryMap(void)
{
    return m_chemical_component_entry_map;
}

bool AtomicModelDataBlock::HasChemicalComponentEntry(ComponentKey comp_key) const
{
    return (m_chemical_component_entry_map.find(comp_key) != m_chemical_component_entry_map.end());
}

bool AtomicModelDataBlock::HasComponentBondEntry(
    ComponentKey comp_key, BondKey bond_key) const
{
    if (HasChemicalComponentEntry(comp_key) == false)
    {
        return false;
    }
    return m_chemical_component_entry_map.at(comp_key)->HasComponentBondEntry(bond_key);
}
