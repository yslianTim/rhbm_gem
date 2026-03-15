#include "internal/object/AtomicModelDataBlock.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <algorithm>
#include <unordered_set>

namespace rhbm_gem {

AtomicModelDataBlock::AtomicModelDataBlock() :
    m_component_key_system{ std::make_unique<ComponentKeySystem>() },
    m_atom_key_system{ std::make_unique<AtomKeySystem>() },
    m_bond_key_system{ std::make_unique<BondKeySystem>() }
{
}

AtomicModelDataBlock::~AtomicModelDataBlock()
{
}

void AtomicModelDataBlock::AddAtomObject(
    int model_number,
    std::unique_ptr<AtomObject> atom_object,
    const std::string & raw_sequence_id_token)
{
    std::string sequence_id_str{ raw_sequence_id_token };
    if (sequence_id_str.empty())
    {
        sequence_id_str = std::to_string(atom_object->GetSequenceID());
        sequence_id_str = (sequence_id_str == "-1") ? "." : sequence_id_str;
    }
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
    const std::string & composite_sheet_id, const SheetRange & range)
{
    m_struct_sheet_range_map[composite_sheet_id] = range;
}

void AtomicModelDataBlock::AddHelixRange(
    const std::string & helix_id, const HelixRange & range)
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
        (void)helix_id;
        if (chain_id == range.chain_id_beg || chain_id == range.chain_id_end)
        {
            if (sequence_id >= range.seq_id_beg && sequence_id <= range.seq_id_end)
            {
                atom_object->SetStructure(ChemicalDataHelper::GetStructureFromString(range.conf_type));
                return;
            }
        }
    }

    for (auto & [composite_sheet_id, range] : m_struct_sheet_range_map)
    {
        (void)composite_sheet_id;
        if (chain_id == range.chain_id_beg || chain_id == range.chain_id_end)
        {
            if (sequence_id >= range.seq_id_beg && sequence_id <= range.seq_id_end)
            {
                atom_object->SetStructure(Structure::SHEET);
                return;
            }
        }
    }

    atom_object->SetStructure(Structure::FREE);
}

std::unique_ptr<ModelObject> AtomicModelDataBlock::TakeModelObject(int preferred_model_number)
{
    auto model_number_list{ GetModelNumberList() };
    if (model_number_list.empty())
    {
        throw std::runtime_error("No atom model found in the input model file.");
    }

    int selected_model_number{ preferred_model_number };
    if (HasModelNumber(selected_model_number) == false)
    {
        selected_model_number = model_number_list.front();
        Logger::Log(
            LogLevel::Warning,
            "Model " + std::to_string(preferred_model_number) + " not found. Fallback to model "
                + std::to_string(selected_model_number) + ".");
    }

    auto model_object{
        std::make_unique<ModelObject>(MoveAtomObjectList(selected_model_number))
    };

    std::unordered_set<const AtomObject *> selected_atom_set;
    selected_atom_set.reserve(model_object->GetAtomList().size());
    for (const auto & atom : model_object->GetAtomList())
    {
        selected_atom_set.insert(atom.get());
    }

    auto bond_list{ MoveBondObjectList() };
    std::vector<std::unique_ptr<BondObject>> filtered_bond_list;
    filtered_bond_list.reserve(bond_list.size());
    for (auto & bond : bond_list)
    {
        if (bond == nullptr) continue;
        auto atom_1{ bond->GetAtomObject1() };
        auto atom_2{ bond->GetAtomObject2() };
        if (selected_atom_set.find(atom_1) == selected_atom_set.end()) continue;
        if (selected_atom_set.find(atom_2) == selected_atom_set.end()) continue;
        filtered_bond_list.emplace_back(std::move(bond));
    }

    model_object->SetPdbID(GetPdbID());
    model_object->SetEmdID(GetEmdID());
    model_object->SetResolution(GetResolution());
    model_object->SetResolutionMethod(GetResolutionMethod());
    model_object->SetChainIDListMap(GetChainIDListMap());
    model_object->SetChemicalComponentEntryMap(GetChemicalComponentEntryMap());
    model_object->SetComponentKeySystem(MoveComponentKeySystem());
    model_object->SetAtomKeySystem(MoveAtomKeySystem());
    model_object->SetBondKeySystem(MoveBondKeySystem());
    model_object->SetBondList(std::move(filtered_bond_list));
    return model_object;
}

std::vector<std::unique_ptr<AtomObject>> AtomicModelDataBlock::MoveAtomObjectList(int model_number)
{
    auto iter{ m_atom_object_list_map.find(model_number) };
    if (iter == m_atom_object_list_map.end())
    {
        Logger::Log(LogLevel::Warning,
            "AtomicModelDataBlock::MoveAtomObjectList() - Model number "
            + std::to_string(model_number) + " not found.");
        return {};
    }
    return std::move(iter->second);
}

std::vector<std::unique_ptr<BondObject>> AtomicModelDataBlock::MoveBondObjectList()
{
    if (m_bond_object_list.empty())
    {
        return {};
    }
    return std::move(m_bond_object_list);
}

const std::string & AtomicModelDataBlock::GetPdbID() const
{
    return m_model_id;
}

const std::string & AtomicModelDataBlock::GetEmdID() const
{
    return m_map_id;
}

double AtomicModelDataBlock::GetResolution() const
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

const std::string & AtomicModelDataBlock::GetResolutionMethod() const
{
    return m_resolution_method;
}

const std::vector<Element> & AtomicModelDataBlock::GetElementTypeList() const
{
    return m_element_type_list;
}

const std::unordered_map<std::string, Entity> &
AtomicModelDataBlock::GetEntityTypeMap() const
{
    return m_entity_type_map;
}

const std::unordered_map<std::string, int> &
AtomicModelDataBlock::GetMoleculesSizeMap() const
{
    return m_molecules_size_map;
}

const std::unordered_map<Entity, std::vector<std::string>> &
AtomicModelDataBlock::GetEntityIDListMap() const
{
    return m_entity_id_list_map;
}

const std::unordered_map<std::string, std::vector<std::string>> &
AtomicModelDataBlock::GetChainIDListMap() const
{
    return m_chain_id_list_map;
}

const std::unordered_map<int, std::vector<std::unique_ptr<AtomObject>>> &
AtomicModelDataBlock::GetAtomObjectMap() const
{
    return m_atom_object_list_map;
}

const std::vector<std::unique_ptr<BondObject>> & AtomicModelDataBlock::GetBondObjectList() const
{
    return m_bond_object_list;
}

std::vector<int> AtomicModelDataBlock::GetModelNumberList() const
{
    std::vector<int> model_number_list;
    model_number_list.reserve(m_atom_object_list_map.size());
    for (const auto & [model_number, _] : m_atom_object_list_map)
    {
        (void)_;
        model_number_list.emplace_back(model_number);
    }
    std::sort(model_number_list.begin(), model_number_list.end());
    return model_number_list;
}

bool AtomicModelDataBlock::HasModelNumber(int model_number) const
{
    return m_atom_object_list_map.find(model_number) != m_atom_object_list_map.end();
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

AtomObject * AtomicModelDataBlock::GetAtomObjectPtrInAnyModel(
    const std::string & chain_id,
    const std::string & comp_id,
    const std::string & seq_id,
    const std::string & atom_id,
    int * model_number) const
{
    auto atom_tuple_key{ std::make_tuple(chain_id, comp_id, seq_id, atom_id) };
    auto model_number_list{ GetModelNumberList() };
    for (const auto current_model_number : model_number_list)
    {
        auto model_iter{ m_atom_object_in_tuple_map.find(current_model_number) };
        if (model_iter == m_atom_object_in_tuple_map.end()) continue;
        auto atom_iter{ model_iter->second.find(atom_tuple_key) };
        if (atom_iter == model_iter->second.end()) continue;
        if (model_number != nullptr)
        {
            *model_number = current_model_number;
        }
        return atom_iter->second;
    }
    return nullptr;
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
AtomicModelDataBlock::GetChemicalComponentEntryMap()
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

} // namespace rhbm_gem
