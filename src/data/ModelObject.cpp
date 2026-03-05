#include "ModelObject.hpp"
#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "DataObjectVisitorBase.hpp"
#include "GroupPotentialEntry.hpp"
#include "ChemicalComponentEntry.hpp"
#include "KDTreeAlgorithm.hpp"
#include "ArrayStats.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"

#include <stdexcept>

namespace rhbm_gem {

ModelObject::ModelObject() :
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" }, m_kd_tree_root{ nullptr },
    m_component_key_system{ std::make_unique<ComponentKeySystem>() },
    m_atom_key_system{ std::make_unique<AtomKeySystem>() },
    m_bond_key_system{ std::make_unique<BondKeySystem>() }
{
}

ModelObject::ModelObject(std::vector<std::unique_ptr<AtomObject>> atom_object_list) :
    m_atom_list{ std::move(atom_object_list) },
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" }, m_kd_tree_root{ nullptr },
    m_component_key_system{ std::make_unique<ComponentKeySystem>() },
    m_atom_key_system{ std::make_unique<AtomKeySystem>() },
    m_bond_key_system{ std::make_unique<BondKeySystem>() }
{
    Update();
}

ModelObject::~ModelObject()
{

}

ModelObject::ModelObject(const ModelObject & other) :
    m_key_tag{ other.m_key_tag }, m_pdb_id{ other.m_pdb_id }, m_emd_id{ other.m_emd_id },
    m_resolution_method{ other.m_resolution_method }, m_resolution{ other.m_resolution },
    m_kd_tree_root{ nullptr },
    m_component_key_system{ std::make_unique<ComponentKeySystem>() },
    m_atom_key_system{ std::make_unique<AtomKeySystem>() },
    m_bond_key_system{ std::make_unique<BondKeySystem>() }
{
    m_atom_list.clear();
    m_atom_list.reserve(other.m_atom_list.size());
    for (auto & atom : other.m_atom_list)
    {
        m_atom_list.emplace_back(atom->AtomObjectClone());
    }
    Update();
}

std::unique_ptr<DataObjectBase> ModelObject::Clone() const
{
    return std::make_unique<ModelObject>(*this);
}

void ModelObject::Display() const
{
    Logger::Log(LogLevel::Info, "ModelObject Display: " + GetKeyTag());
    Logger::Log(LogLevel::Info,
        " - PDB ID = " + m_pdb_id + "\n"
        " - EMD ID = " + m_emd_id + "\n"
        " - Map Resolution = " + StringHelper::ToStringWithPrecision(m_resolution, 2)
            + " A (" + m_resolution_method + ")\n"
        " - #Unique Entities = " + std::to_string(m_chain_id_list_map.size()) + "\n"
        " - #Atoms = "+ std::to_string(m_atom_list.size()) +"\n"
        " - #Bonds = "+ std::to_string(m_bond_list.size()) +"\n"
        " - #Unique Components = " + std::to_string(m_chemical_component_entry_map.size()) + "\n"
    );

    for (auto & [component_key, entry]: m_chemical_component_entry_map)
    {
        Logger::Log(LogLevel::Debug,
            "   - Component ID = " + entry->GetComponentId() + " "
            " (#Atoms = "+ std::to_string(entry->GetComponentAtomEntryMap().size()) +
            ", #Bonds = "+ std::to_string(entry->GetComponentBondEntryMap().size()) +")"
        );
    }
}

void ModelObject::Update()
{
    for (auto & atom : m_atom_list)
    {
        m_serial_id_atom_map[atom->GetSerialID()] = atom.get();
    }
    BuildSelectedAtomList();
    BuildSelectedBondList();
}

void ModelObject::Accept(DataObjectVisitorBase * visitor)
{
    Accept(visitor, ModelVisitMode::AtomsThenSelf);
}

void ModelObject::Accept(DataObjectVisitorBase * visitor, ModelVisitMode mode)
{
    if (visitor == nullptr)
    {
        throw std::invalid_argument("ModelObject::Accept(): visitor is null.");
    }

    switch (mode)
    {
    case ModelVisitMode::AtomsThenSelf:
        for (const auto & atom : m_atom_list) atom->Accept(visitor);
        visitor->VisitModelObject(this);
        break;
    case ModelVisitMode::BondsThenSelf:
        for (const auto & bond : m_bond_list) bond->Accept(visitor);
        visitor->VisitModelObject(this);
        break;
    case ModelVisitMode::AtomsAndBondsThenSelf:
        for (const auto & atom : m_atom_list) atom->Accept(visitor);
        for (const auto & bond : m_bond_list) bond->Accept(visitor);
        visitor->VisitModelObject(this);
        break;
    case ModelVisitMode::SelfOnly:
        visitor->VisitModelObject(this);
        break;
    }
}

void ModelObject::AddAtom(std::unique_ptr<AtomObject> atom)
{
    m_atom_list.emplace_back(std::move(atom));
    m_serial_id_atom_map[atom->GetSerialID()] = atom.get();
}

void ModelObject::AddBond(std::unique_ptr<BondObject> bond)
{
    m_bond_list.emplace_back(std::move(bond));
}

void ModelObject::SetAtomList(std::vector<std::unique_ptr<AtomObject>> atom_list)
{
    m_atom_list = std::move(atom_list);
    Update();
}

void ModelObject::SetBondList(std::vector<std::unique_ptr<BondObject>> bond_list)
{
    m_bond_list = std::move(bond_list);
    BuildSelectedBondList();
}

void ModelObject::AddAtomGroupPotentialEntry(
    const std::string & class_key, std::unique_ptr<GroupPotentialEntry> & entry)
{
    m_atom_group_potential_entry_map[class_key] = std::move(entry);
}

void ModelObject::AddBondGroupPotentialEntry(
    const std::string & class_key, std::unique_ptr<GroupPotentialEntry> & entry)
{
    m_bond_group_potential_entry_map[class_key] = std::move(entry);
}

void ModelObject::AddChemicalComponentEntry(
    ComponentKey key, std::unique_ptr<ChemicalComponentEntry> entry)
{
    m_chemical_component_entry_map[key] = std::move(entry);
}

void ModelObject::SetChemicalComponentEntryMap(
    std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> & entry_map)
{
    m_chemical_component_entry_map = std::move(entry_map);
}

void ModelObject::SetComponentKeySystem(std::unique_ptr<ComponentKeySystem> component_key_system)
{
    if (m_component_key_system != nullptr) m_component_key_system = nullptr;
    m_component_key_system = std::move(component_key_system);
}

void ModelObject::SetAtomKeySystem(std::unique_ptr<AtomKeySystem> atom_key_system)
{
    if (m_atom_key_system != nullptr) m_atom_key_system = nullptr;
    m_atom_key_system = std::move(atom_key_system);
}

void ModelObject::SetBondKeySystem(std::unique_ptr<BondKeySystem> bond_key_system)
{
    if (m_bond_key_system != nullptr) m_bond_key_system = nullptr;
    m_bond_key_system = std::move(bond_key_system);
}

void ModelObject::BuildKDTreeRoot()
{
    if (m_kd_tree_root != nullptr) return;
    std::vector<AtomObject *> atom_ptr_list;
    atom_ptr_list.reserve(m_atom_list.size());
    for (auto & atom : m_atom_list)
    {
        atom_ptr_list.emplace_back(atom.get());
    }
    m_kd_tree_root = KDTreeAlgorithm<AtomObject>::BuildKDTree(atom_ptr_list, 0);
}

std::array<float, 3> ModelObject::GetCenterOfMassPosition()
{
    if (m_center_of_mass_position != nullptr)
    {
        return *(m_center_of_mass_position);
    }
    std::vector<float> pos_x, pos_y, pos_z;
    pos_x.reserve(m_atom_list.size());
    pos_y.reserve(m_atom_list.size());
    pos_z.reserve(m_atom_list.size());
    for (auto & atom : m_atom_list)
    {
        const auto & pos{ atom->GetPositionRef() };
        pos_x.emplace_back(pos.at(0));
        pos_y.emplace_back(pos.at(1));
        pos_z.emplace_back(pos.at(2));
    }
    m_center_of_mass_position = std::make_unique<std::array<float, 3>>();
    m_center_of_mass_position->at(0) = ArrayStats<float>::ComputeMean(pos_x.data(), pos_x.size());
    m_center_of_mass_position->at(1) = ArrayStats<float>::ComputeMean(pos_y.data(), pos_y.size());
    m_center_of_mass_position->at(2) = ArrayStats<float>::ComputeMean(pos_z.data(), pos_z.size());
    return *(m_center_of_mass_position);
}

std::tuple<double, double> ModelObject::GetModelPositionRange(int axis)
{
    if (m_model_position_range[axis] != nullptr)
    {
        return *(m_model_position_range[axis]);
    }
    std::vector<double> position_list;
    position_list.reserve(m_atom_list.size());
    for (auto & atom : m_atom_list)
    {
        position_list.emplace_back(atom->GetPosition().at(static_cast<size_t>(axis)));
    }
    m_model_position_range[axis] = std::make_unique<std::tuple<double, double>>(
        ArrayStats<double>::ComputeScalingRangeTuple(position_list, 0.0));
    return *(m_model_position_range[axis]);
}

double ModelObject::GetModelPosition(int axis, double normalized_pos)
{
    auto range_tuple{ GetModelPositionRange(axis) };
    auto pos_min{ std::get<0>(range_tuple) };
    auto pos_max{ std::get<1>(range_tuple) };
    return pos_min + normalized_pos * (pos_max - pos_min);
}

double ModelObject::GetModelLength(int axis)
{
    auto range_tuple{ GetModelPositionRange(axis) };
    return std::get<1>(range_tuple) - std::get<0>(range_tuple);
}

GroupPotentialEntry * ModelObject::GetAtomGroupPotentialEntry(const std::string & class_key) const
{
    return m_atom_group_potential_entry_map.at(class_key).get();
}

GroupPotentialEntry * ModelObject::GetBondGroupPotentialEntry(const std::string & class_key) const
{
    return m_bond_group_potential_entry_map.at(class_key).get();
}

ChemicalComponentEntry * ModelObject::GetChemicalComponentEntry(ComponentKey key) const
{
    return m_chemical_component_entry_map.at(key).get();
}

const std::map<int, AtomObject *> &
ModelObject::GetSerialIDAtomMap() const
{
    return m_serial_id_atom_map;
}

const std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> &
ModelObject::GetAtomGroupPotentialEntryMap() const
{
    return m_atom_group_potential_entry_map;
}

const std::unordered_map<std::string, std::unique_ptr<GroupPotentialEntry>> &
ModelObject::GetBondGroupPotentialEntryMap() const
{
    return m_bond_group_potential_entry_map;
}

const std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> &
ModelObject::GetChemicalComponentEntryMap() const
{
    return m_chemical_component_entry_map;
}

void ModelObject::BuildSelectedAtomList()
{
    m_selected_atom_list.clear();
    m_selected_atom_list.reserve(m_atom_list.size());
    for (auto & atom : m_atom_list)
    {
        if (atom->GetSelectedFlag() == false) continue;
        m_selected_atom_list.emplace_back(atom.get());
    }
}

void ModelObject::BuildSelectedBondList()
{
    m_selected_bond_list.clear();
    m_selected_bond_list.reserve(m_bond_list.size());
    for (auto & bond : m_bond_list)
    {
        if (bond->GetSelectedFlag() == false) continue;
        m_selected_bond_list.emplace_back(bond.get());
    }
}

void ModelObject::FilterAtomFromSymmetry(bool is_asymmetry)
{
    if (is_asymmetry == true)
    {
        return;
    }
    if (m_chain_id_list_map.empty())
    {
        Logger::Log(LogLevel::Warning,
            "FilterAtomFromSymmetry(): chain metadata is empty. "
            "Skip symmetry filtering and keep current atom selection.");
        return;
    }

    for (auto & atom : m_atom_list)
    {
        auto original_selection_flag{ atom->GetSelectedFlag() };
        auto chain_id{ atom->GetChainID() };
        bool in_candidate_chain{ false };
        for (auto & [entity_id, chain_id_list] : m_chain_id_list_map)
        {
            (void)entity_id;
            if (chain_id_list.empty()) continue;
            if (chain_id == chain_id_list.front())
            {
                atom->SetSelectedFlag(original_selection_flag);
                in_candidate_chain = true;
                break;
            }
        }
        if (in_candidate_chain == false) atom->SetSelectedFlag(false);
    }
}

void ModelObject::FilterBondFromSymmetry(bool is_asymmetry)
{
    if (is_asymmetry == true)
    {
        return;
    }
    if (m_chain_id_list_map.empty())
    {
        Logger::Log(LogLevel::Warning,
            "FilterBondFromSymmetry(): chain metadata is empty. "
            "Skip symmetry filtering and keep current bond selection.");
        return;
    }

    for (auto & bond : m_bond_list)
    {
        auto original_selection_flag{ bond->GetSelectedFlag() };
        auto chain_id{ bond->GetAtomObject1()->GetChainID() };
        bool in_candidate_chain{ false };
        for (auto & [entity_id, chain_id_list] : m_chain_id_list_map)
        {
            (void)entity_id;
            if (chain_id_list.empty()) continue;
            if (chain_id == chain_id_list.front())
            {
                bond->SetSelectedFlag(original_selection_flag);
                in_candidate_chain = true;
                break;
            }
        }
        if (in_candidate_chain == false) bond->SetSelectedFlag(false);
    }
}

std::vector<ComponentKey> ModelObject::GetComponentKeyList() const
{
    std::vector<ComponentKey> component_key_list;
    component_key_list.reserve(m_chemical_component_entry_map.size());
    for (const auto & [component_key, entry] : m_chemical_component_entry_map)
    {
        component_key_list.emplace_back(component_key);
    }
    return component_key_list;
}

bool ModelObject::HasStandardRNAComponent() const
{
    for (const auto & [component_key, entry] : m_chemical_component_entry_map)
    {
        const auto component_id{ entry->GetComponentId() };
        if (component_id == "A" ||
            component_id == "C" ||
            component_id == "G" ||
            component_id == "U")
        {
            return true;
        }
    }
    return false;
}

bool ModelObject::HasStandardDNAComponent() const
{
    for (const auto & [component_key, entry] : m_chemical_component_entry_map)
    {
        const auto component_id{ entry->GetComponentId() };
        if (component_id == "DA" ||
            component_id == "DC" ||
            component_id == "DG" ||
            component_id == "DT")
        {
            return true;
        }
    }
    return false;
}

} // namespace rhbm_gem
