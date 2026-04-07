#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include "data/detail/ModelDerivedState.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/LocalPotentialFitState.hpp"

namespace rhbm_gem {

ModelObject::ModelObject() :
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" },
    m_derived_state{ std::make_unique<ModelDerivedState>() },
    m_component_key_system{ std::make_unique<ComponentKeySystem>() },
    m_atom_key_system{ std::make_unique<AtomKeySystem>() },
    m_bond_key_system{ std::make_unique<BondKeySystem>() },
    m_analysis_data{ std::make_unique<ModelAnalysisData>() }
{
}

ModelObject::ModelObject(std::vector<std::unique_ptr<AtomObject>> atom_object_list) :
    m_atom_list{ std::move(atom_object_list) },
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" },
    m_derived_state{ std::make_unique<ModelDerivedState>() },
    m_component_key_system{ std::make_unique<ComponentKeySystem>() },
    m_atom_key_system{ std::make_unique<AtomKeySystem>() },
    m_bond_key_system{ std::make_unique<BondKeySystem>() },
    m_analysis_data{ std::make_unique<ModelAnalysisData>() }
{
    AttachOwnedObjects();
    InvalidateDerivedState();
    RebuildObjectIndex();
    RebuildSelection();
}

ModelObject::~ModelObject()
{

}

ModelObject::ModelObject(ModelObject && other) noexcept :
    m_atom_list{ std::move(other.m_atom_list) },
    m_bond_list{ std::move(other.m_bond_list) },
    m_selected_atom_list{ std::move(other.m_selected_atom_list) },
    m_selected_bond_list{ std::move(other.m_selected_bond_list) },
    m_key_tag{ std::move(other.m_key_tag) },
    m_pdb_id{ std::move(other.m_pdb_id) },
    m_emd_id{ std::move(other.m_emd_id) },
    m_resolution_method{ std::move(other.m_resolution_method) },
    m_resolution{ other.m_resolution },
    m_serial_id_atom_map{ std::move(other.m_serial_id_atom_map) },
    m_chain_id_list_map{ std::move(other.m_chain_id_list_map) },
    m_chemical_component_entry_map{ std::move(other.m_chemical_component_entry_map) },
    m_derived_state{ std::move(other.m_derived_state) },
    m_component_key_system{ std::move(other.m_component_key_system) },
    m_atom_key_system{ std::move(other.m_atom_key_system) },
    m_bond_key_system{ std::move(other.m_bond_key_system) },
    m_analysis_data{ std::move(other.m_analysis_data) }
{
    if (m_analysis_data == nullptr)
    {
        m_analysis_data = std::make_unique<ModelAnalysisData>();
    }
    if (m_component_key_system == nullptr)
    {
        m_component_key_system = std::make_unique<ComponentKeySystem>();
    }
    if (m_atom_key_system == nullptr)
    {
        m_atom_key_system = std::make_unique<AtomKeySystem>();
    }
    if (m_bond_key_system == nullptr)
    {
        m_bond_key_system = std::make_unique<BondKeySystem>();
    }
    AttachOwnedObjects();
    InvalidateDerivedState();
    RebuildObjectIndex();
    RebuildSelection();
}

ModelObject & ModelObject::operator=(ModelObject && other) noexcept
{
    if (this == &other)
    {
        return *this;
    }

    m_atom_list = std::move(other.m_atom_list);
    m_bond_list = std::move(other.m_bond_list);
    m_selected_atom_list = std::move(other.m_selected_atom_list);
    m_selected_bond_list = std::move(other.m_selected_bond_list);
    m_key_tag = std::move(other.m_key_tag);
    m_pdb_id = std::move(other.m_pdb_id);
    m_emd_id = std::move(other.m_emd_id);
    m_resolution_method = std::move(other.m_resolution_method);
    m_resolution = other.m_resolution;
    m_serial_id_atom_map = std::move(other.m_serial_id_atom_map);
    m_chain_id_list_map = std::move(other.m_chain_id_list_map);
    m_chemical_component_entry_map = std::move(other.m_chemical_component_entry_map);
    m_derived_state = std::move(other.m_derived_state);
    m_component_key_system = std::move(other.m_component_key_system);
    m_atom_key_system = std::move(other.m_atom_key_system);
    m_bond_key_system = std::move(other.m_bond_key_system);
    m_analysis_data = std::move(other.m_analysis_data);

    if (m_analysis_data == nullptr)
    {
        m_analysis_data = std::make_unique<ModelAnalysisData>();
    }
    if (m_component_key_system == nullptr)
    {
        m_component_key_system = std::make_unique<ComponentKeySystem>();
    }
    if (m_atom_key_system == nullptr)
    {
        m_atom_key_system = std::make_unique<AtomKeySystem>();
    }
    if (m_bond_key_system == nullptr)
    {
        m_bond_key_system = std::make_unique<BondKeySystem>();
    }
    AttachOwnedObjects();
    InvalidateDerivedState();
    RebuildObjectIndex();
    RebuildSelection();
    return *this;
}

ModelObject::ModelObject(const ModelObject & other) :
    m_key_tag{ other.m_key_tag }, m_pdb_id{ other.m_pdb_id }, m_emd_id{ other.m_emd_id },
    m_resolution_method{ other.m_resolution_method }, m_resolution{ other.m_resolution },
    m_chain_id_list_map{ other.m_chain_id_list_map },
    m_derived_state{ std::make_unique<ModelDerivedState>() },
    m_component_key_system{
        other.m_component_key_system != nullptr ?
            std::make_unique<ComponentKeySystem>(*other.m_component_key_system) : nullptr },
    m_atom_key_system{
        other.m_atom_key_system != nullptr ?
            std::make_unique<AtomKeySystem>(*other.m_atom_key_system) : nullptr },
    m_bond_key_system{
        other.m_bond_key_system != nullptr ?
            std::make_unique<BondKeySystem>(*other.m_bond_key_system) : nullptr },
    m_analysis_data{ std::make_unique<ModelAnalysisData>() }
{
    m_atom_list.reserve(other.m_atom_list.size());
    std::unordered_map<const AtomObject *, AtomObject *> atom_ptr_map;
    atom_ptr_map.reserve(other.m_atom_list.size());
    for (const auto & atom : other.m_atom_list)
    {
        auto cloned_atom{ std::make_unique<AtomObject>(*atom) };
        atom_ptr_map[atom.get()] = cloned_atom.get();
        m_atom_list.emplace_back(std::move(cloned_atom));
    }

    AttachOwnedObjects();

    m_bond_list.reserve(other.m_bond_list.size());
    std::unordered_map<const BondObject *, BondObject *> bond_ptr_map;
    bond_ptr_map.reserve(other.m_bond_list.size());
    for (const auto & bond : other.m_bond_list)
    {
        auto * atom_1{ atom_ptr_map.at(bond->GetAtomObject1()) };
        auto * atom_2{ atom_ptr_map.at(bond->GetAtomObject2()) };
        auto cloned_bond{ std::make_unique<BondObject>(atom_1, atom_2) };
        cloned_bond->SetSelectedFlag(bond->m_is_selected);
        cloned_bond->SetSpecialBondFlag(bond->GetSpecialBondFlag());
        cloned_bond->SetBondKey(bond->GetBondKey());
        cloned_bond->SetBondType(bond->GetBondType());
        cloned_bond->SetBondOrder(bond->GetBondOrder());
        bond_ptr_map[bond.get()] = cloned_bond.get();
        m_bond_list.emplace_back(std::move(cloned_bond));
    }

    AttachOwnedObjects();

    for (const auto & [component_key, entry] : other.m_chemical_component_entry_map)
    {
        m_chemical_component_entry_map[component_key] =
            std::make_unique<ChemicalComponentEntry>(*entry);
    }

    auto copy_group_entries =
        [&atom_ptr_map, &bond_ptr_map](
            const auto & source_map,
            auto set_entry)
        {
            for (const auto & [class_key, entry] : source_map)
            {
                auto cloned_entry{ std::make_unique<GroupPotentialEntry>() };
                for (const auto & [group_key, bucket] : entry.Entries())
                {
                    auto & cloned_bucket{ cloned_entry->EnsureGroup(group_key) };
                    cloned_bucket.mean = bucket.mean;
                    cloned_bucket.mdpde = bucket.mdpde;
                    cloned_bucket.prior = bucket.prior;
                    cloned_bucket.prior_variance = bucket.prior_variance;
                    cloned_bucket.alpha_g = bucket.alpha_g;
                    cloned_bucket.atom_members.reserve(bucket.atom_members.size());
                    for (auto * atom : bucket.atom_members)
                    {
                        cloned_bucket.atom_members.emplace_back(atom_ptr_map.at(atom));
                    }
                    cloned_bucket.bond_members.reserve(bucket.bond_members.size());
                    for (auto * bond : bucket.bond_members)
                    {
                        cloned_bucket.bond_members.emplace_back(bond_ptr_map.at(bond));
                    }
                }
                set_entry(class_key, std::move(cloned_entry));
            }
        };

    const auto & source_analysis_data{ ModelAnalysisData::Of(other) };
    copy_group_entries(
        source_analysis_data.Atoms().Entries(),
        [this](const std::string & class_key, std::unique_ptr<GroupPotentialEntry> entry)
        {
            m_analysis_data->Atoms().EnsureGroupEntry(class_key) = std::move(*entry);
        });
    copy_group_entries(
        source_analysis_data.Bonds().Entries(),
        [this](const std::string & class_key, std::unique_ptr<GroupPotentialEntry> entry)
        {
            m_analysis_data->Bonds().EnsureGroupEntry(class_key) = std::move(*entry);
        });

    for (const auto & atom : m_atom_list)
    {
        if (const auto * local_entry{ source_analysis_data.Atoms().FindLocalEntry(*atom) };
            local_entry != nullptr)
        {
            m_analysis_data->Atoms().SetLocalEntry(
                *atom,
                std::make_unique<LocalPotentialEntry>(*local_entry));
        }
        if (const auto * fit_state{ source_analysis_data.Atoms().FindFitState(*atom) };
            fit_state != nullptr)
        {
            m_analysis_data->Atoms().EnsureFitState(*atom) = *fit_state;
        }
    }
    for (const auto & bond : m_bond_list)
    {
        if (const auto * local_entry{ source_analysis_data.Bonds().FindLocalEntry(*bond) };
            local_entry != nullptr)
        {
            m_analysis_data->Bonds().SetLocalEntry(
                *bond,
                std::make_unique<LocalPotentialEntry>(*local_entry));
        }
        if (const auto * fit_state{ source_analysis_data.Bonds().FindFitState(*bond) };
            fit_state != nullptr)
        {
            m_analysis_data->Bonds().EnsureFitState(*bond) = *fit_state;
        }
    }

    RebuildObjectIndex();
    RebuildSelection();
}

void ModelObject::RebuildObjectIndex()
{
    m_serial_id_atom_map.clear();
    for (auto & atom : m_atom_list)
    {
        m_serial_id_atom_map[atom->GetSerialID()] = atom.get();
    }
}

void ModelObject::RebuildSelection()
{
    BuildSelectedAtomList();
    BuildSelectedBondList();
}

void ModelObject::AttachOwnedObjects()
{
    for (auto & atom : m_atom_list)
    {
        atom->SetOwnerModel(this);
    }
    for (auto & bond : m_bond_list)
    {
        bond->SetOwnerModel(this);
    }
}

void ModelObject::ApplySymmetrySelection(bool is_asymmetry)
{
    FilterSelectionFromSymmetry(is_asymmetry);
    RebuildSelection();
}

void ModelObject::InvalidateDerivedState()
{
    if (m_derived_state == nullptr)
    {
        m_derived_state = std::make_unique<ModelDerivedState>();
    }
    else
    {
        m_derived_state->Clear();
    }
}

AtomObject * ModelObject::FindAtomPtr(int serial_id) const
{
    return m_serial_id_atom_map.at(serial_id);
}

std::string ModelObject::FindComponentID(ComponentKey component_key) const
{
    return m_component_key_system->GetComponentId(component_key);
}

std::string ModelObject::FindAtomID(AtomKey atom_key) const
{
    return m_atom_key_system->GetAtomId(atom_key);
}

std::string ModelObject::FindBondID(BondKey bond_key) const
{
    return m_bond_key_system->GetBondId(bond_key);
}

std::array<float, 3> ModelObject::GetCenterOfMassPosition()
{
    return ModelDerivedState::Of(*this).GetCenterOfMassPosition(*this);
}

std::tuple<double, double> ModelObject::GetModelPositionRange(int axis)
{
    return ModelDerivedState::Of(*this).GetModelPositionRange(*this, axis);
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

bool ModelObject::HasChemicalComponentEntry(ComponentKey component_key) const
{
    return m_chemical_component_entry_map.find(component_key) != m_chemical_component_entry_map.end();
}

const ChemicalComponentEntry * ModelObject::FindChemicalComponentEntry(ComponentKey component_key) const
{
    const auto iter{ m_chemical_component_entry_map.find(component_key) };
    return iter == m_chemical_component_entry_map.end() ? nullptr : iter->second.get();
}

void ModelObject::SelectAllAtoms(bool selected)
{
    for (auto & atom : m_atom_list)
    {
        atom->SetSelectedFlag(selected);
    }
    RebuildSelection();
}

void ModelObject::SelectAllBonds(bool selected)
{
    for (auto & bond : m_bond_list)
    {
        bond->SetSelectedFlag(selected);
    }
    RebuildSelection();
}

void ModelObject::SelectAtoms(const std::function<bool(const AtomObject &)> & predicate)
{
    for (auto & atom : m_atom_list)
    {
        atom->SetSelectedFlag(predicate(*atom));
    }
    RebuildSelection();
}

void ModelObject::SelectBonds(const std::function<bool(const BondObject &)> & predicate)
{
    for (auto & bond : m_bond_list)
    {
        bond->SetSelectedFlag(predicate(*bond));
    }
    RebuildSelection();
}

void ModelObject::SetAtomSelected(int serial_id, bool selected)
{
    FindAtomPtr(serial_id)->SetSelectedFlag(selected);
    RebuildSelection();
}

void ModelObject::SetBondSelected(int atom_serial_id_1, int atom_serial_id_2, bool selected)
{
    for (auto & bond : m_bond_list)
    {
        if (bond->GetAtomSerialID1() == atom_serial_id_1 &&
            bond->GetAtomSerialID2() == atom_serial_id_2)
        {
            bond->SetSelectedFlag(selected);
            RebuildSelection();
            return;
        }
    }
    throw std::out_of_range("Bond serial pair is not available.");
}

void ModelObject::RestoreAtomSelectionBulk(const std::unordered_set<int> & selected_serial_ids)
{
    for (auto & atom : m_atom_list)
    {
        atom->SetSelectedFlag(
            selected_serial_ids.find(atom->GetSerialID()) != selected_serial_ids.end());
    }
    RebuildSelection();
}

void ModelObject::RestoreBondSelectionBulk(
    const std::set<std::pair<int, int>> & selected_serial_pairs)
{
    for (auto & bond : m_bond_list)
    {
        const auto serial_id_pair{
            std::make_pair(bond->GetAtomSerialID1(), bond->GetAtomSerialID2()) };
        bond->SetSelectedFlag(selected_serial_pairs.find(serial_id_pair) != selected_serial_pairs.end());
    }
    RebuildSelection();
}

void ModelObject::BuildSelectedAtomList()
{
    m_selected_atom_list.clear();
    m_selected_atom_list.reserve(m_atom_list.size());
    for (auto & atom : m_atom_list)
    {
        if (atom->m_is_selected == false) continue;
        m_selected_atom_list.emplace_back(atom.get());
    }
}

void ModelObject::BuildSelectedBondList()
{
    m_selected_bond_list.clear();
    m_selected_bond_list.reserve(m_bond_list.size());
    for (auto & bond : m_bond_list)
    {
        if (bond->m_is_selected == false) continue;
        m_selected_bond_list.emplace_back(bond.get());
    }
}

void ModelObject::FilterSelectionFromSymmetry(bool is_asymmetry)
{
    if (is_asymmetry == true)
    {
        return;
    }
    if (m_chain_id_list_map.empty())
    {
        Logger::Log(LogLevel::Warning,
            "ApplySymmetrySelection(): chain metadata is empty. "
            "Skip symmetry filtering and keep current atom selection.");
        return;
    }

    for (auto & atom : m_atom_list)
    {
        auto original_selection_flag{ atom->m_is_selected };
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

    for (auto & bond : m_bond_list)
    {
        auto original_selection_flag{ bond->m_is_selected };
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
