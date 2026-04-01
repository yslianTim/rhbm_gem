#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include "data/object/LocalPotentialFitState.hpp"
#include "data/object/ModelAnalysisState.hpp"
#include <rhbm_gem/data/object/GroupPotentialEntry.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

namespace rhbm_gem {

ModelObject::ModelObject() :
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" }, m_kd_tree_root{ nullptr },
    m_component_key_system{ std::make_unique<ComponentKeySystem>() },
    m_atom_key_system{ std::make_unique<AtomKeySystem>() },
    m_bond_key_system{ std::make_unique<BondKeySystem>() },
    m_analysis_state{ std::make_unique<ModelAnalysisState>() }
{
}

ModelObject::ModelObject(std::vector<std::unique_ptr<AtomObject>> atom_object_list) :
    m_atom_list{ std::move(atom_object_list) },
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" }, m_kd_tree_root{ nullptr },
    m_component_key_system{ std::make_unique<ComponentKeySystem>() },
    m_atom_key_system{ std::make_unique<AtomKeySystem>() },
    m_bond_key_system{ std::make_unique<BondKeySystem>() },
    m_analysis_state{ std::make_unique<ModelAnalysisState>() }
{
    SyncDerivedState();
}

ModelObject::~ModelObject()
{

}

ModelObject::ModelObject(const ModelObject & other) :
    m_key_tag{ other.m_key_tag }, m_pdb_id{ other.m_pdb_id }, m_emd_id{ other.m_emd_id },
    m_resolution_method{ other.m_resolution_method }, m_resolution{ other.m_resolution },
    m_chain_id_list_map{ other.m_chain_id_list_map },
    m_kd_tree_root{ nullptr },
    m_component_key_system{
        other.m_component_key_system != nullptr ?
            std::make_unique<ComponentKeySystem>(*other.m_component_key_system) : nullptr },
    m_atom_key_system{
        other.m_atom_key_system != nullptr ?
            std::make_unique<AtomKeySystem>(*other.m_atom_key_system) : nullptr },
    m_bond_key_system{
        other.m_bond_key_system != nullptr ?
            std::make_unique<BondKeySystem>(*other.m_bond_key_system) : nullptr },
    m_analysis_state{ std::make_unique<ModelAnalysisState>() }
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

    m_bond_list.reserve(other.m_bond_list.size());
    std::unordered_map<const BondObject *, BondObject *> bond_ptr_map;
    bond_ptr_map.reserve(other.m_bond_list.size());
    for (const auto & bond : other.m_bond_list)
    {
        auto * atom_1{ atom_ptr_map.at(bond->GetAtomObject1()) };
        auto * atom_2{ atom_ptr_map.at(bond->GetAtomObject2()) };
        auto cloned_bond{ std::make_unique<BondObject>(atom_1, atom_2) };
        cloned_bond->SetSelectedFlag(bond->GetSelectedFlag());
        cloned_bond->SetSpecialBondFlag(bond->GetSpecialBondFlag());
        cloned_bond->SetBondKey(bond->GetBondKey());
        cloned_bond->SetBondType(bond->GetBondType());
        cloned_bond->SetBondOrder(bond->GetBondOrder());
        if (bond->GetLocalPotentialEntry() != nullptr)
        {
            cloned_bond->SetLocalPotentialEntry(
                std::make_unique<LocalPotentialEntry>(*bond->GetLocalPotentialEntry()));
        }
        bond_ptr_map[bond.get()] = cloned_bond.get();
        m_bond_list.emplace_back(std::move(cloned_bond));
    }

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
                if (entry == nullptr)
                {
                    set_entry(class_key, nullptr);
                    continue;
                }

                auto cloned_entry{ std::make_unique<GroupPotentialEntry>() };
                for (const auto & [group_key, bucket] : entry->GetGroups())
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

    const auto & source_analysis_state{ other.GetAnalysisState() };
    copy_group_entries(
        source_analysis_state.GetAtomGroupPotentialEntryMap(),
        [this](const std::string & class_key, std::unique_ptr<GroupPotentialEntry> entry)
        {
            m_analysis_state->SetAtomGroupPotentialEntry(class_key, std::move(entry));
        });
    copy_group_entries(
        source_analysis_state.GetBondGroupPotentialEntryMap(),
        [this](const std::string & class_key, std::unique_ptr<GroupPotentialEntry> entry)
        {
            m_analysis_state->SetBondGroupPotentialEntry(class_key, std::move(entry));
        });

    for (const auto & atom : m_atom_list)
    {
        if (const auto * fit_state{ source_analysis_state.FindAtomFitState(*atom) };
            fit_state != nullptr)
        {
            m_analysis_state->EnsureAtomFitState(*atom) = *fit_state;
        }
    }
    for (const auto & bond : m_bond_list)
    {
        if (const auto * fit_state{ source_analysis_state.FindBondFitState(*bond) };
            fit_state != nullptr)
        {
            m_analysis_state->EnsureBondFitState(*bond) = *fit_state;
        }
    }

    RebuildSelectionIndex();
}

void ModelObject::RebuildSelectionIndex()
{
    SyncDerivedState();
}

void ModelObject::FinalizeLoad()
{
    SyncDerivedState();
}

void ModelObject::ApplySymmetrySelection(bool is_asymmetry)
{
    FilterSelectionFromSymmetry(is_asymmetry);
    SyncDerivedState();
}

void ModelObject::RebuildSelectionState()
{
    m_serial_id_atom_map.clear();
    for (auto & atom : m_atom_list)
    {
        m_serial_id_atom_map[atom->GetSerialID()] = atom.get();
    }
    BuildSelectedAtomList();
    BuildSelectedBondList();
}

void ModelObject::InvalidateDerivedCaches()
{
    m_center_of_mass_position.reset();
    for (auto & axis_range : m_model_position_range)
    {
        axis_range.reset();
    }
    m_kd_tree_root.reset();
}

void ModelObject::SyncDerivedState()
{
    InvalidateDerivedCaches();
    RebuildSelectionState();
}

void ModelObject::AddAtom(std::unique_ptr<AtomObject> atom)
{
    m_atom_list.emplace_back(std::move(atom));
    m_analysis_state->Clear();
    SyncDerivedState();
}

void ModelObject::AddBond(std::unique_ptr<BondObject> bond)
{
    m_bond_list.emplace_back(std::move(bond));
    m_analysis_state->Clear();
    SyncDerivedState();
}

void ModelObject::SetAtomList(std::vector<std::unique_ptr<AtomObject>> atom_list)
{
    m_atom_list = std::move(atom_list);
    m_analysis_state->Clear();
    SyncDerivedState();
}

void ModelObject::SetBondList(std::vector<std::unique_ptr<BondObject>> bond_list)
{
    m_bond_list = std::move(bond_list);
    m_analysis_state->Clear();
    SyncDerivedState();
}

void ModelObject::AddChemicalComponentEntry(
    ComponentKey key, std::unique_ptr<ChemicalComponentEntry> entry)
{
    m_chemical_component_entry_map[key] = std::move(entry);
}

void ModelObject::SetChemicalComponentEntryMap(
    std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> entry_map)
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

const std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>> &
ModelObject::GetChemicalComponentEntryMap() const
{
    return m_chemical_component_entry_map;
}

ModelAnalysisState & ModelObject::GetAnalysisState()
{
    return *m_analysis_state;
}

const ModelAnalysisState & ModelObject::GetAnalysisState() const
{
    return *m_analysis_state;
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
