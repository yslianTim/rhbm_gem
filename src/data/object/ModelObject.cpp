#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/NumericValidation.hpp>

#include "data/detail/ModelDerivedState.hpp"
#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"

#include <algorithm>
#include <stdexcept>

namespace {

constexpr double kInitialLocalAlpha{ 0.0 };
constexpr double kInitialGroupAlpha{ 0.0 };

double ComputeDistanceSquare(
    const rhbm_gem::AtomObject & atom_1,
    const rhbm_gem::AtomObject & atom_2)
{
    const auto & position_1{ atom_1.GetPositionRef() };
    const auto & position_2{ atom_2.GetPositionRef() };
    const auto dx{ static_cast<double>(position_1.at(0)) - static_cast<double>(position_2.at(0)) };
    const auto dy{ static_cast<double>(position_1.at(1)) - static_cast<double>(position_2.at(1)) };
    const auto dz{ static_cast<double>(position_1.at(2)) - static_cast<double>(position_2.at(2)) };
    return dx * dx + dy * dy + dz * dz;
}

bool IsBackboneSelectionSpot(Spot spot)
{
    switch (spot)
    {
    case Spot::C:
    case Spot::CA:
    case Spot::N:
    case Spot::O:
    case Spot::H:
    case Spot::HA:
        return true;
    default:
        return false;
    }
}

} // namespace

namespace rhbm_gem {

ModelObject::ModelObject() :
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" },
    m_derived_state{ std::make_unique<ModelDerivedState>() },
    m_analysis_data{ std::make_unique<ModelAnalysisData>() }
{
}

ModelObject::ModelObject(std::vector<std::unique_ptr<AtomObject>> atom_object_list) :
    m_atom_list{ std::move(atom_object_list) },
    m_key_tag{ "" }, m_pdb_id{ "" }, m_emd_id{ "" },
    m_derived_state{ std::make_unique<ModelDerivedState>() },
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

size_t ModelObject::GetNumberOfAtom() const
{
    return m_atom_list.size();
}

size_t ModelObject::GetNumberOfBond() const
{
    return m_bond_list.size();
}

const std::vector<std::unique_ptr<AtomObject>> & ModelObject::GetAtomList() const
{
    return m_atom_list;
}

const std::vector<std::unique_ptr<BondObject>> & ModelObject::GetBondList() const
{
    return m_bond_list;
}

const std::vector<int> & ModelObject::GetSequenceIDList() const
{
    m_sequence_id_list.clear();
    m_sequence_id_list.reserve(m_atom_list.size());
    for (const auto & atom : m_atom_list)
    {
        m_sequence_id_list.emplace_back(atom->GetSequenceID());
    }

    std::sort(m_sequence_id_list.begin(), m_sequence_id_list.end());
    const auto unique_end{ std::unique(m_sequence_id_list.begin(), m_sequence_id_list.end()) };
    m_sequence_id_list.erase(unique_end, m_sequence_id_list.end());
    return m_sequence_id_list;
}

const std::vector<AtomObject *> & ModelObject::GetSelectedAtomList(int residue_id) const
{
    static const std::vector<AtomObject *> empty_atom_list{};

    const auto iter{ m_selected_residue_id_atom_list_map.find(residue_id) };
    return iter == m_selected_residue_id_atom_list_map.end() ?
        empty_atom_list :
        iter->second;
}

ModelObject::ModelObject(ModelObject && other) noexcept :
    m_atom_list{ std::move(other.m_atom_list) },
    m_bond_list{ std::move(other.m_bond_list) },
    m_chain_id_list_map{ std::move(other.m_chain_id_list_map) },
    m_chemical_component_entry_map{ std::move(other.m_chemical_component_entry_map) },
    m_component_key_system{ std::move(other.m_component_key_system) },
    m_atom_key_system{ std::move(other.m_atom_key_system) },
    m_bond_key_system{ std::move(other.m_bond_key_system) },
    m_selected_atom_list{ std::move(other.m_selected_atom_list) },
    m_selected_bond_list{ std::move(other.m_selected_bond_list) },
    m_key_tag{ std::move(other.m_key_tag) },
    m_pdb_id{ std::move(other.m_pdb_id) },
    m_emd_id{ std::move(other.m_emd_id) },
    m_resolution_method{ std::move(other.m_resolution_method) },
    m_resolution{ other.m_resolution },
    m_serial_id_atom_map{ std::move(other.m_serial_id_atom_map) },
    m_derived_state{ std::move(other.m_derived_state) },
    m_analysis_data{ std::move(other.m_analysis_data) }
{
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
    if (m_analysis_data == nullptr)
    {
        m_analysis_data = std::make_unique<ModelAnalysisData>();
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
    m_chain_id_list_map = std::move(other.m_chain_id_list_map);
    m_chemical_component_entry_map = std::move(other.m_chemical_component_entry_map);
    m_component_key_system = std::move(other.m_component_key_system);
    m_atom_key_system = std::move(other.m_atom_key_system);
    m_bond_key_system = std::move(other.m_bond_key_system);
    m_selected_atom_list = std::move(other.m_selected_atom_list);
    m_selected_bond_list = std::move(other.m_selected_bond_list);
    m_key_tag = std::move(other.m_key_tag);
    m_pdb_id = std::move(other.m_pdb_id);
    m_emd_id = std::move(other.m_emd_id);
    m_resolution_method = std::move(other.m_resolution_method);
    m_resolution = other.m_resolution;
    m_serial_id_atom_map = std::move(other.m_serial_id_atom_map);
    m_derived_state = std::move(other.m_derived_state);
    m_analysis_data = std::move(other.m_analysis_data);

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
    if (m_analysis_data == nullptr)
    {
        m_analysis_data = std::make_unique<ModelAnalysisData>();
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
    m_derived_state{ std::make_unique<ModelDerivedState>() },
    m_analysis_data{ std::make_unique<ModelAnalysisData>() }
{
    std::unordered_map<const AtomObject *, AtomObject *> atom_ptr_map;
    m_chain_id_list_map = other.m_chain_id_list_map;
    m_component_key_system =
        other.m_component_key_system != nullptr ?
            std::make_unique<ComponentKeySystem>(*other.m_component_key_system) :
            std::make_unique<ComponentKeySystem>();
    m_atom_key_system =
        other.m_atom_key_system != nullptr ?
            std::make_unique<AtomKeySystem>(*other.m_atom_key_system) :
            std::make_unique<AtomKeySystem>();
    m_bond_key_system =
        other.m_bond_key_system != nullptr ?
            std::make_unique<BondKeySystem>(*other.m_bond_key_system) :
            std::make_unique<BondKeySystem>();

    m_atom_list.reserve(other.m_atom_list.size());
    atom_ptr_map.reserve(other.m_atom_list.size());
    for (const auto & atom : other.m_atom_list)
    {
        auto cloned_atom{ std::make_unique<AtomObject>(*atom) };
        atom_ptr_map[atom.get()] = cloned_atom.get();
        m_atom_list.emplace_back(std::move(cloned_atom));
    }

    m_bond_list.reserve(other.m_bond_list.size());
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
        m_bond_list.emplace_back(std::move(cloned_bond));
    }

    for (const auto & [component_key, entry] : other.m_chemical_component_entry_map)
    {
        m_chemical_component_entry_map[component_key] =
            std::make_unique<ChemicalComponentEntry>(*entry);
    }

    AttachOwnedObjects();

    const auto & source_analysis_data{ ModelAnalysisData::Of(other) };
    {
        const auto & source_entry{ source_analysis_data.AtomGroupEntry() };
        auto & cloned_entry{ m_analysis_data->AtomGroupEntry() };
        for (const auto group_key : source_entry.CollectGroupKeys())
        {
            GroupGaussianResult result;
            result.mean = source_entry.GetMean(group_key);
            result.mdpde = source_entry.GetMDPDE(group_key);
            result.prior = source_entry.GetPriorWithUncertainty(group_key);
            result.alpha_g = source_entry.GetAlphaG(group_key);
            cloned_entry.SetGaussianResult(group_key, result);
            cloned_entry.ReserveMembers(group_key, source_entry.GetMemberCount(group_key));
            for (auto * atom : source_entry.GetMembers(group_key))
            {
                cloned_entry.AddMember(group_key, *atom_ptr_map.at(atom));
            }
        }
    }

    for (const auto & atom : m_atom_list)
    {
        if (const auto * local_entry{ source_analysis_data.FindAtomLocalEntry(*atom) };
            local_entry != nullptr)
        {
            m_analysis_data->SetAtomLocalEntry(
                *atom,
                std::make_unique<LocalPotentialEntry>(*local_entry));
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

ModelAnalysisEditor ModelObject::EditAnalysis()
{
    return ModelAnalysisEditor(*this);
}

ModelAnalysisView ModelObject::GetAnalysisView() const
{
    return ModelAnalysisView(*this);
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

std::vector<AtomObject *> ModelObject::FindNeighborAtoms(
    const AtomObject & center_atom,
    double radius,
    bool include_center) const
{
    numeric_validation::RequireFiniteNonNegative(radius, "ModelObject::FindNeighborAtoms radius");
    if (center_atom.m_owner_model != this)
    {
        throw std::invalid_argument(
            "ModelObject::FindNeighborAtoms center atom does not belong to this model.");
    }

    auto & mutable_model{ const_cast<ModelObject &>(*this) };
    auto results{
        ModelDerivedState::Of(mutable_model).FindAtomsInRange(mutable_model, center_atom, radius) };

    if (!include_center)
    {
        results.erase(
            std::remove(results.begin(), results.end(), const_cast<AtomObject *>(&center_atom)),
            results.end());
    }

    std::sort(
        results.begin(),
        results.end(),
        [&center_atom](const AtomObject * lhs, const AtomObject * rhs)
        {
            const auto lhs_distance{ ComputeDistanceSquare(center_atom, *lhs) };
            const auto rhs_distance{ ComputeDistanceSquare(center_atom, *rhs) };
            if (lhs_distance != rhs_distance)
            {
                return lhs_distance < rhs_distance;
            }
            return lhs->GetSerialID() < rhs->GetSerialID();
        });

    return results;
}

bool ModelObject::HasChemicalComponentEntry(ComponentKey component_key) const
{
    return m_chemical_component_entry_map.find(component_key) !=
           m_chemical_component_entry_map.end();
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

void ModelObject::ApplyElementSelection(Element element, bool is_exclusion)
{
    if (!is_exclusion)
    {
        return;
    }

    for (auto & atom : m_atom_list)
    {
        if (atom->GetElement() == element)
        {
            atom->SetSelectedFlag(false);
        }
    }
    RebuildSelection();
}

void ModelObject::ApplySpotSelection(Spot spot, bool is_exclusion)
{
    if (!is_exclusion)
    {
        return;
    }

    for (auto & atom : m_atom_list)
    {
        if (atom->GetSpot() == spot)
        {
            atom->SetSelectedFlag(false);
        }
    }
    RebuildSelection();
}

void ModelObject::ApplyBackboneSelection(bool is_exclusion)
{
    if (!is_exclusion)
    {
        return;
    }

    for (auto * atom : m_selected_atom_list)
    {
        if (!IsBackboneSelectionSpot(atom->GetSpot()))
        {
            atom->SetSelectedFlag(false);
        }
    }
    RebuildSelection();
}

void ModelObject::ApplyComponentIDSelection(std::string component_id, bool is_exclusion)
{
    if (!is_exclusion)
    {
        return;
    }

    for (auto & atom : m_atom_list)
    {
        if (atom->GetComponentID() == component_id)
        {
            atom->SetSelectedFlag(false);
        }
    }
    RebuildSelection();
}

void ModelObject::ApplySimulationMetadata(double simulated_map_resolution)
{
    SetEmdID("Simulation");
    SetResolution(simulated_map_resolution);
    SetResolutionMethod("Blurring Width");
}

void ModelObject::LocalPotentialInitialization()
{
    auto analysis{ EditAnalysis() };
    analysis.Clear();
    analysis.RebuildAtomGroupsFromSelection();
    analysis.InitializeLocalAlpha(kInitialLocalAlpha);
    analysis.InitializeGroupAlpha(kInitialGroupAlpha);
}

void ModelObject::ClearTransientFitStates()
{
    EditAnalysis().ClearTransientFitStates();
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
    m_selected_residue_id_atom_list_map.clear();
    m_selected_atom_list.reserve(m_atom_list.size());
    m_selected_residue_id_atom_list_map.reserve(m_atom_list.size());
    for (auto & atom : m_atom_list)
    {
        if (atom->m_is_selected == false) continue;
        m_selected_atom_list.emplace_back(atom.get());
        m_selected_residue_id_atom_list_map[atom->GetSequenceID()].emplace_back(atom.get());
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
        (void)component_key;
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
        (void)component_key;
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

void ModelObject::PrintSummary() const
{
    Logger::Log(LogLevel::Info, this->GetAnalysisView().GetAtomCountingSummary());
    Logger::Log(LogLevel::Info, this->GetAnalysisView().GetAtomGroupingSummary());
}

} // namespace rhbm_gem
