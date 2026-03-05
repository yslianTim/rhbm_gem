#include "DataObjectWorkflowVisitors.hpp"

#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <stdexcept>

#include "AtomObject.hpp"
#include "AtomSelector.hpp"
#include "BondObject.hpp"
#include "ComponentHelper.hpp"
#include "LocalPotentialEntry.hpp"
#include "Logger.hpp"
#include "MapObject.hpp"
#include "ModelObject.hpp"

namespace rhbm_gem {

namespace {

[[noreturn]] void ThrowUnsupportedType(const char * visitor_name, const char * supported_type)
{
    throw std::logic_error(
        std::string(visitor_name) + " supports " + supported_type + " only.");
}

double CalculateAtomChargeForSimulation(
    const AtomObject & atom,
    PartialCharge partial_charge_choice)
{
    switch (partial_charge_choice)
    {
    case PartialCharge::NEUTRAL:
        return 0.0;
    case PartialCharge::PARTIAL:
        return ComponentHelper::GetPartialCharge(
            atom.GetResidue(),
            atom.GetSpot(),
            atom.GetStructure());
    case PartialCharge::AMBER:
        return ComponentHelper::GetPartialCharge(
            atom.GetResidue(),
            atom.GetSpot(),
            atom.GetStructure(),
            true);
    default:
        Logger::Log(
            LogLevel::Error,
            "SimulationAtomPreparationVisitor reached invalid partial-charge choice: "
            + std::to_string(static_cast<int>(partial_charge_choice)));
        return 0.0;
    }
}

} // namespace

void MapNormalizeVisitor::VisitAtomObject(AtomObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("MapNormalizeVisitor", "MapObject");
}

void MapNormalizeVisitor::VisitBondObject(BondObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("MapNormalizeVisitor", "MapObject");
}

void MapNormalizeVisitor::VisitModelObject(ModelObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("MapNormalizeVisitor", "MapObject");
}

void MapNormalizeVisitor::VisitMapObject(MapObject & data_object)
{
    data_object.MapValueArrayNormalization();
}

ModelPreparationVisitor::ModelPreparationVisitor(ModelPreparationOptions options) :
    m_options{ options }
{
}

void ModelPreparationVisitor::VisitAtomObject(AtomObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelPreparationVisitor", "ModelObject");
}

void ModelPreparationVisitor::VisitBondObject(BondObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelPreparationVisitor", "ModelObject");
}

void ModelPreparationVisitor::VisitModelObject(ModelObject & data_object)
{
    if (m_options.select_all_atoms)
    {
        for (auto & atom : data_object.GetAtomList())
        {
            atom->SetSelectedFlag(true);
        }
    }
    if (m_options.select_all_bonds)
    {
        for (auto & bond : data_object.GetBondList())
        {
            bond->SetSelectedFlag(true);
        }
    }

    if (m_options.apply_atom_symmetry_filter)
    {
        data_object.FilterAtomFromSymmetry(m_options.asymmetry_flag);
    }
    if (m_options.apply_bond_symmetry_filter)
    {
        data_object.FilterBondFromSymmetry(m_options.asymmetry_flag);
    }

    bool has_updated{ false };
    if (m_options.update_model)
    {
        data_object.Update();
        has_updated = true;
    }

    if ((m_options.initialize_atom_local_entries || m_options.initialize_bond_local_entries)
        && !has_updated)
    {
        data_object.Update();
    }

    if (m_options.initialize_atom_local_entries)
    {
        for (auto * atom : data_object.GetSelectedAtomList())
        {
            atom->AddLocalPotentialEntry(std::make_unique<LocalPotentialEntry>());
        }
    }
    if (m_options.initialize_bond_local_entries)
    {
        for (auto * bond : data_object.GetSelectedBondList())
        {
            bond->AddLocalPotentialEntry(std::make_unique<LocalPotentialEntry>());
        }
    }
}

void ModelPreparationVisitor::VisitMapObject(MapObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelPreparationVisitor", "ModelObject");
}

ModelSelectionVisitor::ModelSelectionVisitor(::AtomSelector * selector) :
    m_selector{ selector }
{
}

void ModelSelectionVisitor::VisitAtomObject(AtomObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelSelectionVisitor", "ModelObject");
}

void ModelSelectionVisitor::VisitBondObject(BondObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelSelectionVisitor", "ModelObject");
}

void ModelSelectionVisitor::VisitModelObject(ModelObject & data_object)
{
    if (m_selector == nullptr)
    {
        throw std::logic_error("ModelSelectionVisitor requires a non-null AtomSelector.");
    }
    for (auto & atom : data_object.GetAtomList())
    {
        atom->SetSelectedFlag(
            m_selector->GetSelectionFlag(
                atom->GetChainID(),
                atom->GetResidue(),
                atom->GetElement()));
    }
}

void ModelSelectionVisitor::VisitMapObject(MapObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelSelectionVisitor", "ModelObject");
}

ModelSelectedAtomCollectorVisitor::ModelSelectedAtomCollectorVisitor(
    ModelAtomCollectorOptions options) :
    m_options{ options }
{
}

void ModelSelectedAtomCollectorVisitor::VisitAtomObject(AtomObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelSelectedAtomCollectorVisitor", "ModelObject");
}

void ModelSelectedAtomCollectorVisitor::VisitBondObject(BondObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelSelectedAtomCollectorVisitor", "ModelObject");
}

void ModelSelectedAtomCollectorVisitor::VisitModelObject(ModelObject & data_object)
{
    m_atom_list.clear();
    m_atom_list.reserve(data_object.GetAtomList().size());
    for (auto & atom : data_object.GetAtomList())
    {
        if (m_options.selected_only && !atom->GetSelectedFlag())
        {
            continue;
        }
        if (m_options.require_local_potential_entry && atom->GetLocalPotentialEntry() == nullptr)
        {
            continue;
        }
        m_atom_list.emplace_back(atom.get());
    }
}

void ModelSelectedAtomCollectorVisitor::VisitMapObject(MapObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelSelectedAtomCollectorVisitor", "ModelObject");
}

SimulationAtomPreparationVisitor::SimulationAtomPreparationVisitor(
    SimulationAtomPreparationOptions options) :
    m_options{ options },
    m_range_minimum{
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max() },
    m_range_maximum{
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest() }
{
}

void SimulationAtomPreparationVisitor::VisitAtomObject(AtomObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("SimulationAtomPreparationVisitor", "ModelObject");
}

void SimulationAtomPreparationVisitor::VisitBondObject(BondObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("SimulationAtomPreparationVisitor", "ModelObject");
}

void SimulationAtomPreparationVisitor::VisitModelObject(ModelObject & data_object)
{
    m_atom_list.clear();
    m_atom_charge_map.clear();
    m_has_atom = false;
    m_range_minimum = {
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max()
    };
    m_range_maximum = {
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest()
    };

    m_atom_list.reserve(data_object.GetNumberOfAtom());
    for (auto & atom : data_object.GetAtomList())
    {
        if (!m_options.include_unknown_atoms && atom->IsUnknownAtom())
        {
            continue;
        }
        m_atom_list.emplace_back(atom.get());
        m_atom_charge_map.emplace(
            atom->GetSerialID(),
            CalculateAtomChargeForSimulation(*atom, m_options.partial_charge_choice));

        const auto & atom_position{ atom->GetPositionRef() };
        m_range_minimum[0] = std::min(m_range_minimum[0], atom_position[0]);
        m_range_minimum[1] = std::min(m_range_minimum[1], atom_position[1]);
        m_range_minimum[2] = std::min(m_range_minimum[2], atom_position[2]);
        m_range_maximum[0] = std::max(m_range_maximum[0], atom_position[0]);
        m_range_maximum[1] = std::max(m_range_maximum[1], atom_position[1]);
        m_range_maximum[2] = std::max(m_range_maximum[2], atom_position[2]);
    }
    m_has_atom = !m_atom_list.empty();
}

void SimulationAtomPreparationVisitor::VisitMapObject(MapObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("SimulationAtomPreparationVisitor", "ModelObject");
}

void ModelAtomBondContextVisitor::VisitAtomObject(AtomObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelAtomBondContextVisitor", "ModelObject");
}

void ModelAtomBondContextVisitor::VisitBondObject(BondObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelAtomBondContextVisitor", "ModelObject");
}

void ModelAtomBondContextVisitor::VisitModelObject(ModelObject & data_object)
{
    m_atom_map.clear();
    m_bond_map.clear();

    const auto & selected_atom_list{ data_object.GetSelectedAtomList() };
    const auto & selected_bond_list{ data_object.GetSelectedBondList() };
    m_atom_map.reserve(selected_atom_list.size());
    m_bond_map.reserve(selected_atom_list.size());

    for (auto * atom : selected_atom_list)
    {
        m_atom_map.emplace(atom->GetSerialID(), atom);
    }
    for (auto * bond : selected_bond_list)
    {
        m_bond_map[bond->GetAtomSerialID1()].emplace_back(bond);
        m_bond_map[bond->GetAtomSerialID2()].emplace_back(bond);
    }
}

void ModelAtomBondContextVisitor::VisitMapObject(MapObject & data_object)
{
    (void)data_object;
    ThrowUnsupportedType("ModelAtomBondContextVisitor", "ModelObject");
}

MapSamplingWorkflow::MapSamplingWorkflow(::SamplerBase * sampler) :
    m_interpolation_visitor{ sampler }
{
}

std::vector<std::tuple<float, float>> MapSamplingWorkflow::Sample(
    const MapObject & map_object,
    const std::array<float, 3> & position,
    const std::array<float, 3> & axis_vector)
{
    m_interpolation_visitor.SetPosition(position);
    m_interpolation_visitor.SetAxisVector(axis_vector);
    map_object.Accept(m_interpolation_visitor);
    return m_interpolation_visitor.ConsumeSamplingDataList();
}

void NormalizeMapObject(MapObject & map_object)
{
    MapNormalizeVisitor visitor;
    map_object.Accept(visitor);
}

void PrepareModelObject(ModelObject & model_object, const ModelPreparationOptions & options)
{
    ModelPreparationVisitor visitor{ options };
    model_object.Accept(visitor, ModelVisitMode::SelfOnly);
}

void ApplyModelSelection(ModelObject & model_object, ::AtomSelector & selector)
{
    ModelSelectionVisitor visitor{ &selector };
    model_object.Accept(visitor, ModelVisitMode::SelfOnly);
}

} // namespace rhbm_gem
