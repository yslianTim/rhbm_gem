#include "workflow/DataObjectWorkflowOps.hpp"

#include <algorithm>
#include <memory>
#include <string>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/utils/domain/AtomSelector.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem {
namespace {

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
            "PrepareSimulationAtoms reached invalid partial-charge choice: "
            + std::to_string(static_cast<int>(partial_charge_choice)));
        return 0.0;
    }
}

} // namespace

void NormalizeMapObject(MapObject & map_object)
{
    map_object.MapValueArrayNormalization();
}

void PrepareModelObject(ModelObject & model_object, const ModelPreparationOptions & options)
{
    if (options.select_all_atoms)
    {
        for (auto & atom : model_object.GetAtomList())
        {
            atom->SetSelectedFlag(true);
        }
    }
    if (options.select_all_bonds)
    {
        for (auto & bond : model_object.GetBondList())
        {
            bond->SetSelectedFlag(true);
        }
    }

    if (options.apply_atom_symmetry_filter)
    {
        model_object.FilterAtomFromSymmetry(options.asymmetry_flag);
    }
    if (options.apply_bond_symmetry_filter)
    {
        model_object.FilterBondFromSymmetry(options.asymmetry_flag);
    }

    bool has_updated{ false };
    if (options.update_model)
    {
        model_object.Update();
        has_updated = true;
    }
    if ((options.initialize_atom_local_entries || options.initialize_bond_local_entries)
        && !has_updated)
    {
        model_object.Update();
    }

    if (options.initialize_atom_local_entries)
    {
        for (auto * atom : model_object.GetSelectedAtomList())
        {
            atom->AddLocalPotentialEntry(std::make_unique<LocalPotentialEntry>());
        }
    }
    if (options.initialize_bond_local_entries)
    {
        for (auto * bond : model_object.GetSelectedBondList())
        {
            bond->AddLocalPotentialEntry(std::make_unique<LocalPotentialEntry>());
        }
    }
}

void ApplyModelSelection(ModelObject & model_object, ::AtomSelector & selector)
{
    for (auto & atom : model_object.GetAtomList())
    {
        atom->SetSelectedFlag(
            selector.GetSelectionFlag(
                atom->GetChainID(),
                atom->GetResidue(),
                atom->GetElement()));
    }
}

std::vector<AtomObject *> CollectModelAtoms(
    ModelObject & model_object,
    ModelAtomCollectorOptions options)
{
    std::vector<AtomObject *> atom_list;
    atom_list.reserve(model_object.GetAtomList().size());
    for (auto & atom : model_object.GetAtomList())
    {
        if (options.selected_only && !atom->GetSelectedFlag())
        {
            continue;
        }
        if (options.require_local_potential_entry && atom->GetLocalPotentialEntry() == nullptr)
        {
            continue;
        }
        atom_list.emplace_back(atom.get());
    }
    return atom_list;
}

SimulationAtomPreparationResult PrepareSimulationAtoms(
    ModelObject & model_object,
    SimulationAtomPreparationOptions options)
{
    SimulationAtomPreparationResult result;
    result.atom_list.reserve(model_object.GetNumberOfAtom());
    for (auto & atom : model_object.GetAtomList())
    {
        if (!options.include_unknown_atoms && atom->IsUnknownAtom())
        {
            continue;
        }
        result.atom_list.emplace_back(atom.get());
        result.atom_charge_map.emplace(
            atom->GetSerialID(),
            CalculateAtomChargeForSimulation(*atom, options.partial_charge_choice));

        const auto & atom_position{ atom->GetPositionRef() };
        result.range_minimum[0] = std::min(result.range_minimum[0], atom_position[0]);
        result.range_minimum[1] = std::min(result.range_minimum[1], atom_position[1]);
        result.range_minimum[2] = std::min(result.range_minimum[2], atom_position[2]);
        result.range_maximum[0] = std::max(result.range_maximum[0], atom_position[0]);
        result.range_maximum[1] = std::max(result.range_maximum[1], atom_position[1]);
        result.range_maximum[2] = std::max(result.range_maximum[2], atom_position[2]);
    }
    result.has_atom = !result.atom_list.empty();
    return result;
}

ModelAtomBondContext BuildModelAtomBondContext(ModelObject & model_object)
{
    ModelAtomBondContext result;
    const auto & selected_atom_list{ model_object.GetSelectedAtomList() };
    const auto & selected_bond_list{ model_object.GetSelectedBondList() };
    result.atom_map.reserve(selected_atom_list.size());
    result.bond_map.reserve(selected_atom_list.size());

    for (auto * atom : selected_atom_list)
    {
        result.atom_map.emplace(atom->GetSerialID(), atom);
    }
    for (auto * bond : selected_bond_list)
    {
        result.bond_map[bond->GetAtomSerialID1()].emplace_back(bond);
        result.bond_map[bond->GetAtomSerialID2()].emplace_back(bond);
    }
    return result;
}

} // namespace rhbm_gem
