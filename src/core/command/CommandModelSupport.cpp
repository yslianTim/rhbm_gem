#include "internal/command/CommandModelSupport.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/AtomSelector.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <algorithm>
#include <memory>
#include <string>

namespace rhbm_gem {
namespace {

void SelectAllAtoms(ModelObject & model_object)
{
    for (auto & atom : model_object.GetAtomList())
    {
        atom->SetSelectedFlag(true);
    }
}

void SelectAllBonds(ModelObject & model_object)
{
    for (auto & bond : model_object.GetBondList())
    {
        bond->SetSelectedFlag(true);
    }
}

void EnsureModelUpdated(ModelObject & model_object)
{
    model_object.Update();
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

void PrepareModelForPotentialAnalysis(ModelObject & model_object, bool asymmetry_flag)
{
    SelectAllAtoms(model_object);
    SelectAllBonds(model_object);
    model_object.FilterAtomFromSymmetry(asymmetry_flag);
    model_object.FilterBondFromSymmetry(asymmetry_flag);
    EnsureModelUpdated(model_object);

    for (auto * atom : model_object.GetSelectedAtomList())
    {
        atom->AddLocalPotentialEntry(std::make_unique<LocalPotentialEntry>());
    }
    for (auto * bond : model_object.GetSelectedBondList())
    {
        bond->AddLocalPotentialEntry(std::make_unique<LocalPotentialEntry>());
    }
}

void PrepareModelForVisualization(ModelObject & model_object)
{
    SelectAllAtoms(model_object);
    SelectAllBonds(model_object);
    EnsureModelUpdated(model_object);
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

std::vector<AtomObject *> CollectSelectedAtoms(ModelObject & model_object)
{
    std::vector<AtomObject *> atom_list;
    atom_list.reserve(model_object.GetAtomList().size());
    for (auto & atom : model_object.GetAtomList())
    {
        if (!atom->GetSelectedFlag())
        {
            continue;
        }
        atom_list.emplace_back(atom.get());
    }
    return atom_list;
}

std::vector<AtomObject *> CollectAtomsWithLocalPotentialEntries(ModelObject & model_object)
{
    std::vector<AtomObject *> atom_list;
    atom_list.reserve(model_object.GetAtomList().size());
    for (auto & atom : model_object.GetAtomList())
    {
        if (atom->GetLocalPotentialEntry() == nullptr)
        {
            continue;
        }
        atom_list.emplace_back(atom.get());
    }
    return atom_list;
}

SimulationAtomPreparationResult PrepareSimulationAtoms(
    ModelObject & model_object,
    PartialCharge partial_charge_choice,
    bool include_unknown_atoms)
{
    SimulationAtomPreparationResult result;
    result.atom_list.reserve(model_object.GetNumberOfAtom());
    for (auto & atom : model_object.GetAtomList())
    {
        if (!include_unknown_atoms && atom->IsUnknownAtom())
        {
            continue;
        }
        result.atom_list.emplace_back(atom.get());
        result.atom_charge_map.emplace(
            atom->GetSerialID(),
            CalculateAtomChargeForSimulation(*atom, partial_charge_choice));

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
