#pragma once

#include <array>
#include <limits>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/core/command/CommandEnumClass.hpp>

class AtomSelector;

namespace rhbm_gem {

class AtomObject;
class BondObject;
class MapObject;
class ModelObject;

struct SimulationAtomPreparationResult
{
    std::vector<AtomObject *> atom_list;
    std::unordered_map<int, double> atom_charge_map;
    std::array<float, 3> range_minimum{
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max() };
    std::array<float, 3> range_maximum{
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest(),
        std::numeric_limits<float>::lowest() };
    bool has_atom{ false };
};

struct ModelAtomBondContext
{
    std::unordered_map<int, AtomObject *> atom_map;
    std::unordered_map<int, std::vector<BondObject *>> bond_map;
};

void NormalizeMapObject(MapObject & map_object);
void PrepareModelForPotentialAnalysis(ModelObject & model_object, bool asymmetry_flag);
void PrepareModelForVisualization(ModelObject & model_object);
void ApplyModelSelection(ModelObject & model_object, ::AtomSelector & selector);
std::vector<AtomObject *> CollectSelectedAtoms(ModelObject & model_object);
std::vector<AtomObject *> CollectAtomsWithLocalPotentialEntries(ModelObject & model_object);
SimulationAtomPreparationResult PrepareSimulationAtoms(
    ModelObject & model_object,
    PartialCharge partial_charge_choice,
    bool include_unknown_atoms = true);
ModelAtomBondContext BuildModelAtomBondContext(ModelObject & model_object);

} // namespace rhbm_gem
