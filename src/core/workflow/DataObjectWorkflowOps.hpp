#pragma once

#include <array>
#include <limits>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/core/OptionEnumClass.hpp>

class AtomSelector;

namespace rhbm_gem {

class MapObject;
class ModelObject;
class AtomObject;
class BondObject;

struct ModelPreparationOptions
{
    bool select_all_atoms{ false };
    bool select_all_bonds{ false };
    bool apply_atom_symmetry_filter{ false };
    bool apply_bond_symmetry_filter{ false };
    bool asymmetry_flag{ false };
    bool update_model{ false };
    bool initialize_atom_local_entries{ false };
    bool initialize_bond_local_entries{ false };
};

struct ModelAtomCollectorOptions
{
    bool selected_only{ true };
    bool require_local_potential_entry{ false };
};

struct SimulationAtomPreparationOptions
{
    PartialCharge partial_charge_choice{ PartialCharge::PARTIAL };
    bool include_unknown_atoms{ true };
};

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
void PrepareModelObject(ModelObject & model_object, const ModelPreparationOptions & options);
void ApplyModelSelection(ModelObject & model_object, ::AtomSelector & selector);
std::vector<AtomObject *> CollectModelAtoms(
    ModelObject & model_object,
    ModelAtomCollectorOptions options = {});
SimulationAtomPreparationResult PrepareSimulationAtoms(
    ModelObject & model_object,
    SimulationAtomPreparationOptions options = {});
ModelAtomBondContext BuildModelAtomBondContext(ModelObject & model_object);

} // namespace rhbm_gem
