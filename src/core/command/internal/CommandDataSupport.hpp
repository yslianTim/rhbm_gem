#pragma once

#include <array>
#include <filesystem>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/core/command/OptionEnumClass.hpp>
#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

class AtomSelector;

namespace rhbm_gem {

class AtomObject;
class BondObject;
class MapObject;
class ModelObject;

namespace command_data_loader {

namespace detail {

template <typename TypedDataObject>
std::shared_ptr<TypedDataObject> ProcessTypedFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label,
    std::string_view object_type_name)
{
    try
    {
        data_manager.ProcessFile(path, std::string(key_tag));
        return data_manager.GetTypedDataObject<TypedDataObject>(std::string(key_tag));
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to process " + std::string(label) + " from '" + path.string()
            + "' as " + std::string(object_type_name) + ": " + ex.what());
    }
}

template <typename TypedDataObject>
std::shared_ptr<TypedDataObject> LoadTypedObject(
    DataObjectManager & data_manager,
    std::string_view key_tag,
    std::string_view label,
    std::string_view object_type_name)
{
    try
    {
        data_manager.LoadDataObject(std::string(key_tag));
        return data_manager.GetTypedDataObject<TypedDataObject>(std::string(key_tag));
    }
    catch (const std::exception & ex)
    {
        throw std::runtime_error(
            "Failed to load " + std::string(label) + " with key tag '"
            + std::string(key_tag) + "' as " + std::string(object_type_name)
            + ": " + ex.what());
    }
}

} // namespace detail

inline std::shared_ptr<ModelObject> ProcessModelFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label)
{
    return detail::ProcessTypedFile<ModelObject>(
        data_manager,
        path,
        key_tag,
        label,
        "ModelObject");
}

inline std::shared_ptr<MapObject> ProcessMapFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label)
{
    return detail::ProcessTypedFile<MapObject>(
        data_manager,
        path,
        key_tag,
        label,
        "MapObject");
}

inline std::shared_ptr<MapObject> OptionalProcessMapFile(
    DataObjectManager & data_manager,
    const std::filesystem::path & path,
    std::string_view key_tag,
    std::string_view label)
{
    if (path.empty())
    {
        return nullptr;
    }
    return ProcessMapFile(data_manager, path, key_tag, label);
}

inline std::shared_ptr<ModelObject> LoadModelObject(
    DataObjectManager & data_manager,
    std::string_view key_tag,
    std::string_view label)
{
    return detail::LoadTypedObject<ModelObject>(
        data_manager,
        key_tag,
        label,
        "ModelObject");
}

} // namespace command_data_loader

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
