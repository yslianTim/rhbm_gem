#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/BondKeySystem.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class ChemicalComponentEntry;
class ModelObject;

// Plain parts payload for internal import and persistence paths.
struct ModelObjectParts
{
    ModelObjectParts();
    ~ModelObjectParts();
    ModelObjectParts(ModelObjectParts && other) noexcept;
    ModelObjectParts & operator=(ModelObjectParts && other) noexcept;
    ModelObjectParts(const ModelObjectParts & other) = delete;
    ModelObjectParts & operator=(const ModelObjectParts & other) = delete;

    // Structure parts.
    std::vector<std::unique_ptr<AtomObject>> atom_list;
    std::vector<std::unique_ptr<BondObject>> bond_list;

    // Metadata and topology parts.
    std::unordered_map<std::string, std::vector<std::string>> chain_id_list_map;
    std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>>
        chemical_component_entry_map;

    // Key systems.
    std::unique_ptr<ComponentKeySystem> component_key_system{
        std::make_unique<ComponentKeySystem>() };
    std::unique_ptr<AtomKeySystem> atom_key_system{
        std::make_unique<AtomKeySystem>() };
    std::unique_ptr<BondKeySystem> bond_key_system{
        std::make_unique<BondKeySystem>() };
};

ModelObject AssembleModelObject(ModelObjectParts parts);

} // namespace rhbm_gem
