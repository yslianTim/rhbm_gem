#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/BondKeySystem.hpp>
#include <rhbm_gem/utils/domain/ComponentKeySystem.hpp>

namespace rhbm_gem {

// Plain assembly data for internal import and persistence paths.
struct ModelObjectAssembly
{
    std::vector<std::unique_ptr<AtomObject>> atom_list;
    std::vector<std::unique_ptr<BondObject>> bond_list;
    std::unordered_map<std::string, std::vector<std::string>> chain_id_list_map;
    std::unordered_map<ComponentKey, std::unique_ptr<ChemicalComponentEntry>>
        chemical_component_entry_map;
    std::unique_ptr<ComponentKeySystem> component_key_system{
        std::make_unique<ComponentKeySystem>() };
    std::unique_ptr<AtomKeySystem> atom_key_system{
        std::make_unique<AtomKeySystem>() };
    std::unique_ptr<BondKeySystem> bond_key_system{
        std::make_unique<BondKeySystem>() };
};

ModelObject AssembleModelObject(ModelObjectAssembly assembly);

} // namespace rhbm_gem
