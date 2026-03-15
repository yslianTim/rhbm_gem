#pragma once

#include <array>
#include <memory>
#include <string>
#include <string_view>

namespace rhbm_gem {

class SQLiteWrapper;
class ModelObject;

class ModelObjectStorage
{
    SQLiteWrapper * m_database;

public:
    explicit ModelObjectStorage(SQLiteWrapper * db_manager);
    ~ModelObjectStorage();

    static void CreateTables(SQLiteWrapper & database);

    void Save(const ModelObject & obj, const std::string & key_tag);
    std::unique_ptr<ModelObject> Load(const std::string & key_tag);
};

namespace model_storage {

inline constexpr std::array<std::string_view, 13> kCanonicalTableNames{
    "model_object",
    "model_chain_map",
    "model_component",
    "model_component_atom",
    "model_component_bond",
    "model_atom",
    "model_bond",
    "model_atom_local_potential",
    "model_bond_local_potential",
    "model_atom_posterior",
    "model_bond_posterior",
    "model_atom_group_potential",
    "model_bond_group_potential"
};

} // namespace model_storage

} // namespace rhbm_gem
