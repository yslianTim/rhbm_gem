#pragma once

#include <array>
#include <cstring>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "data/detail/ModelObjectAssembly.hpp"
#include "io/sqlite/SQLiteWrapper.hpp"
#include "support/CommandTestHelpers.hpp"

namespace data_test {

namespace rg = rhbm_gem;

inline std::shared_ptr<rg::ModelObject> LoadFixtureModel(
    const std::filesystem::path & model_path,
    const std::string & key_tag = "model")
{
    auto model{ rg::ReadModel(model_path) };
    model->SetKeyTag(key_tag);
    return std::shared_ptr<rg::ModelObject>{ std::move(model) };
}

inline std::unique_ptr<rg::ModelObject> MakeModelWithBond()
{
    rg::ModelObjectAssembly assembly;
    assembly.atom_list.reserve(2);

    auto atom_1{ std::make_unique<rg::AtomObject>() };
    atom_1->SetSerialID(1);
    atom_1->SetPosition(0.0f, 0.0f, 0.0f);

    auto atom_2{ std::make_unique<rg::AtomObject>() };
    atom_2->SetSerialID(2);
    atom_2->SetPosition(1.0f, 0.0f, 0.0f);

    auto * atom_1_ptr{ atom_1.get() };
    auto * atom_2_ptr{ atom_2.get() };
    assembly.atom_list.emplace_back(std::move(atom_1));
    assembly.atom_list.emplace_back(std::move(atom_2));

    assembly.bond_list.emplace_back(std::make_unique<rg::BondObject>(
        atom_1_ptr,
        atom_2_ptr));
    return std::make_unique<rg::ModelObject>(rg::AssembleModelObject(std::move(assembly)));
}

inline rg::MapObject MakeMapObject()
{
    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; ++i)
    {
        values[i] = static_cast<float>(i + 1);
    }
    return rg::MapObject{ grid_size, grid_spacing, origin, std::move(values) };
}

inline std::filesystem::path CopyFixtureWithNewName(
    const std::filesystem::path & source_path,
    const std::filesystem::path & output_path)
{
    std::filesystem::copy_file(
        source_path,
        output_path,
        std::filesystem::copy_options::overwrite_existing);
    return output_path;
}

inline rg::MapObject MakeTinyMapObject(float scale = 1.0f)
{
    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; ++i)
    {
        values[i] = static_cast<float>(i) * scale;
    }
    return rg::MapObject{ grid_size, grid_spacing, origin, std::move(values) };
}

inline void SaveTinyMapThroughRepository(
    rg::DataRepository & repository,
    const std::string & key_tag,
    float scale = 1.0f)
{
    auto map_object{ MakeTinyMapObject(scale) };
    repository.SaveMap(map_object, key_tag);
}

constexpr const char * kUpsertObjectMetadataSql =
    "INSERT INTO object_metadata(key_tag, object_type) VALUES (?, ?) "
    "ON CONFLICT(key_tag) DO UPDATE SET object_type = excluded.object_type;";

inline void UpsertObjectMetadata(
    rg::SQLiteWrapper & database,
    const std::string & key_tag,
    const std::string & object_type)
{
    database.Execute(
        "CREATE TABLE IF NOT EXISTS object_metadata (key_tag TEXT PRIMARY KEY, object_type TEXT);");
    database.Prepare(kUpsertObjectMetadataSql);
    rg::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.Bind<std::string>(2, object_type);
    database.StepOnce();
}

inline void CreateLegacyMapListWithoutForeignKeyTable(
    const std::filesystem::path & database_path,
    const rg::MapObject & map_object,
    const std::string & key_tag)
{
    rg::SQLiteWrapper database{ database_path };
    database.Execute("DROP TABLE IF EXISTS map_list;");
    database.Execute(
        "CREATE TABLE map_list ("
        "key_tag TEXT PRIMARY KEY, "
        "grid_size_x INTEGER, grid_size_y INTEGER, grid_size_z INTEGER, "
        "grid_spacing_x DOUBLE, grid_spacing_y DOUBLE, grid_spacing_z DOUBLE, "
        "origin_x DOUBLE, origin_y DOUBLE, origin_z DOUBLE, "
        "map_value_array BLOB"
        ");");

    const auto grid_size{ map_object.GetGridSize() };
    const auto grid_spacing{ map_object.GetGridSpacing() };
    const auto origin{ map_object.GetOrigin() };
    std::vector<float> values(map_object.GetMapValueArraySize());
    std::memcpy(values.data(), map_object.GetMapValueArray(), values.size() * sizeof(float));

    database.Prepare(
        "INSERT INTO map_list("
        "key_tag, grid_size_x, grid_size_y, grid_size_z, "
        "grid_spacing_x, grid_spacing_y, grid_spacing_z, "
        "origin_x, origin_y, origin_z, map_value_array"
        ") VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);");
    rg::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.Bind<int>(2, grid_size[0]);
    database.Bind<int>(3, grid_size[1]);
    database.Bind<int>(4, grid_size[2]);
    database.Bind<double>(5, grid_spacing[0]);
    database.Bind<double>(6, grid_spacing[1]);
    database.Bind<double>(7, grid_spacing[2]);
    database.Bind<double>(8, origin[0]);
    database.Bind<double>(9, origin[1]);
    database.Bind<double>(10, origin[2]);
    database.Bind<std::vector<float>>(11, values);
    database.StepOnce();
}

inline void CreateVersion2MetadataBasedMapShapeDatabase(const std::filesystem::path & database_path)
{
    auto map_object{ MakeTinyMapObject() };

    CreateLegacyMapListWithoutForeignKeyTable(database_path, map_object, "map_only");
    rg::SQLiteWrapper database{ database_path };
    UpsertObjectMetadata(database, "map_only", "map");
    database.Execute("PRAGMA user_version = 2;");
}

inline int GetUserVersion(const std::filesystem::path & database_path)
{
    rg::SQLiteWrapper database{ database_path };
    database.Prepare("PRAGMA user_version;");
    rg::SQLiteWrapper::StatementGuard guard(database);
    const auto rc{ database.StepNext() };
    if (rc != rg::SQLiteWrapper::StepRow())
    {
        return 0;
    }
    return database.GetColumn<int>(0);
}

inline void SetUserVersion(const std::filesystem::path & database_path, int user_version)
{
    rg::SQLiteWrapper database{ database_path };
    database.Execute("PRAGMA user_version = " + std::to_string(user_version) + ";");
}

inline void ExecuteSql(const std::filesystem::path & database_path, const std::string & sql)
{
    rg::SQLiteWrapper database{ database_path };
    database.Execute(sql);
}

inline void ExecuteSqlWithForeignKeysOff(
    const std::filesystem::path & database_path,
    const std::string & sql)
{
    rg::SQLiteWrapper database{ database_path };
    database.Execute("PRAGMA foreign_keys = OFF;");
    database.Execute(sql);
}

inline bool HasTable(const std::filesystem::path & database_path, const std::string & table_name)
{
    rg::SQLiteWrapper database{ database_path };
    database.Prepare(
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;");
    rg::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, table_name);
    return database.StepNext() == rg::SQLiteWrapper::StepRow();
}

inline int CountRows(
    const std::filesystem::path & database_path,
    const std::string & table_name,
    const std::string & key_tag = "")
{
    rg::SQLiteWrapper database{ database_path };
    const auto sql{
        key_tag.empty()
            ? "SELECT COUNT(*) FROM " + table_name + ";"
            : "SELECT COUNT(*) FROM " + table_name + " WHERE key_tag = ?;" };
    database.Prepare(sql);
    rg::SQLiteWrapper::StatementGuard guard(database);
    if (!key_tag.empty())
    {
        database.Bind<std::string>(1, key_tag);
    }
    if (database.StepNext() != rg::SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Failed to count rows in " + table_name);
    }
    return database.GetColumn<int>(0);
}

inline bool HasForeignKey(
    const std::filesystem::path & database_path,
    const std::string & table_name,
    const std::string & from_column,
    const std::string & referenced_table,
    const std::string & referenced_column,
    const std::string & on_delete)
{
    rg::SQLiteWrapper database{ database_path };
    database.Prepare("PRAGMA foreign_key_list(" + table_name + ");");
    rg::SQLiteWrapper::StatementGuard guard(database);
    while (true)
    {
        const auto rc{ database.StepNext() };
        if (rc == rg::SQLiteWrapper::StepDone())
        {
            return false;
        }
        if (rc != rg::SQLiteWrapper::StepRow())
        {
            throw std::runtime_error("Failed to inspect foreign_key_list for " + table_name);
        }
        if (database.GetColumn<std::string>(3) == from_column
            && database.GetColumn<std::string>(2) == referenced_table
            && database.GetColumn<std::string>(4) == referenced_column
            && database.GetColumn<std::string>(6) == on_delete)
        {
            return true;
        }
    }
}

} // namespace data_test
