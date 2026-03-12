#pragma once


#include <array>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include "CommandTestHelpers.hpp"
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/io/FileIO.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "support/DataInternalTestSeam.hpp"

namespace rg = rhbm_gem;

namespace data_object_io_schema_test {

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

    auto grid_size{ map_object.GetGridSize() };
    auto grid_spacing{ map_object.GetGridSpacing() };
    auto origin{ map_object.GetOrigin() };
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
    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; i++) values[i] = static_cast<float>(i);
    rg::MapObject map_object{ grid_size, grid_spacing, origin, std::move(values) };

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
    if (rc != rg::SQLiteWrapper::StepRow()) return 0;
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
            : "SELECT COUNT(*) FROM " + table_name + " WHERE key_tag = ?;"
    };
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
        if (database.GetColumn<std::string>(3) == from_column &&
            database.GetColumn<std::string>(2) == referenced_table &&
            database.GetColumn<std::string>(4) == referenced_column &&
            database.GetColumn<std::string>(6) == on_delete)
        {
            return true;
        }
    }
}

inline std::string SanitizeLegacyKey(const std::string & key_tag)
{
    std::string sanitized_key_tag;
    sanitized_key_tag.reserve(key_tag.size());
    for (char ch : key_tag)
    {
        if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_')
        {
            sanitized_key_tag.push_back(ch);
        }
        else
        {
            sanitized_key_tag.push_back('_');
        }
    }
    return sanitized_key_tag;
}

inline std::filesystem::path LegacyFixturePath()
{
    return command_test::TestDataPath("legacy/legacy_model_v1.sqlite");
}

inline void CopyLegacyFixtureDatabase(const std::filesystem::path & database_path)
{
    std::filesystem::copy_file(
        LegacyFixturePath(),
        database_path,
        std::filesystem::copy_options::overwrite_existing);
}

inline void RenameLegacyModel(
    const std::filesystem::path & database_path,
    const std::string & old_key_tag,
    const std::string & new_key_tag,
    const std::string & new_pdb_id)
{
    const auto old_sanitized{ SanitizeLegacyKey(old_key_tag) };
    const auto new_sanitized{ SanitizeLegacyKey(new_key_tag) };
    ExecuteSql(
        database_path,
        "ALTER TABLE atom_list_in_" + old_sanitized + " RENAME TO atom_list_in_" + new_sanitized +
            ";");
    ExecuteSql(
        database_path,
        "ALTER TABLE bond_list_in_" + old_sanitized + " RENAME TO bond_list_in_" + new_sanitized +
            ";");
    ExecuteSql(
        database_path,
        "UPDATE model_list SET key_tag = '" + new_key_tag + "', pdb_id = '" + new_pdb_id +
            "' WHERE key_tag = '" + old_key_tag + "';");
    ExecuteSql(
        database_path,
        "UPDATE object_metadata SET key_tag = '" + new_key_tag +
            "' WHERE key_tag = '" + old_key_tag + "' AND object_type = 'model';");
}

inline void DuplicateLegacyModel(
    const std::filesystem::path & database_path,
    const std::string & source_key_tag,
    const std::string & new_key_tag,
    const std::string & new_pdb_id)
{
    const auto source_sanitized{ SanitizeLegacyKey(source_key_tag) };
    const auto new_sanitized{ SanitizeLegacyKey(new_key_tag) };
    ExecuteSql(
        database_path,
        "CREATE TABLE atom_list_in_" + new_sanitized +
            " AS SELECT * FROM atom_list_in_" + source_sanitized + ";");
    ExecuteSql(
        database_path,
        "CREATE TABLE bond_list_in_" + new_sanitized +
            " AS SELECT * FROM bond_list_in_" + source_sanitized + ";");
    ExecuteSql(
        database_path,
        "INSERT INTO model_list(key_tag, atom_size, pdb_id, emd_id, map_resolution, resolution_method) "
        "SELECT '" + new_key_tag + "', atom_size, '" + new_pdb_id +
            "', emd_id, map_resolution, resolution_method "
        "FROM model_list WHERE key_tag = '" + source_key_tag + "';");
    ExecuteSql(
        database_path,
        "INSERT INTO object_metadata(key_tag, object_type) VALUES ('" + new_key_tag + "', 'model');");
}

inline std::shared_ptr<rg::ModelObject> LoadFixtureModel(
    const std::filesystem::path & model_path,
    const std::string & key_tag = "model")
{
    rg::DataObjectManager manager{};
    manager.ProcessFile(model_path, key_tag);
    return manager.GetTypedDataObject<rg::ModelObject>(key_tag);
}

inline rg::MapObject MakeTinyMapObject(float scale = 1.0f)
{
    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; i++) values[i] = static_cast<float>(i) * scale;
    return rg::MapObject{ grid_size, grid_spacing, origin, std::move(values) };
}

inline void SaveTinyMapThroughManager(
    rg::DataObjectManager & manager,
    const std::filesystem::path & map_path,
    const std::string & key_tag,
    float scale = 1.0f)
{
    auto map_object{ MakeTinyMapObject(scale) };
    rg::WriteMap(map_path, map_object);
    manager.ProcessFile(map_path, key_tag);
    manager.SaveDataObject(key_tag);
}

class FailingDataObject final : public rg::DataObjectBase
{
    std::string m_key_tag;

public:
    std::unique_ptr<rg::DataObjectBase> Clone() const override
    {
        return std::make_unique<FailingDataObject>(*this);
    }

    void Display() const override {}
    void Update() override {}
    void SetKeyTag(const std::string & label) override { m_key_tag = label; }
    std::string GetKeyTag() const override { return m_key_tag; }
};

} // namespace
