#include <gtest/gtest.h>

#include <array>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "CommandTestHelpers.hpp"
#include "DataObjectBase.hpp"
#include "DataObjectDAOBase.hpp"
#include "DataObjectDAOFactoryRegistry.hpp"
#include "DataObjectManager.hpp"
#include "DatabaseManager.hpp"
#include "FileFormatRegistry.hpp"
#include "FileProcessFactoryBase.hpp"
#include "FileProcessFactoryRegistry.hpp"
#include "FileProcessFactoryResolver.hpp"
#include "LocalPotentialEntry.hpp"
#include "MapFileReader.hpp"
#include "MapFileWriter.hpp"
#include "MapObject.hpp"
#include "MapObjectDAO.hpp"
#include "ModelFileReader.hpp"
#include "ModelObject.hpp"
#include "ModelObjectDAO.hpp"
#include "SQLiteWrapper.hpp"

#include "../src/data/model_io/ModelSchemaSql.hpp"
#include "../src/data/persistence/MapStoreSql.hpp"

namespace rg = rhbm_gem;

namespace {

constexpr const char * kUpsertObjectMetadataSql =
    "INSERT INTO object_metadata(key_tag, object_type) VALUES (?, ?) "
    "ON CONFLICT(key_tag) DO UPDATE SET object_type = excluded.object_type;";

std::string ReplaceAll(std::string value, std::string_view from, std::string_view to)
{
    std::size_t search_pos{ 0 };
    while ((search_pos = value.find(from, search_pos)) != std::string::npos)
    {
        value.replace(search_pos, from.size(), to);
        search_pos += to.size();
    }
    return value;
}

std::string BuildLegacyMetadataBasedCreateSql(std::string sql)
{
    sql = ReplaceAll(
        std::move(sql),
        " REFERENCES object_catalog(key_tag) ON DELETE CASCADE",
        "");
    sql = ReplaceAll(
        std::move(sql),
        ",\n        FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE",
        "");
    return sql;
}

void CreateLegacyMetadataBasedV2Tables(
    rg::SQLiteWrapper & database,
    bool include_model_tables,
    bool include_map_table)
{
    database.Execute(
        "CREATE TABLE IF NOT EXISTS object_metadata (key_tag TEXT PRIMARY KEY, object_type TEXT);");
    if (include_model_tables)
    {
        for (const auto create_sql : rg::model_io::kCreateModelTableSqlList)
        {
            database.Execute(BuildLegacyMetadataBasedCreateSql(std::string(create_sql)));
        }
    }
    if (include_map_table)
    {
        database.Execute(
            BuildLegacyMetadataBasedCreateSql(std::string(rg::persistence::kCreateMapTableSql)));
    }
}

void UpsertObjectMetadata(
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

void CreateLegacyMetadataBasedV2MapDatabase(
    const std::filesystem::path & database_path,
    const rg::MapObject & map_object,
    int user_version)
{
    rg::SQLiteWrapper database{ database_path };
    CreateLegacyMetadataBasedV2Tables(database, false, true);
    rg::MapObjectDAO map_dao{ &database };
    map_dao.Save(&map_object, "map_only");
    UpsertObjectMetadata(database, "map_only", "map");
    database.Execute("PRAGMA user_version = " + std::to_string(user_version) + ";");
}

void CreateLegacyMetadataBasedV2ModelDatabase(
    const std::filesystem::path & database_path,
    const rg::ModelObject & model_object,
    const std::string & key_tag,
    int user_version)
{
    rg::SQLiteWrapper database{ database_path };
    CreateLegacyMetadataBasedV2Tables(database, true, false);
    rg::ModelObjectDAO model_dao{ &database };
    model_dao.Save(&model_object, key_tag);
    UpsertObjectMetadata(database, key_tag, "model");
    database.Execute("PRAGMA user_version = " + std::to_string(user_version) + ";");
}

int GetUserVersion(const std::filesystem::path & database_path)
{
    rg::SQLiteWrapper database{ database_path };
    database.Prepare("PRAGMA user_version;");
    rg::SQLiteWrapper::StatementGuard guard(database);
    const auto rc{ database.StepNext() };
    if (rc != rg::SQLiteWrapper::StepRow()) return 0;
    return database.GetColumn<int>(0);
}

void SetUserVersion(const std::filesystem::path & database_path, int user_version)
{
    rg::SQLiteWrapper database{ database_path };
    database.Execute("PRAGMA user_version = " + std::to_string(user_version) + ";");
}

void ExecuteSql(const std::filesystem::path & database_path, const std::string & sql)
{
    rg::SQLiteWrapper database{ database_path };
    database.Execute(sql);
}

bool HasTable(const std::filesystem::path & database_path, const std::string & table_name)
{
    rg::SQLiteWrapper database{ database_path };
    database.Prepare(
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;");
    rg::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, table_name);
    return database.StepNext() == rg::SQLiteWrapper::StepRow();
}

int CountRows(
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

std::string SanitizeLegacyKey(const std::string & key_tag)
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

std::filesystem::path LegacyFixturePath()
{
    return command_test::TestDataPath("legacy/legacy_model_v1.sqlite");
}

void CopyLegacyFixtureDatabase(const std::filesystem::path & database_path)
{
    std::filesystem::copy_file(
        LegacyFixturePath(),
        database_path,
        std::filesystem::copy_options::overwrite_existing);
}

void RenameLegacyModel(
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

void DuplicateLegacyModel(
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

std::shared_ptr<rg::ModelObject> LoadFixtureModel(
    const std::filesystem::path & model_path,
    const std::string & key_tag = "model")
{
    rg::DataObjectManager manager;
    manager.ProcessFile(model_path, key_tag);
    return manager.GetTypedDataObject<rg::ModelObject>(key_tag);
}

rg::MapObject MakeTinyMapObject(float scale = 1.0f)
{
    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; i++) values[i] = static_cast<float>(i) * scale;
    return rg::MapObject{ grid_size, grid_spacing, origin, std::move(values) };
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
    void Accept(rg::DataObjectVisitorBase * /*visitor*/) override {}
    void SetKeyTag(const std::string & label) override { m_key_tag = label; }
    std::string GetKeyTag() const override { return m_key_tag; }
};

class FailingDataObjectDAO final : public rg::DataObjectDAOBase
{
    rg::SQLiteWrapper * m_database;

public:
    explicit FailingDataObjectDAO(rg::SQLiteWrapper * database) : m_database{ database } {}

    void Save(const rg::DataObjectBase * /*obj*/, const std::string & key_tag) override
    {
        m_database->Execute("CREATE TABLE IF NOT EXISTS failing_payload (key_tag TEXT PRIMARY KEY);");
        m_database->Prepare("INSERT INTO failing_payload(key_tag) VALUES (?);");
        rg::SQLiteWrapper::StatementGuard guard(*m_database);
        m_database->Bind<std::string>(1, key_tag);
        m_database->StepOnce();
        throw std::runtime_error("Injected DAO failure");
    }

    std::unique_ptr<rg::DataObjectBase> Load(const std::string & /*key_tag*/) override
    {
        throw std::runtime_error("Load not implemented for failing DAO");
    }
};

rg::DataObjectDAORegistrar<FailingDataObject, FailingDataObjectDAO>
    registrar_failing_dao("failing");

class OverrideFileFactory final : public rg::FileProcessFactoryBase
{
public:
    std::unique_ptr<rg::DataObjectBase> CreateDataObject(const std::string & /*filename*/) override
    {
        return std::make_unique<FailingDataObject>();
    }

    void OutputDataObject(
        const std::string & /*filename*/, rg::DataObjectBase * /*data_object*/) override
    {
    }
};

class ScopedFactoryOverride
{
    std::string m_extension;

public:
    ScopedFactoryOverride(
        const std::string & extension,
        std::function<std::unique_ptr<rg::FileProcessFactoryBase>()> creator) :
        m_extension{ extension }
    {
        rg::FileProcessFactoryRegistry::Instance().RegisterFactory(m_extension, std::move(creator));
    }

    ~ScopedFactoryOverride()
    {
        rg::FileProcessFactoryRegistry::Instance().UnregisterFactory(m_extension);
    }
};

} // namespace

TEST(DataObjectIOSchemaTest, EmptyDatabaseBootstrapsNormalizedSchema)
{
    const command_test::ScopedTempDir temp_dir{ "schema_bootstrap" };
    const auto database_path{ temp_dir.path() / "bootstrap.sqlite" };

    rg::DatabaseManager database_manager{ database_path };

    EXPECT_EQ(database_manager.GetSchemaVersion(), rg::DatabaseSchemaVersion::NormalizedV2);
    EXPECT_TRUE(HasTable(database_path, "object_catalog"));
    EXPECT_FALSE(HasTable(database_path, "object_metadata"));
    EXPECT_TRUE(HasTable(database_path, "model_object"));
    EXPECT_TRUE(HasTable(database_path, "map_list"));
    EXPECT_EQ(GetUserVersion(database_path), 2);
}

TEST(DataObjectIOSchemaTest, SaveRollbackLeavesNoPayloadOrMetadata)
{
    const command_test::ScopedTempDir temp_dir{ "schema_atomic_save" };
    const auto database_path{ temp_dir.path() / "atomic.sqlite" };

    rg::DatabaseManager database_manager{ database_path };
    FailingDataObject data_object;
    data_object.SetKeyTag("failing");

    EXPECT_THROW(database_manager.SaveDataObject(&data_object, "failing"), std::runtime_error);
    EXPECT_FALSE(HasTable(database_path, "failing_payload"));
    EXPECT_EQ(CountRows(database_path, "object_catalog"), 0);
}

TEST(DataObjectIOSchemaTest, LegacyFixtureMigratesToNormalizedSchema)
{
    const command_test::ScopedTempDir temp_dir{ "schema_migrate_legacy" };
    const auto database_path{ temp_dir.path() / "legacy.sqlite" };

    CopyLegacyFixtureDatabase(database_path);
    ASSERT_EQ(GetUserVersion(database_path), 0);
    ASSERT_TRUE(HasTable(database_path, "model_list"));

    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    EXPECT_EQ(
        manager.GetDatabaseManager()->GetSchemaVersion(),
        rg::DatabaseSchemaVersion::NormalizedV2);

    ASSERT_NO_THROW(manager.LoadDataObject("legacy_model"));
    auto model{ manager.GetTypedDataObject<rg::ModelObject>("legacy_model") };
    EXPECT_EQ(model->GetPdbID(), "LEGACY");
    EXPECT_TRUE(HasTable(database_path, "model_object"));
    EXPECT_FALSE(HasTable(database_path, "model_list"));
    EXPECT_EQ(GetUserVersion(database_path), 2);
}

TEST(DataObjectIOSchemaTest, FailedLegacyMigrationRollsBackToLegacyState)
{
    const command_test::ScopedTempDir temp_dir{ "schema_migrate_rollback" };
    const auto database_path{ temp_dir.path() / "legacy_rollback.sqlite" };

    CopyLegacyFixtureDatabase(database_path);
    RenameLegacyModel(database_path, "legacy_model", "key_a", "MODEL_A");
    DuplicateLegacyModel(database_path, "key_a", "key_b", "MODEL_B");
    ExecuteSql(
        database_path,
        "UPDATE model_list SET atom_size = 999 WHERE key_tag = 'key_b';");

    rg::DataObjectManager manager;
    EXPECT_THROW(manager.SetDatabaseManager(database_path), std::runtime_error);
    EXPECT_TRUE(HasTable(database_path, "model_list"));
    EXPECT_FALSE(HasTable(database_path, "model_object"));
    EXPECT_EQ(GetUserVersion(database_path), 0);
}

TEST(DataObjectIOSchemaTest, MapOnlyLegacyMetadataBasedV2DatabaseRemainsLoadable)
{
    const command_test::ScopedTempDir temp_dir{ "schema_map_only" };
    const auto database_path{ temp_dir.path() / "map_only.sqlite" };
    const auto expected_map{ MakeTinyMapObject() };
    const auto expected_grid_size{ expected_map.GetGridSize() };

    CreateLegacyMetadataBasedV2MapDatabase(database_path, expected_map, 2);

    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    EXPECT_EQ(
        manager.GetDatabaseManager()->GetSchemaVersion(),
        rg::DatabaseSchemaVersion::NormalizedV2);
    ASSERT_NO_THROW(manager.LoadDataObject("map_only"));
    auto loaded_map{ manager.GetTypedDataObject<rg::MapObject>("map_only") };
    EXPECT_EQ(loaded_map->GetGridSize(), expected_grid_size);
    EXPECT_TRUE(HasTable(database_path, "object_catalog"));
    EXPECT_FALSE(HasTable(database_path, "object_metadata"));
    EXPECT_EQ(GetUserVersion(database_path), 2);
}

TEST(DataObjectIOSchemaTest, UnknownSchemaVersionThrows)
{
    const command_test::ScopedTempDir temp_dir{ "schema_unknown_version" };
    const auto database_path{ temp_dir.path() / "unknown.sqlite" };

    SetUserVersion(database_path, 99);
    EXPECT_THROW((void)rg::DatabaseManager{ database_path }, std::runtime_error);
}

TEST(DataObjectIOSchemaTest, ManagedButUnversionedMapOnlyDatabaseRepairsAndUpgradesToV2)
{
    const command_test::ScopedTempDir temp_dir{ "schema_unversioned_map" };
    const auto database_path{ temp_dir.path() / "managed_unversioned.sqlite" };
    const auto expected_map{ MakeTinyMapObject(2.0f) };
    const auto expected_grid_size{ expected_map.GetGridSize() };

    CreateLegacyMetadataBasedV2MapDatabase(database_path, expected_map, 0);

    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    EXPECT_EQ(
        manager.GetDatabaseManager()->GetSchemaVersion(),
        rg::DatabaseSchemaVersion::NormalizedV2);
    ASSERT_NO_THROW(manager.LoadDataObject("map_only"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::MapObject>("map_only")->GetGridSize(), expected_grid_size);
    EXPECT_EQ(CountRows(database_path, "object_catalog"), 1);
    EXPECT_FALSE(HasTable(database_path, "object_metadata"));
}

TEST(DataObjectIOSchemaTest, LegacyMetadataBasedV2DatabaseUpgradesToFinalV2)
{
    const command_test::ScopedTempDir temp_dir{ "schema_legacy_v2_upgrade" };
    const auto database_path{ temp_dir.path() / "legacy_v2.sqlite" };
    const auto model{ LoadFixtureModel(command_test::TestDataPath("test_model.cif")) };

    CreateLegacyMetadataBasedV2ModelDatabase(database_path, *model, "legacy_v2_model", 2);

    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    ASSERT_NO_THROW(manager.LoadDataObject("legacy_v2_model"));
    EXPECT_EQ(
        manager.GetTypedDataObject<rg::ModelObject>("legacy_v2_model")->GetPdbID(),
        model->GetPdbID());
    EXPECT_TRUE(HasTable(database_path, "object_catalog"));
    EXPECT_FALSE(HasTable(database_path, "object_metadata"));
    EXPECT_EQ(CountRows(database_path, "object_catalog"), 1);
}

TEST(DataObjectIOSchemaTest, CatalogLookupLoadsCorrectDaoType)
{
    const command_test::ScopedTempDir temp_dir{ "schema_catalog_lookup" };
    const auto database_path{ temp_dir.path() / "catalog_lookup.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto map_object{ MakeTinyMapObject() };

    rg::DataObjectManager manager;
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    manager.SaveDataObject("model");
    manager.GetDatabaseManager()->SaveDataObject(&map_object, "map");
    manager.ClearDataObjects();

    ASSERT_NO_THROW(manager.LoadDataObject("model"));
    ASSERT_NO_THROW(manager.LoadDataObject("map"));
    EXPECT_NE(dynamic_cast<rg::ModelObject *>(manager.GetDataObject("model").get()), nullptr);
    EXPECT_NE(dynamic_cast<rg::MapObject *>(manager.GetDataObject("map").get()), nullptr);
}

TEST(DataObjectIOSchemaTest, GhostMetadataRowsDoNotSurviveFinalV2Upgrade)
{
    const command_test::ScopedTempDir temp_dir{ "schema_legacy_v2_ghost" };
    const auto database_path{ temp_dir.path() / "legacy_v2_ghost.sqlite" };
    const auto model{ LoadFixtureModel(command_test::TestDataPath("test_model.cif")) };

    CreateLegacyMetadataBasedV2ModelDatabase(database_path, *model, "real_model", 2);
    {
        rg::SQLiteWrapper database{ database_path };
        UpsertObjectMetadata(database, "ghost_model", "model");
    }

    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    ASSERT_NO_THROW(manager.LoadDataObject("real_model"));
    EXPECT_THROW(manager.LoadDataObject("ghost_model"), std::runtime_error);
    EXPECT_EQ(CountRows(database_path, "object_catalog"), 1);
    EXPECT_FALSE(HasTable(database_path, "object_metadata"));
}

TEST(DataObjectIOSchemaTest, NormalizedV2DatabaseMissingRequiredTableThrows)
{
    const command_test::ScopedTempDir temp_dir{ "schema_missing_v2_table" };
    const auto database_path{ temp_dir.path() / "missing_v2.sqlite" };

    {
        rg::DatabaseManager database_manager{ database_path };
        EXPECT_EQ(database_manager.GetSchemaVersion(), rg::DatabaseSchemaVersion::NormalizedV2);
    }
    ExecuteSql(database_path, "DROP TABLE model_bond_group_potential;");
    SetUserVersion(database_path, 2);

    EXPECT_THROW((void)rg::DatabaseManager{ database_path }, std::runtime_error);
}

TEST(DataObjectIOSchemaTest, DistinctUnsanitizedKeysDoNotCollideInV2Schema)
{
    const command_test::ScopedTempDir temp_dir{ "schema_key_collision" };
    const auto database_path{ temp_dir.path() / "collision.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataObjectManager manager;
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "mem_a");
    manager.ProcessFile(model_path, "mem_b");
    manager.GetTypedDataObject<rg::ModelObject>("mem_a")->SetPdbID("MODEL_A");
    manager.GetTypedDataObject<rg::ModelObject>("mem_b")->SetPdbID("MODEL_B");

    manager.SaveDataObject("mem_a", "a-b");
    manager.SaveDataObject("mem_b", "a_b");
    manager.ClearDataObjects();

    manager.LoadDataObject("a-b");
    manager.LoadDataObject("a_b");
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("a-b")->GetPdbID(), "MODEL_A");
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("a_b")->GetPdbID(), "MODEL_B");
    EXPECT_EQ(CountRows(database_path, "model_object"), 2);
}

TEST(DataObjectIOSchemaTest, LegacyMigrationUsesModelListWhenMetadataIncomplete)
{
    const command_test::ScopedTempDir temp_dir{ "schema_partial_metadata" };
    const auto database_path{ temp_dir.path() / "partial_metadata.sqlite" };

    CopyLegacyFixtureDatabase(database_path);
    RenameLegacyModel(database_path, "legacy_model", "key_a", "MODEL_A");
    DuplicateLegacyModel(database_path, "key_a", "key_b", "MODEL_B");
    ExecuteSql(
        database_path,
        "DELETE FROM object_metadata WHERE key_tag = 'key_b' AND object_type = 'model';");

    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    ASSERT_NO_THROW(manager.LoadDataObject("key_a"));
    ASSERT_NO_THROW(manager.LoadDataObject("key_b"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("key_a")->GetPdbID(), "MODEL_A");
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("key_b")->GetPdbID(), "MODEL_B");
    EXPECT_EQ(CountRows(database_path, "model_object"), 2);
    EXPECT_EQ(CountRows(database_path, "object_catalog"), 2);
    EXPECT_FALSE(HasTable(database_path, "object_metadata"));
}

TEST(DataObjectIOSchemaTest, LegacyMigrationIgnoresMetadataOnlyGhostKeys)
{
    const command_test::ScopedTempDir temp_dir{ "schema_ghost_metadata" };
    const auto database_path{ temp_dir.path() / "ghost_metadata.sqlite" };

    CopyLegacyFixtureDatabase(database_path);
    RenameLegacyModel(database_path, "legacy_model", "real_model", "REAL");
    {
        rg::SQLiteWrapper database{ database_path };
        UpsertObjectMetadata(database, "ghost_model", "model");
    }

    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    ASSERT_NO_THROW(manager.LoadDataObject("real_model"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("real_model")->GetPdbID(), "REAL");
    EXPECT_THROW(manager.LoadDataObject("ghost_model"), std::runtime_error);
    EXPECT_EQ(CountRows(database_path, "model_object"), 1);
    EXPECT_EQ(CountRows(database_path, "object_catalog"), 1);
    EXPECT_FALSE(HasTable(database_path, "object_metadata"));
}

TEST(DataObjectIOSchemaTest, LegacyMigrationDropsOnlyOwnedTables)
{
    const command_test::ScopedTempDir temp_dir{ "schema_owned_table_drop" };
    const auto database_path{ temp_dir.path() / "owned_drop.sqlite" };

    CopyLegacyFixtureDatabase(database_path);
    ExecuteSql(database_path, "CREATE TABLE custom_atom_list_in_notes (value INTEGER);");

    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    EXPECT_TRUE(HasTable(database_path, "custom_atom_list_in_notes"));
    EXPECT_FALSE(HasTable(database_path, "atom_list_in_legacy_model"));
}

TEST(DataObjectIOSchemaTest, MixedUnknownSchemaFailsFast)
{
    const command_test::ScopedTempDir temp_dir{ "schema_mixed_unknown" };
    const auto database_path{ temp_dir.path() / "mixed_unknown.sqlite" };

    ExecuteSql(database_path, "CREATE TABLE foreign_table (id INTEGER PRIMARY KEY);");
    EXPECT_THROW((void)rg::DatabaseManager{ database_path }, std::runtime_error);
}

TEST(DataObjectIOSchemaTest, ForeignKeyRejectsOrphanModelChildRows)
{
    const command_test::ScopedTempDir temp_dir{ "schema_fk_orphan" };
    const auto database_path{ temp_dir.path() / "orphan.sqlite" };
    rg::DatabaseManager database_manager{ database_path };

    EXPECT_THROW(
        database_manager.GetDatabase()->Execute(
            "INSERT INTO model_chain_map(key_tag, entity_id, chain_ordinal, chain_id) "
            "VALUES ('missing', '1', 0, 'A');"),
        std::runtime_error);
}

TEST(DataObjectIOSchemaTest, DeletingCatalogRootCascadesPayloadRows)
{
    const command_test::ScopedTempDir temp_dir{ "schema_catalog_cascade" };
    const auto database_path{ temp_dir.path() / "cascade.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataObjectManager manager;
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    manager.SaveDataObject("model");

    manager.GetDatabaseManager()->GetDatabase()->Execute(
        "DELETE FROM object_catalog WHERE key_tag = 'model';");

    EXPECT_EQ(CountRows(database_path, "object_catalog"), 0);
    EXPECT_EQ(CountRows(database_path, "model_object"), 0);
    EXPECT_EQ(CountRows(database_path, "model_atom"), 0);
}

TEST(DataObjectIOSchemaTest, ChainMetadataPersistsAcrossDatabaseRoundTrip)
{
    const command_test::ScopedTempDir temp_dir{ "schema_chain_roundtrip" };
    const auto database_path{ temp_dir.path() / "chain_roundtrip.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model_keyvalue_entity.cif") };

    rg::DataObjectManager manager;
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    auto original_model{ manager.GetTypedDataObject<rg::ModelObject>("model") };
    const auto original_chain_map{ original_model->GetChainIDListMap() };

    manager.SaveDataObject("model");
    manager.ClearDataObjects();
    manager.LoadDataObject("model");

    auto loaded_model{ manager.GetTypedDataObject<rg::ModelObject>("model") };
    EXPECT_EQ(loaded_model->GetChainIDListMap(), original_chain_map);
    EXPECT_GT(CountRows(database_path, "model_chain_map", "model"), 0);
}

TEST(DataObjectIOSchemaTest, SymmetryFilteringMatchesAfterDatabaseReload)
{
    const command_test::ScopedTempDir temp_dir{ "schema_symmetry_roundtrip" };
    const auto database_path{ temp_dir.path() / "symmetry.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model_keyvalue_entity.cif") };

    auto original_model{ LoadFixtureModel(model_path) };
    for (const auto & atom : original_model->GetAtomList()) atom->SetSelectedFlag(true);
    original_model->FilterAtomFromSymmetry(false);
    original_model->Update();
    const auto original_selected_count{ original_model->GetNumberOfSelectedAtom() };

    rg::DataObjectManager manager;
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    for (const auto & atom : manager.GetTypedDataObject<rg::ModelObject>("model")->GetAtomList())
    {
        atom->SetSelectedFlag(true);
    }
    manager.GetTypedDataObject<rg::ModelObject>("model")->Update();
    manager.SaveDataObject("model");
    manager.ClearDataObjects();
    manager.LoadDataObject("model");

    auto loaded_model{ manager.GetTypedDataObject<rg::ModelObject>("model") };
    for (const auto & atom : loaded_model->GetAtomList()) atom->SetSelectedFlag(true);
    loaded_model->FilterAtomFromSymmetry(false);
    loaded_model->Update();
    EXPECT_EQ(loaded_model->GetNumberOfSelectedAtom(), original_selected_count);
}

TEST(DataObjectIOSchemaTest, FileImportTransfersBondKeySystem)
{
    rg::DataObjectManager manager;
    const auto model_path{
        command_test::TestDataPath("test_model_auth_seq_alnum_struct_conn.cif") };

    ASSERT_NO_THROW(manager.ProcessFile(model_path, "model"));
    auto model{ manager.GetTypedDataObject<rg::ModelObject>("model") };
    ASSERT_GE(model->GetNumberOfBond(), 1);
    ASSERT_NE(model->GetBondKeySystemPtr(), nullptr);
    EXPECT_NE(model->GetBondKeySystemPtr()->GetBondId(model->GetBondList().front()->GetBondKey()), "UNK");
}

TEST(DataObjectIOSchemaTest, UppercaseExtensionsDispatchCorrectly)
{
    const command_test::ScopedTempDir temp_dir{ "schema_uppercase_ext" };
    const auto source_model_path{ command_test::TestDataPath("test_model.cif") };
    const auto uppercase_model_path{ temp_dir.path() / "TEST_MODEL.MMCIF" };
    std::filesystem::copy_file(source_model_path, uppercase_model_path);

    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.ProcessFile(uppercase_model_path, "model"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("model")->GetNumberOfAtom(), 1);

    const auto map_object{ MakeTinyMapObject() };
    const auto uppercase_map_path{ temp_dir.path() / "TEST_MAP.MAP" };
    rg::MapFileWriter writer{ uppercase_map_path.string(), &map_object };
    ASSERT_NO_THROW(writer.Write());
    rg::MapFileReader reader{ uppercase_map_path.string() };
    ASSERT_NO_THROW(reader.Read());
    EXPECT_EQ(reader.GetGridSizeArray(), map_object.GetGridSize());
}

TEST(DataObjectIOSchemaTest, BuiltInFormatsStillWorkWithoutAnyRegistrationStep)
{
    auto factory{ rg::FileProcessFactoryRegistry::Instance().CreateFactory(".cif") };
    auto data_object{ factory->CreateDataObject(command_test::TestDataPath("test_model.cif").string()) };
    EXPECT_NE(dynamic_cast<rg::ModelObject *>(data_object.get()), nullptr);
}

TEST(DataObjectIOSchemaTest, DataObjectManagerDoesNotRequireDefaultFactoryRegistration)
{
    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.ProcessFile(command_test::TestDataPath("test_model.cif"), "model"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("model")->GetNumberOfAtom(), 1);
}

TEST(DataObjectIOSchemaTest, InjectedOverrideTakesPrecedenceOverDescriptorFallback)
{
    auto resolver{ std::make_shared<rg::OverrideableFileProcessFactoryResolver>() };
    resolver->RegisterFactory(".cif", []() { return std::make_unique<OverrideFileFactory>(); });

    rg::DataObjectManager manager{ resolver };
    ASSERT_NO_THROW(manager.ProcessFile(command_test::TestDataPath("test_model.cif"), "override"));
    EXPECT_NE(dynamic_cast<FailingDataObject *>(manager.GetDataObject("override").get()), nullptr);
}

TEST(DataObjectIOSchemaTest, FileProcessFactoryRegistryOverrideCanBeRemoved)
{
    auto & registry{ rg::FileProcessFactoryRegistry::Instance() };
    registry.RegisterFactory(".cif", []() { return std::make_unique<OverrideFileFactory>(); });

    auto override_factory{ registry.CreateFactory(".cif") };
    auto override_object{
        override_factory->CreateDataObject(command_test::TestDataPath("test_model.cif").string()) };
    EXPECT_NE(dynamic_cast<FailingDataObject *>(override_object.get()), nullptr);

    registry.UnregisterFactory(".cif");

    auto default_factory{ registry.CreateFactory(".cif") };
    auto default_object{
        default_factory->CreateDataObject(command_test::TestDataPath("test_model.cif").string()) };
    EXPECT_NE(dynamic_cast<rg::ModelObject *>(default_object.get()), nullptr);
}

TEST(DataObjectIOSchemaTest, CustomFactoryOverrideTakesPrecedenceOverDescriptorFallback)
{
    ScopedFactoryOverride override{
        ".cif",
        []() { return std::make_unique<OverrideFileFactory>(); }
    };

    auto factory{ rg::FileProcessFactoryRegistry::Instance().CreateFactory(".cif") };
    auto data_object{
        factory->CreateDataObject(command_test::TestDataPath("test_model.cif").string()) };
    EXPECT_NE(dynamic_cast<FailingDataObject *>(data_object.get()), nullptr);
}

TEST(DataObjectIOSchemaTest, CustomFactoryOverrideDoesNotLeakAcrossTests)
{
    {
        ScopedFactoryOverride override{
            ".cif",
            []() { return std::make_unique<OverrideFileFactory>(); }
        };
        auto factory{ rg::FileProcessFactoryRegistry::Instance().CreateFactory(".cif") };
        auto data_object{
            factory->CreateDataObject(command_test::TestDataPath("test_model.cif").string()) };
        EXPECT_NE(dynamic_cast<FailingDataObject *>(data_object.get()), nullptr);
    }

    auto factory{ rg::FileProcessFactoryRegistry::Instance().CreateFactory(".cif") };
    auto data_object{ factory->CreateDataObject(command_test::TestDataPath("test_model.cif").string()) };
    EXPECT_NE(dynamic_cast<rg::ModelObject *>(data_object.get()), nullptr);
}

TEST(DataObjectIOSchemaTest, FileFormatDescriptorsAreUniqueAndConsistent)
{
    std::unordered_set<std::string> extension_set;
    for (const auto & descriptor : rg::FileFormatRegistry::Instance().GetAllDescriptors())
    {
        EXPECT_TRUE(extension_set.insert(descriptor.extension).second);
        if (descriptor.kind == rg::DataObjectKind::Model)
        {
            EXPECT_TRUE(descriptor.model_backend.has_value());
            EXPECT_FALSE(descriptor.map_backend.has_value());
        }
        else
        {
            EXPECT_TRUE(descriptor.map_backend.has_value());
            EXPECT_FALSE(descriptor.model_backend.has_value());
        }

        if (descriptor.extension == ".mmcif" || descriptor.extension == ".mcif")
        {
            EXPECT_TRUE(descriptor.supports_read);
            EXPECT_FALSE(descriptor.supports_write);
        }
    }
}

TEST(DataObjectIOSchemaTest, ModelWriteSupportMatrixRejectsMmcifOutput)
{
    const command_test::ScopedTempDir temp_dir{ "schema_output_matrix" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto output_path{ temp_dir.path() / "unsupported_output.mmcif" };

    rg::DataObjectManager manager;
    manager.ProcessFile(model_path, "model");
    EXPECT_THROW(manager.ProduceFile(output_path, "model"), std::runtime_error);
}

TEST(DataObjectIOSchemaTest, SaveRenamedKeyDoesNotRenameInMemoryObject)
{
    const command_test::ScopedTempDir temp_dir{ "schema_rename_semantics" };
    const auto database_path{ temp_dir.path() / "rename.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataObjectManager manager;
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    manager.SaveDataObject("model", "saved_model");

    EXPECT_TRUE(manager.HasDataObject("model"));
    EXPECT_FALSE(manager.HasDataObject("saved_model"));

    ASSERT_NO_THROW(manager.LoadDataObject("saved_model"));
    EXPECT_TRUE(manager.HasDataObject("saved_model"));
}

TEST(DataObjectIOSchemaTest, DatabaseManagerConstructorEnsuresSchemaAndHotPathDoesNotRecheck)
{
    const command_test::ScopedTempDir temp_dir{ "schema_hot_path" };
    const auto database_path{ temp_dir.path() / "hot_path.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataObjectManager manager;
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    EXPECT_TRUE(HasTable(database_path, "model_object"));

    manager.GetDatabaseManager()->GetDatabase()->Execute("DROP TABLE model_object;");
    EXPECT_FALSE(HasTable(database_path, "model_object"));
    EXPECT_THROW(manager.SaveDataObject("model"), std::runtime_error);
    EXPECT_FALSE(HasTable(database_path, "model_object"));
}

TEST(DataObjectIOSchemaTest, ModelDaoSaveDoesNotCreateSchemaImplicitly)
{
    const command_test::ScopedTempDir temp_dir{ "schema_model_dao_contract" };
    const auto database_path{ temp_dir.path() / "model_dao.sqlite" };
    const auto model{ LoadFixtureModel(command_test::TestDataPath("test_model.cif")) };

    rg::SQLiteWrapper database{ database_path };
    rg::ModelObjectDAO dao{ &database };

    EXPECT_THROW(dao.Save(model.get(), "model"), std::runtime_error);
    EXPECT_FALSE(HasTable(database_path, "model_object"));
}

TEST(DataObjectIOSchemaTest, MapDaoSaveDoesNotCreateSchemaImplicitly)
{
    const command_test::ScopedTempDir temp_dir{ "schema_map_dao_contract" };
    const auto database_path{ temp_dir.path() / "map_dao.sqlite" };
    const auto map_object{ MakeTinyMapObject() };

    rg::SQLiteWrapper database{ database_path };
    rg::MapObjectDAO dao{ &database };

    EXPECT_THROW(dao.Save(&map_object, "map"), std::runtime_error);
    EXPECT_FALSE(HasTable(database_path, "map_list"));
}

TEST(DataObjectIOSchemaTest, ProcessFileThrowsOnMalformedModelInput)
{
    const command_test::ScopedTempDir temp_dir{ "bad_model_file" };
    const auto malformed_path{ temp_dir.path() / "bad_model.cif" };
    {
        std::ofstream output{ malformed_path };
        output << "data_bad\nloop_\n_atom_site.id\n";
    }

    rg::DataObjectManager manager;
    try
    {
        manager.ProcessFile(malformed_path, "broken");
        FAIL() << "Expected malformed model input to throw.";
    }
    catch (const std::runtime_error & ex)
    {
        const std::string message{ ex.what() };
        EXPECT_NE(message.find(malformed_path.string()), std::string::npos);
        EXPECT_NE(message.find("broken"), std::string::npos);
    }
}

TEST(DataObjectIOSchemaTest, ProcessFileThrowsOnMalformedMapInput)
{
    const command_test::ScopedTempDir temp_dir{ "bad_map_file" };
    const auto malformed_path{ temp_dir.path() / "bad_map.map" };
    {
        std::ofstream output{ malformed_path, std::ios::binary };
        output << "bad";
    }

    rg::DataObjectManager manager;
    try
    {
        manager.ProcessFile(malformed_path, "broken_map");
        FAIL() << "Expected malformed map input to throw.";
    }
    catch (const std::runtime_error & ex)
    {
        const std::string message{ ex.what() };
        EXPECT_NE(message.find(malformed_path.string()), std::string::npos);
        EXPECT_NE(message.find("broken_map"), std::string::npos);
    }
}

TEST(DataObjectIOSchemaTest, ProduceFileThrowsWhenWriterCannotOpenTarget)
{
    const command_test::ScopedTempDir temp_dir{ "bad_output_target" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto output_path{ temp_dir.path() / "missing_dir" / "output.cif" };

    rg::DataObjectManager manager;
    manager.ProcessFile(model_path, "model");

    try
    {
        manager.ProduceFile(output_path, "model");
        FAIL() << "Expected writer open failure to throw.";
    }
    catch (const std::runtime_error & ex)
    {
        const std::string message{ ex.what() };
        EXPECT_NE(message.find(output_path.string()), std::string::npos);
        EXPECT_NE(message.find("model"), std::string::npos);
    }
}

TEST(DataObjectIOSchemaTest, ReaderGettersThrowIfReadDidNotSucceed)
{
    const command_test::ScopedTempDir temp_dir{ "reader_failure_contract" };
    const auto missing_model_path{ temp_dir.path() / "missing_model.cif" };
    const auto malformed_map_path{ temp_dir.path() / "bad_map.map" };
    {
        std::ofstream map_output{ malformed_map_path, std::ios::binary };
        map_output << "bad";
    }

    rg::ModelFileReader model_reader{ missing_model_path.string() };
    EXPECT_THROW(model_reader.Read(), std::runtime_error);
    EXPECT_THROW(model_reader.GetDataBlockPtr(), std::runtime_error);

    rg::MapFileReader map_reader{ malformed_map_path.string() };
    EXPECT_THROW(map_reader.Read(), std::runtime_error);
    EXPECT_THROW(map_reader.GetGridSizeArray(), std::runtime_error);
    EXPECT_THROW(map_reader.GetGridSpacingArray(), std::runtime_error);
    EXPECT_THROW(map_reader.GetOriginArray(), std::runtime_error);
    EXPECT_THROW((void)map_reader.GetMapValueArray(), std::runtime_error);
}
