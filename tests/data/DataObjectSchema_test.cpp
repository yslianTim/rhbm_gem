#include <gtest/gtest.h>

#include <array>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include <rhbm_gem/data/io/FileIO.hpp>
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "internal/sqlite/DatabaseManager.hpp"
#include "internal/sqlite/SQLiteWrapper.hpp"
#include "CommandTestHelpers.hpp"

namespace rg = rhbm_gem;

namespace {

constexpr const char* kUpsertObjectMetadataSql =
    "INSERT INTO object_metadata(key_tag, object_type) VALUES (?, ?) "
    "ON CONFLICT(key_tag) DO UPDATE SET object_type = excluded.object_type;";

void UpsertObjectMetadata(
    rg::SQLiteWrapper& database,
    const std::string& key_tag,
    const std::string& object_type) {
    database.Execute(
        "CREATE TABLE IF NOT EXISTS object_metadata (key_tag TEXT PRIMARY KEY, object_type TEXT);");
    database.Prepare(kUpsertObjectMetadataSql);
    rg::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.Bind<std::string>(2, object_type);
    database.StepOnce();
}

void CreateLegacyMapListWithoutForeignKeyTable(
    const std::filesystem::path& database_path,
    const rg::MapObject& map_object,
    const std::string& key_tag) {
    rg::SQLiteWrapper database{database_path};
    database.Execute("DROP TABLE IF EXISTS map_list;");
    database.Execute(
        "CREATE TABLE map_list ("
        "key_tag TEXT PRIMARY KEY, "
        "grid_size_x INTEGER, grid_size_y INTEGER, grid_size_z INTEGER, "
        "grid_spacing_x DOUBLE, grid_spacing_y DOUBLE, grid_spacing_z DOUBLE, "
        "origin_x DOUBLE, origin_y DOUBLE, origin_z DOUBLE, "
        "map_value_array BLOB"
        ");");

    const auto grid_size{map_object.GetGridSize()};
    const auto grid_spacing{map_object.GetGridSpacing()};
    const auto origin{map_object.GetOrigin()};
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

void CreateVersion2MetadataBasedMapShapeDatabase(const std::filesystem::path& database_path) {
    std::array<int, 3> grid_size{2, 2, 2};
    std::array<float, 3> grid_spacing{1.0f, 1.0f, 1.0f};
    std::array<float, 3> origin{0.0f, 0.0f, 0.0f};
    auto values{std::make_unique<float[]>(8)};
    for (size_t i = 0; i < 8; ++i) {
        values[i] = static_cast<float>(i);
    }
    rg::MapObject map_object{grid_size, grid_spacing, origin, std::move(values)};

    CreateLegacyMapListWithoutForeignKeyTable(database_path, map_object, "map_only");
    rg::SQLiteWrapper database{database_path};
    UpsertObjectMetadata(database, "map_only", "map");
    database.Execute("PRAGMA user_version = 2;");
}

int GetUserVersion(const std::filesystem::path& database_path) {
    rg::SQLiteWrapper database{database_path};
    database.Prepare("PRAGMA user_version;");
    rg::SQLiteWrapper::StatementGuard guard(database);
    const auto rc{database.StepNext()};
    if (rc != rg::SQLiteWrapper::StepRow()) {
        return 0;
    }
    return database.GetColumn<int>(0);
}

void SetUserVersion(const std::filesystem::path& database_path, int user_version) {
    rg::SQLiteWrapper database{database_path};
    database.Execute("PRAGMA user_version = " + std::to_string(user_version) + ";");
}

void ExecuteSql(const std::filesystem::path& database_path, const std::string& sql) {
    rg::SQLiteWrapper database{database_path};
    database.Execute(sql);
}

void ExecuteSqlWithForeignKeysOff(
    const std::filesystem::path& database_path,
    const std::string& sql) {
    rg::SQLiteWrapper database{database_path};
    database.Execute("PRAGMA foreign_keys = OFF;");
    database.Execute(sql);
}

bool HasTable(const std::filesystem::path& database_path, const std::string& table_name) {
    rg::SQLiteWrapper database{database_path};
    database.Prepare(
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;");
    rg::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, table_name);
    return database.StepNext() == rg::SQLiteWrapper::StepRow();
}

int CountRows(
    const std::filesystem::path& database_path,
    const std::string& table_name,
    const std::string& key_tag = "") {
    rg::SQLiteWrapper database{database_path};
    const auto sql{
        key_tag.empty()
            ? "SELECT COUNT(*) FROM " + table_name + ";"
            : "SELECT COUNT(*) FROM " + table_name + " WHERE key_tag = ?;"};
    database.Prepare(sql);
    rg::SQLiteWrapper::StatementGuard guard(database);
    if (!key_tag.empty()) {
        database.Bind<std::string>(1, key_tag);
    }
    if (database.StepNext() != rg::SQLiteWrapper::StepRow()) {
        throw std::runtime_error("Failed to count rows in " + table_name);
    }
    return database.GetColumn<int>(0);
}

bool HasForeignKey(
    const std::filesystem::path& database_path,
    const std::string& table_name,
    const std::string& from_column,
    const std::string& referenced_table,
    const std::string& referenced_column,
    const std::string& on_delete) {
    rg::SQLiteWrapper database{database_path};
    database.Prepare("PRAGMA foreign_key_list(" + table_name + ");");
    rg::SQLiteWrapper::StatementGuard guard(database);
    while (true) {
        const auto rc{database.StepNext()};
        if (rc == rg::SQLiteWrapper::StepDone()) {
            return false;
        }
        if (rc != rg::SQLiteWrapper::StepRow()) {
            throw std::runtime_error("Failed to inspect foreign_key_list for " + table_name);
        }
        if (database.GetColumn<std::string>(3) == from_column &&
            database.GetColumn<std::string>(2) == referenced_table &&
            database.GetColumn<std::string>(4) == referenced_column &&
            database.GetColumn<std::string>(6) == on_delete) {
            return true;
        }
    }
}

std::string SanitizeLegacyKey(const std::string& key_tag) {
    std::string sanitized_key_tag;
    sanitized_key_tag.reserve(key_tag.size());
    for (char ch : key_tag) {
        if (std::isalnum(static_cast<unsigned char>(ch)) || ch == '_') {
            sanitized_key_tag.push_back(ch);
        } else {
            sanitized_key_tag.push_back('_');
        }
    }
    return sanitized_key_tag;
}

std::filesystem::path LegacyFixturePath() {
    return command_test::TestDataPath("legacy/legacy_model_v1.sqlite");
}

void CopyLegacyFixtureDatabase(const std::filesystem::path& database_path) {
    std::filesystem::copy_file(
        LegacyFixturePath(),
        database_path,
        std::filesystem::copy_options::overwrite_existing);
}

void RenameLegacyModel(
    const std::filesystem::path& database_path,
    const std::string& old_key_tag,
    const std::string& new_key_tag,
    const std::string& new_pdb_id) {
    const auto old_sanitized{SanitizeLegacyKey(old_key_tag)};
    const auto new_sanitized{SanitizeLegacyKey(new_key_tag)};
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
    const std::filesystem::path& database_path,
    const std::string& source_key_tag,
    const std::string& new_key_tag,
    const std::string& new_pdb_id) {
    const auto source_sanitized{SanitizeLegacyKey(source_key_tag)};
    const auto new_sanitized{SanitizeLegacyKey(new_key_tag)};
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
        "SELECT '" +
            new_key_tag + "', atom_size, '" + new_pdb_id +
            "', emd_id, map_resolution, resolution_method "
            "FROM model_list WHERE key_tag = '" +
            source_key_tag + "';");
    ExecuteSql(
        database_path,
        "INSERT INTO object_metadata(key_tag, object_type) VALUES ('" + new_key_tag + "', 'model');");
}

std::shared_ptr<rg::ModelObject> LoadFixtureModel(
    const std::filesystem::path& model_path,
    const std::string& key_tag = "model") {
    rg::DataObjectManager manager{};
    manager.ProcessFile(model_path, key_tag);
    return manager.GetTypedDataObject<rg::ModelObject>(key_tag);
}

rg::MapObject MakeTinyMapObject(float scale = 1.0f) {
    std::array<int, 3> grid_size{2, 2, 2};
    std::array<float, 3> grid_spacing{1.0f, 1.0f, 1.0f};
    std::array<float, 3> origin{0.0f, 0.0f, 0.0f};
    auto values{std::make_unique<float[]>(8)};
    for (size_t i = 0; i < 8; ++i) {
        values[i] = static_cast<float>(i) * scale;
    }
    return rg::MapObject{grid_size, grid_spacing, origin, std::move(values)};
}

void SaveTinyMapThroughManager(
    rg::DataObjectManager& manager,
    const std::filesystem::path& map_path,
    const std::string& key_tag,
    float scale = 1.0f) {
    auto map_object{MakeTinyMapObject(scale)};
    rg::WriteMap(map_path, map_object);
    manager.ProcessFile(map_path, key_tag);
    manager.SaveDataObject(key_tag);
}

class FailingDataObject final : public rg::DataObjectBase {
  public:
    std::unique_ptr<rg::DataObjectBase> Clone() const override {
        return std::make_unique<FailingDataObject>(*this);
    }

    void Display() const override {}
    void Update() override {}
    void SetKeyTag(const std::string& label) override { m_key_tag = label; }
    std::string GetKeyTag() const override { return m_key_tag; }

  private:
    std::string m_key_tag;
};

} // namespace

TEST(DataObjectSchemaBootstrapTest, EmptyDatabaseBootstrapsNormalizedSchema) {
    const command_test::ScopedTempDir temp_dir{"data_schema_bootstrap"};
    const auto database_path{temp_dir.path() / "bootstrap.sqlite"};

    rg::DatabaseManager database_manager{database_path};

    EXPECT_EQ(GetUserVersion(database_path), 2);
    EXPECT_TRUE(HasTable(database_path, "object_catalog"));
    EXPECT_FALSE(HasTable(database_path, "object_metadata"));
    EXPECT_TRUE(HasTable(database_path, "model_object"));
    EXPECT_TRUE(HasTable(database_path, "map_list"));
    EXPECT_EQ(GetUserVersion(database_path), 2);
}

TEST(DataObjectSchemaBootstrapTest, SaveUnsupportedTypeRollsBackCatalogMutation) {
    const command_test::ScopedTempDir temp_dir{"data_schema_atomic_save"};
    const auto database_path{temp_dir.path() / "atomic.sqlite"};

    rg::DatabaseManager database_manager{database_path};
    FailingDataObject data_object;
    data_object.SetKeyTag("failing");

    EXPECT_THROW(database_manager.SaveDataObject(&data_object, "failing"), std::runtime_error);
    EXPECT_EQ(CountRows(database_path, "object_catalog"), 0);
}

TEST(DataObjectSchemaBootstrapTest, UnknownSchemaVersionThrows) {
    const command_test::ScopedTempDir temp_dir{"data_schema_unknown_version"};
    const auto database_path{temp_dir.path() / "unknown.sqlite"};

    SetUserVersion(database_path, 99);
    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaPersistenceTest, FinalV2CatalogDatabaseRemainsLoadable) {
    const command_test::ScopedTempDir temp_dir{"data_schema_final_v2_loadable"};
    const auto database_path{temp_dir.path() / "final_v2.sqlite"};
    const auto model_path{command_test::TestDataPath("test_model.cif")};
    const auto map_path{temp_dir.path() / "map_only.map"};

    {
        rg::DataObjectManager manager{};
        manager.SetDatabaseManager(database_path);
        manager.ProcessFile(model_path, "model");
        manager.SaveDataObject("model");
        SaveTinyMapThroughManager(manager, map_path, "map", 3.0f);
    }

    SetUserVersion(database_path, 2);

    rg::DataObjectManager manager{};
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    ASSERT_NO_THROW(manager.LoadDataObject("model"));
    ASSERT_NO_THROW(manager.LoadDataObject("map"));
    EXPECT_NE(manager.GetTypedDataObject<rg::ModelObject>("model"), nullptr);
    EXPECT_NE(manager.GetTypedDataObject<rg::MapObject>("map"), nullptr);
}

TEST(DataObjectSchemaPersistenceTest, DatabaseRoundTripPreservesChainMetadataAndSymmetryFiltering) {
    const command_test::ScopedTempDir temp_dir{"data_schema_roundtrip"};
    const auto database_path{temp_dir.path() / "roundtrip.sqlite"};
    const auto model_path{command_test::TestDataPath("test_model_keyvalue_entity.cif")};

    auto original_model{LoadFixtureModel(model_path)};
    const auto original_chain_map{original_model->GetChainIDListMap()};
    for (const auto& atom : original_model->GetAtomList()) {
        atom->SetSelectedFlag(true);
    }
    original_model->FilterAtomFromSymmetry(false);
    original_model->Update();
    const auto original_selected_count{original_model->GetNumberOfSelectedAtom()};

    rg::DataObjectManager manager{};
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    manager.SaveDataObject("model");
    manager.ClearDataObjects();
    manager.LoadDataObject("model");

    auto loaded_model{manager.GetTypedDataObject<rg::ModelObject>("model")};
    EXPECT_EQ(loaded_model->GetChainIDListMap(), original_chain_map);
    EXPECT_GT(CountRows(database_path, "model_chain_map", "model"), 0);

    for (const auto& atom : loaded_model->GetAtomList()) {
        atom->SetSelectedFlag(true);
    }
    loaded_model->FilterAtomFromSymmetry(false);
    loaded_model->Update();
    EXPECT_EQ(loaded_model->GetNumberOfSelectedAtom(), original_selected_count);
}

TEST(DataObjectSchemaPersistenceTest, DistinctUnsanitizedKeysDoNotCollideInV2Schema) {
    const command_test::ScopedTempDir temp_dir{"data_schema_key_collision"};
    const auto database_path{temp_dir.path() / "collision.sqlite"};
    const auto model_path{command_test::TestDataPath("test_model.cif")};

    rg::DataObjectManager manager{};
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

TEST(DataObjectSchemaPersistenceTest, SaveRenamedKeyDoesNotRenameInMemoryObject) {
    const command_test::ScopedTempDir temp_dir{"data_schema_rename_semantics"};
    const auto database_path{temp_dir.path() / "rename.sqlite"};
    const auto model_path{command_test::TestDataPath("test_model.cif")};

    rg::DataObjectManager manager{};
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    manager.SaveDataObject("model", "saved_model");

    EXPECT_TRUE(manager.HasDataObject("model"));
    EXPECT_FALSE(manager.HasDataObject("saved_model"));

    ASSERT_NO_THROW(manager.LoadDataObject("saved_model"));
    EXPECT_TRUE(manager.HasDataObject("saved_model"));
}

TEST(DataObjectSchemaCompatibilityTest, Version2WithObjectMetadataFailsFast) {
    const command_test::ScopedTempDir temp_dir{"data_schema_v2_metadata_row"};
    const auto database_path{temp_dir.path() / "metadata_row.sqlite"};
    const auto map_object{MakeTinyMapObject()};

    {
        rg::DatabaseManager database_manager{database_path};
        database_manager.SaveDataObject(&map_object, "map_only");
    }

    ExecuteSql(
        database_path,
        "CREATE TABLE IF NOT EXISTS object_metadata (key_tag TEXT PRIMARY KEY, object_type TEXT);");
    ExecuteSql(
        database_path,
        "INSERT INTO object_metadata(key_tag, object_type) VALUES ('map_only', 'map');");
    SetUserVersion(database_path, 2);

    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaCompatibilityTest, Version2MetadataBasedShapeFailsFast) {
    const command_test::ScopedTempDir temp_dir{"data_schema_v2_metadata_shape"};
    const auto database_path{temp_dir.path() / "metadata_shape.sqlite"};

    CreateVersion2MetadataBasedMapShapeDatabase(database_path);
    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaCompatibilityTest, ManagedButUnversionedDatabaseFailsFast) {
    const command_test::ScopedTempDir temp_dir{"data_schema_unversioned_nonlegacy"};
    const auto database_path{temp_dir.path() / "managed_unversioned.sqlite"};
    const auto model_path{command_test::TestDataPath("test_model.cif")};

    {
        rg::DataObjectManager manager{};
        manager.SetDatabaseManager(database_path);
        manager.ProcessFile(model_path, "model");
        manager.SaveDataObject("model");
    }

    SetUserVersion(database_path, 0);
    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaCompatibilityTest, MixedUnknownSchemaFailsFast) {
    const command_test::ScopedTempDir temp_dir{"data_schema_mixed_unknown"};
    const auto database_path{temp_dir.path() / "mixed_unknown.sqlite"};

    ExecuteSql(database_path, "CREATE TABLE foreign_table (id INTEGER PRIMARY KEY);");
    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, NormalizedV2DatabaseMissingRequiredTableThrows) {
    const command_test::ScopedTempDir temp_dir{"data_schema_missing_v2_table"};
    const auto database_path{temp_dir.path() / "missing_v2.sqlite"};

    {
        rg::DatabaseManager database_manager{database_path};
        EXPECT_EQ(GetUserVersion(database_path), 2);
    }
    ExecuteSql(database_path, "DROP TABLE model_bond_group_potential;");
    SetUserVersion(database_path, 2);

    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, FinalV2SchemaValidationRejectsMissingForeignKeys) {
    const command_test::ScopedTempDir temp_dir{"data_schema_missing_fk"};
    const auto database_path{temp_dir.path() / "missing_fk.sqlite"};

    {
        rg::DatabaseManager database_manager{database_path};
        EXPECT_EQ(GetUserVersion(database_path), 2);
    }

    ExecuteSql(database_path, "DROP TABLE map_list;");
    ExecuteSql(
        database_path,
        "CREATE TABLE map_list ("
        "key_tag TEXT PRIMARY KEY, "
        "grid_size_x INTEGER, grid_size_y INTEGER, grid_size_z INTEGER, "
        "grid_spacing_x DOUBLE, grid_spacing_y DOUBLE, grid_spacing_z DOUBLE, "
        "origin_x DOUBLE, origin_y DOUBLE, origin_z DOUBLE, "
        "map_value_array BLOB"
        ");");
    SetUserVersion(database_path, 2);

    EXPECT_FALSE(
        HasForeignKey(
            database_path,
            "map_list",
            "key_tag",
            "object_catalog",
            "key_tag",
            "CASCADE"));
    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, FinalV2SchemaValidationRejectsMissingRequiredCatalogColumns) {
    const command_test::ScopedTempDir temp_dir{"data_schema_bad_catalog_columns"};
    const auto database_path{temp_dir.path() / "bad_catalog_columns.sqlite"};

    {
        rg::DatabaseManager database_manager{database_path};
        EXPECT_EQ(GetUserVersion(database_path), 2);
    }

    ExecuteSqlWithForeignKeysOff(database_path, "DROP TABLE object_catalog;");
    ExecuteSqlWithForeignKeysOff(
        database_path,
        "CREATE TABLE object_catalog (key_tag TEXT PRIMARY KEY);");
    SetUserVersion(database_path, 2);

    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, FinalV2SchemaValidationRejectsUnknownObjectTypeValue) {
    const command_test::ScopedTempDir temp_dir{"data_schema_bad_catalog_type"};
    const auto database_path{temp_dir.path() / "bad_catalog_type.sqlite"};

    {
        rg::DatabaseManager database_manager{database_path};
        EXPECT_EQ(GetUserVersion(database_path), 2);
    }

    ExecuteSqlWithForeignKeysOff(database_path, "DROP TABLE object_catalog;");
    ExecuteSqlWithForeignKeysOff(
        database_path,
        "CREATE TABLE object_catalog (key_tag TEXT PRIMARY KEY, object_type TEXT NOT NULL);");
    ExecuteSqlWithForeignKeysOff(
        database_path,
        "INSERT INTO object_catalog(key_tag, object_type) VALUES ('unknown_root', 'unknown');");
    SetUserVersion(database_path, 2);

    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, ForeignKeyRejectsOrphanModelChildRows) {
    const command_test::ScopedTempDir temp_dir{"data_schema_fk_orphan"};
    const auto database_path{temp_dir.path() / "orphan.sqlite"};
    rg::DatabaseManager database_manager{database_path};

    EXPECT_THROW(
        ExecuteSql(
            database_path,
            "INSERT INTO model_chain_map(key_tag, entity_id, chain_ordinal, chain_id) "
            "VALUES ('missing', '1', 0, 'A');"),
        std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, DeletingCatalogRootCascadesPayloadRows) {
    const command_test::ScopedTempDir temp_dir{"data_schema_catalog_cascade"};
    const auto database_path{temp_dir.path() / "cascade.sqlite"};
    const auto model_path{command_test::TestDataPath("test_model.cif")};

    rg::DataObjectManager manager{};
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    manager.SaveDataObject("model");

    ExecuteSql(database_path, "DELETE FROM object_catalog WHERE key_tag = 'model';");

    EXPECT_EQ(CountRows(database_path, "object_catalog"), 0);
    EXPECT_EQ(CountRows(database_path, "model_object"), 0);
    EXPECT_EQ(CountRows(database_path, "model_atom"), 0);
}

#ifdef RHBM_GEM_LEGACY_V1_SUPPORT

TEST(DataObjectSchemaMigrationTest, LegacyFixtureMigratesToNormalizedSchema) {
    const command_test::ScopedTempDir temp_dir{"data_schema_migrate_legacy"};
    const auto database_path{temp_dir.path() / "legacy.sqlite"};

    CopyLegacyFixtureDatabase(database_path);
    ASSERT_EQ(GetUserVersion(database_path), 0);
    ASSERT_TRUE(HasTable(database_path, "model_list"));

    rg::DataObjectManager manager{};
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    rg::DatabaseManager verifier{database_path};
    EXPECT_EQ(GetUserVersion(database_path), 2);

    ASSERT_NO_THROW(manager.LoadDataObject("legacy_model"));
    auto model{manager.GetTypedDataObject<rg::ModelObject>("legacy_model")};
    EXPECT_EQ(model->GetPdbID(), "LEGACY");
    EXPECT_TRUE(HasTable(database_path, "model_object"));
    EXPECT_FALSE(HasTable(database_path, "model_list"));
    EXPECT_EQ(GetUserVersion(database_path), 2);
}

TEST(DataObjectSchemaMigrationTest, FailedLegacyMigrationRollsBackToLegacyState) {
    const command_test::ScopedTempDir temp_dir{"data_schema_migrate_rollback"};
    const auto database_path{temp_dir.path() / "legacy_rollback.sqlite"};

    CopyLegacyFixtureDatabase(database_path);
    RenameLegacyModel(database_path, "legacy_model", "key_a", "MODEL_A");
    DuplicateLegacyModel(database_path, "key_a", "key_b", "MODEL_B");
    ExecuteSql(
        database_path,
        "UPDATE model_list SET atom_size = 999 WHERE key_tag = 'key_b';");

    rg::DataObjectManager manager{};
    EXPECT_THROW(manager.SetDatabaseManager(database_path), std::runtime_error);
    EXPECT_TRUE(HasTable(database_path, "model_list"));
    EXPECT_FALSE(HasTable(database_path, "model_object"));
    EXPECT_EQ(GetUserVersion(database_path), 0);
}

TEST(DataObjectSchemaMigrationTest, LegacyMigrationUsesModelListWhenMetadataIncomplete) {
    const command_test::ScopedTempDir temp_dir{"data_schema_partial_metadata"};
    const auto database_path{temp_dir.path() / "partial_metadata.sqlite"};

    CopyLegacyFixtureDatabase(database_path);
    RenameLegacyModel(database_path, "legacy_model", "key_a", "MODEL_A");
    DuplicateLegacyModel(database_path, "key_a", "key_b", "MODEL_B");
    ExecuteSql(
        database_path,
        "DELETE FROM object_metadata WHERE key_tag = 'key_b' AND object_type = 'model';");

    rg::DataObjectManager manager{};
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    ASSERT_NO_THROW(manager.LoadDataObject("key_a"));
    ASSERT_NO_THROW(manager.LoadDataObject("key_b"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("key_a")->GetPdbID(), "MODEL_A");
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("key_b")->GetPdbID(), "MODEL_B");
    EXPECT_EQ(CountRows(database_path, "model_object"), 2);
    EXPECT_EQ(CountRows(database_path, "object_catalog"), 2);
    EXPECT_FALSE(HasTable(database_path, "object_metadata"));
}

TEST(DataObjectSchemaMigrationTest, LegacyMigrationIgnoresMetadataOnlyGhostKeys) {
    const command_test::ScopedTempDir temp_dir{"data_schema_ghost_metadata"};
    const auto database_path{temp_dir.path() / "ghost_metadata.sqlite"};

    CopyLegacyFixtureDatabase(database_path);
    RenameLegacyModel(database_path, "legacy_model", "real_model", "REAL");
    {
        rg::SQLiteWrapper database{database_path};
        UpsertObjectMetadata(database, "ghost_model", "model");
    }

    rg::DataObjectManager manager{};
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    ASSERT_NO_THROW(manager.LoadDataObject("real_model"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("real_model")->GetPdbID(), "REAL");
    EXPECT_THROW(manager.LoadDataObject("ghost_model"), std::runtime_error);
    EXPECT_EQ(CountRows(database_path, "model_object"), 1);
    EXPECT_EQ(CountRows(database_path, "object_catalog"), 1);
    EXPECT_FALSE(HasTable(database_path, "object_metadata"));
}

TEST(DataObjectSchemaMigrationTest, LegacyMigrationDropsOnlyOwnedTables) {
    const command_test::ScopedTempDir temp_dir{"data_schema_owned_table_drop"};
    const auto database_path{temp_dir.path() / "owned_drop.sqlite"};

    CopyLegacyFixtureDatabase(database_path);
    ExecuteSql(database_path, "CREATE TABLE custom_atom_list_in_notes (value INTEGER);");

    rg::DataObjectManager manager{};
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    EXPECT_TRUE(HasTable(database_path, "custom_atom_list_in_notes"));
    EXPECT_FALSE(HasTable(database_path, "atom_list_in_legacy_model"));
}

TEST(DataObjectSchemaMigrationTest, LegacyV1MapListWithoutFkMigratesToFinalV2) {
    const command_test::ScopedTempDir temp_dir{"data_schema_legacy_map_fk_rebuild"};
    const auto database_path{temp_dir.path() / "legacy_map_fk.sqlite"};
    const auto map_object{MakeTinyMapObject(4.0f)};

    CopyLegacyFixtureDatabase(database_path);
    CreateLegacyMapListWithoutForeignKeyTable(database_path, map_object, "legacy_map");
    SetUserVersion(database_path, 1);

    EXPECT_FALSE(
        HasForeignKey(
            database_path,
            "map_list",
            "key_tag",
            "object_catalog",
            "key_tag",
            "CASCADE"));

    rg::DataObjectManager manager{};
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    EXPECT_EQ(GetUserVersion(database_path), 2);
    EXPECT_TRUE(
        HasForeignKey(
            database_path,
            "map_list",
            "key_tag",
            "object_catalog",
            "key_tag",
            "CASCADE"));
    ASSERT_NO_THROW(manager.LoadDataObject("legacy_map"));
    EXPECT_EQ(
        manager.GetTypedDataObject<rg::MapObject>("legacy_map")->GetGridSize(),
        map_object.GetGridSize());
    EXPECT_EQ(GetUserVersion(database_path), 2);
}

#endif
