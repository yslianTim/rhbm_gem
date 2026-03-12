#include <gtest/gtest.h>

#include "data/schema/DataObjectIOSchemaTestSupport.hpp"

using namespace data_object_io_schema_test;

TEST(DataObjectSchemaValidationTest, NormalizedV2DatabaseMissingRequiredTableThrows)
{
    const command_test::ScopedTempDir temp_dir{ "schema_missing_v2_table" };
    const auto database_path{ temp_dir.path() / "missing_v2.sqlite" };

    {
        rg::DatabaseManager database_manager{ database_path };
        EXPECT_EQ(database_manager.GetSchemaVersion(), rg::DatabaseSchemaVersion::NormalizedV2);
    }
    ExecuteSql(database_path, "DROP TABLE model_bond_group_potential;");
    SetUserVersion(database_path, 2);

    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, FinalV2SchemaValidationRejectsMissingForeignKeys)
{
    const command_test::ScopedTempDir temp_dir{ "schema_missing_fk" };
    const auto database_path{ temp_dir.path() / "missing_fk.sqlite" };

    {
        rg::DatabaseManager database_manager{ database_path };
        EXPECT_EQ(database_manager.GetSchemaVersion(), rg::DatabaseSchemaVersion::NormalizedV2);
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

TEST(DataObjectSchemaValidationTest, FinalV2SchemaValidationRejectsMissingRequiredCatalogColumns)
{
    const command_test::ScopedTempDir temp_dir{ "schema_bad_catalog_columns" };
    const auto database_path{ temp_dir.path() / "bad_catalog_columns.sqlite" };

    {
        rg::DatabaseManager database_manager{ database_path };
        EXPECT_EQ(database_manager.GetSchemaVersion(), rg::DatabaseSchemaVersion::NormalizedV2);
    }

    ExecuteSqlWithForeignKeysOff(database_path, "DROP TABLE object_catalog;");
    ExecuteSqlWithForeignKeysOff(
        database_path,
        "CREATE TABLE object_catalog (key_tag TEXT PRIMARY KEY);");
    SetUserVersion(database_path, 2);

    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, FinalV2SchemaValidationRejectsUnknownObjectTypeValue)
{
    const command_test::ScopedTempDir temp_dir{ "schema_bad_catalog_type" };
    const auto database_path{ temp_dir.path() / "bad_catalog_type.sqlite" };

    {
        rg::DatabaseManager database_manager{ database_path };
        EXPECT_EQ(database_manager.GetSchemaVersion(), rg::DatabaseSchemaVersion::NormalizedV2);
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

TEST(DataObjectSchemaValidationTest, DistinctUnsanitizedKeysDoNotCollideInV2Schema)
{
    const command_test::ScopedTempDir temp_dir{ "schema_key_collision" };
    const auto database_path{ temp_dir.path() / "collision.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

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

TEST(DataObjectSchemaValidationTest, MixedUnknownSchemaFailsFast)
{
    const command_test::ScopedTempDir temp_dir{ "schema_mixed_unknown" };
    const auto database_path{ temp_dir.path() / "mixed_unknown.sqlite" };

    ExecuteSql(database_path, "CREATE TABLE foreign_table (id INTEGER PRIMARY KEY);");
    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, ForeignKeyRejectsOrphanModelChildRows)
{
    const command_test::ScopedTempDir temp_dir{ "schema_fk_orphan" };
    const auto database_path{ temp_dir.path() / "orphan.sqlite" };
    rg::DatabaseManager database_manager{ database_path };

    EXPECT_THROW(
        ExecuteSql(
            database_path,
            "INSERT INTO model_chain_map(key_tag, entity_id, chain_ordinal, chain_id) "
            "VALUES ('missing', '1', 0, 'A');"),
        std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, DeletingCatalogRootCascadesPayloadRows)
{
    const command_test::ScopedTempDir temp_dir{ "schema_catalog_cascade" };
    const auto database_path{ temp_dir.path() / "cascade.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataObjectManager manager{};
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    manager.SaveDataObject("model");

    ExecuteSql(database_path, "DELETE FROM object_catalog WHERE key_tag = 'model';");

    EXPECT_EQ(CountRows(database_path, "object_catalog"), 0);
    EXPECT_EQ(CountRows(database_path, "model_object"), 0);
    EXPECT_EQ(CountRows(database_path, "model_atom"), 0);
}

TEST(DataObjectSchemaValidationTest, SaveRenamedKeyDoesNotRenameInMemoryObject)
{
    const command_test::ScopedTempDir temp_dir{ "schema_rename_semantics" };
    const auto database_path{ temp_dir.path() / "rename.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataObjectManager manager{};
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    manager.SaveDataObject("model", "saved_model");

    EXPECT_TRUE(manager.HasDataObject("model"));
    EXPECT_FALSE(manager.HasDataObject("saved_model"));

    ASSERT_NO_THROW(manager.LoadDataObject("saved_model"));
    EXPECT_TRUE(manager.HasDataObject("saved_model"));
}
