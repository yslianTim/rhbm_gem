#include <gtest/gtest.h>

#include "data/schema/DataObjectIOSchemaTestSupport.hpp"

using namespace data_object_io_schema_test;

TEST(DataObjectSchemaMigrationTest, LegacyFixtureMigratesToNormalizedSchema)
{
    const command_test::ScopedTempDir temp_dir{ "schema_migrate_legacy" };
    const auto database_path{ temp_dir.path() / "legacy.sqlite" };

    CopyLegacyFixtureDatabase(database_path);
    ASSERT_EQ(GetUserVersion(database_path), 0);
    ASSERT_TRUE(HasTable(database_path, "model_list"));

    rg::DataObjectManager manager{ command_test::BuildDataIoServices() };
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    const auto dao_registry{ BuildDefaultDaoFactoryRegistry() };
    rg::DatabaseManager verifier{ database_path, dao_registry };
    EXPECT_EQ(verifier.GetSchemaVersion(), rg::DatabaseSchemaVersion::NormalizedV2);

    ASSERT_NO_THROW(manager.LoadDataObject("legacy_model"));
    auto model{ manager.GetTypedDataObject<rg::ModelObject>("legacy_model") };
    EXPECT_EQ(model->GetPdbID(), "LEGACY");
    EXPECT_TRUE(HasTable(database_path, "model_object"));
    EXPECT_FALSE(HasTable(database_path, "model_list"));
    EXPECT_EQ(GetUserVersion(database_path), 2);
}

TEST(DataObjectSchemaMigrationTest, FailedLegacyMigrationRollsBackToLegacyState)
{
    const command_test::ScopedTempDir temp_dir{ "schema_migrate_rollback" };
    const auto database_path{ temp_dir.path() / "legacy_rollback.sqlite" };

    CopyLegacyFixtureDatabase(database_path);
    RenameLegacyModel(database_path, "legacy_model", "key_a", "MODEL_A");
    DuplicateLegacyModel(database_path, "key_a", "key_b", "MODEL_B");
    ExecuteSql(
        database_path,
        "UPDATE model_list SET atom_size = 999 WHERE key_tag = 'key_b';");

    rg::DataObjectManager manager{ command_test::BuildDataIoServices() };
    EXPECT_THROW(manager.SetDatabaseManager(database_path), std::runtime_error);
    EXPECT_TRUE(HasTable(database_path, "model_list"));
    EXPECT_FALSE(HasTable(database_path, "model_object"));
    EXPECT_EQ(GetUserVersion(database_path), 0);
}

TEST(DataObjectSchemaMigrationTest, LegacyMigrationUsesModelListWhenMetadataIncomplete)
{
    const command_test::ScopedTempDir temp_dir{ "schema_partial_metadata" };
    const auto database_path{ temp_dir.path() / "partial_metadata.sqlite" };

    CopyLegacyFixtureDatabase(database_path);
    RenameLegacyModel(database_path, "legacy_model", "key_a", "MODEL_A");
    DuplicateLegacyModel(database_path, "key_a", "key_b", "MODEL_B");
    ExecuteSql(
        database_path,
        "DELETE FROM object_metadata WHERE key_tag = 'key_b' AND object_type = 'model';");

    rg::DataObjectManager manager{ command_test::BuildDataIoServices() };
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    ASSERT_NO_THROW(manager.LoadDataObject("key_a"));
    ASSERT_NO_THROW(manager.LoadDataObject("key_b"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("key_a")->GetPdbID(), "MODEL_A");
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("key_b")->GetPdbID(), "MODEL_B");
    EXPECT_EQ(CountRows(database_path, "model_object"), 2);
    EXPECT_EQ(CountRows(database_path, "object_catalog"), 2);
    EXPECT_FALSE(HasTable(database_path, "object_metadata"));
}

TEST(DataObjectSchemaMigrationTest, LegacyMigrationIgnoresMetadataOnlyGhostKeys)
{
    const command_test::ScopedTempDir temp_dir{ "schema_ghost_metadata" };
    const auto database_path{ temp_dir.path() / "ghost_metadata.sqlite" };

    CopyLegacyFixtureDatabase(database_path);
    RenameLegacyModel(database_path, "legacy_model", "real_model", "REAL");
    {
        rg::SQLiteWrapper database{ database_path };
        UpsertObjectMetadata(database, "ghost_model", "model");
    }

    rg::DataObjectManager manager{ command_test::BuildDataIoServices() };
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    ASSERT_NO_THROW(manager.LoadDataObject("real_model"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("real_model")->GetPdbID(), "REAL");
    EXPECT_THROW(manager.LoadDataObject("ghost_model"), std::runtime_error);
    EXPECT_EQ(CountRows(database_path, "model_object"), 1);
    EXPECT_EQ(CountRows(database_path, "object_catalog"), 1);
    EXPECT_FALSE(HasTable(database_path, "object_metadata"));
}

TEST(DataObjectSchemaMigrationTest, LegacyMigrationDropsOnlyOwnedTables)
{
    const command_test::ScopedTempDir temp_dir{ "schema_owned_table_drop" };
    const auto database_path{ temp_dir.path() / "owned_drop.sqlite" };

    CopyLegacyFixtureDatabase(database_path);
    ExecuteSql(database_path, "CREATE TABLE custom_atom_list_in_notes (value INTEGER);");

    rg::DataObjectManager manager{ command_test::BuildDataIoServices() };
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    EXPECT_TRUE(HasTable(database_path, "custom_atom_list_in_notes"));
    EXPECT_FALSE(HasTable(database_path, "atom_list_in_legacy_model"));
}

TEST(DataObjectSchemaMigrationTest, LegacyV1MapListWithoutFkMigratesToFinalV2)
{
    const command_test::ScopedTempDir temp_dir{ "schema_legacy_map_fk_rebuild" };
    const auto database_path{ temp_dir.path() / "legacy_map_fk.sqlite" };
    const auto map_object{ MakeTinyMapObject(4.0f) };

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

    rg::DataObjectManager manager{ command_test::BuildDataIoServices() };
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
    EXPECT_EQ(manager.GetTypedDataObject<rg::MapObject>("legacy_map")->GetGridSize(), map_object.GetGridSize());
    EXPECT_EQ(GetUserVersion(database_path), 2);
}
