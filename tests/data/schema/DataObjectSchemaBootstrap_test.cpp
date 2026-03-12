#include <gtest/gtest.h>

#include "data/schema/DataObjectIOSchemaTestSupport.hpp"

using namespace data_object_io_schema_test;

TEST(DataObjectSchemaBootstrapTest, EmptyDatabaseBootstrapsNormalizedSchema)
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

TEST(DataObjectSchemaBootstrapTest, SaveUnsupportedTypeRollsBackCatalogMutation)
{
    const command_test::ScopedTempDir temp_dir{ "schema_atomic_save" };
    const auto database_path{ temp_dir.path() / "atomic.sqlite" };

    rg::DatabaseManager database_manager{ database_path };
    FailingDataObject data_object;
    data_object.SetKeyTag("failing");

    EXPECT_THROW(database_manager.SaveDataObject(&data_object, "failing"), std::runtime_error);
    EXPECT_EQ(CountRows(database_path, "object_catalog"), 0);
}

TEST(DataObjectSchemaBootstrapTest, UnknownSchemaVersionThrows)
{
    const command_test::ScopedTempDir temp_dir{ "schema_unknown_version" };
    const auto database_path{ temp_dir.path() / "unknown.sqlite" };

    SetUserVersion(database_path, 99);
    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaBootstrapTest, Version2WithObjectMetadataFailsFast)
{
    const command_test::ScopedTempDir temp_dir{ "schema_v2_metadata_row" };
    const auto database_path{ temp_dir.path() / "metadata_row.sqlite" };
    const auto map_object{ MakeTinyMapObject() };

    {
        rg::DatabaseManager database_manager{ database_path };
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

TEST(DataObjectSchemaBootstrapTest, Version2MetadataBasedShapeFailsFast)
{
    const command_test::ScopedTempDir temp_dir{ "schema_v2_metadata_shape" };
    const auto database_path{ temp_dir.path() / "metadata_shape.sqlite" };

    CreateVersion2MetadataBasedMapShapeDatabase(database_path);
    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaBootstrapTest, ManagedButUnversionedDatabaseFailsFast)
{
    const command_test::ScopedTempDir temp_dir{ "schema_unversioned_nonlegacy" };
    const auto database_path{ temp_dir.path() / "managed_unversioned.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    {
        rg::DataObjectManager manager{};
        manager.SetDatabaseManager(database_path);
        manager.ProcessFile(model_path, "model");
        manager.SaveDataObject("model");
    }

    SetUserVersion(database_path, 0);
    EXPECT_THROW((void)rg::DatabaseManager(database_path), std::runtime_error);
}

TEST(DataObjectSchemaBootstrapTest, FinalV2CatalogDatabaseRemainsLoadable)
{
    const command_test::ScopedTempDir temp_dir{ "schema_final_v2_loadable" };
    const auto database_path{ temp_dir.path() / "final_v2.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto map_path{ temp_dir.path() / "map_only.map" };

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
