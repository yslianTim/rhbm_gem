#include <gtest/gtest.h>

#include <filesystem>
#include <stdexcept>

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include "io/sqlite/SQLitePersistence.hpp"
#include "support/CommandTestHelpers.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

namespace {

template <typename SetupFn>
void ExpectSchemaOpenFails(
    const char * temp_dir_name,
    const char * database_name,
    SetupFn setup_database)
{
    const command_test::ScopedTempDir temp_dir{ temp_dir_name };
    const auto database_path{ temp_dir.path() / database_name };
    setup_database(database_path);
    EXPECT_THROW((void)rg::SQLitePersistence(database_path), std::runtime_error);
}

} // namespace

TEST(DataObjectSchemaLifecycleTest, EmptyDatabaseBootstrapsNormalizedSchema)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_bootstrap" };
    const auto database_path{ temp_dir.path() / "bootstrap.sqlite" };

    rg::SQLitePersistence database_manager{ database_path };

    EXPECT_EQ(data_test::GetUserVersion(database_path), 2);
    EXPECT_TRUE(data_test::HasTable(database_path, "object_catalog"));
    EXPECT_FALSE(data_test::HasTable(database_path, "object_metadata"));
    EXPECT_TRUE(data_test::HasTable(database_path, "model_object"));
    EXPECT_TRUE(data_test::HasTable(database_path, "map_list"));
}

TEST(DataObjectSchemaLifecycleTest, UnknownSchemaVersionThrows)
{
    ExpectSchemaOpenFails(
        "data_schema_unknown_version",
        "unknown.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::SetUserVersion(database_path, 99);
        });
}

TEST(DataObjectSchemaLifecycleTest, VersionOneSchemaFailsFast)
{
    ExpectSchemaOpenFails(
        "data_schema_version_one",
        "version_one.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::SetUserVersion(database_path, 1);
        });
}

TEST(DataObjectSchemaLifecycleTest, Version2WithObjectMetadataFailsFast)
{
    ExpectSchemaOpenFails(
        "data_schema_v2_metadata_row",
        "metadata_row.sqlite",
        [](const std::filesystem::path & database_path)
        {
            const auto map_object{ data_test::MakeTinyMapObject() };
            {
                rg::SQLitePersistence database_manager{ database_path };
                database_manager.SaveMap(map_object, "map_only");
            }

            data_test::ExecuteSql(
                database_path,
                "CREATE TABLE IF NOT EXISTS object_metadata (key_tag TEXT PRIMARY KEY, object_type TEXT);");
            data_test::ExecuteSql(
                database_path,
                "INSERT INTO object_metadata(key_tag, object_type) VALUES ('map_only', 'map');");
            data_test::SetUserVersion(database_path, 2);
        });
}

TEST(DataObjectSchemaLifecycleTest, Version2MetadataBasedShapeFailsFast)
{
    ExpectSchemaOpenFails(
        "data_schema_v2_metadata_shape",
        "metadata_shape.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::CreateVersion2MetadataBasedMapShapeDatabase(database_path);
        });
}

TEST(DataObjectSchemaLifecycleTest, ManagedButUnversionedDatabaseFailsFast)
{
    ExpectSchemaOpenFails(
        "data_schema_unversioned_nonlegacy",
        "managed_unversioned.sqlite",
        [](const std::filesystem::path & database_path)
        {
            const auto model_path{ command_test::TestDataPath("test_model.cif") };

            {
                rg::DataRepository repository{ database_path };
                auto model{ rg::ReadModel(model_path) };
                model->SetKeyTag("model");
                repository.SaveModel(*model, "model");
            }

            data_test::SetUserVersion(database_path, 0);
        });
}

TEST(DataObjectSchemaLifecycleTest, MixedUnknownSchemaFailsFast)
{
    ExpectSchemaOpenFails(
        "data_schema_mixed_unknown",
        "mixed_unknown.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::ExecuteSql(database_path, "CREATE TABLE foreign_table (id INTEGER PRIMARY KEY);");
        });
}
