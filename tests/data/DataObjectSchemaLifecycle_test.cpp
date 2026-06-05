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

    EXPECT_EQ(data_test::GetUserVersion(database_path), 4);
    EXPECT_TRUE(data_test::HasTable(database_path, "object_catalog"));
    EXPECT_FALSE(data_test::HasTable(database_path, "object_metadata"));
    EXPECT_TRUE(data_test::HasTable(database_path, "model_object"));
    EXPECT_TRUE(data_test::HasTable(database_path, "map_list"));
}

TEST(DataObjectSchemaLifecycleTest, EmptyDatabaseBootstrapsGaussianInterceptColumns)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_gaussian_intercept_columns" };
    const auto database_path{ temp_dir.path() / "gaussian_intercept_columns.sqlite" };

    rg::SQLitePersistence database_manager{ database_path };

    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_local_potential", "intercept_estimate_ols"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_local_potential", "intercept_estimate_mdpde"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_bond_local_potential", "intercept_estimate_ols"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_bond_local_potential", "intercept_estimate_mdpde"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_posterior", "intercept_estimate_posterior"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_posterior", "intercept_variance_posterior"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_bond_posterior", "intercept_estimate_posterior"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_bond_posterior", "intercept_variance_posterior"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_group_potential", "intercept_estimate_mean"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_group_potential", "intercept_estimate_mdpde"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_group_potential", "intercept_estimate_prior"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_group_potential", "intercept_variance_prior"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_bond_group_potential", "intercept_estimate_mean"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_bond_group_potential", "intercept_estimate_mdpde"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_bond_group_potential", "intercept_estimate_prior"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_bond_group_potential", "intercept_variance_prior"));
}

TEST(DataObjectSchemaLifecycleTest, VersionThreeSchemaMigratesGaussianInterceptColumns)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_v3_gaussian_intercept_migration" };
    const auto database_path{ temp_dir.path() / "v3_gaussian_intercept.sqlite" };

    {
        rg::SQLitePersistence database_manager{ database_path };
    }
    data_test::ExecuteSqlWithForeignKeysOff(database_path, "DROP TABLE model_atom_local_potential;");
    data_test::ExecuteSqlWithForeignKeysOff(
        database_path,
        "CREATE TABLE model_atom_local_potential ("
        "key_tag TEXT, "
        "serial_id INTEGER, "
        "sampling_size INTEGER, "
        "distance_and_map_value_list BLOB, "
        "amplitude_estimate_ols DOUBLE, "
        "width_estimate_ols DOUBLE, "
        "amplitude_estimate_mdpde DOUBLE, "
        "width_estimate_mdpde DOUBLE, "
        "alpha_r DOUBLE, "
        "PRIMARY KEY (key_tag, serial_id), "
        "FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE"
        ");");
    data_test::SetUserVersion(database_path, 3);

    rg::SQLitePersistence database_manager{ database_path };

    EXPECT_EQ(data_test::GetUserVersion(database_path), 4);
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_local_potential", "intercept_estimate_ols"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_local_potential", "intercept_estimate_mdpde"));
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
