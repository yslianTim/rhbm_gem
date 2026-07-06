#include <gtest/gtest.h>

#include <filesystem>
#include <stdexcept>

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include "io/sqlite/SQLitePersistence.hpp"
#include "io/sqlite/SQLiteWrapper.hpp"
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

int QuerySingleInt(const std::filesystem::path & database_path, const std::string & sql)
{
    rg::SQLiteWrapper database{ database_path };
    database.Prepare(sql);
    rg::SQLiteWrapper::StatementGuard guard(database);
    if (database.StepNext() != rg::SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Expected one integer row.");
    }
    return database.GetColumn<int>(0);
}

void RecreateVersion4AtomResultTables(const std::filesystem::path & database_path)
{
    data_test::ExecuteSqlWithForeignKeysOff(database_path, "DROP TABLE model_atom_posterior;");
    data_test::ExecuteSqlWithForeignKeysOff(
        database_path,
        "CREATE TABLE model_atom_posterior ("
        "key_tag TEXT, "
        "class_key TEXT, "
        "serial_id INTEGER, "
        "amplitude_estimate_posterior DOUBLE, "
        "width_estimate_posterior DOUBLE, "
        "intercept_estimate_posterior DOUBLE, "
        "amplitude_variance_posterior DOUBLE, "
        "width_variance_posterior DOUBLE, "
        "intercept_variance_posterior DOUBLE, "
        "outlier_tag INTEGER, "
        "statistical_distance DOUBLE, "
        "PRIMARY KEY (key_tag, class_key, serial_id), "
        "FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE"
        ");");

    data_test::ExecuteSqlWithForeignKeysOff(database_path, "DROP TABLE model_atom_group_potential;");
    data_test::ExecuteSqlWithForeignKeysOff(
        database_path,
        "CREATE TABLE model_atom_group_potential ("
        "key_tag TEXT, "
        "class_key TEXT, "
        "group_key INTEGER, "
        "member_size INTEGER, "
        "amplitude_estimate_mean DOUBLE, "
        "width_estimate_mean DOUBLE, "
        "intercept_estimate_mean DOUBLE, "
        "amplitude_estimate_mdpde DOUBLE, "
        "width_estimate_mdpde DOUBLE, "
        "intercept_estimate_mdpde DOUBLE, "
        "amplitude_estimate_prior DOUBLE, "
        "width_estimate_prior DOUBLE, "
        "intercept_estimate_prior DOUBLE, "
        "amplitude_variance_prior DOUBLE, "
        "width_variance_prior DOUBLE, "
        "intercept_variance_prior DOUBLE, "
        "alpha_g DOUBLE, "
        "PRIMARY KEY (key_tag, class_key, group_key), "
        "FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE"
        ");");
}

} // namespace

TEST(DataObjectSchemaLifecycleTest, EmptyDatabaseBootstrapsNormalizedSchema)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_bootstrap" };
    const auto database_path{ temp_dir.path() / "bootstrap.sqlite" };

    rg::SQLitePersistence database_manager{ database_path };

    EXPECT_EQ(data_test::GetUserVersion(database_path), 6);
    EXPECT_TRUE(data_test::HasTable(database_path, "object_catalog"));
    EXPECT_FALSE(data_test::HasTable(database_path, "object_metadata"));
    EXPECT_TRUE(data_test::HasTable(database_path, "model_object"));
    EXPECT_TRUE(data_test::HasTable(database_path, "map_list"));
    EXPECT_FALSE(data_test::HasColumn(database_path, "model_atom_posterior", "class_key"));
    EXPECT_FALSE(data_test::HasColumn(database_path, "model_atom_group_potential", "class_key"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_bond_posterior", "class_key"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_bond_group_potential", "class_key"));
}

TEST(DataObjectSchemaLifecycleTest, EmptyDatabaseBootstrapsUpdatedSamplingEntryColumns)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_updated_sampling_columns" };
    const auto database_path{ temp_dir.path() / "updated_sampling_columns.sqlite" };

    rg::SQLitePersistence database_manager{ database_path };

    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_local_potential", "updated_sampling_size"));
    EXPECT_TRUE(data_test::HasColumn(
        database_path,
        "model_atom_local_potential",
        "updated_distance_and_map_value_list"));
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
    RecreateVersion4AtomResultTables(database_path);
    data_test::SetUserVersion(database_path, 3);

    rg::SQLitePersistence database_manager{ database_path };

    EXPECT_EQ(data_test::GetUserVersion(database_path), 6);
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_local_potential", "intercept_estimate_ols"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_local_potential", "intercept_estimate_mdpde"));
    EXPECT_FALSE(data_test::HasColumn(database_path, "model_atom_posterior", "class_key"));
    EXPECT_FALSE(data_test::HasColumn(database_path, "model_atom_group_potential", "class_key"));
}

TEST(DataObjectSchemaLifecycleTest, VersionFourSchemaMigratesAtomTablesToSingleComponentClass)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_v4_atom_class_migration" };
    const auto database_path{ temp_dir.path() / "v4_atom_class.sqlite" };

    {
        rg::SQLitePersistence database_manager{ database_path };
    }
    data_test::ExecuteSql(
        database_path,
        "INSERT OR REPLACE INTO object_catalog(key_tag, object_type) VALUES ('model', 'model');");
    data_test::ExecuteSql(
        database_path,
        "INSERT OR REPLACE INTO model_object("
        "key_tag, atom_size, pdb_id, emd_id, map_resolution, resolution_method"
        ") VALUES ('model', 0, 'model', '', 0.0, '');");
    RecreateVersion4AtomResultTables(database_path);
    data_test::ExecuteSql(
        database_path,
        "INSERT INTO model_atom_posterior VALUES "
        "('model', 'component_atom_class', 1, 1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0, 1.5),"
        "('model', 'simple_atom_class', 2, 4.0, 5.0, 6.0, 0.4, 0.5, 0.6, 1, 2.5);");
    data_test::ExecuteSql(
        database_path,
        "INSERT INTO model_atom_group_potential VALUES "
        "('model', 'component_atom_class', 10, 2, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 0.1, 0.2, 0.3, 0.4),"
        "('model', 'structure_atom_class', 20, 3, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 1.1, 1.2, 1.3, 1.4);");
    data_test::SetUserVersion(database_path, 4);

    rg::SQLitePersistence database_manager{ database_path };

    EXPECT_EQ(data_test::GetUserVersion(database_path), 6);
    EXPECT_FALSE(data_test::HasColumn(database_path, "model_atom_posterior", "class_key"));
    EXPECT_FALSE(data_test::HasColumn(database_path, "model_atom_group_potential", "class_key"));
    EXPECT_EQ(data_test::CountRows(database_path, "model_atom_posterior", "model"), 1);
    EXPECT_EQ(data_test::CountRows(database_path, "model_atom_group_potential", "model"), 1);
    EXPECT_EQ(
        QuerySingleInt(database_path, "SELECT serial_id FROM model_atom_posterior WHERE key_tag = 'model';"),
        1);
    EXPECT_EQ(
        QuerySingleInt(database_path, "SELECT group_key FROM model_atom_group_potential WHERE key_tag = 'model';"),
        10);
}

TEST(DataObjectSchemaLifecycleTest, VersionFiveSchemaMigratesUpdatedSamplingEntryColumns)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_v5_updated_sampling_migration" };
    const auto database_path{ temp_dir.path() / "v5_updated_sampling.sqlite" };

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
        "intercept_estimate_ols DOUBLE, "
        "amplitude_estimate_mdpde DOUBLE, "
        "width_estimate_mdpde DOUBLE, "
        "intercept_estimate_mdpde DOUBLE, "
        "alpha_r DOUBLE, "
        "PRIMARY KEY (key_tag, serial_id), "
        "FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE"
        ");");
    data_test::SetUserVersion(database_path, 5);

    rg::SQLitePersistence database_manager{ database_path };

    EXPECT_EQ(data_test::GetUserVersion(database_path), 6);
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_local_potential", "updated_sampling_size"));
    EXPECT_TRUE(data_test::HasColumn(
        database_path,
        "model_atom_local_potential",
        "updated_distance_and_map_value_list"));
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
