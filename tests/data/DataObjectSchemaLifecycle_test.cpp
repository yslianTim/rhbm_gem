#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <stdexcept>
#include <string_view>

#include <rhbm_gem/data/io/DataRepository.hpp>
#include "io/sqlite/SQLiteWrapper.hpp"
#include "support/CommandTestHelpers.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

namespace {

void CreateVersionedMarkerDatabase(
    const std::filesystem::path & database_path,
    int user_version)
{
    rg::SQLiteWrapper database{ database_path };
    database.Execute("CREATE TABLE legacy_marker (value INTEGER PRIMARY KEY);");
    database.Execute("INSERT INTO legacy_marker(value) VALUES (1);");
    database.Execute("PRAGMA user_version = " + std::to_string(user_version) + ";");
}

void ExpectVersionedDatabaseRejectedWithoutMutation(int user_version)
{
    const command_test::ScopedTempDir temp_dir{
        "data_schema_reject_" + std::to_string(user_version) };
    const auto database_path{ temp_dir.path() / "legacy.sqlite" };
    CreateVersionedMarkerDatabase(database_path, user_version);

    EXPECT_THROW((void)rg::DataRepository(database_path), std::runtime_error);
    EXPECT_EQ(data_test::GetUserVersion(database_path), user_version);
    EXPECT_TRUE(data_test::HasTable(database_path, "legacy_marker"));
    EXPECT_EQ(data_test::CountRows(database_path, "legacy_marker"), 1);
}

} // namespace

TEST(DataObjectSchemaLifecycleTest, EmptyDatabaseBootstrapsNormalizedSchema)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_bootstrap" };
    const auto database_path{ temp_dir.path() / "bootstrap.sqlite" };

    { rg::DataRepository repository{ database_path }; }

    EXPECT_EQ(data_test::GetUserVersion(database_path), 15);
    for (const auto table_name : std::array<std::string_view, 10>{
             "model_object",
             "model_chain_map",
             "model_component",
             "model_component_atom",
             "model_component_bond",
             "model_atom",
             "model_bond",
             "model_atom_local_potential",
             "model_atom_posterior",
             "model_atom_group_potential" })
    {
        EXPECT_TRUE(data_test::HasTable(database_path, std::string(table_name)));
    }
    EXPECT_FALSE(data_test::HasTable(database_path, "object_catalog"));
    EXPECT_FALSE(data_test::HasTable(database_path, "map_list"));
    EXPECT_FALSE(data_test::HasTable(database_path, "model_bond_local_potential"));
    EXPECT_FALSE(data_test::HasTable(database_path, "model_bond_posterior"));
    EXPECT_FALSE(data_test::HasTable(database_path, "model_bond_group_potential"));

    EXPECT_NO_THROW((void)rg::DataRepository(database_path));
}

TEST(DataObjectSchemaLifecycleTest, EmptyDatabaseBootstrapsRawAndPeelingSamplingEntryColumns)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_sampling_columns" };
    const auto database_path{ temp_dir.path() / "sampling.sqlite" };
    { rg::DataRepository repository{ database_path }; }

    EXPECT_TRUE(data_test::HasColumn(
        database_path,
        "model_atom_local_potential",
        "raw_distance_and_map_value_list"));
    EXPECT_TRUE(data_test::HasColumn(
        database_path,
        "model_atom_local_potential",
        "peeling_distance_and_map_value_list"));
    EXPECT_FALSE(data_test::HasColumn(
        database_path, "model_atom_local_potential", "raw_sampling_size"));
    EXPECT_FALSE(data_test::HasColumn(
        database_path, "model_atom_local_potential", "peeling_sampling_size"));
    EXPECT_FALSE(data_test::HasColumn(
        database_path, "model_atom_group_potential", "member_size"));
}

TEST(DataObjectSchemaLifecycleTest, EmptyDatabaseBootstrapsGaussianInterceptColumns)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_gaussian_columns" };
    const auto database_path{ temp_dir.path() / "gaussian.sqlite" };
    { rg::DataRepository repository{ database_path }; }

    for (const auto suffix : { "1st", "2nd", "3rd" })
    {
        EXPECT_TRUE(data_test::HasColumn(
            database_path,
            "model_atom_local_potential",
            "intercept_estimate_ols_" + std::string(suffix)));
        EXPECT_FALSE(data_test::HasColumn(
            database_path,
            "model_atom_group_potential",
            "intercept_estimate_prior_" + std::string(suffix)));
    }
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom", "is_selected"));
    EXPECT_TRUE(data_test::HasColumn(database_path, "model_bond", "is_selected"));
}

TEST(DataObjectSchemaLifecycleTest, VersionNineSchemaIsRejectedWithoutModification)
{
    ExpectVersionedDatabaseRejectedWithoutMutation(9);
}

TEST(DataObjectSchemaLifecycleTest, VersionTenSchemaIsRejectedWithoutModification)
{
    ExpectVersionedDatabaseRejectedWithoutMutation(10);
}

TEST(DataObjectSchemaLifecycleTest, VersionElevenSchemaIsRejectedWithoutModification)
{
    ExpectVersionedDatabaseRejectedWithoutMutation(11);
}

TEST(DataObjectSchemaLifecycleTest, VersionTwelveSchemaIsRejectedWithoutModification)
{
    ExpectVersionedDatabaseRejectedWithoutMutation(12);
}

TEST(DataObjectSchemaLifecycleTest, VersionThirteenSchemaIsRejectedWithoutModification)
{
    ExpectVersionedDatabaseRejectedWithoutMutation(13);
}

TEST(DataObjectSchemaLifecycleTest, UnknownSchemaVersionThrows)
{
    ExpectVersionedDatabaseRejectedWithoutMutation(99);
}

TEST(DataObjectSchemaLifecycleTest, VersionOneSchemaFailsFast)
{
    ExpectVersionedDatabaseRejectedWithoutMutation(1);
}

TEST(DataObjectSchemaLifecycleTest, Version2MetadataBasedShapeFailsFast)
{
    ExpectVersionedDatabaseRejectedWithoutMutation(2);
}

TEST(DataObjectSchemaLifecycleTest, ManagedButUnversionedDatabaseFailsFast)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_unversioned" };
    const auto database_path{ temp_dir.path() / "unversioned.sqlite" };
    CreateVersionedMarkerDatabase(database_path, 0);

    EXPECT_THROW((void)rg::DataRepository(database_path), std::runtime_error);
    EXPECT_EQ(data_test::GetUserVersion(database_path), 0);
    EXPECT_EQ(data_test::CountRows(database_path, "legacy_marker"), 1);
}

TEST(DataObjectSchemaLifecycleTest, MixedUnknownSchemaFailsFast)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_mixed_unknown" };
    const auto database_path{ temp_dir.path() / "mixed.sqlite" };
    CreateVersionedMarkerDatabase(database_path, 15);

    EXPECT_THROW((void)rg::DataRepository(database_path), std::runtime_error);
    EXPECT_EQ(data_test::GetUserVersion(database_path), 15);
    EXPECT_EQ(data_test::CountRows(database_path, "legacy_marker"), 1);
}

TEST(DataObjectSchemaLifecycleTest, VersionFourteenSchemaIsRejectedWithoutModification)
{
    ExpectVersionedDatabaseRejectedWithoutMutation(14);
}

TEST(DataObjectSchemaLifecycleTest, EmptyDatabaseBootstrapsSingleGroupGaussianResult)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_single_group_result" };
    const auto database_path{ temp_dir.path() / "group.sqlite" };
    { rg::DataRepository repository{ database_path }; }

    for (const auto column : {
             "amplitude_estimate_mean", "width_estimate_mean", "intercept_estimate_mean",
             "amplitude_estimate_mdpde", "width_estimate_mdpde", "intercept_estimate_mdpde",
             "amplitude_estimate_prior", "width_estimate_prior", "intercept_estimate_prior",
             "amplitude_variance_prior", "width_variance_prior", "intercept_variance_prior",
             "alpha_g" })
    {
        EXPECT_TRUE(data_test::HasColumn(database_path, "model_atom_group_potential", column));
        for (const auto suffix : { "_1st", "_2nd", "_3rd" })
        {
            EXPECT_FALSE(data_test::HasColumn(database_path, "model_atom_group_potential",
                std::string(column) + suffix));
        }
    }
}
