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

    EXPECT_EQ(data_test::GetUserVersion(database_path), 19);
    for (const auto table_name : std::array<std::string_view, 9>{
             "model_object",
             "model_joint_result",
             "model_chain_map",
             "model_component",
             "model_component_atom",
             "model_component_bond",
             "model_atom",
             "model_bond",
             "model_stage_result" })
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
    CreateVersionedMarkerDatabase(database_path, 17);

    EXPECT_THROW((void)rg::DataRepository(database_path), std::runtime_error);
    EXPECT_EQ(data_test::GetUserVersion(database_path), 17);
    EXPECT_EQ(data_test::CountRows(database_path, "legacy_marker"), 1);
}

TEST(DataObjectSchemaLifecycleTest, VersionFourteenAndFifteenSchemasAreRejectedWithoutModification)
{
    ExpectVersionedDatabaseRejectedWithoutMutation(14);
    ExpectVersionedDatabaseRejectedWithoutMutation(15);
}

TEST(DataObjectSchemaLifecycleTest, VersionSixteenSchemaIsRejectedWithoutModification)
{
    ExpectVersionedDatabaseRejectedWithoutMutation(16);
}

TEST(DataObjectSchemaLifecycleTest, CanonicalAnalysisReplacesLegacyTables)
{
    const command_test::ScopedTempDir dir{"canonical_analysis_tables"}; const auto path=dir.path()/"model.sqlite";
    {rg::DataRepository repo(path);}
    EXPECT_TRUE(data_test::HasColumn(path,"model_stage_result","result_json"));
    for(const auto table:{"model_atom_local_potential","model_atom_posterior","model_atom_group_potential"})
        EXPECT_FALSE(data_test::HasTable(path,table));
    EXPECT_TRUE(data_test::HasColumn(path,"model_atom","is_selected"));
    EXPECT_TRUE(data_test::HasColumn(path,"model_bond","is_selected"));
}
