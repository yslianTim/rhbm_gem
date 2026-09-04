#include <gtest/gtest.h>

#include <filesystem>
#include <stdexcept>

#include <rhbm_gem/data/io/DataRepository.hpp>
#include "support/CommandTestHelpers.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

namespace {

template <typename MutateFn>
void ExpectCurrentSchemaValidationFailure(
    const char * temp_dir_name,
    const char * database_name,
    MutateFn mutate_database)
{
    const command_test::ScopedTempDir temp_dir{ temp_dir_name };
    const auto database_path{ temp_dir.path() / database_name };
    { rg::DataRepository repository{ database_path }; }
    ASSERT_EQ(data_test::GetUserVersion(database_path), 15);

    mutate_database(database_path);
    EXPECT_THROW((void)rg::DataRepository(database_path), std::runtime_error);
    EXPECT_EQ(data_test::GetUserVersion(database_path), 15);
}

void RecreateChainMapTable(
    const std::filesystem::path & database_path,
    const std::string & table_constraints)
{
    data_test::ExecuteSqlWithForeignKeysOff(
        database_path,
        "DROP TABLE model_chain_map;");
    data_test::ExecuteSqlWithForeignKeysOff(
        database_path,
        "CREATE TABLE model_chain_map ("
        "key_tag TEXT, entity_id TEXT, chain_ordinal INTEGER, chain_id TEXT, "
        + table_constraints + ");");
}

} // namespace

TEST(DataObjectSchemaValidationTest, CurrentSchemaMissingPeelingNeighborCountColumnThrows)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_missing_neighbor_count",
        "missing_neighbor_count.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::ExecuteSqlWithForeignKeysOff(
                database_path,
                "ALTER TABLE model_atom_local_potential "
                "DROP COLUMN neighbor_count_for_peeling;");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaMissingLocalGaussianStageColumnThrows)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_missing_local_stage",
        "missing_local_stage.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::ExecuteSqlWithForeignKeysOff(
                database_path,
                "ALTER TABLE model_atom_local_potential DROP COLUMN alpha_r_2nd;");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaRejectsLegacyLocalGaussianColumn)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_legacy_local_column",
        "legacy_local_column.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::ExecuteSql(
                database_path,
                "ALTER TABLE model_atom_local_potential "
                "ADD COLUMN alpha_r DOUBLE DEFAULT 0.0;");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaMissingGroupGaussianColumnThrows)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_missing_group_result",
        "missing_group_result.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::ExecuteSqlWithForeignKeysOff(
                database_path,
                "ALTER TABLE model_atom_group_potential DROP COLUMN alpha_g;");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaRejectsLegacyGroupGaussianColumn)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_legacy_group_column",
        "legacy_group_column.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::ExecuteSql(
                database_path,
                "ALTER TABLE model_atom_group_potential "
                "ADD COLUMN member_size INTEGER;");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaMissingRequiredTableThrows)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_missing_table",
        "missing_table.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::ExecuteSqlWithForeignKeysOff(
                database_path, "DROP TABLE model_atom_posterior;");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaMissingStandardQScoreColumnThrows)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_missing_qscore",
        "missing_qscore.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::ExecuteSqlWithForeignKeysOff(
                database_path,
                "ALTER TABLE model_atom DROP COLUMN standard_qscore;");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaMissingReferenceHeightColumnThrows)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_missing_reference_height",
        "missing_reference_height.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::ExecuteSqlWithForeignKeysOff(
                database_path,
                "ALTER TABLE model_object DROP COLUMN reference_height;");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaMissingReferenceOffsetColumnThrows)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_missing_reference_offset",
        "missing_reference_offset.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::ExecuteSqlWithForeignKeysOff(
                database_path,
                "ALTER TABLE model_object DROP COLUMN reference_offset;");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaRejectsUnexpectedTable)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_unexpected_table",
        "unexpected_table.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::ExecuteSql(
                database_path,
                "CREATE TABLE object_catalog (key_tag TEXT PRIMARY KEY);");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaValidationRejectsMissingForeignKeys)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_missing_foreign_key",
        "missing_foreign_key.sqlite",
        [](const std::filesystem::path & database_path)
        {
            RecreateChainMapTable(
                database_path,
                "PRIMARY KEY (key_tag, entity_id, chain_ordinal)");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaValidationRejectsWrongPrimaryKey)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_wrong_primary_key",
        "wrong_primary_key.sqlite",
        [](const std::filesystem::path & database_path)
        {
            RecreateChainMapTable(
                database_path,
                "PRIMARY KEY (key_tag), "
                "FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE");
        });
}

TEST(DataObjectSchemaValidationTest, ForeignKeyRejectsOrphanModelChildRows)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_orphan_child" };
    const auto database_path{ temp_dir.path() / "orphan.sqlite" };
    { rg::DataRepository repository{ database_path }; }

    EXPECT_THROW(
        data_test::ExecuteSql(
            database_path,
            "INSERT INTO model_atom(key_tag, serial_id, is_selected) "
            "VALUES ('missing', 1, 0);"),
        std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, DeletingModelRootCascadesPayloadRows)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_model_cascade" };
    const auto database_path{ temp_dir.path() / "cascade.sqlite" };
    auto model{ data_test::MakeModelWithBond() };

    {
        rg::DataRepository repository{ database_path };
        repository.SaveModel(*model, "model");
    }
    ASSERT_GT(data_test::CountRows(database_path, "model_atom", "model"), 0);

    data_test::ExecuteSql(
        database_path,
        "DELETE FROM model_object WHERE key_tag = 'model';");

    EXPECT_EQ(data_test::CountRows(database_path, "model_object"), 0);
    EXPECT_EQ(data_test::CountRows(database_path, "model_atom"), 0);
    EXPECT_EQ(data_test::CountRows(database_path, "model_bond"), 0);
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaRejectsGroupStageColumn)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_group_stage_column",
        "group_stage_column.sqlite",
        [](const std::filesystem::path & database_path)
        {
            data_test::ExecuteSql(database_path,
                "ALTER TABLE model_atom_group_potential ADD COLUMN alpha_g_3rd DOUBLE;");
        });
}
