#include <gtest/gtest.h>

#include <filesystem>
#include <stdexcept>
#include <string>

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
    const command_test::ScopedTempDir temp_dir{temp_dir_name};
    const auto path{temp_dir.path() / database_name};
    { rg::DataRepository repository{path}; }
    ASSERT_EQ(data_test::GetUserVersion(path), 20);

    mutate_database(path);
    EXPECT_THROW((void)rg::DataRepository(path), std::runtime_error);
    EXPECT_EQ(data_test::GetUserVersion(path), 20);
}

void RecreateChainMapTable(
    const std::filesystem::path & path,
    const std::string & constraints)
{
    data_test::ExecuteSqlWithForeignKeysOff(path, "DROP TABLE model_chain_map;");
    data_test::ExecuteSqlWithForeignKeysOff(
        path,
        "CREATE TABLE model_chain_map (key_tag TEXT, entity_id TEXT, chain_ordinal INTEGER, "
        "chain_id TEXT, " + constraints + ");");
}

} // namespace

TEST(DataObjectSchemaValidationTest, CurrentSchemaMissingStageResultTableThrows)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_missing_stage_table", "missing_stage_table.sqlite",
        [](const auto & path)
        {
            data_test::ExecuteSqlWithForeignKeysOff(path, "DROP TABLE model_stage_result;");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaMissingStandardQScoreColumnThrows)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_missing_qscore", "missing_qscore.sqlite",
        [](const auto & path)
        {
            data_test::ExecuteSqlWithForeignKeysOff(path,
                "ALTER TABLE model_atom DROP COLUMN standard_qscore;");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaMissingReferenceHeightColumnThrows)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_missing_reference_height", "missing_reference_height.sqlite",
        [](const auto & path)
        {
            data_test::ExecuteSqlWithForeignKeysOff(path,
                "ALTER TABLE model_object DROP COLUMN reference_height;");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaMissingReferenceOffsetColumnThrows)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_missing_reference_offset", "missing_reference_offset.sqlite",
        [](const auto & path)
        {
            data_test::ExecuteSqlWithForeignKeysOff(path,
                "ALTER TABLE model_object DROP COLUMN reference_offset;");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaRejectsUnexpectedTable)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_unexpected_table", "unexpected_table.sqlite",
        [](const auto & path)
        {
            data_test::ExecuteSql(path, "CREATE TABLE object_catalog (key_tag TEXT PRIMARY KEY);");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaRejectsMissingForeignKeys)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_missing_foreign_key", "missing_foreign_key.sqlite",
        [](const auto & path)
        {
            RecreateChainMapTable(path, "PRIMARY KEY (key_tag, entity_id, chain_ordinal)");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaRejectsWrongPrimaryKey)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_wrong_primary_key", "wrong_primary_key.sqlite",
        [](const auto & path)
        {
            RecreateChainMapTable(path,
                "PRIMARY KEY (key_tag), FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE");
        });
}

TEST(DataObjectSchemaValidationTest, CurrentSchemaRequiresNotNullSelectionFlag)
{
    ExpectCurrentSchemaValidationFailure(
        "data_schema_nullable_selection", "nullable_selection.sqlite",
        [](const auto & path)
        {
            data_test::ExecuteSqlWithForeignKeysOff(path, "DROP TABLE model_atom;");
            data_test::ExecuteSqlWithForeignKeysOff(path,
                "CREATE TABLE model_atom (key_tag TEXT, serial_id INTEGER, sequence_id INTEGER, "
                "component_id TEXT, atom_id TEXT, chain_id TEXT, indicator TEXT, occupancy DOUBLE, "
                "temperature DOUBLE, element INTEGER, structure INTEGER, is_special_atom INTEGER, "
                "is_selected INTEGER DEFAULT 0, position_x DOUBLE, position_y DOUBLE, position_z DOUBLE, "
                "component_key INTEGER, atom_key INTEGER, standard_qscore DOUBLE DEFAULT 0.0, "
                "PRIMARY KEY (key_tag, serial_id), "
                "FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE);");
        });
}

TEST(DataObjectSchemaValidationTest, JointResultPayloadMustBeNotNull)
{
    ExpectCurrentSchemaValidationFailure(
        "joint_nullable_payload", "nullable_joint.sqlite",
        [](const auto & path)
        {
            data_test::ExecuteSql(path, "DROP TABLE model_joint_result;");
            data_test::ExecuteSql(path,
                "CREATE TABLE model_joint_result (key_tag TEXT PRIMARY KEY, result_json TEXT, "
                "FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE);");
        });
}

TEST(DataObjectSchemaValidationTest, StageResultPayloadMustBeNotNull)
{
    ExpectCurrentSchemaValidationFailure(
        "stage_nullable_payload", "nullable_stage.sqlite",
        [](const auto & path)
        {
            data_test::ExecuteSql(path, "DROP TABLE model_stage_result;");
            data_test::ExecuteSql(path,
                "CREATE TABLE model_stage_result (key_tag TEXT PRIMARY KEY, result_json TEXT, "
                "FOREIGN KEY(key_tag) REFERENCES model_object(key_tag) ON DELETE CASCADE);");
        });
}

TEST(DataObjectSchemaValidationTest, ForeignKeyRejectsOrphanModelChildRows)
{
    const command_test::ScopedTempDir temp_dir{"data_schema_orphan_child"};
    const auto path{temp_dir.path() / "orphan.sqlite"};
    { rg::DataRepository repository{path}; }

    EXPECT_THROW(
        data_test::ExecuteSql(path,
            "INSERT INTO model_atom(key_tag, serial_id, is_selected) VALUES ('missing', 1, 0);"),
        std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, DeletingModelRootCascadesPayloadRows)
{
    const command_test::ScopedTempDir temp_dir{"data_schema_model_cascade"};
    const auto path{temp_dir.path() / "cascade.sqlite"};
    auto model{data_test::MakeModelWithBond()};

    {
        rg::DataRepository repository{path};
        repository.SaveModel(*model, "model");
    }
    ASSERT_GT(data_test::CountRows(path, "model_atom", "model"), 0);
    data_test::ExecuteSql(path, "DELETE FROM model_object WHERE key_tag = 'model';");
    EXPECT_EQ(data_test::CountRows(path, "model_object"), 0);
    EXPECT_EQ(data_test::CountRows(path, "model_atom"), 0);
    EXPECT_EQ(data_test::CountRows(path, "model_bond"), 0);
}
