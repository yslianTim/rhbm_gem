#include <gtest/gtest.h>

#include <filesystem>
#include <stdexcept>

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/FileIO.hpp>
#include "io/sqlite/SQLitePersistence.hpp"
#include "support/CommandTestHelpers.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

TEST(DataObjectSchemaValidationTest, NormalizedV2DatabaseMissingRequiredTableThrows)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_missing_v2_table" };
    const auto database_path{ temp_dir.path() / "missing_v2.sqlite" };

    {
        rg::SQLitePersistence database_manager{ database_path };
        EXPECT_EQ(data_test::GetUserVersion(database_path), 2);
    }
    data_test::ExecuteSql(database_path, "DROP TABLE model_bond_group_potential;");
    data_test::SetUserVersion(database_path, 2);

    EXPECT_THROW((void)rg::SQLitePersistence(database_path), std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, FinalV2SchemaValidationRejectsMissingForeignKeys)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_missing_fk" };
    const auto database_path{ temp_dir.path() / "missing_fk.sqlite" };

    {
        rg::SQLitePersistence database_manager{ database_path };
        EXPECT_EQ(data_test::GetUserVersion(database_path), 2);
    }

    data_test::ExecuteSql(database_path, "DROP TABLE map_list;");
    data_test::ExecuteSql(
        database_path,
        "CREATE TABLE map_list ("
        "key_tag TEXT PRIMARY KEY, "
        "grid_size_x INTEGER, grid_size_y INTEGER, grid_size_z INTEGER, "
        "grid_spacing_x DOUBLE, grid_spacing_y DOUBLE, grid_spacing_z DOUBLE, "
        "origin_x DOUBLE, origin_y DOUBLE, origin_z DOUBLE, "
        "map_value_array BLOB"
        ");");
    data_test::SetUserVersion(database_path, 2);

    EXPECT_FALSE(
        data_test::HasForeignKey(
            database_path,
            "map_list",
            "key_tag",
            "object_catalog",
            "key_tag",
            "CASCADE"));
    EXPECT_THROW((void)rg::SQLitePersistence(database_path), std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, FinalV2SchemaValidationRejectsMissingRequiredCatalogColumns)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_bad_catalog_columns" };
    const auto database_path{ temp_dir.path() / "bad_catalog_columns.sqlite" };

    {
        rg::SQLitePersistence database_manager{ database_path };
        EXPECT_EQ(data_test::GetUserVersion(database_path), 2);
    }

    data_test::ExecuteSqlWithForeignKeysOff(database_path, "DROP TABLE object_catalog;");
    data_test::ExecuteSqlWithForeignKeysOff(
        database_path,
        "CREATE TABLE object_catalog (key_tag TEXT PRIMARY KEY);");
    data_test::SetUserVersion(database_path, 2);

    EXPECT_THROW((void)rg::SQLitePersistence(database_path), std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, FinalV2SchemaValidationRejectsUnknownObjectTypeValue)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_bad_catalog_type" };
    const auto database_path{ temp_dir.path() / "bad_catalog_type.sqlite" };

    {
        rg::SQLitePersistence database_manager{ database_path };
        EXPECT_EQ(data_test::GetUserVersion(database_path), 2);
    }

    data_test::ExecuteSqlWithForeignKeysOff(database_path, "DROP TABLE object_catalog;");
    data_test::ExecuteSqlWithForeignKeysOff(
        database_path,
        "CREATE TABLE object_catalog (key_tag TEXT PRIMARY KEY, object_type TEXT NOT NULL);");
    data_test::ExecuteSqlWithForeignKeysOff(
        database_path,
        "INSERT INTO object_catalog(key_tag, object_type) VALUES ('unknown_root', 'unknown');");
    data_test::SetUserVersion(database_path, 2);

    EXPECT_THROW((void)rg::SQLitePersistence(database_path), std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, ForeignKeyRejectsOrphanModelChildRows)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_fk_orphan" };
    const auto database_path{ temp_dir.path() / "orphan.sqlite" };
    rg::SQLitePersistence database_manager{ database_path };

    EXPECT_THROW(
        data_test::ExecuteSql(
            database_path,
            "INSERT INTO model_chain_map(key_tag, entity_id, chain_ordinal, chain_id) "
            "VALUES ('missing', '1', 0, 'A');"),
        std::runtime_error);
}

TEST(DataObjectSchemaValidationTest, DeletingCatalogRootCascadesPayloadRows)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_catalog_cascade" };
    const auto database_path{ temp_dir.path() / "cascade.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataRepository repository{ database_path };
    auto model{ rg::ReadModel(model_path) };
    model->SetKeyTag("model");
    repository.SaveModel(*model, "model");

    data_test::ExecuteSql(database_path, "DELETE FROM object_catalog WHERE key_tag = 'model';");

    EXPECT_EQ(data_test::CountRows(database_path, "object_catalog"), 0);
    EXPECT_EQ(data_test::CountRows(database_path, "model_object"), 0);
    EXPECT_EQ(data_test::CountRows(database_path, "model_atom"), 0);
}
