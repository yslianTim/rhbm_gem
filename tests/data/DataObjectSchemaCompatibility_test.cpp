#include <gtest/gtest.h>

#include <filesystem>
#include <stdexcept>

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include "io/sqlite/SQLitePersistence.hpp"
#include "support/CommandTestHelpers.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

TEST(DataObjectSchemaCompatibilityTest, Version2WithObjectMetadataFailsFast)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_v2_metadata_row" };
    const auto database_path{ temp_dir.path() / "metadata_row.sqlite" };
    const auto map_object{ data_test::MakeTinyMapObject() };

    {
        rg::SQLitePersistence database_manager{ database_path };
        database_manager.Save(map_object, "map_only");
    }

    data_test::ExecuteSql(
        database_path,
        "CREATE TABLE IF NOT EXISTS object_metadata (key_tag TEXT PRIMARY KEY, object_type TEXT);");
    data_test::ExecuteSql(
        database_path,
        "INSERT INTO object_metadata(key_tag, object_type) VALUES ('map_only', 'map');");
    data_test::SetUserVersion(database_path, 2);

    EXPECT_THROW((void)rg::SQLitePersistence(database_path), std::runtime_error);
}

TEST(DataObjectSchemaCompatibilityTest, Version2MetadataBasedShapeFailsFast)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_v2_metadata_shape" };
    const auto database_path{ temp_dir.path() / "metadata_shape.sqlite" };

    data_test::CreateVersion2MetadataBasedMapShapeDatabase(database_path);
    EXPECT_THROW((void)rg::SQLitePersistence(database_path), std::runtime_error);
}

TEST(DataObjectSchemaCompatibilityTest, ManagedButUnversionedDatabaseFailsFast)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_unversioned_nonlegacy" };
    const auto database_path{ temp_dir.path() / "managed_unversioned.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    {
        rg::DataRepository repository{ database_path };
        auto model{ rg::ReadModel(model_path) };
        model->SetKeyTag("model");
        repository.SaveModel(*model, "model");
    }

    data_test::SetUserVersion(database_path, 0);
    EXPECT_THROW((void)rg::SQLitePersistence(database_path), std::runtime_error);
}

TEST(DataObjectSchemaCompatibilityTest, MixedUnknownSchemaFailsFast)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_mixed_unknown" };
    const auto database_path{ temp_dir.path() / "mixed_unknown.sqlite" };

    data_test::ExecuteSql(database_path, "CREATE TABLE foreign_table (id INTEGER PRIMARY KEY);");
    EXPECT_THROW((void)rg::SQLitePersistence(database_path), std::runtime_error);
}
