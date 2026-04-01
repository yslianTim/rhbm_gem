#include <gtest/gtest.h>

#include <filesystem>
#include <stdexcept>

#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include "io/sqlite/SQLitePersistence.hpp"
#include "support/CommandTestHelpers.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

TEST(SQLitePersistenceTypedApiTest, SaveAndLoadModelRoundTrip)
{
    const command_test::ScopedTempDir temp_dir{ "sqlite_persistence_typed_model" };
    const auto database_path{ temp_dir.path() / "typed_model.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::SQLitePersistence persistence{ database_path };
    auto model{ rg::ReadModel(model_path) };
    model->SetKeyTag("model");
    model->SetPdbID("SQLITE_TYPED_MODEL");

    persistence.SaveModel(*model, "model");

    auto loaded_model{ persistence.LoadModel("model") };
    ASSERT_NE(loaded_model, nullptr);
    EXPECT_EQ(loaded_model->GetPdbID(), "SQLITE_TYPED_MODEL");
    EXPECT_EQ(loaded_model->GetNumberOfAtom(), model->GetNumberOfAtom());
}

TEST(SQLitePersistenceTypedApiTest, SaveAndLoadMapRoundTrip)
{
    const command_test::ScopedTempDir temp_dir{ "sqlite_persistence_typed_map" };
    const auto database_path{ temp_dir.path() / "typed_map.sqlite" };

    rg::SQLitePersistence persistence{ database_path };
    auto map{ data_test::MakeMapObject() };
    map.SetKeyTag("map");

    persistence.SaveMap(map, "map");

    auto loaded_map{ persistence.LoadMap("map") };
    ASSERT_NE(loaded_map, nullptr);
    EXPECT_EQ(loaded_map->GetGridSize(), map.GetGridSize());
    EXPECT_FLOAT_EQ(loaded_map->GetMapValue(0), map.GetMapValue(0));
}

TEST(SQLitePersistenceTypedApiTest, LoadModelThrowsWhenCatalogRowIsMap)
{
    const command_test::ScopedTempDir temp_dir{ "sqlite_persistence_typed_model_mismatch" };
    const auto database_path{ temp_dir.path() / "typed_model_mismatch.sqlite" };

    rg::SQLitePersistence persistence{ database_path };
    auto map{ data_test::MakeMapObject() };

    persistence.SaveMap(map, "shared_key");

    EXPECT_THROW((void)persistence.LoadModel("shared_key"), std::runtime_error);
}

TEST(SQLitePersistenceTypedApiTest, LoadMapThrowsWhenCatalogRowIsModel)
{
    const command_test::ScopedTempDir temp_dir{ "sqlite_persistence_typed_map_mismatch" };
    const auto database_path{ temp_dir.path() / "typed_map_mismatch.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::SQLitePersistence persistence{ database_path };
    auto model{ rg::ReadModel(model_path) };

    persistence.SaveModel(*model, "shared_key");

    EXPECT_THROW((void)persistence.LoadMap("shared_key"), std::runtime_error);
}
