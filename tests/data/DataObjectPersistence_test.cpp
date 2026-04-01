#include <gtest/gtest.h>

#include <filesystem>
#include <stdexcept>

#include <rhbm_gem/data/io/DataRepository.hpp>
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include "support/CommandTestHelpers.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

TEST(DataObjectPersistenceTest, SaveAndLoadModelWithoutRegistryState)
{
    const command_test::ScopedTempDir temp_dir{ "data_repository_schema_roundtrip" };
    const auto database_path{ temp_dir.path() / "repository.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataRepository repository{ database_path };
    auto model{ rg::ReadModel(model_path) };
    model->SetKeyTag("model");
    model->SetPdbID("MODEL_REPOSITORY");
    repository.SaveModel(*model, "model");

    auto loaded_model{ repository.LoadModel("model") };
    ASSERT_NE(loaded_model, nullptr);
    EXPECT_EQ(loaded_model->GetPdbID(), "MODEL_REPOSITORY");
    EXPECT_EQ(loaded_model->GetNumberOfAtom(), model->GetNumberOfAtom());
}

TEST(DataObjectPersistenceTest, FinalV2CatalogDatabaseRemainsLoadable)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_final_v2_loadable" };
    const auto database_path{ temp_dir.path() / "final_v2.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    {
        rg::DataRepository repository{ database_path };
        auto model{ rg::ReadModel(model_path) };
        model->SetKeyTag("model");
        repository.SaveModel(*model, "model");
        data_test::SaveTinyMapThroughRepository(repository, "map", 3.0f);
    }

    data_test::SetUserVersion(database_path, 2);

    rg::DataRepository repository{ database_path };
    EXPECT_NE(repository.LoadModel("model"), nullptr);
    EXPECT_NE(repository.LoadMap("map"), nullptr);
}

TEST(DataObjectPersistenceTest, DatabaseRoundTripPreservesChainMetadataAndSymmetryFiltering)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_roundtrip" };
    const auto database_path{ temp_dir.path() / "roundtrip.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model_keyvalue_entity.cif") };

    auto original_model{ data_test::LoadFixtureModel(model_path) };
    const auto original_chain_map{ original_model->GetChainIDListMap() };
    for (const auto & atom : original_model->GetAtomList())
    {
        atom->SetSelectedFlag(true);
    }
    original_model->ApplySymmetrySelection(false);
    const auto original_selected_count{ original_model->GetNumberOfSelectedAtom() };

    rg::DataRepository repository{ database_path };
    auto stored_model{ rg::ReadModel(model_path) };
    stored_model->SetKeyTag("model");
    repository.SaveModel(*stored_model, "model");

    auto loaded_model{ repository.LoadModel("model") };
    EXPECT_EQ(loaded_model->GetChainIDListMap(), original_chain_map);
    EXPECT_GT(data_test::CountRows(database_path, "model_chain_map", "model"), 0);

    for (const auto & atom : loaded_model->GetAtomList())
    {
        atom->SetSelectedFlag(true);
    }
    loaded_model->ApplySymmetrySelection(false);
    EXPECT_EQ(loaded_model->GetNumberOfSelectedAtom(), original_selected_count);
}

TEST(DataObjectPersistenceTest, DistinctUnsanitizedKeysDoNotCollideInV2Schema)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_key_collision" };
    const auto database_path{ temp_dir.path() / "collision.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataRepository repository{ database_path };
    auto model_a{ rg::ReadModel(model_path) };
    auto model_b{ rg::ReadModel(model_path) };
    model_a->SetKeyTag("mem_a");
    model_b->SetKeyTag("mem_b");
    model_a->SetPdbID("MODEL_A");
    model_b->SetPdbID("MODEL_B");

    repository.SaveModel(*model_a, "a-b");
    repository.SaveModel(*model_b, "a_b");

    auto loaded_model_a{ repository.LoadModel("a-b") };
    auto loaded_model_b{ repository.LoadModel("a_b") };
    EXPECT_EQ(loaded_model_a->GetPdbID(), "MODEL_A");
    EXPECT_EQ(loaded_model_b->GetPdbID(), "MODEL_B");
    EXPECT_EQ(data_test::CountRows(database_path, "model_object"), 2);
}

TEST(DataObjectPersistenceTest, SaveRenamedKeyDoesNotRenameInMemoryObject)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_rename_semantics" };
    const auto database_path{ temp_dir.path() / "rename.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataRepository repository{ database_path };
    auto model{ rg::ReadModel(model_path) };
    model->SetKeyTag("model");

    repository.SaveModel(*model, "saved_model");

    EXPECT_EQ(model->GetKeyTag(), "model");
    auto loaded_model{ repository.LoadModel("saved_model") };
    ASSERT_NE(loaded_model, nullptr);
    EXPECT_EQ(loaded_model->GetNumberOfAtom(), model->GetNumberOfAtom());
}

TEST(DataObjectPersistenceTest, LoadModelThrowsWhenDatabaseKeyIsMissing)
{
    const command_test::ScopedTempDir temp_dir{ "data_schema_missing_database_key" };
    const auto database_path{ temp_dir.path() / "missing.sqlite" };

    rg::DataRepository repository{ database_path };
    EXPECT_THROW((void)repository.LoadModel("missing"), std::runtime_error);
}
