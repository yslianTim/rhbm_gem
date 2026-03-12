#include <gtest/gtest.h>

#include "data/schema/DataObjectIOSchemaTestSupport.hpp"

using namespace data_object_io_schema_test;

TEST(DataObjectIOContractTest, ChainMetadataPersistsAcrossDatabaseRoundTrip)
{
    const command_test::ScopedTempDir temp_dir{ "schema_chain_roundtrip" };
    const auto database_path{ temp_dir.path() / "chain_roundtrip.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model_keyvalue_entity.cif") };

    rg::DataObjectManager manager{};
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    auto original_model{ manager.GetTypedDataObject<rg::ModelObject>("model") };
    const auto original_chain_map{ original_model->GetChainIDListMap() };

    manager.SaveDataObject("model");
    manager.ClearDataObjects();
    manager.LoadDataObject("model");

    auto loaded_model{ manager.GetTypedDataObject<rg::ModelObject>("model") };
    EXPECT_EQ(loaded_model->GetChainIDListMap(), original_chain_map);
    EXPECT_GT(CountRows(database_path, "model_chain_map", "model"), 0);
}

TEST(DataObjectIOContractTest, SymmetryFilteringMatchesAfterDatabaseReload)
{
    const command_test::ScopedTempDir temp_dir{ "schema_symmetry_roundtrip" };
    const auto database_path{ temp_dir.path() / "symmetry.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model_keyvalue_entity.cif") };

    auto original_model{ LoadFixtureModel(model_path) };
    for (const auto & atom : original_model->GetAtomList()) atom->SetSelectedFlag(true);
    original_model->FilterAtomFromSymmetry(false);
    original_model->Update();
    const auto original_selected_count{ original_model->GetNumberOfSelectedAtom() };

    rg::DataObjectManager manager{};
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    for (const auto & atom : manager.GetTypedDataObject<rg::ModelObject>("model")->GetAtomList())
    {
        atom->SetSelectedFlag(true);
    }
    manager.GetTypedDataObject<rg::ModelObject>("model")->Update();
    manager.SaveDataObject("model");
    manager.ClearDataObjects();
    manager.LoadDataObject("model");

    auto loaded_model{ manager.GetTypedDataObject<rg::ModelObject>("model") };
    for (const auto & atom : loaded_model->GetAtomList()) atom->SetSelectedFlag(true);
    loaded_model->FilterAtomFromSymmetry(false);
    loaded_model->Update();
    EXPECT_EQ(loaded_model->GetNumberOfSelectedAtom(), original_selected_count);
}

TEST(DataObjectIOContractTest, FileImportTransfersBondKeySystem)
{
    rg::DataObjectManager manager{};
    const auto model_path{
        command_test::TestDataPath("test_model_auth_seq_alnum_struct_conn.cif") };

    ASSERT_NO_THROW(manager.ProcessFile(model_path, "model"));
    auto model{ manager.GetTypedDataObject<rg::ModelObject>("model") };
    ASSERT_GE(model->GetNumberOfBond(), 1);
    ASSERT_NE(model->GetBondKeySystemPtr(), nullptr);
    EXPECT_NE(model->GetBondKeySystemPtr()->GetBondId(model->GetBondList().front()->GetBondKey()), "UNK");
}

TEST(DataObjectIOContractTest, UppercaseExtensionsDispatchCorrectly)
{
    const command_test::ScopedTempDir temp_dir{ "schema_uppercase_ext" };
    const auto source_model_path{ command_test::TestDataPath("test_model.cif") };
    const auto uppercase_model_path{ temp_dir.path() / "TEST_MODEL.MMCIF" };
    std::filesystem::copy_file(source_model_path, uppercase_model_path);

    rg::DataObjectManager manager{};
    ASSERT_NO_THROW(manager.ProcessFile(uppercase_model_path, "model"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("model")->GetNumberOfAtom(), 1);

    const auto map_object{ MakeTinyMapObject() };
    const auto uppercase_map_path{ temp_dir.path() / "TEST_MAP.MAP" };
    ASSERT_NO_THROW(rg::WriteMap(uppercase_map_path, map_object));
    auto loaded_map{ rg::ReadMap(uppercase_map_path) };
    ASSERT_NE(loaded_map, nullptr);
    EXPECT_EQ(loaded_map->GetGridSize(), map_object.GetGridSize());
}

TEST(DataObjectIOContractTest, FileFormatDescriptorsAreUniqueAndConsistent)
{
    std::unordered_set<std::string> extension_set;
    const auto registry{ rg::BuildDefaultFileFormatRegistry() };
    for (const auto & descriptor : registry.GetAllDescriptors())
    {
        EXPECT_TRUE(extension_set.insert(descriptor.extension).second);
        if (descriptor.kind == rg::DataObjectKind::Model)
        {
            EXPECT_TRUE(descriptor.model_backend.has_value());
            EXPECT_FALSE(descriptor.map_backend.has_value());
        }
        else
        {
            EXPECT_TRUE(descriptor.map_backend.has_value());
            EXPECT_FALSE(descriptor.model_backend.has_value());
        }

        if (descriptor.extension == ".mmcif" || descriptor.extension == ".mcif")
        {
            EXPECT_TRUE(descriptor.supports_read);
            EXPECT_FALSE(descriptor.supports_write);
        }
    }
}

TEST(DataObjectIOContractTest, PdbWriteRoundTripBasicFields)
{
    const command_test::ScopedTempDir temp_dir{ "pdb_write_roundtrip" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto output_path{ temp_dir.path() / "roundtrip.pdb" };

    rg::DataObjectManager manager{};
    manager.ProcessFile(model_path, "source");
    auto source_model{ manager.GetTypedDataObject<rg::ModelObject>("source") };
    ASSERT_NO_THROW(manager.ProduceFile(output_path, "source"));
    ASSERT_NO_THROW(manager.ProcessFile(output_path, "roundtrip"));

    auto roundtrip_model{ manager.GetTypedDataObject<rg::ModelObject>("roundtrip") };
    ASSERT_GT(roundtrip_model->GetNumberOfAtom(), 0);
    EXPECT_EQ(roundtrip_model->GetNumberOfAtom(), source_model->GetNumberOfAtom());
    EXPECT_EQ(
        roundtrip_model->GetAtomList().front()->GetSerialID(),
        source_model->GetAtomList().front()->GetSerialID());
    EXPECT_EQ(
        roundtrip_model->GetAtomList().front()->GetChainID(),
        source_model->GetAtomList().front()->GetChainID());
}

TEST(DataObjectIOContractTest, ModelWriteSupportMatrixAllowsPdbAndCifAndRejectsMmcif)
{
    const command_test::ScopedTempDir temp_dir{ "schema_output_matrix" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto pdb_output_path{ temp_dir.path() / "supported_output.pdb" };
    const auto cif_output_path{ temp_dir.path() / "supported_output.cif" };
    const auto output_path{ temp_dir.path() / "unsupported_output.mmcif" };

    rg::DataObjectManager manager{};
    manager.ProcessFile(model_path, "model");
    EXPECT_NO_THROW(manager.ProduceFile(pdb_output_path, "model"));
    EXPECT_NO_THROW(manager.ProduceFile(cif_output_path, "model"));
    EXPECT_THROW(manager.ProduceFile(output_path, "model"), std::runtime_error);
}

TEST(DataObjectIOContractTest, ProcessFileThrowsOnMalformedModelInput)
{
    const command_test::ScopedTempDir temp_dir{ "bad_model_file" };
    const auto malformed_path{ temp_dir.path() / "bad_model.cif" };
    {
        std::ofstream output{ malformed_path };
        output << "data_bad\nloop_\n_atom_site.id\n";
    }

    rg::DataObjectManager manager{};
    EXPECT_THROW(manager.ProcessFile(malformed_path, "broken"), std::runtime_error);
}

TEST(DataObjectIOContractTest, ProcessFileThrowsOnMalformedMapInput)
{
    const command_test::ScopedTempDir temp_dir{ "bad_map_file" };
    const auto malformed_path{ temp_dir.path() / "bad_map.map" };
    {
        std::ofstream output{ malformed_path, std::ios::binary };
        output << "bad";
    }

    rg::DataObjectManager manager{};
    EXPECT_THROW(manager.ProcessFile(malformed_path, "broken_map"), std::runtime_error);
}

TEST(DataObjectIOContractTest, ProduceFileThrowsWhenWriterCannotOpenTarget)
{
    const command_test::ScopedTempDir temp_dir{ "bad_output_target" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto output_path{ temp_dir.path() / "missing_dir" / "output.cif" };

    rg::DataObjectManager manager{};
    manager.ProcessFile(model_path, "model");

    EXPECT_THROW(manager.ProduceFile(output_path, "model"), std::runtime_error);
}

TEST(DataObjectIOContractTest, FunctionFileIoThrowsWhenReadFails)
{
    const command_test::ScopedTempDir temp_dir{ "file_io_failure_contract" };
    const auto missing_model_path{ temp_dir.path() / "missing_model.cif" };
    const auto malformed_map_path{ temp_dir.path() / "bad_map.map" };
    {
        std::ofstream map_output{ malformed_map_path, std::ios::binary };
        map_output << "bad";
    }

    EXPECT_THROW((void)rg::ReadModel(missing_model_path), std::runtime_error);
    EXPECT_THROW((void)rg::ReadMap(malformed_map_path), std::runtime_error);
}
