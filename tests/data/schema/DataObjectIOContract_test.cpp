#include <gtest/gtest.h>

#include "data/schema/DataObjectIOSchemaTestSupport.hpp"

using namespace data_object_io_schema_test;

TEST(DataObjectIOContractTest, CatalogLookupLoadsCorrectDaoType)
{
    const command_test::ScopedTempDir temp_dir{ "schema_catalog_lookup" };
    const auto database_path{ temp_dir.path() / "catalog_lookup.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto map_path{ temp_dir.path() / "map_only.map" };

    rg::DataObjectManager manager;
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    manager.SaveDataObject("model");
    SaveTinyMapThroughManager(manager, map_path, "map");
    manager.ClearDataObjects();

    ASSERT_NO_THROW(manager.LoadDataObject("model"));
    ASSERT_NO_THROW(manager.LoadDataObject("map"));
    EXPECT_NE(dynamic_cast<rg::ModelObject *>(manager.GetDataObject("model").get()), nullptr);
    EXPECT_NE(dynamic_cast<rg::MapObject *>(manager.GetDataObject("map").get()), nullptr);
}

TEST(DataObjectIOContractTest, ChainMetadataPersistsAcrossDatabaseRoundTrip)
{
    const command_test::ScopedTempDir temp_dir{ "schema_chain_roundtrip" };
    const auto database_path{ temp_dir.path() / "chain_roundtrip.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model_keyvalue_entity.cif") };

    rg::DataObjectManager manager;
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

    rg::DataObjectManager manager;
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
    rg::DataObjectManager manager;
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

    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.ProcessFile(uppercase_model_path, "model"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("model")->GetNumberOfAtom(), 1);

    const auto map_object{ MakeTinyMapObject() };
    const auto uppercase_map_path{ temp_dir.path() / "TEST_MAP.MAP" };
    rg::MapFileWriter writer{ uppercase_map_path.string(), &map_object };
    ASSERT_NO_THROW(writer.Write());
    rg::MapFileReader reader{ uppercase_map_path.string() };
    ASSERT_NO_THROW(reader.Read());
    EXPECT_EQ(reader.GetGridSizeArray(), map_object.GetGridSize());
}

TEST(DataObjectIOContractTest, BuiltInFormatsStillWorkWithoutAnyRegistrationStep)
{
    rg::DefaultFileProcessFactoryResolver resolver;
    auto factory{ resolver.CreateFactory(".cif") };
    auto data_object{ factory->CreateDataObject(command_test::TestDataPath("test_model.cif").string()) };
    EXPECT_NE(dynamic_cast<rg::ModelObject *>(data_object.get()), nullptr);
}

TEST(DataObjectIOContractTest, DataObjectManagerDoesNotRequireDefaultFactoryRegistration)
{
    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.ProcessFile(command_test::TestDataPath("test_model.cif"), "model"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("model")->GetNumberOfAtom(), 1);
}

TEST(DataObjectIOContractTest, OverrideResolverTakesPrecedenceOverDescriptorFallback)
{
    auto resolver{ std::make_shared<rg::OverrideableFileProcessFactoryResolver>() };
    resolver->RegisterFactory(".cif", []() { return std::make_unique<OverrideFileFactory>(); });

    auto override_factory{ resolver->CreateFactory(".cif") };
    auto override_object{
        override_factory->CreateDataObject(command_test::TestDataPath("test_model.cif").string()) };
    EXPECT_NE(dynamic_cast<FailingDataObject *>(override_object.get()), nullptr);
}

TEST(DataObjectIOContractTest, ResolverInjectionSupportsOverrideWithoutGlobalRegistry)
{
    rg::OverrideableFileProcessFactoryResolver resolver;
    resolver.RegisterFactory(".cif", []() { return std::make_unique<OverrideFileFactory>(); });

    auto override_factory{ resolver.CreateFactory(".cif") };
    auto override_object{
        override_factory->CreateDataObject(command_test::TestDataPath("test_model.cif").string()) };
    EXPECT_NE(dynamic_cast<FailingDataObject *>(override_object.get()), nullptr);

    resolver.UnregisterFactory(".cif");

    auto default_factory{ resolver.CreateFactory(".cif") };
    auto default_object{
        default_factory->CreateDataObject(command_test::TestDataPath("test_model.cif").string()) };
    EXPECT_NE(dynamic_cast<rg::ModelObject *>(default_object.get()), nullptr);
}

TEST(DataObjectIOContractTest, CustomFactoryOverrideTakesPrecedenceOverDescriptorFallback)
{
    rg::OverrideableFileProcessFactoryResolver resolver;
    ScopedFactoryOverride override{
        resolver,
        ".cif",
        []() { return std::make_unique<OverrideFileFactory>(); }
    };

    auto factory{ resolver.CreateFactory(".cif") };
    auto data_object{
        factory->CreateDataObject(command_test::TestDataPath("test_model.cif").string()) };
    EXPECT_NE(dynamic_cast<FailingDataObject *>(data_object.get()), nullptr);
}

TEST(DataObjectIOContractTest, CustomFactoryOverrideDoesNotLeakAcrossTests)
{
    rg::OverrideableFileProcessFactoryResolver resolver;
    {
        ScopedFactoryOverride override{
            resolver,
            ".cif",
            []() { return std::make_unique<OverrideFileFactory>(); }
        };
        auto factory{ resolver.CreateFactory(".cif") };
        auto data_object{
            factory->CreateDataObject(command_test::TestDataPath("test_model.cif").string()) };
        EXPECT_NE(dynamic_cast<FailingDataObject *>(data_object.get()), nullptr);
    }

    auto factory{ resolver.CreateFactory(".cif") };
    auto data_object{ factory->CreateDataObject(command_test::TestDataPath("test_model.cif").string()) };
    EXPECT_NE(dynamic_cast<rg::ModelObject *>(data_object.get()), nullptr);
}

TEST(DataObjectIOContractTest, FileFormatDescriptorsAreUniqueAndConsistent)
{
    std::unordered_set<std::string> extension_set;
    for (const auto & descriptor : rg::FileFormatRegistry::Instance().GetAllDescriptors())
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

TEST(DataObjectIOContractTest, FileFormatRegistryIndexMatchesDescriptorSource)
{
    const auto & registry{ rg::FileFormatRegistry::Instance() };
    for (const auto & descriptor : registry.GetAllDescriptors())
    {
        const auto & looked_up_lower{ registry.Lookup(descriptor.extension) };
        EXPECT_EQ(looked_up_lower.extension, descriptor.extension);
        EXPECT_EQ(looked_up_lower.kind, descriptor.kind);
        EXPECT_EQ(looked_up_lower.supports_read, descriptor.supports_read);
        EXPECT_EQ(looked_up_lower.supports_write, descriptor.supports_write);
        EXPECT_EQ(looked_up_lower.model_backend, descriptor.model_backend);
        EXPECT_EQ(looked_up_lower.map_backend, descriptor.map_backend);

        auto upper_extension{ descriptor.extension };
        for (auto & ch : upper_extension)
        {
            ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        }
        const auto & looked_up_upper{ registry.Lookup(upper_extension) };
        EXPECT_EQ(looked_up_upper.extension, descriptor.extension);
    }
}

TEST(DataObjectIOContractTest, ResolverConcurrentOverrideAndLookupIsSafe)
{
    rg::OverrideableFileProcessFactoryResolver resolver;
    constexpr int kIterationCount{ 800 };
    std::atomic<bool> writer_done{ false };
    std::atomic<int> failure_count{ 0 };

    std::thread writer{
        [&]()
        {
            for (int i = 0; i < kIterationCount; i++)
            {
                resolver.RegisterFactory(".cif", []() { return std::make_unique<OverrideFileFactory>(); });
                resolver.UnregisterFactory(".cif");
            }
            writer_done = true;
        }
    };

    auto reader_task{
        [&]()
        {
            while (!writer_done.load())
            {
                try
                {
                    auto factory{ resolver.CreateFactory(".cif") };
                    if (!factory)
                    {
                        ++failure_count;
                    }
                }
                catch (const std::exception &)
                {
                    ++failure_count;
                }
            }
        }
    };

    std::thread reader_a{ reader_task };
    std::thread reader_b{ reader_task };

    writer.join();
    reader_a.join();
    reader_b.join();
    EXPECT_EQ(failure_count.load(), 0);
}

TEST(DataObjectIOContractTest, PdbWriteProducesNonEmptyOutput)
{
    const command_test::ScopedTempDir temp_dir{ "pdb_write_non_empty" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto output_path{ temp_dir.path() / "output.pdb" };

    rg::DataObjectManager manager;
    manager.ProcessFile(model_path, "model");
    ASSERT_NO_THROW(manager.ProduceFile(output_path, "model"));

    std::ifstream output_stream{ output_path };
    ASSERT_TRUE(output_stream.is_open());
    const std::string output_content{
        std::istreambuf_iterator<char>(output_stream),
        std::istreambuf_iterator<char>()
    };
    EXPECT_FALSE(output_content.empty());
    EXPECT_NE(output_content.find("ATOM"), std::string::npos);
}

TEST(DataObjectIOContractTest, PdbWriteRoundTripBasicFields)
{
    const command_test::ScopedTempDir temp_dir{ "pdb_write_roundtrip" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto output_path{ temp_dir.path() / "roundtrip.pdb" };

    rg::DataObjectManager manager;
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

TEST(DataObjectIOContractTest, ModelWriteSupportMatrixStillAcceptsPdbAndCifOnly)
{
    const command_test::ScopedTempDir temp_dir{ "schema_output_matrix" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto pdb_output_path{ temp_dir.path() / "supported_output.pdb" };
    const auto cif_output_path{ temp_dir.path() / "supported_output.cif" };
    const auto output_path{ temp_dir.path() / "unsupported_output.mmcif" };

    rg::DataObjectManager manager;
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

    rg::DataObjectManager manager;
    try
    {
        manager.ProcessFile(malformed_path, "broken");
        FAIL() << "Expected malformed model input to throw.";
    }
    catch (const std::runtime_error & ex)
    {
        const std::string message{ ex.what() };
        EXPECT_NE(message.find(malformed_path.string()), std::string::npos);
        EXPECT_NE(message.find("broken"), std::string::npos);
    }
}

TEST(DataObjectIOContractTest, ProcessFileThrowsOnMalformedMapInput)
{
    const command_test::ScopedTempDir temp_dir{ "bad_map_file" };
    const auto malformed_path{ temp_dir.path() / "bad_map.map" };
    {
        std::ofstream output{ malformed_path, std::ios::binary };
        output << "bad";
    }

    rg::DataObjectManager manager;
    try
    {
        manager.ProcessFile(malformed_path, "broken_map");
        FAIL() << "Expected malformed map input to throw.";
    }
    catch (const std::runtime_error & ex)
    {
        const std::string message{ ex.what() };
        EXPECT_NE(message.find(malformed_path.string()), std::string::npos);
        EXPECT_NE(message.find("broken_map"), std::string::npos);
    }
}

TEST(DataObjectIOContractTest, ProduceFileThrowsWhenWriterCannotOpenTarget)
{
    const command_test::ScopedTempDir temp_dir{ "bad_output_target" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto output_path{ temp_dir.path() / "missing_dir" / "output.cif" };

    rg::DataObjectManager manager;
    manager.ProcessFile(model_path, "model");

    try
    {
        manager.ProduceFile(output_path, "model");
        FAIL() << "Expected writer open failure to throw.";
    }
    catch (const std::runtime_error & ex)
    {
        const std::string message{ ex.what() };
        EXPECT_NE(message.find(output_path.string()), std::string::npos);
        EXPECT_NE(message.find("model"), std::string::npos);
    }
}

TEST(DataObjectIOContractTest, ReaderGettersThrowIfReadDidNotSucceed)
{
    const command_test::ScopedTempDir temp_dir{ "reader_failure_contract" };
    const auto missing_model_path{ temp_dir.path() / "missing_model.cif" };
    const auto malformed_map_path{ temp_dir.path() / "bad_map.map" };
    {
        std::ofstream map_output{ malformed_map_path, std::ios::binary };
        map_output << "bad";
    }

    rg::ModelFileReader model_reader{ missing_model_path.string() };
    EXPECT_THROW(model_reader.Read(), std::runtime_error);
    EXPECT_THROW(model_reader.GetDataBlockPtr(), std::runtime_error);

    rg::MapFileReader map_reader{ malformed_map_path.string() };
    EXPECT_THROW(map_reader.Read(), std::runtime_error);
    EXPECT_THROW(map_reader.GetGridSizeArray(), std::runtime_error);
    EXPECT_THROW(map_reader.GetGridSpacingArray(), std::runtime_error);
    EXPECT_THROW(map_reader.GetOriginArray(), std::runtime_error);
    EXPECT_THROW((void)map_reader.GetMapValueArray(), std::runtime_error);
}
