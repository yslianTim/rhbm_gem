#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "CommandTestHelpers.hpp"
#include "BondObject.hpp"
#include "DataObjectBase.hpp"
#include "DataObjectDAOBase.hpp"
#include "DataObjectDAOFactoryRegistry.hpp"
#include "DataObjectManager.hpp"
#include "DatabaseManager.hpp"
#include "LegacyModelObjectDAO.hpp"
#include "LocalPotentialEntry.hpp"
#include "MapFileReader.hpp"
#include "MapFileWriter.hpp"
#include "MapObject.hpp"
#include "MapObjectDAO.hpp"
#include "ModelObject.hpp"
#include "SQLiteWrapper.hpp"

namespace rg = rhbm_gem;

namespace {

constexpr const char * kUpsertObjectMetadataSql =
    "INSERT INTO object_metadata(key_tag, object_type) VALUES (?, ?) "
    "ON CONFLICT(key_tag) DO UPDATE SET object_type = excluded.object_type;";

void UpsertObjectMetadata(
    rg::SQLiteWrapper & database,
    const std::string & key_tag,
    const std::string & object_type)
{
    database.Execute(
        "CREATE TABLE IF NOT EXISTS object_metadata (key_tag TEXT PRIMARY KEY, object_type TEXT);");
    database.Prepare(kUpsertObjectMetadataSql);
    rg::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, key_tag);
    database.Bind<std::string>(2, object_type);
    database.StepOnce();
}

int GetUserVersion(const std::filesystem::path & database_path)
{
    rg::SQLiteWrapper database{ database_path };
    database.Prepare("PRAGMA user_version;");
    rg::SQLiteWrapper::StatementGuard guard(database);
    const auto rc{ database.StepNext() };
    if (rc != rg::SQLiteWrapper::StepRow()) return 0;
    return database.GetColumn<int>(0);
}

bool HasTable(const std::filesystem::path & database_path, const std::string & table_name)
{
    rg::SQLiteWrapper database{ database_path };
    database.Prepare(
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name = ? LIMIT 1;");
    rg::SQLiteWrapper::StatementGuard guard(database);
    database.Bind<std::string>(1, table_name);
    return database.StepNext() == rg::SQLiteWrapper::StepRow();
}

int CountRows(
    const std::filesystem::path & database_path,
    const std::string & table_name,
    const std::string & key_tag = "")
{
    rg::SQLiteWrapper database{ database_path };
    const auto sql{
        key_tag.empty()
            ? "SELECT COUNT(*) FROM " + table_name + ";"
            : "SELECT COUNT(*) FROM " + table_name + " WHERE key_tag = ?;"
    };
    database.Prepare(sql);
    rg::SQLiteWrapper::StatementGuard guard(database);
    if (!key_tag.empty())
    {
        database.Bind<std::string>(1, key_tag);
    }
    const auto rc{ database.StepNext() };
    if (rc != rg::SQLiteWrapper::StepRow())
    {
        throw std::runtime_error("Failed to count rows in " + table_name);
    }
    return database.GetColumn<int>(0);
}

std::shared_ptr<rg::ModelObject> LoadFixtureModel(
    const std::filesystem::path & model_path,
    const std::string & key_tag = "model")
{
    rg::DataObjectManager manager;
    manager.ProcessFile(model_path, key_tag);
    return manager.GetTypedDataObject<rg::ModelObject>(key_tag);
}

void SeedLegacyModel(
    const std::filesystem::path & database_path,
    const std::filesystem::path & model_path,
    const std::string & key_tag,
    const std::string & pdb_id)
{
    auto model{ LoadFixtureModel(model_path, key_tag) };
    model->SetPdbID(pdb_id);
    if (!model->GetAtomList().empty())
    {
        auto entry{ std::make_unique<rg::LocalPotentialEntry>() };
        entry->AddDistanceAndMapValueList({ std::make_tuple(0.0f, 1.0f) });
        entry->AddGausEstimateOLS(1.0, 2.0);
        entry->AddGausEstimateMDPDE(1.5, 2.5);
        entry->SetAlphaR(0.25);
        model->GetAtomList().front()->SetSelectedFlag(true);
        model->GetAtomList().front()->AddLocalPotentialEntry(std::move(entry));
        model->Update();
    }

    rg::SQLiteWrapper database{ database_path };
    rg::LegacyModelObjectDAO legacy_dao{ &database };
    legacy_dao.Save(model.get(), key_tag);
    UpsertObjectMetadata(database, key_tag, "model");
}

class FailingDataObject final : public rg::DataObjectBase
{
    std::string m_key_tag;

public:
    std::unique_ptr<rg::DataObjectBase> Clone() const override
    {
        return std::make_unique<FailingDataObject>(*this);
    }

    void Display() const override {}
    void Update() override {}
    void Accept(rg::DataObjectVisitorBase * /*visitor*/) override {}
    void SetKeyTag(const std::string & label) override { m_key_tag = label; }
    std::string GetKeyTag() const override { return m_key_tag; }
};

class FailingDataObjectDAO final : public rg::DataObjectDAOBase
{
    rg::SQLiteWrapper * m_database;

public:
    explicit FailingDataObjectDAO(rg::SQLiteWrapper * database) : m_database{ database } {}

    void Save(const rg::DataObjectBase * /*obj*/, const std::string & key_tag) override
    {
        m_database->Execute("CREATE TABLE IF NOT EXISTS failing_payload (key_tag TEXT PRIMARY KEY);");
        m_database->Prepare("INSERT INTO failing_payload(key_tag) VALUES (?);");
        rg::SQLiteWrapper::StatementGuard guard(*m_database);
        m_database->Bind<std::string>(1, key_tag);
        m_database->StepOnce();
        throw std::runtime_error("Injected DAO failure");
    }

    std::unique_ptr<rg::DataObjectBase> Load(const std::string & /*key_tag*/) override
    {
        throw std::runtime_error("Load not implemented for failing DAO");
    }
};

rg::DataObjectDAORegistrar<FailingDataObject, FailingDataObjectDAO>
    registrar_failing_dao("failing");

} // namespace

TEST(DataObjectIOSchemaTest, EmptyDatabaseBootstrapsNormalizedSchema)
{
    const command_test::ScopedTempDir temp_dir{ "schema_bootstrap" };
    const auto database_path{ temp_dir.path() / "bootstrap.sqlite" };

    rg::DatabaseManager database_manager{ database_path };

    EXPECT_EQ(database_manager.GetSchemaVersion(), rg::DatabaseSchemaVersion::NormalizedV2);
    EXPECT_TRUE(HasTable(database_path, "object_metadata"));
    EXPECT_TRUE(HasTable(database_path, "model_object"));
    EXPECT_TRUE(HasTable(database_path, "map_list"));
    EXPECT_EQ(GetUserVersion(database_path), 2);
}

TEST(DataObjectIOSchemaTest, SaveRollbackLeavesNoPayloadOrMetadata)
{
    const command_test::ScopedTempDir temp_dir{ "schema_atomic_save" };
    const auto database_path{ temp_dir.path() / "atomic.sqlite" };

    rg::DatabaseManager database_manager{ database_path };
    FailingDataObject data_object;
    data_object.SetKeyTag("failing");

    EXPECT_THROW(database_manager.SaveDataObject(&data_object, "failing"), std::runtime_error);
    EXPECT_FALSE(HasTable(database_path, "failing_payload"));
    EXPECT_EQ(CountRows(database_path, "object_metadata"), 0);
}

TEST(DataObjectIOSchemaTest, LegacyModelDatabaseMigratesToNormalizedSchema)
{
    const command_test::ScopedTempDir temp_dir{ "schema_migrate_legacy" };
    const auto database_path{ temp_dir.path() / "legacy.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model_keyvalue_entity.cif") };

    SeedLegacyModel(database_path, model_path, "legacy_model", "LEGACY");
    ASSERT_EQ(GetUserVersion(database_path), 0);
    ASSERT_TRUE(HasTable(database_path, "model_list"));

    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    EXPECT_EQ(
        manager.GetDatabaseManager()->GetSchemaVersion(),
        rg::DatabaseSchemaVersion::NormalizedV2);

    ASSERT_NO_THROW(manager.LoadDataObject("legacy_model"));
    auto model{ manager.GetTypedDataObject<rg::ModelObject>("legacy_model") };
    EXPECT_EQ(model->GetPdbID(), "LEGACY");
    EXPECT_TRUE(HasTable(database_path, "model_object"));
    EXPECT_FALSE(HasTable(database_path, "model_list"));
    EXPECT_EQ(GetUserVersion(database_path), 2);
}

TEST(DataObjectIOSchemaTest, FailedLegacyMigrationRollsBackToLegacyState)
{
    const command_test::ScopedTempDir temp_dir{ "schema_migrate_rollback" };
    const auto database_path{ temp_dir.path() / "legacy_rollback.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    SeedLegacyModel(database_path, model_path, "key_a", "MODEL_A");
    SeedLegacyModel(database_path, model_path, "key_b", "MODEL_B");

    {
        rg::SQLiteWrapper database{ database_path };
        database.Prepare("UPDATE model_list SET atom_size = ? WHERE key_tag = ?;");
        rg::SQLiteWrapper::StatementGuard guard(database);
        database.Bind<int>(1, 999);
        database.Bind<std::string>(2, "key_b");
        database.StepOnce();
    }

    rg::DataObjectManager manager;
    EXPECT_THROW(manager.SetDatabaseManager(database_path), std::runtime_error);
    EXPECT_TRUE(HasTable(database_path, "model_list"));
    EXPECT_FALSE(HasTable(database_path, "model_object"));
    EXPECT_EQ(GetUserVersion(database_path), 0);
}

TEST(DataObjectIOSchemaTest, MapOnlyLegacyDatabaseRemainsLoadable)
{
    const command_test::ScopedTempDir temp_dir{ "schema_map_only" };
    const auto database_path{ temp_dir.path() / "map_only.sqlite" };

    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; i++) values[i] = static_cast<float>(i);
    rg::MapObject map_object{ grid_size, grid_spacing, origin, std::move(values) };

    {
        rg::SQLiteWrapper database{ database_path };
        rg::MapObjectDAO map_dao{ &database };
        map_dao.Save(&map_object, "map_only");
        UpsertObjectMetadata(database, "map_only", "map");
    }

    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.SetDatabaseManager(database_path));
    EXPECT_EQ(manager.GetDatabaseManager()->GetSchemaVersion(), rg::DatabaseSchemaVersion::NormalizedV2);
    ASSERT_NO_THROW(manager.LoadDataObject("map_only"));
    auto loaded_map{ manager.GetTypedDataObject<rg::MapObject>("map_only") };
    EXPECT_EQ(loaded_map->GetGridSize(), grid_size);
    EXPECT_EQ(GetUserVersion(database_path), 2);
}

TEST(DataObjectIOSchemaTest, DistinctUnsanitizedKeysDoNotCollideInV2Schema)
{
    const command_test::ScopedTempDir temp_dir{ "schema_key_collision" };
    const auto database_path{ temp_dir.path() / "collision.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataObjectManager manager;
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "mem_a");
    manager.ProcessFile(model_path, "mem_b");
    manager.GetTypedDataObject<rg::ModelObject>("mem_a")->SetPdbID("MODEL_A");
    manager.GetTypedDataObject<rg::ModelObject>("mem_b")->SetPdbID("MODEL_B");

    manager.SaveDataObject("mem_a", "a-b");
    manager.SaveDataObject("mem_b", "a_b");
    manager.ClearDataObjects();

    manager.LoadDataObject("a-b");
    manager.LoadDataObject("a_b");
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("a-b")->GetPdbID(), "MODEL_A");
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("a_b")->GetPdbID(), "MODEL_B");
    EXPECT_EQ(CountRows(database_path, "model_object"), 2);
}

TEST(DataObjectIOSchemaTest, ChainMetadataPersistsAcrossDatabaseRoundTrip)
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

TEST(DataObjectIOSchemaTest, SymmetryFilteringMatchesAfterDatabaseReload)
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

TEST(DataObjectIOSchemaTest, FileImportTransfersBondKeySystem)
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

TEST(DataObjectIOSchemaTest, UppercaseExtensionsDispatchCorrectly)
{
    const command_test::ScopedTempDir temp_dir{ "schema_uppercase_ext" };
    const auto source_model_path{ command_test::TestDataPath("test_model.cif") };
    const auto uppercase_model_path{ temp_dir.path() / "TEST_MODEL.MMCIF" };
    std::filesystem::copy_file(source_model_path, uppercase_model_path);

    rg::DataObjectManager manager;
    ASSERT_NO_THROW(manager.ProcessFile(uppercase_model_path, "model"));
    EXPECT_EQ(manager.GetTypedDataObject<rg::ModelObject>("model")->GetNumberOfAtom(), 1);

    std::array<int, 3> grid_size{ 2, 2, 2 };
    std::array<float, 3> grid_spacing{ 1.0f, 1.0f, 1.0f };
    std::array<float, 3> origin{ 0.0f, 0.0f, 0.0f };
    auto values{ std::make_unique<float[]>(8) };
    for (size_t i = 0; i < 8; i++) values[i] = static_cast<float>(i);
    rg::MapObject map_object{ grid_size, grid_spacing, origin, std::move(values) };

    const auto uppercase_map_path{ temp_dir.path() / "TEST_MAP.MAP" };
    rg::MapFileWriter writer{ uppercase_map_path.string(), &map_object };
    ASSERT_NO_THROW(writer.Write());
    rg::MapFileReader reader{ uppercase_map_path.string() };
    ASSERT_NO_THROW(reader.Read());
    EXPECT_EQ(reader.GetGridSizeArray(), grid_size);
}

TEST(DataObjectIOSchemaTest, ModelWriteSupportMatrixRejectsMmcifOutput)
{
    const command_test::ScopedTempDir temp_dir{ "schema_output_matrix" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };
    const auto output_path{ temp_dir.path() / "unsupported_output.mmcif" };

    rg::DataObjectManager manager;
    manager.ProcessFile(model_path, "model");
    EXPECT_THROW(manager.ProduceFile(output_path, "model"), std::runtime_error);
}

TEST(DataObjectIOSchemaTest, SaveRenamedKeyDoesNotRenameInMemoryObject)
{
    const command_test::ScopedTempDir temp_dir{ "schema_rename_semantics" };
    const auto database_path{ temp_dir.path() / "rename.sqlite" };
    const auto model_path{ command_test::TestDataPath("test_model.cif") };

    rg::DataObjectManager manager;
    manager.SetDatabaseManager(database_path);
    manager.ProcessFile(model_path, "model");
    manager.SaveDataObject("model", "saved_model");

    EXPECT_TRUE(manager.HasDataObject("model"));
    EXPECT_FALSE(manager.HasDataObject("saved_model"));

    ASSERT_NO_THROW(manager.LoadDataObject("saved_model"));
    EXPECT_TRUE(manager.HasDataObject("saved_model"));
}
