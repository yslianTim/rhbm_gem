#include <gtest/gtest.h>
#include <memory>
#include <string>
#include <unordered_map>

#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include "CommandTestHelpers.hpp"
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rg = rhbm_gem;
using DataObjectManager = rg::DataObjectManager;
using ModelObject = rg::ModelObject;

namespace {

std::shared_ptr<ModelObject> LoadModelFixture(
    DataObjectManager& manager,
    const std::string& fixture_name) {
    manager.ProcessFile(command_test::TestDataPath(fixture_name), "model");
    return manager.GetTypedDataObject<ModelObject>("model");
}

struct FixtureAtomCountCase {
    const char* name;
    const char* fixture_name;
    size_t expected_atom_count;
};

class DataObjectManagerAtomCountFixtureTest
    : public ::testing::TestWithParam<FixtureAtomCountCase> {};

TEST_P(DataObjectManagerAtomCountFixtureTest, LoadsFixtureWithExpectedAtomCount) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto model{LoadModelFixture(manager, GetParam().fixture_name)};
    EXPECT_EQ(model->GetNumberOfAtom(), GetParam().expected_atom_count);
}

INSTANTIATE_TEST_SUITE_P(
    AtomCountFixtures,
    DataObjectManagerAtomCountFixtureTest,
    ::testing::Values(
        FixtureAtomCountCase{ "CifExtension", "test_model.cif", 1 },
        FixtureAtomCountCase{ "MmcifExtension", "test_model.mmcif", 1 }),
    [](const ::testing::TestParamInfo<FixtureAtomCountCase>& info) {
        return info.param.name;
    });

struct FixtureSerialIdCase {
    const char* name;
    const char* fixture_name;
    size_t expected_serial_id;
};

class DataObjectManagerSerialIdFixtureTest
    : public ::testing::TestWithParam<FixtureSerialIdCase> {};

TEST_P(DataObjectManagerSerialIdFixtureTest, UsesExpectedSerialIdFallback) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto model{LoadModelFixture(manager, GetParam().fixture_name)};
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_EQ(model->GetAtomList().front()->GetSerialID(), GetParam().expected_serial_id);
}

INSTANTIATE_TEST_SUITE_P(
    SerialIdFixtures,
    DataObjectManagerSerialIdFixtureTest,
    ::testing::Values(
        FixtureSerialIdCase{ "MissingModelNumberDefaultsToOne", "test_model_missing_model_num.cif", 1 },
        FixtureSerialIdCase{ "ModelTwoFallback", "test_model_model2_only.cif", 1 }),
    [](const ::testing::TestParamInfo<FixtureSerialIdCase>& info) {
        return info.param.name;
    });

struct FixtureAtomIdCase {
    const char* name;
    const char* fixture_name;
    const char* expected_atom_id;
};

class DataObjectManagerAtomIdFixtureTest
    : public ::testing::TestWithParam<FixtureAtomIdCase> {};

TEST_P(DataObjectManagerAtomIdFixtureTest, KeepsExpectedAtomIdentifier) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto model{LoadModelFixture(manager, GetParam().fixture_name)};
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_EQ(model->GetAtomList().front()->GetAtomID(), GetParam().expected_atom_id);
}

INSTANTIATE_TEST_SUITE_P(
    AtomIdFixtures,
    DataObjectManagerAtomIdFixtureTest,
    ::testing::Values(
        FixtureAtomIdCase{ "DoubleQuotedAtomId", "test_model_atom_id_double_quote.cif", "CA A" },
        FixtureAtomIdCase{ "LoopMultilineQuotedToken", "test_model_loop_multiline.cif", "CA A" }),
    [](const ::testing::TestParamInfo<FixtureAtomIdCase>& info) {
        return info.param.name;
    });

} // namespace

TEST(DataObjectManagerTest, LoadCifFileWithKeyValueEntityMetadata) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{
        command_test::TestDataPath("test_model_keyvalue_entity.cif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    ASSERT_EQ(model->GetNumberOfAtom(), 1);

    for (const auto& atom : model->GetAtomList())
        atom->SetSelectedFlag(true);
    model->FilterAtomFromSymmetry(false);
    model->Update();
    EXPECT_EQ(model->GetNumberOfSelectedAtom(), 1);
}

TEST(DataObjectManagerTest, FilterAtomFromSymmetrySkipsWhenChainMapMissing) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{command_test::TestDataPath("test_model.cif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    ASSERT_EQ(model->GetNumberOfAtom(), 1);

    for (const auto& atom : model->GetAtomList())
        atom->SetSelectedFlag(true);
    model->SetChainIDListMap(std::unordered_map<std::string, std::vector<std::string>>{});
    model->FilterAtomFromSymmetry(false);
    model->Update();
    EXPECT_EQ(model->GetNumberOfSelectedAtom(), 1);
}

TEST(DataObjectManagerTest, LoadCifFileAuthOnlyAtomSiteColumns) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{
        command_test::TestDataPath("test_model_auth_only.cif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_EQ(model->GetAtomList().front()->GetComponentID(), "ALA");
    EXPECT_EQ(model->GetAtomList().front()->GetChainID(), "A");
    EXPECT_EQ(model->GetAtomList().front()->GetSequenceID(), 1);
}

TEST(DataObjectManagerTest, LoadCifFileMissingNumericUsesDefaults) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{
        command_test::TestDataPath("test_model_missing_numeric.cif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_FLOAT_EQ(model->GetAtomList().front()->GetOccupancy(), 1.0f);
    EXPECT_FLOAT_EQ(model->GetAtomList().front()->GetTemperature(), 0.0f);
}

TEST(DataObjectManagerTest, LoadCifFileDatabaseOrderKeepsEmdb) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{
        command_test::TestDataPath("test_model_database_order.cif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    EXPECT_EQ(model->GetEmdID(), "EMD-1234");
}

TEST(DataObjectManagerTest, LoadCifFileAltBOnlyKeepsAltIndicator) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{
        command_test::TestDataPath("test_model_alt_b_only.cif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_EQ(model->GetAtomList().front()->GetIndicator(), "B");
}

TEST(DataObjectManagerTest, LoadCifFileInvalidSecondaryRangeDoesNotDropAtoms) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{
        command_test::TestDataPath("test_model_invalid_secondary_range.cif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
}

TEST(DataObjectManagerTest, LoadCifFileAuthSeqAlnumStructConnBuildsBond) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{
        command_test::TestDataPath("test_model_auth_seq_alnum_struct_conn.cif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    ASSERT_EQ(model->GetNumberOfAtom(), 2);
    EXPECT_GE(model->GetNumberOfBond(), 1);
}
