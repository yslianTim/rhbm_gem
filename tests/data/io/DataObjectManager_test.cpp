#include <gtest/gtest.h>
#include <unordered_map>

#include <rhbm_gem/data/io/DataObjectManager.hpp>
#include "CommandTestHelpers.hpp"
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rg = rhbm_gem;
using DataObjectManager = rg::DataObjectManager;
using ModelObject = rg::ModelObject;

TEST(DataObjectManagerTest, LoadCifFile) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{command_test::TestDataPath("test_model.cif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    EXPECT_EQ(model->GetNumberOfAtom(), 1);
}

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

TEST(DataObjectManagerTest, LoadCifFileMissingModelNumberDefaultsToOne) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{
        command_test::TestDataPath("test_model_missing_model_num.cif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_EQ(model->GetAtomList().front()->GetSerialID(), 1);
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

TEST(DataObjectManagerTest, LoadCifFileModelTwoFallback) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{
        command_test::TestDataPath("test_model_model2_only.cif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_EQ(model->GetAtomList().front()->GetSerialID(), 1);
}

TEST(DataObjectManagerTest, LoadCifFileWithDoubleQuotedAtomId) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{
        command_test::TestDataPath("test_model_atom_id_double_quote.cif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_EQ(model->GetAtomList().front()->GetAtomID(), "CA A");
}

TEST(DataObjectManagerTest, LoadMmcifExtension) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{command_test::TestDataPath("test_model.mmcif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    EXPECT_EQ(model->GetNumberOfAtom(), 1);
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

TEST(DataObjectManagerTest, LoadCifFileLoopMultilineAndQuotedToken) {
    DataObjectManager manager{command_test::BuildDataIoServices()};
    auto cif_path{
        command_test::TestDataPath("test_model_loop_multiline.cif")};
    manager.ProcessFile(cif_path, "model");
    auto model{manager.GetTypedDataObject<ModelObject>("model")};
    ASSERT_EQ(model->GetNumberOfAtom(), 1);
    EXPECT_EQ(model->GetAtomList().front()->GetAtomID(), "CA A");
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
