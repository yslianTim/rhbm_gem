#include <gtest/gtest.h>

#include <filesystem>
#include <functional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include "support/CommandTestHelpers.hpp"
#include "support/DataObjectTestSupport.hpp"

namespace rg = rhbm_gem;

namespace {

struct ImportRegressionCase
{
    const char * name;
    std::filesystem::path path;
    std::function<void(rg::ModelObject &)> verify;
};

} // namespace

TEST(DataObjectImportRegressionTest, CifEdgeCaseMatrix)
{
    const command_test::ScopedTempDir temp_dir{ "data_runtime_mmcif_fixture" };
    const auto mmcif_path{
        data_test::CopyFixtureWithNewName(
            command_test::TestDataPath("test_model.cif"),
            temp_dir.path() / "test_model_runtime.mmcif") };

    const std::vector<ImportRegressionCase> cases{
        { "CifExtensionLoadsOneAtom",
            command_test::TestDataPath("test_model.cif"),
            [](rg::ModelObject & model)
            {
                EXPECT_EQ(model.GetNumberOfAtom(), 1);
            } },
        { "MmcifExtensionLoadsOneAtom",
            mmcif_path,
            [](rg::ModelObject & model)
            {
                EXPECT_EQ(model.GetNumberOfAtom(), 1);
            } },
        { "MissingModelNumberDefaultsSerialId",
            command_test::TestDataPath("test_model_missing_model_num.cif"),
            [](rg::ModelObject & model)
            {
                ASSERT_EQ(model.GetNumberOfAtom(), 1);
                EXPECT_EQ(model.GetAtomList().front()->GetSerialID(), 1);
            } },
        { "ModelTwoUsesSerialFallback",
            command_test::TestDataPath("test_model_model2_only.cif"),
            [](rg::ModelObject & model)
            {
                ASSERT_EQ(model.GetNumberOfAtom(), 1);
                EXPECT_EQ(model.GetAtomList().front()->GetSerialID(), 1);
            } },
        { "DoubleQuotedAtomIdPreserved",
            command_test::TestDataPath("test_model_atom_id_double_quote.cif"),
            [](rg::ModelObject & model)
            {
                ASSERT_EQ(model.GetNumberOfAtom(), 1);
                EXPECT_EQ(model.GetAtomList().front()->GetAtomID(), "CA A");
            } },
        { "LoopMultilineQuotedTokenPreserved",
            command_test::TestDataPath("test_model_loop_multiline.cif"),
            [](rg::ModelObject & model)
            {
                ASSERT_EQ(model.GetNumberOfAtom(), 1);
                EXPECT_EQ(model.GetAtomList().front()->GetAtomID(), "CA A");
            } },
        { "AuthOnlyColumnsPopulateFallbackFields",
            command_test::TestDataPath("test_model_auth_only.cif"),
            [](rg::ModelObject & model)
            {
                ASSERT_EQ(model.GetNumberOfAtom(), 1);
                EXPECT_EQ(model.GetAtomList().front()->GetComponentID(), "ALA");
                EXPECT_EQ(model.GetAtomList().front()->GetChainID(), "A");
                EXPECT_EQ(model.GetAtomList().front()->GetSequenceID(), 1);
            } },
        { "MissingNumericFieldsUseDefaults",
            command_test::TestDataPath("test_model_missing_numeric.cif"),
            [](rg::ModelObject & model)
            {
                ASSERT_EQ(model.GetNumberOfAtom(), 1);
                EXPECT_FLOAT_EQ(model.GetAtomList().front()->GetOccupancy(), 1.0f);
                EXPECT_FLOAT_EQ(model.GetAtomList().front()->GetTemperature(), 0.0f);
            } },
        { "DatabaseOrderKeepsEmdbIdentifier",
            command_test::TestDataPath("test_model_database_order.cif"),
            [](rg::ModelObject & model)
            {
                EXPECT_EQ(model.GetEmdID(), "EMD-1234");
            } },
        { "AltBIndicatorPreserved",
            command_test::TestDataPath("test_model_alt_b_only.cif"),
            [](rg::ModelObject & model)
            {
                ASSERT_EQ(model.GetNumberOfAtom(), 1);
                EXPECT_EQ(model.GetAtomList().front()->GetIndicator(), "B");
            } },
        { "InvalidSecondaryRangeDoesNotDropAtoms",
            command_test::TestDataPath("test_model_invalid_secondary_range.cif"),
            [](rg::ModelObject & model)
            {
                EXPECT_EQ(model.GetNumberOfAtom(), 1);
            } },
        { "StructConnImportBuildsBondAndBondKeySystem",
            command_test::TestDataPath("test_model_auth_seq_alnum_struct_conn.cif"),
            [](rg::ModelObject & model)
            {
                ASSERT_EQ(model.GetNumberOfAtom(), 2);
                ASSERT_GE(model.GetNumberOfBond(), 1);
                EXPECT_NE(
                    model.FindBondID(model.GetBondList().front()->GetBondKey()),
                    "UNK");
            } },
        { "ModifiedPeptideFallbackBuildsObservedComponentBonds",
            command_test::TestDataPath("test_model_modified_peptide_fallback.cif"),
            [](rg::ModelObject & model)
            {
                ASSERT_EQ(model.GetNumberOfAtom(), 6);
                std::unordered_set<std::string> bond_id_set;
                for (const auto & bond : model.GetBondList())
                {
                    bond_id_set.insert(model.FindBondID(bond->GetBondKey()));
                }
                EXPECT_TRUE(bond_id_set.find("N_CA") != bond_id_set.end());
                EXPECT_TRUE(bond_id_set.find("CA_C") != bond_id_set.end());
                EXPECT_TRUE(bond_id_set.find("C_O") != bond_id_set.end());
                EXPECT_TRUE(bond_id_set.find("C_N") != bond_id_set.end());
                EXPECT_GE(model.GetNumberOfBond(), 5);
            } },
    };

    for (const auto & case_data : cases)
    {
        SCOPED_TRACE(case_data.name);
        auto model{ data_test::LoadFixtureModel(case_data.path) };
        case_data.verify(*model);
    }
}

TEST(DataObjectImportRegressionTest, ModifiedPeptideFallbackDoesNotEmitComponentBondMissWarnings)
{
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    auto model{
        data_test::LoadFixtureModel(
            command_test::TestDataPath("test_model_modified_peptide_fallback.cif")) };
    const auto stdout_output{ testing::internal::GetCapturedStdout() };
    const auto stderr_output{ testing::internal::GetCapturedStderr() };
    ASSERT_NE(model, nullptr);
    EXPECT_TRUE(stdout_output.empty());
    EXPECT_EQ(
        stderr_output.find("Component bond entry (bond_key: 62) not found in chemical component map."),
        std::string::npos);
    EXPECT_EQ(stderr_output.find("Component bond entry"), std::string::npos);
    Logger::SetLogLevel(LogLevel::Info);
}
