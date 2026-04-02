#include <gtest/gtest.h>

#include <filesystem>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/data/object/ModelObject.hpp>
#include "data/detail/ModelObjectAccess.hpp"
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
                auto & bond_key_system{ rg::ModelObjectAccess::BondKeySystemRef(model) };
                EXPECT_NE(
                    bond_key_system.GetBondId(model.GetBondList().front()->GetBondKey()),
                    "UNK");
            } },
    };

    for (const auto & case_data : cases)
    {
        SCOPED_TRACE(case_data.name);
        auto model{ data_test::LoadFixtureModel(case_data.path) };
        case_data.verify(*model);
    }
}
