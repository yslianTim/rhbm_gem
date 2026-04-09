#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <tuple>

#include "support/CommandTestHelpers.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include <rhbm_gem/data/io/ModelMapFileIO.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/math/GausLinearTransformHelper.hpp>

namespace rg = rhbm_gem;

namespace {

std::shared_ptr<rg::ModelObject> LoadModelFixture(const std::string & fixture_name)
{
    auto model{ rg::ReadModel(command_test::TestDataPath(fixture_name)) };
    model->SetKeyTag("model");
    return std::shared_ptr<rg::ModelObject>{ std::move(model) };
}

} // namespace

TEST(LocalPotentialSeriesTest, EntryComputesRangeAndBinningForDistanceMapSeries)
{
    rg::LocalPotentialEntry entry;
    entry.SetDistanceAndMapValueList({
        {0.0f, 2.0f},
        {0.2f, 4.0f},
        {0.7f, 8.0f},
        {1.0f, -1.0f},
    });

    const auto distance_range{ entry.GetDistanceRange(0.0) };
    const auto map_value_range{ entry.GetMapValueRange(0.0) };
    const auto binned{ entry.GetBinnedDistanceAndMapValueList(2, 0.0, 1.0) };

    EXPECT_FLOAT_EQ(std::get<0>(distance_range), 0.0f);
    EXPECT_FLOAT_EQ(std::get<1>(distance_range), 1.0f);
    EXPECT_FLOAT_EQ(std::get<0>(map_value_range), -1.0f);
    EXPECT_FLOAT_EQ(std::get<1>(map_value_range), 8.0f);
    ASSERT_EQ(binned.size(), 2U);
    EXPECT_FLOAT_EQ(std::get<1>(binned.at(0)), 3.0f);
    EXPECT_FLOAT_EQ(std::get<1>(binned.at(1)), 8.0f);
}

TEST(LocalPotentialSeriesTest, EntryLinearModelTransformKeepsPositiveSamplesOnly)
{
    rg::LocalPotentialEntry entry;
    entry.SetDistanceAndMapValueList({
        {0.1f, 4.0f},
        {0.2f, -2.0f},
        {0.3f, 8.0f},
    });

    const auto transformed{ entry.GetLinearModelDistanceAndMapValueList() };
    ASSERT_EQ(transformed.size(), 2U);

    Eigen::VectorXd init{ Eigen::VectorXd::Zero(3) };
    const auto expected0{
        GausLinearTransformHelper::BuildLinearModelDataVector(0.1f, 4.0f, init) };
    const auto expected1{
        GausLinearTransformHelper::BuildLinearModelDataVector(0.3f, 8.0f, init) };

    EXPECT_NEAR(std::get<0>(transformed.at(0)), expected0(1), 1e-6);
    EXPECT_NEAR(std::get<1>(transformed.at(0)), expected0(2), 1e-6);
    EXPECT_NEAR(std::get<0>(transformed.at(1)), expected1(1), 1e-6);
    EXPECT_NEAR(std::get<1>(transformed.at(1)), expected1(2), 1e-6);
}

TEST(LocalPotentialSeriesTest, EntryBinningRespectsNonZeroMinimum)
{
    rg::LocalPotentialEntry entry;
    entry.SetDistanceAndMapValueList({
        {0.20f, 1.0f},
        {0.60f, 4.0f},
        {0.90f, 6.0f},
        {1.20f, 8.0f},
    });

    const auto binned{ entry.GetBinnedDistanceAndMapValueList(2, 0.5, 1.5) };

    ASSERT_EQ(binned.size(), 2U);
    EXPECT_FLOAT_EQ(std::get<0>(binned.at(0)), 0.75f);
    EXPECT_FLOAT_EQ(std::get<0>(binned.at(1)), 1.25f);
    EXPECT_FLOAT_EQ(std::get<1>(binned.at(0)), 5.0f);
    EXPECT_FLOAT_EQ(std::get<1>(binned.at(1)), 8.0f);
}

TEST(LocalPotentialSeriesTest, ViewForwardsSeriesDerivationsFromResolvedEntry)
{
    auto model{ LoadModelFixture("test_model_auth_seq_alnum_struct_conn.cif") };
    ASSERT_NE(model, nullptr);
    ASSERT_FALSE(model->GetAtomList().empty());

    auto * atom{ model->GetAtomList().front().get() };
    auto analysis{ model->EditAnalysis() };
    auto mutable_view{ analysis.EnsureAtomLocalPotential(*atom) };
    mutable_view.SetDistanceAndMapValueList({
        {0.1f, 2.0f},
        {0.4f, 4.0f},
        {0.8f, 6.0f},
    });

    const auto view{ rg::LocalPotentialView::RequireFor(*atom) };
    const auto map_value_range{ view.GetMapValueRange(0.0) };
    const auto binned{ view.GetBinnedDistanceAndMapValueList(2, 0.0, 1.0) };

    EXPECT_FLOAT_EQ(std::get<0>(map_value_range), 2.0f);
    EXPECT_FLOAT_EQ(std::get<1>(map_value_range), 6.0f);
    ASSERT_EQ(binned.size(), 2U);
    EXPECT_FLOAT_EQ(std::get<1>(binned.at(0)), 3.0f);
    EXPECT_FLOAT_EQ(std::get<1>(binned.at(1)), 6.0f);
}
