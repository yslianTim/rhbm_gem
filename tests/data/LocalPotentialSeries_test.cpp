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
    entry.SetSamplingEntries({
        {0.0f, 2.0f, 1.0f},
        {0.2f, 4.0f, 3.0f},
        {0.7f, 8.0f, 2.0f},
        {1.0f, -1.0f, 9.0f},
    });

    const auto distance_range{ entry.GetDistanceRange(0.0) };
    const auto map_value_range{ entry.GetResponseRange(0.0) };
    const auto binned{ entry.GetBinnedDistanceResponseSeries(2, 0.0, 1.0) };

    EXPECT_FLOAT_EQ(std::get<0>(distance_range), 0.0f);
    EXPECT_FLOAT_EQ(std::get<1>(distance_range), 1.0f);
    EXPECT_FLOAT_EQ(std::get<0>(map_value_range), -1.0f);
    EXPECT_FLOAT_EQ(std::get<1>(map_value_range), 8.0f);
    ASSERT_EQ(binned.size(), 2U);
    EXPECT_FLOAT_EQ(binned.at(0).response, 3.0f);
    EXPECT_FLOAT_EQ(binned.at(1).response, 8.0f);
    EXPECT_FLOAT_EQ(binned.at(0).weight, 2.0f);
    EXPECT_FLOAT_EQ(binned.at(1).weight, 2.0f);
}

TEST(LocalPotentialSeriesTest, EntryLinearModelTransformKeepsPositiveSamplesOnly)
{
    rg::LocalPotentialEntry entry;
    entry.SetSamplingEntries({
        {0.1f, 4.0f, 0.5f},
        {0.2f, -2.0f, 7.0f},
        {0.3f, 8.0f, 2.5f},
    });

    const auto transformed{ entry.GetLinearModelSeries() };
    ASSERT_EQ(transformed.size(), 2U);

    Eigen::VectorXd init{ Eigen::VectorXd::Zero(3) };
    const auto expected0{
        GausLinearTransformHelper::BuildLinearModelDataVector(0.1f, 4.0f, init) };
    const auto expected1{
        GausLinearTransformHelper::BuildLinearModelDataVector(0.3f, 8.0f, init) };

    EXPECT_NEAR(transformed.at(0).GetBasisValue(0), expected0(1), 1e-6);
    EXPECT_NEAR(transformed.at(0).response, expected0(2), 1e-6);
    EXPECT_NEAR(transformed.at(1).GetBasisValue(0), expected1(1), 1e-6);
    EXPECT_NEAR(transformed.at(1).response, expected1(2), 1e-6);
    EXPECT_FLOAT_EQ(transformed.at(0).weight, 0.5f);
    EXPECT_FLOAT_EQ(transformed.at(1).weight, 2.5f);
}

TEST(LocalPotentialSeriesTest, EntryFitDatasetSeriesKeepsFullBasisWithinFitRange)
{
    rg::LocalPotentialEntry entry;
    entry.SetSamplingEntries({
        {0.1f, 4.0f, 0.5f},
        {0.2f, -2.0f, 7.0f},
        {0.3f, 8.0f, 2.5f},
        {1.1f, 16.0f, 3.0f},
    });

    const auto transformed{ entry.GetFitDatasetSeries(0.1, 0.3) };
    ASSERT_EQ(transformed.size(), 2U);

    Eigen::VectorXd init{ Eigen::VectorXd::Zero(3) };
    const auto expected0{
        GausLinearTransformHelper::BuildLinearModelDataVector(0.1f, 4.0f, init) };
    const auto expected1{
        GausLinearTransformHelper::BuildLinearModelDataVector(0.3f, 8.0f, init) };

    ASSERT_EQ(transformed.at(0).GetBasisSize(), 2U);
    ASSERT_EQ(transformed.at(1).GetBasisSize(), 2U);
    EXPECT_NEAR(transformed.at(0).GetBasisValue(0), expected0(0), 1e-6);
    EXPECT_NEAR(transformed.at(0).GetBasisValue(1), expected0(1), 1e-6);
    EXPECT_NEAR(transformed.at(0).response, expected0(2), 1e-6);
    EXPECT_NEAR(transformed.at(1).GetBasisValue(0), expected1(0), 1e-6);
    EXPECT_NEAR(transformed.at(1).GetBasisValue(1), expected1(1), 1e-6);
    EXPECT_NEAR(transformed.at(1).response, expected1(2), 1e-6);
    EXPECT_FLOAT_EQ(transformed.at(0).weight, 0.5f);
    EXPECT_FLOAT_EQ(transformed.at(1).weight, 2.5f);
}

TEST(LocalPotentialSeriesTest, EntryBinningRespectsNonZeroMinimum)
{
    rg::LocalPotentialEntry entry;
    entry.SetSamplingEntries({
        {0.20f, 1.0f, 8.0f},
        {0.60f, 4.0f, 2.0f},
        {0.90f, 6.0f, 4.0f},
        {1.20f, 8.0f, 6.0f},
    });

    const auto binned{ entry.GetBinnedDistanceResponseSeries(2, 0.5, 1.5) };

    ASSERT_EQ(binned.size(), 2U);
    EXPECT_FLOAT_EQ(binned.at(0).GetBasisValue(0), 0.75f);
    EXPECT_FLOAT_EQ(binned.at(1).GetBasisValue(0), 1.25f);
    EXPECT_FLOAT_EQ(binned.at(0).response, 5.0f);
    EXPECT_FLOAT_EQ(binned.at(1).response, 8.0f);
    EXPECT_FLOAT_EQ(binned.at(0).weight, 3.0f);
    EXPECT_FLOAT_EQ(binned.at(1).weight, 6.0f);
}

TEST(LocalPotentialSeriesTest, ViewForwardsSeriesDerivationsFromResolvedEntry)
{
    auto model{ LoadModelFixture("test_model_auth_seq_alnum_struct_conn.cif") };
    ASSERT_NE(model, nullptr);
    ASSERT_FALSE(model->GetAtomList().empty());

    auto * atom{ model->GetAtomList().front().get() };
    auto analysis{ model->EditAnalysis() };
    auto mutable_view{ analysis.EnsureAtomLocalPotential(*atom) };
    mutable_view.SetSamplingEntries({
        {0.1f, 2.0f, 1.0f},
        {0.4f, 4.0f, 2.0f},
        {0.8f, 6.0f, 5.0f},
    });

    const auto view{ rg::LocalPotentialView::RequireFor(*atom) };
    const auto map_value_range{ view.GetResponseRange(0.0) };
    const auto binned{ view.GetBinnedDistanceResponseSeries(2, 0.0, 1.0) };

    EXPECT_FLOAT_EQ(std::get<0>(map_value_range), 2.0f);
    EXPECT_FLOAT_EQ(std::get<1>(map_value_range), 6.0f);
    ASSERT_EQ(binned.size(), 2U);
    EXPECT_FLOAT_EQ(binned.at(0).response, 3.0f);
    EXPECT_FLOAT_EQ(binned.at(1).response, 6.0f);
    EXPECT_FLOAT_EQ(binned.at(0).weight, 1.5f);
    EXPECT_FLOAT_EQ(binned.at(1).weight, 5.0f);

    const auto fit_dataset_series{ view.GetFitDatasetSeries(0.0, 0.5) };
    ASSERT_EQ(fit_dataset_series.size(), 2U);
    EXPECT_EQ(fit_dataset_series.at(0).GetBasisSize(), 2U);
}

TEST(LocalPotentialSeriesTest, EntryBinningReturnsZeroWeightForEmptyBins)
{
    rg::LocalPotentialEntry entry;
    entry.SetSamplingEntries({
        {0.10f, 2.0f, 4.0f},
    });

    const auto binned{ entry.GetBinnedDistanceResponseSeries(2, 0.0, 1.0) };

    ASSERT_EQ(binned.size(), 2U);
    EXPECT_FLOAT_EQ(binned.at(0).weight, 4.0f);
    EXPECT_FLOAT_EQ(binned.at(1).weight, 0.0f);
}
