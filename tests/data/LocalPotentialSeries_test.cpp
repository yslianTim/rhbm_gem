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
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>

namespace rg = rhbm_gem;
namespace ls = rhbm_gem::linearization_service;

namespace {

std::shared_ptr<rg::ModelObject> LoadModelFixture(const std::string & fixture_name)
{
    auto model{ rg::ReadModel(command_test::TestDataPath(fixture_name)) };
    model->SetKeyTag("model");
    return std::shared_ptr<rg::ModelObject>{ std::move(model) };
}

const ls::LinearizationSpec & DatasetLinearizationSpec()
{
    static const auto spec{ ls::LinearizationSpec::DefaultDataset() };
    return spec;
}

SeriesPointList BuildExpectedDatasetSeries(
    const LocalPotentialSampleList & sampling_entries,
    double x_min,
    double x_max)
{
    return ls::BuildDatasetSeries(DatasetLinearizationSpec(), sampling_entries, x_min, x_max);
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
    EXPECT_FLOAT_EQ(binned.at(0).score, 2.0f);
    EXPECT_FLOAT_EQ(binned.at(1).score, 2.0f);
}

TEST(LocalPotentialSeriesTest, EntryFitDatasetSeriesKeepsFullBasisWithinFitRange)
{
    rg::LocalPotentialEntry entry;
    const LocalPotentialSampleList sampling_entries{
        {0.1f, 4.0f, 0.5f},
        {0.2f, -2.0f, 7.0f},
        {0.25f, 9.0f, 0.0f},
        {0.3f, 8.0f, 2.5f},
        {1.1f, 16.0f, 3.0f},
    };
    entry.SetSamplingEntries(sampling_entries);

    const auto transformed{ BuildExpectedDatasetSeries(sampling_entries, 0.1, 0.3) };
    ASSERT_EQ(transformed.size(), 2U);

    ASSERT_EQ(transformed.at(0).GetBasisSize(), 2U);
    ASSERT_EQ(transformed.at(1).GetBasisSize(), 2U);
    EXPECT_NEAR(transformed.at(0).GetBasisValue(0), 1.0, 1e-6);
    EXPECT_NEAR(transformed.at(0).GetBasisValue(1), -0.5 * 0.1 * 0.1, 1e-6);
    EXPECT_NEAR(transformed.at(0).response, std::log(4.0), 1e-6);
    EXPECT_NEAR(transformed.at(1).GetBasisValue(0), 1.0, 1e-6);
    EXPECT_NEAR(transformed.at(1).GetBasisValue(1), -0.5 * 0.3 * 0.3, 1e-6);
    EXPECT_NEAR(transformed.at(1).response, std::log(8.0), 1e-6);
    EXPECT_FLOAT_EQ(transformed.at(0).score, 0.5f);
    EXPECT_FLOAT_EQ(transformed.at(1).score, 2.5f);
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
    EXPECT_FLOAT_EQ(binned.at(0).score, 3.0f);
    EXPECT_FLOAT_EQ(binned.at(1).score, 6.0f);
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
    EXPECT_FLOAT_EQ(binned.at(0).score, 1.5f);
    EXPECT_FLOAT_EQ(binned.at(1).score, 5.0f);

    const auto fit_dataset_series{
        ls::BuildDatasetSeries(DatasetLinearizationSpec(), view.GetSamplingEntries(), 0.0, 0.5)
    };
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
    EXPECT_FLOAT_EQ(binned.at(0).score, 4.0f);
    EXPECT_FLOAT_EQ(binned.at(1).score, 0.0f);
}
