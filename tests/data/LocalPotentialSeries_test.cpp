#include <gtest/gtest.h>

#include <memory>
#include <string>
#include <tuple>
#include <utility>

#include "support/CommandTestHelpers.hpp"
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

SeriesPointList BuildExpectedDatasetSeries(
    const LocalPotentialSampleList & sampling_entries,
    double x_min,
    double x_max)
{
    return ls::BuildDatasetSeries(sampling_entries, x_min, x_max);
}

struct LocalPotentialViewFixture
{
    std::shared_ptr<rg::ModelObject> model;
    rg::AtomObject * atom;
};

LocalPotentialViewFixture MakeLocalPotentialViewFixture(
    LocalPotentialSampleList sampling_entries)
{
    auto model{ LoadModelFixture("test_model_auth_seq_alnum_struct_conn.cif") };
    auto * atom{ model->GetAtomList().front().get() };
    auto analysis{ model->EditAnalysis() };
    auto editor{ analysis.EnsureAtomLocalPotential(*atom) };
    editor.SetSamplingEntries(std::move(sampling_entries));
    return LocalPotentialViewFixture{ std::move(model), atom };
}

} // namespace

TEST(LocalPotentialSeriesTest, ViewComputesRangeAndBinningForDistanceMapSeries)
{
    auto fixture{ MakeLocalPotentialViewFixture({
        {2.0f, SamplingPoint{ 0.0f }},
        {4.0f, SamplingPoint{ 0.2f }},
        {8.0f, SamplingPoint{ 0.7f }},
        {-1.0f, SamplingPoint{ 1.0f }},
    }) };
    const auto view{ rg::AtomLocalPotentialView::RequireFor(*fixture.atom) };

    const auto distance_range{ view.GetDistanceRange(0.0) };
    const auto map_value_range{ view.GetResponseRange(0.0) };
    const auto binned{ view.GetBinnedDistanceResponseSeries(2, 0.0, 1.0) };

    EXPECT_FLOAT_EQ(std::get<0>(distance_range), 0.0f);
    EXPECT_FLOAT_EQ(std::get<1>(distance_range), 1.0f);
    EXPECT_FLOAT_EQ(std::get<0>(map_value_range), -1.0f);
    EXPECT_FLOAT_EQ(std::get<1>(map_value_range), 8.0f);
    ASSERT_EQ(binned.size(), 2U);
    EXPECT_FLOAT_EQ(binned.at(0).response, 3.0f);
    EXPECT_FLOAT_EQ(binned.at(1).response, 8.0f);
}

TEST(LocalPotentialSeriesTest, FitDatasetSeriesKeepsFullBasisWithinFitRange)
{
    const LocalPotentialSampleList sampling_entries{
        {4.0f, SamplingPoint{ 0.1f }},
        {-2.0f, SamplingPoint{ 0.2f }},
        {9.0f, SamplingPoint{ 0.25f }},
        {8.0f, SamplingPoint{ 0.3f }},
        {16.0f, SamplingPoint{ 1.1f }},
    };

    const auto transformed{ BuildExpectedDatasetSeries(sampling_entries, 0.1, 0.3) };
    ASSERT_EQ(transformed.size(), 3U);

    ASSERT_EQ(transformed.at(0).GetBasisSize(), 2U);
    ASSERT_EQ(transformed.at(1).GetBasisSize(), 2U);
    ASSERT_EQ(transformed.at(2).GetBasisSize(), 2U);
    EXPECT_NEAR(transformed.at(0).GetBasisValue(0), 1.0, 1e-6);
    EXPECT_NEAR(transformed.at(0).GetBasisValue(1), -0.5 * 0.1 * 0.1, 1e-6);
    EXPECT_NEAR(transformed.at(0).response, std::log(4.0), 1e-6);
    EXPECT_NEAR(transformed.at(1).GetBasisValue(0), 1.0, 1e-6);
    EXPECT_NEAR(transformed.at(1).GetBasisValue(1), -0.5 * 0.25 * 0.25, 1e-6);
    EXPECT_NEAR(transformed.at(1).response, std::log(9.0), 1e-6);
    EXPECT_NEAR(transformed.at(2).GetBasisValue(0), 1.0, 1e-6);
    EXPECT_NEAR(transformed.at(2).GetBasisValue(1), -0.5 * 0.3 * 0.3, 1e-6);
    EXPECT_NEAR(transformed.at(2).response, std::log(8.0), 1e-6);
}

TEST(LocalPotentialSeriesTest, ViewBinningRespectsNonZeroMinimum)
{
    auto fixture{ MakeLocalPotentialViewFixture({
        {1.0f, SamplingPoint{ 0.20f }},
        {4.0f, SamplingPoint{ 0.60f }},
        {6.0f, SamplingPoint{ 0.90f }},
        {8.0f, SamplingPoint{ 1.20f }},
    }) };
    const auto view{ rg::AtomLocalPotentialView::RequireFor(*fixture.atom) };

    const auto binned{ view.GetBinnedDistanceResponseSeries(2, 0.5, 1.5) };

    ASSERT_EQ(binned.size(), 2U);
    EXPECT_FLOAT_EQ(binned.at(0).GetBasisValue(0), 0.75f);
    EXPECT_FLOAT_EQ(binned.at(1).GetBasisValue(0), 1.25f);
    EXPECT_FLOAT_EQ(binned.at(0).response, 5.0f);
    EXPECT_FLOAT_EQ(binned.at(1).response, 8.0f);
}

TEST(LocalPotentialSeriesTest, ViewForwardsSeriesDerivationsFromResolvedEntry)
{
    auto model{ LoadModelFixture("test_model_auth_seq_alnum_struct_conn.cif") };
    ASSERT_NE(model, nullptr);
    ASSERT_FALSE(model->GetAtomList().empty());

    auto * atom{ model->GetAtomList().front().get() };
    auto analysis{ model->EditAnalysis() };
    auto editor{ analysis.EnsureAtomLocalPotential(*atom) };
    editor.SetSamplingEntries({
        {2.0f, SamplingPoint{ 0.1f }},
        {4.0f, SamplingPoint{ 0.4f }},
        {6.0f, SamplingPoint{ 0.8f }},
    });

    const auto view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
    const auto map_value_range{ view.GetResponseRange(0.0) };
    const auto binned{ view.GetBinnedDistanceResponseSeries(2, 0.0, 1.0) };

    EXPECT_FLOAT_EQ(std::get<0>(map_value_range), 2.0f);
    EXPECT_FLOAT_EQ(std::get<1>(map_value_range), 6.0f);
    ASSERT_EQ(binned.size(), 2U);
    EXPECT_FLOAT_EQ(binned.at(0).response, 3.0f);
    EXPECT_FLOAT_EQ(binned.at(1).response, 6.0f);

    const auto fit_dataset_series{
        ls::BuildDatasetSeries(view.GetSamplingEntries(), 0.0, 0.5)
    };
    ASSERT_EQ(fit_dataset_series.size(), 2U);
    EXPECT_EQ(fit_dataset_series.at(0).GetBasisSize(), 2U);
}

TEST(LocalPotentialSeriesTest, ViewBinningReturnsZeroResponseForEmptyBins)
{
    auto fixture{ MakeLocalPotentialViewFixture({
        {2.0f, SamplingPoint{ 0.10f }},
    }) };
    const auto view{ rg::AtomLocalPotentialView::RequireFor(*fixture.atom) };

    const auto binned{ view.GetBinnedDistanceResponseSeries(2, 0.0, 1.0) };

    ASSERT_EQ(binned.size(), 2U);
    EXPECT_FLOAT_EQ(binned.at(0).response, 2.0f);
    EXPECT_FLOAT_EQ(binned.at(1).response, 0.0f);
}
