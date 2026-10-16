#include <gtest/gtest.h>
#include "support/JointHaloObservationExclusion.hpp"
#include <cmath>

using namespace joint_test::halo_exclusion;

TEST(JointHaloObservationExclusionTest, RadiusZeroRetainsEveryObservationAndNoHaloIsUnmasked)
{
    const std::array<Position, 2> observations{{{0,0,0},{1,0,0}}};
    const auto distances = NearestHaloDistances(observations, std::span<const Position>{});
    const auto mask = RetainedRows(distances, 0.0);
    EXPECT_EQ(mask, (RowMask{1,1}));
    EXPECT_TRUE(std::isinf(distances[0]));
}

TEST(JointHaloObservationExclusionTest, NearestDistanceUsesAllHalosAndStrictRadiusBoundary)
{
    const std::array<Position, 3> observations{{{0,0,0},{1,0,0},{2,0,0}}};
    const std::array<Position, 2> halos{{{0.25,0,0},{1.5,0,0}}};
    const auto distances = NearestHaloDistances(observations, halos);
    EXPECT_DOUBLE_EQ(distances[0], 0.25);
    EXPECT_DOUBLE_EQ(distances[1], 0.5);
    EXPECT_DOUBLE_EQ(distances[2], 0.5);
    EXPECT_EQ(RetainedRows(distances, 0.5), (RowMask{0,1,1}));
}

TEST(JointHaloObservationExclusionTest, FarHalosLeaveRowsAndNearHaloCentersAreExcluded)
{
    const std::array<Position, 2> observations{{{0,0,0},{1,0,0}}};
    const std::array<Position, 1> far_halo{{{10,0,0}}};
    EXPECT_EQ(RetainedRows(NearestHaloDistances(observations, far_halo), 2.5), (RowMask{1,1}));
    const std::array<Position, 1> near_halo{{{0,0,0}}};
    EXPECT_EQ(RetainedRows(NearestHaloDistances(observations, near_halo), 0.5), (RowMask{0,1}));
}

TEST(JointHaloObservationExclusionTest, CoverageCountsSharedRowsMixedRowsRadialShellsAndFragmentation)
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids = {"a", "b", "halo"};
    input.row_ids = {"0", "1", "2", "3"};
    input.observations.assign(4, 0.0);
    input.support = {{{0,0.04},{1,0.36}}, {{1,0.64},{2,1.44}}, {{1,0.64},{3,2.25}}};
    input.selection_domain.emplace();
    input.selection_domain->target_indices = {0,1};
    const auto summary = AnalyzeCoverage(input, RowMask{1,0,1,1}, std::array<double,2>{0.1,0.1},
        std::array<double,4>{0.1,0.1,0.1,0.1});
    EXPECT_EQ(summary.original_rows, 4u);
    EXPECT_EQ(summary.retained_rows, 3u);
    EXPECT_EQ(summary.mixed_target_halo_rows, 1u);
    EXPECT_EQ(summary.excluded_mixed_target_halo_rows, 1u);
    EXPECT_EQ(summary.rows_shared_by_targets, 1u);
    EXPECT_EQ(summary.retained_rows_shared_by_targets, 0u);
    EXPECT_EQ(summary.target_components, 2u);
    EXPECT_EQ(summary.targets[0].original_rows, 2u);
    EXPECT_EQ(summary.targets[0].retained_rows, 1u);
    EXPECT_EQ(summary.radial[0].original, 1u);
    EXPECT_EQ(summary.radial[1].original, 2u);
}

TEST(JointHaloObservationExclusionTest, EmptyRetainedDomainIsReportedWithoutInventingCoverage)
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids = {"a", "halo"};
    input.row_ids = {"0"};
    input.observations = {0.0};
    input.support = {{{0,0.0}}, {{0,0.0}}};
    input.selection_domain.emplace();
    input.selection_domain->target_indices = {0};
    const auto summary = AnalyzeCoverage(input, RowMask{0}, std::array<double,1>{0.0}, std::array<double,1>{0.0});
    EXPECT_EQ(summary.targets_without_retained_rows, 1u);
    EXPECT_EQ(summary.target_components, 1u);
    EXPECT_EQ(summary.radial[0].retained, 0u);
}

TEST(JointHaloObservationExclusionTest, LocalDesignReportsRankDeficiencyAndProjectedWidthLoss)
{
    Eigen::MatrixXd ac(4, 4), width(4, 2);
    ac << 1,0,1,0,
          0,1,0,1,
          1,1,1,1,
          2,1,2,1;
    width << 1,1,
             0,0,
             1,1,
             2,2;
    const auto result = AnalyzeLocalDesign(ac, width);
    ASSERT_TRUE(result.ac.available);
    ASSERT_TRUE(result.abc.available);
    ASSERT_TRUE(result.projected_width.available);
    EXPECT_EQ(result.ac.rank, 2u);
    EXPECT_LT(result.ac.rank, result.ac.columns);
    EXPECT_LT(result.abc.rank, result.abc.columns);
    EXPECT_LT(result.projected_width.rank, result.projected_width.columns);
    EXPECT_NEAR(result.minimum_projected_width_fraction, 0.0, 1e-14);
}
