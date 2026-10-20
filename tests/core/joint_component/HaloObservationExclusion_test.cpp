#include <gtest/gtest.h>
#include "support/JointHaloObservationExclusion.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include <algorithm>
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

TEST(JointHaloObservationExclusionTest, MultipleHalosAroundATargetUseTheNearestHalo)
{
    const std::array<Position, 3> target_observations{{{0,0,0},{0,1.5,0},{0,3,0}}};
    const std::array<Position, 2> surrounding_halos{{{-0.5,0,0},{0.5,0,0}}};
    const auto nearest = NearestHaloDistances(target_observations, surrounding_halos);
    EXPECT_DOUBLE_EQ(nearest[0], 0.5);
    EXPECT_NEAR(nearest[1], std::sqrt(2.5), 1e-14);
    EXPECT_DOUBLE_EQ(nearest[2], std::sqrt(9.25));
    EXPECT_EQ(RetainedRows(nearest, 1.0), (RowMask{0,1,1}));
}

TEST(JointHaloObservationExclusionTest, PositiveRadiusCanExcludeEveryObservation)
{
    EXPECT_EQ(RetainedRows(std::array<double, 4>{0.0, 0.25, 1.0, 2.5}, 3.0), (RowMask{0,0,0,0}));
}

TEST(JointHaloObservationExclusionTest, ExactSupportMaskHasNoEffectWithoutSupportedHalos)
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids = {"t0", "unsupported-halo"};
    input.row_ids = {"10", "11"};
    input.observations = {1.0, 2.0};
    input.support = {{{0,0.0},{1,6.25}}, {}};
    input.selection_domain.emplace();
    input.selection_domain->target_indices = {0};
    EXPECT_EQ(HaloSupportExclusionRows(input), (RowMask{1,1}));
}

TEST(JointHaloObservationExclusionTest, ExactSupportMaskExcludesEveryRowSupportedByOneHalo)
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids = {"target", "halo"};
    input.row_ids = {"10", "11", "12"};
    input.observations = {1.0, 2.0, 3.0};
    input.support = {{{0,0.0},{1,0.25},{2,0.5}}, {{1,0.25}}};
    input.selection_domain.emplace();
    input.selection_domain->target_indices = {0};
    EXPECT_EQ(HaloSupportExclusionRows(input), (RowMask{1,0,1}));
}

TEST(JointHaloObservationExclusionTest, ExactSupportMaskExcludesSingleAndMultipleHaloIncidencesAtBoundary)
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids = {"t0", "t1", "h0", "h1"};
    input.row_ids = {"10", "11", "12", "13"};
    input.observations = {1.0, 2.0, 3.0, 4.0};
    input.support = {{{0,0.0},{1,0.25}}, {{1,0.25},{2,0.5},{3,1.0}},
        {{1,6.25},{3,6.25}}, {{2,0.0}}};
    input.selection_domain.emplace();
    input.selection_domain->target_indices = {0,1};

    const auto retained = HaloSupportExclusionRows(input);
    EXPECT_EQ(retained, (RowMask{1,0,0,0}));
    for (std::size_t row = 0; row < retained.size(); ++row)
        if (retained[row])
            for (std::size_t atom = 0; atom < input.support.size(); ++atom)
                if (!std::binary_search(input.selection_domain->target_indices.begin(),
                        input.selection_domain->target_indices.end(), atom))
                    EXPECT_FALSE(std::any_of(input.support[atom].begin(), input.support[atom].end(),
                        [&](const auto & edge) { return edge.row == row; }));
}

TEST(JointHaloObservationExclusionTest, ExactSupportMaskCanFragmentSharedTargetRows)
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids = {"t0", "t1", "halo"};
    input.row_ids = {"10", "11", "12"};
    input.observations = {1.0, 2.0, 3.0};
    input.support = {{{0,0.0},{1,0.25}}, {{1,0.25},{2,0.5}}, {{1,6.25}}};
    input.selection_domain.emplace();
    input.selection_domain->target_indices = {0,1};

    const auto retained = HaloSupportExclusionRows(input);
    const auto summary = AnalyzeCoverage(input, retained, std::array<double,2>{3.0,3.0},
        std::array<double,3>{3.0,3.0,3.0});
    EXPECT_EQ(retained, (RowMask{1,0,1}));
    EXPECT_EQ(summary.target_components, 2u);
    EXPECT_EQ(summary.targets[0].retained_rows, 1u);
    EXPECT_EQ(summary.targets[1].retained_rows, 1u);
}

TEST(JointHaloObservationExclusionTest, GaussianNoiseIsPairedAndReproducibleForAFixedSeed)
{
    const std::vector<double> clean(10000, 0.0);
    std::vector<double> first(clean.size()), repeated(clean.size()), other_seed(clean.size());
    const auto summary = FillGaussianNoise(clean, first, 0.05, 17);
    EXPECT_EQ(FillGaussianNoise(clean, repeated, 0.05, 17).standard_deviation, summary.standard_deviation);
    EXPECT_EQ(first, repeated);
    FillGaussianNoise(clean, other_seed, 0.05, 19);
    EXPECT_NE(first, other_seed);
    EXPECT_NEAR(summary.mean, 0.0, 0.001);
    EXPECT_NEAR(summary.standard_deviation, 0.05, 0.001);
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

TEST(JointHaloObservationExclusionTest, MaskedInputsKeepOriginalRowsAndDropZeroSupportHalos)
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids = {"t0", "h0", "t1", "h1"};
    input.row_ids = {"10", "11", "12", "13"};
    input.observations = {1, 2, 3, 4};
    input.support = {{{0,0.1},{1,0.2}}, {{0,0.3},{1,0.4}},
        {{1,0.5},{2,0.6}}, {{3,0.7}}};
    input.selection_domain.emplace();
    input.selection_domain->target_indices = {0,2};

    const auto targets_only = BuildMaskedInput(input, RowMask{1,0,1,0});
    EXPECT_EQ(targets_only.input.row_ids, (std::vector<std::string>{"10", "12"}));
    EXPECT_EQ(targets_only.source_rows, (std::vector<std::size_t>{0, 2}));
    EXPECT_EQ(targets_only.input.observations, (std::vector<double>{1, 3}));
    EXPECT_EQ(targets_only.input.atom_ids, (std::vector<std::string>{"t0", "t1"}));
    EXPECT_EQ(targets_only.input.support[0].size(), 1u);
    EXPECT_EQ(targets_only.input.support[0][0].row, 0u);
    EXPECT_EQ(targets_only.input.support[1].size(), 1u);
    EXPECT_EQ(targets_only.input.support[1][0].row, 1u);
    EXPECT_EQ(targets_only.source_contributors, (std::vector<std::size_t>{0, 2}));

}

TEST(JointHaloObservationExclusionTest, EmptyDomainAndTargetWithoutRowsAreRejected)
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids = {"t0", "t1", "halo"};
    input.row_ids = {"10", "11"};
    input.observations = {1.0, 2.0};
    input.support = {{{0,0.0}}, {{1,0.25}}, {{0,0.25},{1,0.25}}};
    input.selection_domain.emplace();
    input.selection_domain->target_indices = {0,1};

    const auto all_halo_supported = HaloSupportExclusionRows(input);
    EXPECT_EQ(all_halo_supported, (RowMask{0,0}));
    EXPECT_THROW(BuildMaskedInput(input, all_halo_supported), std::invalid_argument);
    EXPECT_THROW(BuildMaskedInput(input, RowMask{1,0}), std::invalid_argument);
}

TEST(JointHaloObservationExclusionTest, MaskedFitRetainsOriginalObservationNormalization)
{
    rhbm_gem::core::JointProblemInput input;
    input.atom_ids = {"t0", "h0"};
    input.row_ids = {"0", "1", "2"};
    input.observations = {2, 3, 4};
    input.support = {{{0,0.1},{1,0.2},{2,0.3}}, {{0,0.4}}};
    input.selection_domain.emplace();
    input.selection_domain->target_indices = {0};
    const rhbm_gem::core::JointProblem reference(input);
    const auto masked_input = BuildMaskedInput(input, RowMask{1,0,0});
    const rhbm_gem::core::JointProblem masked(masked_input.input);
    const auto normalized = rhbm_gem::core::JointProblemAccess::WithReferenceNormalization(masked, reference);
    EXPECT_DOUBLE_EQ(normalized.ObservationScale(), reference.ObservationScale());
    EXPECT_DOUBLE_EQ(rhbm_gem::core::JointProblemAccess::Get(normalized).context.linear.release_response_norm,
        rhbm_gem::core::JointProblemAccess::Get(reference).context.linear.release_response_norm);
}
