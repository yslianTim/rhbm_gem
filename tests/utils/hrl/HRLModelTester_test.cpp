#include <gtest/gtest.h>

#include <initializer_list>

#include <rhbm_gem/utils/hrl/HRLModelTester.hpp>

class HRLModelTesterTestPeer
{
public:
    static LocalPotentialSampleList BuildRandomGausSamplingEntryWithNeighborhood(
        HRLModelTester & tester,
        size_t sampling_entry_size,
        const Eigen::VectorXd & gaus_par,
        double neighbor_distance,
        size_t neighbor_count,
        double angle)
    {
        return tester.BuildRandomGausSamplingEntryWithNeighborhood(
            sampling_entry_size,
            0.0, 1.0,
            gaus_par,
            neighbor_distance,
            neighbor_count,
            angle);
    }
};

namespace {

Eigen::VectorXd MakeVector(std::initializer_list<double> values)
{
    Eigen::VectorXd result(static_cast<Eigen::Index>(values.size()));
    Eigen::Index index{ 0 };
    for (double value : values)
    {
        result(index++) = value;
    }
    return result;
}

} // namespace

TEST(HRLModelTesterTest, NeighborhoodSamplingKeepsPointCountWhenAngleIsZero)
{
    HRLModelTester tester(3, 2, 2);
    tester.SetFittingRange(0.0, 1.0);
    const Eigen::VectorXd gaus_par{ MakeVector({ 1.0, 0.5, 0.0 }) };

    const auto sampling_entries{
        HRLModelTesterTestPeer::BuildRandomGausSamplingEntryWithNeighborhood(
            tester,
            4,
            gaus_par,
            2.0,
            1,
            0.0)
    };

    EXPECT_EQ(sampling_entries.size(), 44u);
}

TEST(HRLModelTesterTest, NeighborhoodSamplingRemovesPointsWhenAngleFilteringEnabled)
{
    HRLModelTester tester(3, 2, 2);
    tester.SetFittingRange(0.0, 1.0);
    const Eigen::VectorXd gaus_par{ MakeVector({ 1.0, 0.5, 0.0 }) };

    const auto unfiltered_sampling_entries{
        HRLModelTesterTestPeer::BuildRandomGausSamplingEntryWithNeighborhood(
            tester,
            4,
            gaus_par,
            2.0,
            1,
            0.0)
    };
    const auto filtered_sampling_entries{
        HRLModelTesterTestPeer::BuildRandomGausSamplingEntryWithNeighborhood(
            tester,
            4,
            gaus_par,
            2.0,
            1,
            120.0)
    };

    ASSERT_FALSE(unfiltered_sampling_entries.empty());
    EXPECT_LT(filtered_sampling_entries.size(), unfiltered_sampling_entries.size());
}

TEST(HRLModelTesterTest, PublicNeighborhoodWorkflowAcceptsExplicitAngle)
{
    HRLModelTester tester(3, 2, 2);
    tester.SetFittingRange(0.0, 1.0);

    std::vector<Eigen::VectorXd> residual_mean_ols_list;
    std::vector<Eigen::VectorXd> residual_mean_mdpde_list;
    std::vector<Eigen::VectorXd> residual_sigma_ols_list;
    std::vector<Eigen::VectorXd> residual_sigma_mdpde_list;
    double training_alpha_r_average{ 0.0 };

    const auto success{
        tester.RunBetaMDPDEWithNeighborhoodTest(
            { 0.1 },
            residual_mean_ols_list,
            residual_mean_mdpde_list,
            residual_sigma_ols_list,
            residual_sigma_mdpde_list,
            MakeVector({ 1.0, 0.5, 0.0 }),
            training_alpha_r_average,
            4,
            0.0,
            2.0,
            1,
            120.0)
    };

    EXPECT_TRUE(success);
    EXPECT_EQ(residual_mean_ols_list.size(), 2u);
    EXPECT_EQ(residual_mean_mdpde_list.size(), 2u);
    EXPECT_EQ(residual_sigma_ols_list.size(), 2u);
    EXPECT_EQ(residual_sigma_mdpde_list.size(), 2u);
}
