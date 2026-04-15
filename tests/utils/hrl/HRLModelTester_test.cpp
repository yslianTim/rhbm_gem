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
