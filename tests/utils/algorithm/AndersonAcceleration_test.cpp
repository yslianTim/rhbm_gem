#include <gtest/gtest.h>

#include <rhbm_gem/utils/algorithm/AndersonAcceleration.hpp>

#include <cmath>
#include <limits>

namespace {
namespace alg = rhbm_gem::algorithm;

std::vector<Eigen::VectorXd> MakeState(double value)
{
    return std::vector<Eigen::VectorXd>{ Eigen::VectorXd::Constant(1, value) };
}

} // namespace

TEST(AndersonAccelerationTest, AcceleratesLinearFixedPointStep)
{
    alg::AndersonAccelerationHistory history{
        alg::AndersonAccelerationOptions{
            5,
            100.0,
            10.0,
            1.0e-12
        }
    };
    const std::vector<std::size_t> active_index_list{ 0 };
    history.Commit(active_index_list, MakeState(0.0), MakeState(1.0));

    const auto candidate{
        history.BuildCandidate(active_index_list, MakeState(1.0), MakeState(1.5))
    };

    ASSERT_TRUE(candidate.has_value());
    ASSERT_EQ(1u, candidate->size());
    EXPECT_LT(std::abs(candidate->front()(0) - 2.0), std::abs(1.5 - 2.0));
    EXPECT_NEAR(candidate->front()(0), 2.0, 1.0e-6);
}

TEST(AndersonAccelerationTest, RequiresCommittedHistory)
{
    alg::AndersonAccelerationHistory history;

    const auto candidate{
        history.BuildCandidate({ 0 }, MakeState(0.0), MakeState(1.0))
    };

    EXPECT_FALSE(candidate.has_value());
}

TEST(AndersonAccelerationTest, RejectsNonFiniteSamples)
{
    alg::AndersonAccelerationHistory history;
    const std::vector<std::size_t> active_index_list{ 0 };
    history.Commit(active_index_list, MakeState(0.0), MakeState(1.0));
    auto bad_state{ MakeState(1.0) };
    bad_state.front()(0) = std::numeric_limits<double>::quiet_NaN();

    const auto candidate{
        history.BuildCandidate(active_index_list, bad_state, MakeState(1.5))
    };

    EXPECT_FALSE(candidate.has_value());
}

TEST(AndersonAccelerationTest, RejectsExcessiveExtrapolation)
{
    alg::AndersonAccelerationHistory history{
        alg::AndersonAccelerationOptions{
            5,
            100.0,
            2.0,
            1.0e-12
        }
    };
    const std::vector<std::size_t> active_index_list{ 0 };
    history.Commit(active_index_list, MakeState(0.0), MakeState(1.0));

    const auto candidate{
        history.BuildCandidate(active_index_list, MakeState(1.0), MakeState(1.5))
    };

    EXPECT_FALSE(candidate.has_value());
}

TEST(AndersonAccelerationTest, InvalidatesWhenActiveIndexesDiffer)
{
    alg::AndersonAccelerationHistory history;
    history.Commit({ 0 }, MakeState(0.0), MakeState(1.0));

    EXPECT_FALSE(history.HasCompatibleActiveIndexList({ 1 }));
    EXPECT_FALSE(history.BuildCandidate({ 1 }, MakeState(1.0), MakeState(1.5)).has_value());

    history.Commit({ 1 }, MakeState(1.0), MakeState(1.5));

    EXPECT_TRUE(history.HasCompatibleActiveIndexList({ 1 }));
    EXPECT_FALSE(history.BuildCandidate({ 1 }, MakeState(1.5), MakeState(1.75)).has_value());
}
