#include <gtest/gtest.h>

#include <rhbm_gem/utils/algorithm/ClusteredAndersonAcceleration.hpp>

#include <algorithm>

namespace {
namespace alg = rhbm_gem::algorithm;

std::vector<Eigen::VectorXd> MakeState(double left, double right)
{
    return std::vector<Eigen::VectorXd>{
        Eigen::VectorXd::Constant(1, left),
        Eigen::VectorXd::Constant(1, right)
    };
}

} // namespace

TEST(ClusteredAndersonAccelerationTest, BuildsCandidateFromIndependentClusters)
{
    alg::ClusteredAndersonAccelerationHistorySet history{
        alg::AndersonAccelerationOptions{
            5,
            100.0,
            10.0,
            1.0e-12
        }
    };
    const std::vector<alg::ClusterKey> key_list{ { 0 }, { 1 } };
    history.Reconcile(key_list);
    history.Commit(key_list, MakeState(0.0, 0.0), MakeState(1.0, 1.0));

    const auto candidate{
        history.BuildCandidate(key_list, MakeState(1.0, 1.0), MakeState(1.5, 1.5))
    };

    ASSERT_TRUE(candidate.has_value());
    EXPECT_EQ(2u, candidate->used_cluster_key_list.size());
    ASSERT_EQ(2u, candidate->state_list.size());
    EXPECT_NEAR(candidate->state_list.at(0)(0), 2.0, 1.0e-6);
    EXPECT_NEAR(candidate->state_list.at(1)(0), 2.0, 1.0e-6);
}

TEST(ClusteredAndersonAccelerationTest, SuppressesOnlyRequestedCluster)
{
    alg::ClusteredAndersonAccelerationHistorySet history{
        alg::AndersonAccelerationOptions{
            5,
            100.0,
            10.0,
            1.0e-12
        }
    };
    const std::vector<alg::ClusterKey> key_list{ { 0 }, { 1 } };
    history.Reconcile(key_list);
    history.Commit(key_list, MakeState(0.0, 0.0), MakeState(1.0, 1.0));
    history.ClearAndSuppress({ { 0 } });

    const auto candidate{
        history.BuildCandidate(key_list, MakeState(1.0, 1.0), MakeState(1.5, 1.5))
    };

    ASSERT_TRUE(candidate.has_value());
    ASSERT_EQ(1u, candidate->used_cluster_key_list.size());
    EXPECT_EQ((alg::ClusterKey{ 1 }), candidate->used_cluster_key_list.front());
    EXPECT_NEAR(candidate->state_list.at(0)(0), 1.5, 1.0e-12);
    EXPECT_NEAR(candidate->state_list.at(1)(0), 2.0, 1.0e-6);
}

TEST(ClusteredAndersonAccelerationTest, SkipsClusterWhenSingleCoefficientLimitIsExceeded)
{
    alg::ClusteredAndersonAccelerationHistorySet history{
        alg::AndersonAccelerationOptions{
            5,
            100.0,
            10.0,
            1.0e-12,
            1.5
        }
    };
    const std::vector<alg::ClusterKey> key_list{ { 0 }, { 1 } };
    history.Reconcile(key_list);
    history.Commit(key_list, MakeState(0.0, 0.0), MakeState(1.0, 1.0));

    const auto candidate{
        history.BuildCandidate(key_list, MakeState(1.0, 1.0), MakeState(1.5, 2.0))
    };

    ASSERT_TRUE(candidate.has_value());
    ASSERT_EQ(1u, candidate->used_cluster_key_list.size());
    EXPECT_EQ((alg::ClusterKey{ 1 }), candidate->used_cluster_key_list.front());
    ASSERT_EQ(2u, candidate->state_list.size());
    EXPECT_NEAR(candidate->state_list.at(0)(0), 1.5, 1.0e-12);
    EXPECT_NEAR(candidate->state_list.at(1)(0), 1.5, 1.0e-6);
}

TEST(ClusteredAndersonAccelerationTest, SuppressedClusterRequiresNewCommitAfterRelease)
{
    alg::ClusteredAndersonAccelerationHistorySet history{
        alg::AndersonAccelerationOptions{
            5,
            100.0,
            10.0,
            1.0e-12
        }
    };
    const std::vector<alg::ClusterKey> key_list{ { 0 }, { 1 } };
    history.Reconcile(key_list);
    history.Commit(key_list, MakeState(0.0, 0.0), MakeState(1.0, 1.0));
    history.ClearAndSuppress({ { 0 } });
    history.ReleaseSuppression({ { 0 } });

    const auto released_candidate{
        history.BuildCandidate(key_list, MakeState(1.0, 1.0), MakeState(1.5, 1.5))
    };

    ASSERT_TRUE(released_candidate.has_value());
    ASSERT_EQ(1u, released_candidate->used_cluster_key_list.size());
    EXPECT_EQ((alg::ClusterKey{ 1 }), released_candidate->used_cluster_key_list.front());

    history.Commit({ { 0 } }, MakeState(1.0, 1.0), MakeState(1.5, 1.5));
    const auto recommitted_candidate{
        history.BuildCandidate(key_list, MakeState(1.5, 1.5), MakeState(1.75, 1.75))
    };

    ASSERT_TRUE(recommitted_candidate.has_value());
    EXPECT_NE(
        std::find(
            recommitted_candidate->used_cluster_key_list.begin(),
            recommitted_candidate->used_cluster_key_list.end(),
            alg::ClusterKey{ 0 }),
        recommitted_candidate->used_cluster_key_list.end());
}

TEST(ClusteredAndersonAccelerationTest, SuppressesContainingAtomWithoutClearingRemoteCluster)
{
    alg::ClusteredAndersonAccelerationHistorySet history{
        alg::AndersonAccelerationOptions{
            5,
            100.0,
            10.0,
            1.0e-12
        }
    };
    const std::vector<alg::ClusterKey> key_list{ { 0 }, { 1 } };
    history.Reconcile(key_list);
    history.Commit(key_list, MakeState(0.0, 0.0), MakeState(1.0, 1.0));
    history.ClearAndSuppressContaining({ 0 });

    const auto candidate{
        history.BuildCandidate(key_list, MakeState(1.0, 1.0), MakeState(1.5, 1.5))
    };

    ASSERT_TRUE(candidate.has_value());
    ASSERT_EQ(1u, candidate->used_cluster_key_list.size());
    EXPECT_EQ((alg::ClusterKey{ 1 }), candidate->used_cluster_key_list.front());
}

TEST(ClusteredAndersonAccelerationTest, ReconcileDropsMissingClusterState)
{
    alg::ClusteredAndersonAccelerationHistorySet history{
        alg::AndersonAccelerationOptions{
            5,
            100.0,
            10.0,
            1.0e-12
        }
    };
    history.Reconcile({ { 0 }, { 1 } });
    history.Commit({ { 0 }, { 1 } }, MakeState(0.0, 0.0), MakeState(1.0, 1.0));
    history.Reconcile({ { 1 } });

    const auto candidate{
        history.BuildCandidate({ { 0 }, { 1 } }, MakeState(1.0, 1.0), MakeState(1.5, 1.5))
    };

    ASSERT_TRUE(candidate.has_value());
    ASSERT_EQ(1u, candidate->used_cluster_key_list.size());
    EXPECT_EQ((alg::ClusterKey{ 1 }), candidate->used_cluster_key_list.front());
}
