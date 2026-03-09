#include <gtest/gtest.h>

#include <rhbm_gem/utils/KDTreeAlgorithm.hpp>

struct TestNode
{
    std::array<double, 3> m_position;
    TestNode(double x, double y, double z) : m_position{ x, y, z } {}
    const std::array<double, 3> & GetPosition() const { return m_position; }
};

std::unique_ptr<KDNode<TestNode>> BuildSimpleTree(
    std::vector<std::unique_ptr<TestNode>> & storage)
{
    storage.clear();
    for (int i = 0; i < 5; i++)
    {
        storage.emplace_back(std::make_unique<TestNode>(i, 0.0, 0.0));
    }
    std::vector<TestNode *> node_ptr_list;
    for (auto & ptr : storage) node_ptr_list.emplace_back(ptr.get());
    return KDTreeAlgorithm<TestNode>::BuildKDTree(node_ptr_list);
}

template <typename NodeType>
std::unique_ptr<KDNode<NodeType>> BuildSampleTree(std::vector<NodeType> & node_list)
{
    std::vector<NodeType *> node_ptr_list;
    node_ptr_list.reserve(node_list.size());
    for (auto & node : node_list) node_ptr_list.emplace_back(&node);
    return KDTreeAlgorithm<NodeType>::BuildKDTree(node_ptr_list, 0);
}

// Building a KD-tree with an empty vector should yield a null root.
TEST(KDTreeAlgorithmTest, BuildKDTreeReturnsNullOnEmpty)
{
    std::vector<TestNode *> node_list;
    auto root{ KDTreeAlgorithm<TestNode>::BuildKDTree(node_list) };
    EXPECT_EQ(root, nullptr);
}

TEST(KDTreeAlgorithmTest, BuildKDTreeSingleNode)
{
    TestNode node{ 1.0, 2.0, 3.0 };
    std::vector<TestNode *> node_ptr_list{ &node };
    auto root{ KDTreeAlgorithm<TestNode>::BuildKDTree(node_ptr_list) };
    ASSERT_NE(root, nullptr);
    EXPECT_EQ(0u, root->m_axis);
    EXPECT_EQ(nullptr, root->m_left);
    EXPECT_EQ(nullptr, root->m_right);
    EXPECT_EQ(&node, root->m_node);
}

TEST(KDTreeAlgorithmTest, KNearestNeighborsAscendingOrder)
{
    std::vector<std::unique_ptr<TestNode>> storage;
    auto root{ BuildSimpleTree(storage) };

    auto knn{
        KDTreeAlgorithm<TestNode>::KNearestNeighbors(root.get(), storage[0].get(), 3)
    };
    ASSERT_EQ(3u, knn.size());
    EXPECT_EQ(storage[0].get(), knn[0]);
    EXPECT_EQ(storage[1].get(), knn[1]);
    EXPECT_EQ(storage[2].get(), knn[2]);
}

TEST(KDTreeAlgorithmTest, KNearestNeighborsWithBuffer)
{
    std::vector<std::unique_ptr<TestNode>> storage;
    auto root{ BuildSimpleTree(storage) };
    std::vector<TestNode *> buffer;
    KDTreeAlgorithm<TestNode>::KNearestNeighbors(root.get(), storage[0].get(), 3, buffer);
    ASSERT_EQ(3u, buffer.size());
    EXPECT_EQ(storage[0].get(), buffer[0]);
    EXPECT_EQ(storage[1].get(), buffer[1]);
    EXPECT_EQ(storage[2].get(), buffer[2]);
}

TEST(KDTreeAlgorithmTest, BuildKDTreeMultiAxisStructure)
{
    std::vector<TestNode> node_list{
        {0.0, 1.0, 2.0},
        {3.0, 4.0, 5.0},
        {6.0, 7.0, 8.0},
        {6.0, 7.0, 8.0}
    };
    auto root{BuildSampleTree(node_list)};

    ASSERT_NE(root, nullptr);
    EXPECT_EQ(0u, root->m_axis);
    EXPECT_EQ(node_list[2].GetPosition(), root->m_node->GetPosition());

    ASSERT_NE(root->m_left, nullptr);
    EXPECT_EQ(1u, root->m_left->m_axis);
    EXPECT_EQ(node_list[1].GetPosition(), root->m_left->m_node->GetPosition());
    ASSERT_NE(root->m_left->m_left, nullptr);
    EXPECT_EQ(2u, root->m_left->m_left->m_axis);
    EXPECT_EQ(node_list[0].GetPosition(), root->m_left->m_left->m_node->GetPosition());
    EXPECT_EQ(nullptr, root->m_left->m_right);

    ASSERT_NE(root->m_right, nullptr);
    EXPECT_EQ(1u, root->m_right->m_axis);
    EXPECT_EQ(root->m_node->GetPosition(), root->m_right->m_node->GetPosition());
    EXPECT_EQ(nullptr, root->m_right->m_left);
    EXPECT_EQ(nullptr, root->m_right->m_right);
}

TEST(KDTreeAlgorithmTest, KSizeLargerThanTreeReturnsAll)
{
    std::vector<std::unique_ptr<TestNode>> storage;
    auto root{ BuildSimpleTree(storage) };

    auto knn{
        KDTreeAlgorithm<TestNode>::KNearestNeighbors(root.get(), storage[0].get(), 10)
    };
    ASSERT_EQ(storage.size(), knn.size());
    for (size_t i = 0; i < storage.size(); i++)
    {
        EXPECT_EQ(storage[i].get(), knn[i]);
    }
}

TEST(KDTreeAlgorithmTest, NullRootReturnsEmpty)
{
    std::vector<std::unique_ptr<TestNode>> storage;
    storage.emplace_back(std::make_unique<TestNode>(0.0, 0.0, 0.0));
    auto result{ KDTreeAlgorithm<TestNode>::KNearestNeighbors(nullptr, storage[0].get(), 3) };
    EXPECT_TRUE(result.empty());
}

TEST(KDTreeAlgorithmTest, NullQueryReturnsEmpty)
{
    std::vector<std::unique_ptr<TestNode>> storage;
    auto root{ BuildSimpleTree(storage) };
    auto result{ KDTreeAlgorithm<TestNode>::KNearestNeighbors(root.get(), nullptr, 3) };
    EXPECT_TRUE(result.empty());
}

TEST(KDTreeAlgorithmTest, KNearestNeighborsZeroKReturnsEmpty)
{
    std::vector<std::unique_ptr<TestNode>> storage;
    auto root{ BuildSimpleTree(storage) };
    auto knn{ KDTreeAlgorithm<TestNode>::KNearestNeighbors(root.get(), storage[0].get(), 0) };
    EXPECT_TRUE(knn.empty());
}

TEST(KDTreeAlgorithmTest, KNearestNeighborsExternalQuery)
{
    std::vector<TestNode> node_list{
        {1.0, 1.0, 1.0},
        {2.0, 2.0, 2.0},
        {3.0, 3.0, 3.0},
        {4.0, 4.0, 4.0}
    };
    auto root{ BuildSampleTree(node_list) };
    TestNode query{ 0.0, 0.0, 0.0 };
    auto knn{ KDTreeAlgorithm<TestNode>::KNearestNeighbors(root.get(), &query, 3) };
    ASSERT_EQ(3u, knn.size());
    EXPECT_EQ(&node_list[0], knn[0]);
    EXPECT_EQ(&node_list[1], knn[1]);
    EXPECT_EQ(&node_list[2], knn[2]);
}

TEST(KDTreeAlgorithmTest, KNearestNeighborsBatchProcessesQueries)
{
    std::vector<std::unique_ptr<TestNode>> storage;
    auto root{ BuildSimpleTree(storage) };

    std::vector<TestNode *> queries{ storage[0].get(), storage[4].get() };
    auto results{ KDTreeAlgorithm<TestNode>::KNearestNeighborsBatch(root.get(), queries, 2) };

    ASSERT_EQ(queries.size(), results.size());
    ASSERT_EQ(2u, results[0].size());
    EXPECT_EQ(storage[0].get(), results[0][0]);
    EXPECT_EQ(storage[1].get(), results[0][1]);
    ASSERT_EQ(2u, results[1].size());
    EXPECT_EQ(storage[4].get(), results[1][0]);
    EXPECT_EQ(storage[3].get(), results[1][1]);
}

TEST(KDTreeAlgorithmTest, KNearestNeighborsTieDistances)
{
    std::vector<TestNode> node_list{
        {1.0, 0.0, 0.0},
        {-1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, -1.0, 0.0}
    };
    auto root{ BuildSampleTree(node_list) };
    TestNode query{ 0.0, 0.0, 0.0 };
    auto one{ KDTreeAlgorithm<TestNode>::KNearestNeighbors(root.get(), &query, 1) };
    ASSERT_EQ(1u, one.size());
    EXPECT_TRUE(one[0] == &node_list[0] || one[0] == &node_list[1] ||
                one[0] == &node_list[2] || one[0] == &node_list[3]);

    auto all{ KDTreeAlgorithm<TestNode>::KNearestNeighbors(root.get(), &query, 4) };
    ASSERT_EQ(4u, all.size());
    EXPECT_NE(all.end(), std::find(all.begin(), all.end(), &node_list[0]));
    EXPECT_NE(all.end(), std::find(all.begin(), all.end(), &node_list[1]));
    EXPECT_NE(all.end(), std::find(all.begin(), all.end(), &node_list[2]));
    EXPECT_NE(all.end(), std::find(all.begin(), all.end(), &node_list[3]));
}

TEST(KDTreeAlgorithmTest, KNearestNeighborsPrunesFarSubtree)
{
    struct TrackingNode
    {
        std::array<double, 3> m_position;
        mutable bool m_visited{ false };
        TrackingNode(double x, double y, double z) : m_position{ x, y, z } {}
        const std::array<double, 3> & GetPosition() const
        {
            m_visited = true;
            return m_position;
        }
        void Reset() { m_visited = false; }
    };

    std::vector<TrackingNode> node_list{
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {100.0, 0.0, 0.0},
        {101.0, 0.0, 0.0},
        {102.0, 0.0, 0.0}
    };
    auto root{ BuildSampleTree(node_list) };
    for (auto & node : node_list) node.Reset();

    TrackingNode query{ 0.0, 0.0, 0.0 };
    auto knn{ KDTreeAlgorithm<TrackingNode>::KNearestNeighbors(root.get(), &query, 1) };
    ASSERT_EQ(1u, knn.size());
    EXPECT_EQ(&node_list[0], knn[0]);

    EXPECT_TRUE(node_list[0].m_visited);
    EXPECT_TRUE(node_list[1].m_visited);
    EXPECT_TRUE(node_list[2].m_visited);
    EXPECT_FALSE(node_list[3].m_visited);
    EXPECT_FALSE(node_list[4].m_visited);
}

TEST(KDTreeAlgorithmTest, RangeSearchIncludesBoundary)
{
    std::vector<TestNode> node_list{
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {2.0, 0.0, 0.0},
        {5.0, 0.0, 0.0}
    };
    auto root{ BuildSampleTree(node_list) };

    auto results{ KDTreeAlgorithm<TestNode>::RangeSearch(root.get(), &node_list[0], 1.0) };
    EXPECT_EQ(2u, results.size());
    EXPECT_NE(results.end(), std::find(results.begin(), results.end(), &node_list[0]));
    EXPECT_NE(results.end(), std::find(results.begin(), results.end(), &node_list[1]));

    auto boundary{ KDTreeAlgorithm<TestNode>::RangeSearch(root.get(), &node_list[0], 2.0) };
    EXPECT_EQ(3u, boundary.size());
    EXPECT_NE(boundary.end(), std::find(boundary.begin(), boundary.end(), &node_list[2]));
}

TEST(KDTreeAlgorithmTest, RangeSearchReturnsEmptyWhenRangeTooSmall)
{
    std::vector<TestNode> node_list{
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {2.0, 0.0, 0.0}
    };
    auto root{ BuildSampleTree(node_list) };

    TestNode query{ 10.0, 10.0, 10.0 };
    auto results{ KDTreeAlgorithm<TestNode>::RangeSearch(root.get(), &query, 1.0) };
    EXPECT_TRUE(results.empty());
}

TEST(KDTreeAlgorithmTest, RangeSearchClearsOutputVector)
{
    std::vector<TestNode> node_list{
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {2.0, 0.0, 0.0}
    };
    auto root{ BuildSampleTree(node_list) };

    std::vector<TestNode *> out{ &node_list[2] };
    KDTreeAlgorithm<TestNode>::RangeSearch(root.get(), &node_list[0], 1.0, out);
    EXPECT_EQ(2u, out.size());
    EXPECT_NE(out.end(), std::find(out.begin(), out.end(), &node_list[0]));
    EXPECT_NE(out.end(), std::find(out.begin(), out.end(), &node_list[1]));
    EXPECT_EQ(out.end(), std::find(out.begin(), out.end(), &node_list[2]));
}

TEST(KDTreeAlgorithmTest, RangeSearchNullRoot)
{
    TestNode query{ 0.0, 0.0, 0.0 };
    auto results{ KDTreeAlgorithm<TestNode>::RangeSearch(nullptr, &query, 1.0) };
    EXPECT_TRUE(results.empty());

    std::vector<TestNode *> out{ &query };
    KDTreeAlgorithm<TestNode>::RangeSearch(nullptr, &query, 1.0, out);
    EXPECT_TRUE(out.empty());
}

TEST(KDTreeAlgorithmTest, RangeSearchNullQuery)
{
    std::vector<TestNode> node_list{
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
    };
    auto root{ BuildSampleTree(node_list) };
    auto results{ KDTreeAlgorithm<TestNode>::RangeSearch(root.get(), nullptr, 1.0) };
    EXPECT_TRUE(results.empty());
    std::vector<TestNode *> out{ &node_list[0] };
    KDTreeAlgorithm<TestNode>::RangeSearch(root.get(), nullptr, 1.0, out);
    EXPECT_TRUE(out.empty());
}

TEST(KDTreeAlgorithmTest, RangeSearchDirectionalPruning)
{
    struct TrackingNode
    {
        std::array<double, 3> m_position;
        mutable bool m_visited{ false };
        TrackingNode(double x, double y, double z) : m_position{ x, y, z } {}
        const std::array<double, 3> & GetPosition() const
        {
            m_visited = true;
            return m_position;
        }
        void Reset() { m_visited = false; }
    };

    std::vector<TrackingNode> node_list{
        {5.0, 0.0, 0.0},
        {2.0, 1.0, 0.0},
        {1.0, -1.0, 0.0},
        {8.0, 1.0, 0.0},
        {9.0, -1.0, 0.0}
    };
    auto root{ BuildSampleTree(node_list) };
    for (auto & node : node_list) node.Reset();

    TrackingNode query_left{1.0, 0.0, 0.0};
    auto left_results{ KDTreeAlgorithm<TrackingNode>::RangeSearch(root.get(), &query_left, 1.0) };
    ASSERT_EQ(1u, left_results.size());
    EXPECT_EQ(&node_list[2], left_results[0]);
    EXPECT_TRUE(node_list[0].m_visited);
    EXPECT_TRUE(node_list[1].m_visited);
    EXPECT_TRUE(node_list[2].m_visited);
    EXPECT_FALSE(node_list[3].m_visited);
    EXPECT_FALSE(node_list[4].m_visited);

    for (auto & node : node_list) node.Reset();

    TrackingNode query_right{9.0, 0.0, 0.0};
    auto right_results{ KDTreeAlgorithm<TrackingNode>::RangeSearch(root.get(), &query_right, 1.0) };
    ASSERT_EQ(1u, right_results.size());
    EXPECT_EQ(&node_list[4], right_results[0]);
    EXPECT_TRUE(node_list[0].m_visited);
    EXPECT_FALSE(node_list[1].m_visited);
    EXPECT_FALSE(node_list[2].m_visited);
    EXPECT_TRUE(node_list[3].m_visited);
    EXPECT_TRUE(node_list[4].m_visited);
}

TEST(KDTreeAlgorithmTest, RangeSearchZeroRange)
{
    std::vector<TestNode> node_list{
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {2.0, 0.0, 0.0}
    };
    auto root{ BuildSampleTree(node_list) };

    TestNode query{1.0, 0.0, 0.0};
    auto results{ KDTreeAlgorithm<TestNode>::RangeSearch(root.get(), &query, 0.0) };
    ASSERT_EQ(1u, results.size());
    EXPECT_EQ(&node_list[1], results[0]);
}

TEST(KDTreeAlgorithmTest, RangeSearchNegativeRangeReturnsEmpty)
{
    std::vector<TestNode> node_list{
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
    };
    auto root{ BuildSampleTree(node_list) };

    TestNode query{ 0.0, 0.0, 0.0 };
    auto results{ KDTreeAlgorithm<TestNode>::RangeSearch(root.get(), &query, -1.0) };
    EXPECT_TRUE(results.empty());

    std::vector<TestNode *> out{ &node_list[1] };
    KDTreeAlgorithm<TestNode>::RangeSearch(root.get(), &query, -1.0, out);
    EXPECT_TRUE(out.empty());
}

TEST(KDTreeAlgorithmTest, RangeSearchLargeRangeReturnsAllNodes)
{
    std::vector<TestNode> node_list{
        {0.0, 0.0, 0.0},
        {1.0, 1.0, 1.0},
        {2.0, 2.0, 2.0},
        {-1.0, -1.0, -1.0},
        {3.0, 0.0, 0.0}
    };
    auto root{ BuildSampleTree(node_list) };

    double max_dist{ 0.0 };
    for (size_t i = 0; i < node_list.size(); i++)
    {
        for (size_t j = i + 1; j < node_list.size(); j++)
        {
            const auto & a{ node_list[i].GetPosition() };
            const auto & b{ node_list[j].GetPosition() };
            double dist{ 0.0 };
            for (size_t axis = 0; axis < a.size(); axis++)
            {
                auto diff{ a[axis] - b[axis] };
                dist += diff * diff;
            }
            dist = std::sqrt(dist);
            max_dist = std::max(max_dist, dist);
        }
    }

    auto results{
        KDTreeAlgorithm<TestNode>::RangeSearch(root.get(), root->m_node, max_dist + 1.0)
    };
    ASSERT_EQ(node_list.size(), results.size());
    for (auto & node : node_list)
    {
        EXPECT_NE(results.end(), std::find(results.begin(), results.end(), &node));
    }
}
