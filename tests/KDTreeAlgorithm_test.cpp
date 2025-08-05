#include <gtest/gtest.h>

#include "KDTreeAlgorithm.hpp"

struct TestNode
{
    std::array<double, 3> m_position;
    TestNode(double x, double y, double z) : m_position{ x, y, z } {}
    const std::array<double, 3> & GetPosition(void) const { return m_position; }
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

std::unique_ptr<KDNode<TestNode>> BuildSampleTree(std::vector<TestNode> & node_list)
{
    std::vector<TestNode *> node_ptr_list;
    node_ptr_list.reserve(node_list.size());
    for (auto & node : node_list) node_ptr_list.emplace_back(&node);
    return KDTreeAlgorithm<TestNode>::BuildKDTree(node_ptr_list, 0);
}

// Building a KD-tree with an empty vector should yield a null root.
TEST(KDTreeAlgorithmTest, BuildKDTreeReturnsNullOnEmpty)
{
    std::vector<TestNode *> node_list;
    auto root{ KDTreeAlgorithm<TestNode>::BuildKDTree(node_list) };
    EXPECT_EQ(root, nullptr);
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

TEST(KDTreeAlgorithmTest, ZeroKReturnsEmpty)
{
    std::vector<std::unique_ptr<TestNode>> storage;
    auto root{ BuildSimpleTree(storage) };

    auto zero_result{
        KDTreeAlgorithm<TestNode>::KNearestNeighbors(root.get(), storage[0].get(), 0)
    };
    EXPECT_TRUE(zero_result.empty());
}

TEST(KDTreeAlgorithmTest, NullRootReturnsEmpty)
{
    std::vector<std::unique_ptr<TestNode>> storage;
    storage.emplace_back(std::make_unique<TestNode>(0.0, 0.0, 0.0));
    auto result{
        KDTreeAlgorithm<TestNode>::KNearestNeighbors(nullptr, storage[0].get(), 3)
    };
    EXPECT_TRUE(result.empty());
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

    TestNode query{10.0, 10.0, 10.0};
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
    TestNode query{0.0, 0.0, 0.0};
    auto results{ KDTreeAlgorithm<TestNode>::RangeSearch(nullptr, &query, 1.0) };
    EXPECT_TRUE(results.empty());

    std::vector<TestNode *> out{ &query };
    KDTreeAlgorithm<TestNode>::RangeSearch(nullptr, &query, 1.0, out);
    EXPECT_TRUE(out.empty());
}
