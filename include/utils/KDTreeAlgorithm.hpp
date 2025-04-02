#pragma once

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <memory>
#include <queue>

template <typename NodeType>
struct KDNode
{
    NodeType * m_node;
    int m_axis;
    std::unique_ptr<KDNode> m_left;
    std::unique_ptr<KDNode> m_right;
    KDNode(NodeType * node, int axis) :
        m_node{ node }, m_axis{ axis }, m_left{ nullptr }, m_right{ nullptr } {}
};

template <typename NodeType>
class KDTreeAlgorithm
{
    struct DistNode
    {
        double m_distance;
        NodeType * m_node;
        DistNode(double dist, NodeType * node) :
            m_distance{ dist }, m_node{ node } {}
    };

    struct DistNodeComparator
    {
        bool operator()(const DistNode & a, const DistNode & b) const
        {
            return a.m_distance < b.m_distance;
        }
    };

public:
    KDTreeAlgorithm(void) = default;
    ~KDTreeAlgorithm() {}
    static std::unique_ptr<KDNode<NodeType>> BuildKDTree(
        const std::vector<NodeType *> & node_list, int depth = 0)
    {
        return BuildKDTree(node_list.begin(), node_list.end(), depth);
    }

    static std::vector<NodeType *> KNearestNeighbors(
        const KDNode<NodeType> * root_kd_node, const NodeType * query_node, int k_size)
    {
        if (k_size <= 0) return {};

        std::priority_queue<DistNode, std::vector<DistNode>, DistNodeComparator> max_heap;
        KNearestNeighborsHelper(root_kd_node, query_node, k_size, max_heap);
        std::vector<NodeType *> knn_list;
        knn_list.reserve(max_heap.size());
        while (!max_heap.empty())
        {
            knn_list.emplace_back(max_heap.top().m_node);
            max_heap.pop();
        }
        std::reverse(knn_list.begin(), knn_list.end());

        return knn_list;
    }

    static std::vector<NodeType *> RangeSearch(
        const KDNode<NodeType> * root_kd_node, const NodeType * query_node, double range)
    {
        std::vector<NodeType *> knn_list;
        knn_list.reserve(64);
        RangeSearchHelper(root_kd_node, query_node, range, knn_list);
        return knn_list;
    }

private:
    static std::unique_ptr<KDNode<NodeType>> BuildKDTree(
        typename std::vector<NodeType *>::iterator begin,
        typename std::vector<NodeType *>::iterator end,
        int depth = 0)
    {
        if (begin == end)
        {
            return nullptr;
        }

        constexpr int dimension{ 3 };
        int axis{ depth % dimension };
        auto count{ std::distance(begin, end) };
        auto mid_iter{ begin + count / 2 };

        std::nth_element(
            begin,
            mid_iter,
            end,
            [axis](const NodeType * a, const NodeType * b)
            {
                return a->GetPosition().at(axis) < b->GetPosition().at(axis);
            }
        );

        std::unique_ptr<KDNode<NodeType>> kd_node{ std::make_unique<KDNode<NodeType>>(*mid_iter, axis) };

        kd_node->m_left  = BuildKDTree(begin, mid_iter, depth + 1);
        kd_node->m_right = BuildKDTree(mid_iter + 1, end, depth + 1);

        return kd_node;
    }

    static double ComputeNodeDistanceSquare(
        const NodeType * node_a, const NodeType * node_b)
    {
        auto position_a{ node_a->GetPosition() };
        auto position_b{ node_b->GetPosition() };
        return (position_a - position_b).squaredNorm();
    }

    static double ComputeNodeDifference(
        const NodeType * node_a, const NodeType * node_b, int axis)
    {
        auto position_a{ node_a->GetPosition() };
        auto position_b{ node_b->GetPosition() };
        return position_a.at(axis) - position_b.at(axis);
    }

    static void KNearestNeighborsHelper(
        const KDNode<NodeType> * kd_node, const NodeType * query_node, int k,
        std::priority_queue<DistNode, std::vector<DistNode>, DistNodeComparator> & max_heap)
    {
        if (kd_node == nullptr) return;

        double dist{ ComputeNodeDistanceSquare(query_node, kd_node->m_node) };
        if (static_cast<int>(max_heap.size()) < k)
        {
            max_heap.push(DistNode(dist, kd_node->m_node));
        }
        else if (dist < max_heap.top().m_distance)
        {
            max_heap.pop();
            max_heap.push(DistNode(dist, kd_node->m_node));
        }

        int axis{ kd_node->m_axis };
        double diff{ ComputeNodeDifference(query_node, kd_node->m_node, axis) };

        const KDNode<NodeType> * next_branch{ (diff < 0.0) ? kd_node->m_left.get() : kd_node->m_right.get() };
        const KDNode<NodeType> * other_branch{ (diff < 0.0) ? kd_node->m_right.get() : kd_node->m_left.get() };

        KNearestNeighborsHelper(next_branch, query_node, k, max_heap);

        double dist_axis{ diff * diff };
        if (static_cast<int>(max_heap.size()) < k || dist_axis < max_heap.top().m_distance)
        {
            KNearestNeighborsHelper(other_branch, query_node, k, max_heap);
        }
    }

    static void RangeSearchHelper(
        const KDNode<NodeType> * kd_node, const NodeType * query_node,
        double range, std::vector<NodeType *> & results)
    {
        if (kd_node == nullptr) return;

        auto dist_square{ ComputeNodeDistanceSquare(query_node, kd_node->m_node) };
        auto range_square{ range * range };
        if (dist_square <= range_square) results.emplace_back(kd_node->m_node);

        auto axis{ kd_node->m_axis };
        auto diff{ ComputeNodeDifference(query_node, kd_node->m_node, axis) };
        if (diff < -range)
        {
            RangeSearchHelper(kd_node->m_left.get(), query_node, range, results);
        }
        else if (diff > range)
        {
            RangeSearchHelper(kd_node->m_right.get(), query_node, range, results);
        }
        else
        {
            RangeSearchHelper(kd_node->m_left.get(), query_node, range, results);
            RangeSearchHelper(kd_node->m_right.get(), query_node, range, results);
        }
    }


};