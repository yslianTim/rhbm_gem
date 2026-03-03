#pragma once

#include <cstddef>
#include <vector>
#include <algorithm>
#include <memory>
#include <queue>
#include <cmath>
#include <atomic>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "Logger.hpp"

template <typename NodeType>
struct KDNode
{
    NodeType * m_node;
    size_t m_axis;
    std::unique_ptr<KDNode> m_left;
    std::unique_ptr<KDNode> m_right;
    KDNode(NodeType * node, size_t axis) :
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
    KDTreeAlgorithm() = default;
    ~KDTreeAlgorithm() {}
    static std::unique_ptr<KDNode<NodeType>> BuildKDTree(
        std::vector<NodeType *> & node_list,
        int depth = 0,
        int thread_size = 1,
        bool show_progress = false)
    {
        std::atomic<size_t> progress{ 0 };
        auto total{ node_list.size() };
        return BuildKDTree(
            node_list.begin(),
            node_list.end(),
            depth,
            thread_size,
            show_progress ? &progress : nullptr,
            total
        );
    }

    static std::unique_ptr<KDNode<NodeType>> BuildKDTree(
        std::vector<NodeType> & node_list,
        int depth = 0,
        int thread_size = 1,
        bool show_progress = false)
    {
        std::vector<NodeType *> node_ptr_list;
        node_ptr_list.reserve(node_list.size());
        for (auto & node : node_list)
        {
            node_ptr_list.emplace_back(&node);
        }
        return BuildKDTree(node_ptr_list, depth, thread_size, show_progress);
    }

    static void KNearestNeighbors(
        const KDNode<NodeType> * root_kd_node,
        const NodeType * query_node,
        size_t k_size,
        std::vector<NodeType *> & knn_list)
    {
        knn_list.clear();
        if (root_kd_node == nullptr || query_node == nullptr || k_size == 0) return;

        if (knn_list.capacity() < k_size) knn_list.reserve(k_size);
        std::priority_queue<DistNode, std::vector<DistNode>, DistNodeComparator> max_heap;
        KNearestNeighborsHelper(root_kd_node, query_node, static_cast<int>(k_size), max_heap);
        while (!max_heap.empty())
        {
            knn_list.emplace_back(max_heap.top().m_node);
            max_heap.pop();
        }
        std::reverse(knn_list.begin(), knn_list.end());
    }

    static std::vector<NodeType *> KNearestNeighbors(
        const KDNode<NodeType> * root_kd_node,
        const NodeType * query_node,
        size_t k_size)
    {
        std::vector<NodeType *> knn_list;
        if (k_size > 0) knn_list.reserve(k_size);
        KNearestNeighbors(root_kd_node, query_node, k_size, knn_list);
        return knn_list;
    }

    static std::vector<std::vector<NodeType *>> KNearestNeighborsBatch(
        const KDNode<NodeType> * root_kd_node,
        const std::vector<NodeType *> & query_nodes,
        size_t k_size, int thread_size = 1)
    {
        std::vector<std::vector<NodeType *>> results(query_nodes.size());

#ifndef USE_OPENMP
        (void)thread_size;
#endif

#ifdef USE_OPENMP
        #pragma omp parallel for num_threads(thread_size)
#endif
        for (size_t i = 0; i < query_nodes.size(); i++)
        {
            results[i] = KNearestNeighbors(root_kd_node, query_nodes[i], k_size);
        }
        return results;
    }

    static std::vector<NodeType *> RangeSearch(
        const KDNode<NodeType> * root_kd_node,
        const NodeType * query_node,
        double range)
    {
        if (range < 0.0 || root_kd_node == nullptr || query_node == nullptr)
        {
            return {};
        }
        std::vector<NodeType *> knn_list;
        knn_list.reserve(64);
        RangeSearchHelper(root_kd_node, query_node, range, knn_list);
        return knn_list;
    }

    static void RangeSearch(
        const KDNode<NodeType> * root_kd_node,
        const NodeType * query_node,
        double range,
        std::vector<NodeType *> & knn_list)
    {
        if (range < 0.0 || root_kd_node == nullptr || query_node == nullptr)
        {
            knn_list = {};
            return;
        }
        knn_list.clear();
        RangeSearchHelper(root_kd_node, query_node, range, knn_list);
    }

private:
    static std::unique_ptr<KDNode<NodeType>> BuildKDTree(
        typename std::vector<NodeType *>::iterator begin,
        typename std::vector<NodeType *>::iterator end,
        int depth = 0,
        int thread_size = 1,
        std::atomic<size_t> * progress = nullptr,
        size_t total = 0)
    {
        if (begin == end)
        {
            return nullptr;
        }

        constexpr int dimension{ 3 };
        auto axis{ static_cast<size_t>(depth % dimension) };
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

        std::unique_ptr<KDNode<NodeType>> kd_node{
            std::make_unique<KDNode<NodeType>>(*mid_iter, axis)
        };

        if (progress)
        {
            auto current{ progress->fetch_add(1) + 1 };
            Logger::ProgressPercent(current, total);
        }

#ifdef USE_OPENMP
        constexpr int parallel_threshold{ 1024 };
        if (thread_size > 1 && count >= parallel_threshold)
        {
            auto left_threads{ thread_size / 2 };
            auto right_threads{ thread_size - left_threads };
            if (depth == 0)
            {
                #pragma omp parallel num_threads(thread_size)
                #pragma omp single nowait
                {
                    #pragma omp taskgroup
                    {
                        #pragma omp task shared(kd_node, progress, total) if(left_threads > 0)
                        kd_node->m_left = BuildKDTree(
                            begin, mid_iter, depth + 1, left_threads, progress, total);

                        #pragma omp task shared(kd_node, progress, total) if(right_threads > 0)
                        kd_node->m_right = BuildKDTree(
                            mid_iter + 1, end, depth + 1, right_threads, progress, total);
                    }
                }
            }
            else
            {
                #pragma omp taskgroup
                {
                    #pragma omp task shared(kd_node, progress, total) if(left_threads > 0)
                    kd_node->m_left = BuildKDTree(
                        begin, mid_iter, depth + 1, left_threads, progress, total);

                    #pragma omp task shared(kd_node, progress, total) if(right_threads > 0)
                    kd_node->m_right = BuildKDTree(
                        mid_iter + 1, end, depth + 1, right_threads, progress, total);
                }
            }
        }
        else
        {
            kd_node->m_left  = BuildKDTree(begin, mid_iter, depth + 1, 1, progress, total);
            kd_node->m_right = BuildKDTree(mid_iter + 1, end, depth + 1, 1, progress, total);
        }
#else
        kd_node->m_left  = BuildKDTree(begin, mid_iter, depth + 1, thread_size, progress, total);
        kd_node->m_right = BuildKDTree(mid_iter + 1, end, depth + 1, thread_size, progress, total);
#endif

        return kd_node;
    }

    static double ComputeNodeDistanceSquare(
        const NodeType * node_a, const NodeType * node_b)
    {
        auto position_a{ node_a->GetPosition() };
        auto position_b{ node_b->GetPosition() };
        auto dimension{ position_a.size() };
        auto squared_norm{ 0.0 };
        for (size_t i = 0; i < dimension; i++)
        {
            squared_norm += std::pow(position_a.at(i) - position_b.at(i), 2.0);
        }
        return squared_norm;
    }

    static double ComputeNodeDifference(
        const NodeType * node_a, const NodeType * node_b, size_t axis)
    {
        auto position_a{ node_a->GetPosition() };
        auto position_b{ node_b->GetPosition() };
        return position_a.at(axis) - position_b.at(axis);
    }

    static void KNearestNeighborsHelper(
        const KDNode<NodeType> * kd_node,
        const NodeType * query_node,
        int k,
        std::priority_queue<DistNode, std::vector<DistNode>, DistNodeComparator> & max_heap)
    {
        if (kd_node == nullptr) return;

        auto dist{ ComputeNodeDistanceSquare(query_node, kd_node->m_node) };
        if (static_cast<int>(max_heap.size()) < k)
        {
            max_heap.push(DistNode(dist, kd_node->m_node));
        }
        else if (dist < max_heap.top().m_distance)
        {
            max_heap.pop();
            max_heap.push(DistNode(dist, kd_node->m_node));
        }

        auto axis{ kd_node->m_axis };
        auto diff{ ComputeNodeDifference(query_node, kd_node->m_node, axis) };

        const KDNode<NodeType> * next_branch{
            (diff < 0.0) ? kd_node->m_left.get() : kd_node->m_right.get()
        };
        const KDNode<NodeType> * other_branch{
            (diff < 0.0) ? kd_node->m_right.get() : kd_node->m_left.get()
        };

        KNearestNeighborsHelper(next_branch, query_node, k, max_heap);

        auto dist_axis{ diff * diff };
        if (static_cast<int>(max_heap.size()) < k || dist_axis < max_heap.top().m_distance)
        {
            KNearestNeighborsHelper(other_branch, query_node, k, max_heap);
        }
    }

    static void RangeSearchHelper(
        const KDNode<NodeType> * kd_node,
        const NodeType * query_node,
        double range,
        std::vector<NodeType *> & results)
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
