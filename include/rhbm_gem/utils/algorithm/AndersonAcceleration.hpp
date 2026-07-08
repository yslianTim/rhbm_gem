#pragma once

#include <cstddef>
#include <optional>
#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem::algorithm {

struct AndersonAccelerationOptions
{
    std::size_t history_depth{ 5 };
    double scale_floor{ 1.0 };
    double coefficient_l1_limit{ 10.0 };
    double regularization{ 1.0e-12 };
};

class AndersonAccelerationHistory
{
public:
    struct Sample
    {
        std::vector<Eigen::VectorXd> input_list{};
        std::vector<Eigen::VectorXd> output_list{};
    };

private:
    AndersonAccelerationOptions m_options{};
    std::vector<std::size_t> m_active_index_list{};
    std::vector<Sample> m_sample_list{};

public:
    explicit AndersonAccelerationHistory(AndersonAccelerationOptions options = {});

    void Clear();
    bool HasCompatibleActiveIndexList(const std::vector<std::size_t> & active_index_list) const;
    void Commit(
        const std::vector<std::size_t> & active_index_list,
        const std::vector<Eigen::VectorXd> & input_list,
        const std::vector<Eigen::VectorXd> & output_list);
    std::optional<std::vector<Eigen::VectorXd>> BuildCandidate(
        const std::vector<std::size_t> & active_index_list,
        const std::vector<Eigen::VectorXd> & input_list,
        const std::vector<Eigen::VectorXd> & output_list) const;
};

} // namespace rhbm_gem::algorithm
