#pragma once

#include <optional>

#include <Eigen/Dense>

#include <rhbm_gem/core/TestDataFactory.hpp>

namespace rhbm_gem::core {

struct BiasStatistics
{
    Eigen::VectorXd mean;
    Eigen::VectorXd sigma;
};

struct BiasStatisticsSeries
{
    BiasStatistics requested_alpha;
    std::optional<BiasStatistics> trained_alpha;
    std::optional<double> trained_alpha_median;
};

struct LocalTestBias
{
    BiasStatistics ols;
    BiasStatisticsSeries mdpde;
};

struct GroupTestBias
{
    BiasStatistics median;
    BiasStatisticsSeries mdpde;
};

struct LocalTestOptions
{
    double requested_alpha_r{ 0.0 };
    bool alpha_training{ true };
    int thread_size{ 1 };
};

struct GroupTestOptions
{
    double requested_alpha_g{ 0.0 };
    bool alpha_training{ true };
    int thread_size{ 1 };
};

LocalTestBias RunLocalEstimationTest(
    const LocalTestData & input,
    const LocalTestOptions & options
);

GroupTestBias RunGroupEstimationTest(
    const GroupTestData & input,
    const GroupTestOptions & options
);

} // namespace rhbm_gem::core
