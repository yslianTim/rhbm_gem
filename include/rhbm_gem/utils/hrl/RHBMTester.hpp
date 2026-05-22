#pragma once

#include <optional>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/TestDataFactory.hpp>

namespace rhbm_gem::rhbm_tester
{

struct BiasStatistics
{
    Eigen::VectorXd mean;
    Eigen::VectorXd sigma;
};

struct BiasStatisticsSeries
{
    std::vector<BiasStatistics> requested_alpha;
    std::optional<BiasStatistics> trained_alpha;
    std::optional<double> trained_alpha_average;
};

struct BetaMDPDETestBias
{
    BiasStatistics ols;
    BiasStatisticsSeries mdpde;
};

struct MuMDPDETestBias
{
    BiasStatistics median;
    BiasStatisticsSeries mdpde;
};

BetaMDPDETestBias RunBetaMDPDETest(
    const test_data_factory::LocalTestData & test_input,
    int thread_size
);

MuMDPDETestBias RunMuMDPDETest(
    const test_data_factory::GroupTestData & test_input,
    int thread_size
);

} // namespace rhbm_gem::rhbm_tester
