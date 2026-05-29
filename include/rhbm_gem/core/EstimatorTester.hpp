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

struct BetaMDPDETestOptions
{
    double requested_alpha_r{ 0.0 };
    bool alpha_training{ true };
    int thread_size{ 1 };
};

struct MuMDPDETestOptions
{
    double requested_alpha_g{ 0.0 };
    bool alpha_training{ true };
    int thread_size{ 1 };
};

BetaMDPDETestBias RunBetaMDPDETest(
    const LocalTestData & input,
    const BetaMDPDETestOptions & options
);

MuMDPDETestBias RunMuMDPDETest(
    const GroupTestData & input,
    const MuMDPDETestOptions & options
);

} // namespace rhbm_gem::core
