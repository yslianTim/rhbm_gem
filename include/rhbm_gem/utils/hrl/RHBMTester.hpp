#pragma once

#include <optional>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/TestDataFactory.hpp>

namespace rhbm_gem::rhbm_tester
{

struct ResidualStatistics
{
    Eigen::VectorXd mean;
    Eigen::VectorXd sigma;
};

struct ResidualStatisticsSeries
{
    std::vector<ResidualStatistics> requested_alpha;
    std::optional<ResidualStatistics> trained_alpha;
    std::optional<double> trained_alpha_average;
};

struct BetaMDPDETestResidual
{
    ResidualStatistics ols;
    ResidualStatisticsSeries mdpde;
};

struct MuMDPDETestResidual
{
    ResidualStatistics median;
    ResidualStatisticsSeries mdpde;
};

bool RunBetaMDPDETest(
    BetaMDPDETestResidual & result,
    const test_data_factory::RHBMBetaTestInput & test_input,
    int thread_size = 1
);

bool RunMuMDPDETest(
    MuMDPDETestResidual & result,
    const test_data_factory::RHBMMuTestInput & test_input,
    int thread_size = 1
);

} // namespace rhbm_gem::rhbm_tester
