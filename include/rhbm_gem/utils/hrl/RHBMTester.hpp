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

struct NeighborhoodMDPDETestResidual
{
    ResidualStatistics no_cut_ols;
    ResidualStatistics no_cut_mdpde;
    ResidualStatistics cut_mdpde;
    double trained_alpha_r_average{ 0.0 };
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

bool RunBetaMDPDEWithNeighborhoodTest(
    NeighborhoodMDPDETestResidual & result,
    const test_data_factory::RHBMNeighborhoodTestInput & test_input,
    int thread_size = 1
);

} // namespace rhbm_gem::rhbm_tester
