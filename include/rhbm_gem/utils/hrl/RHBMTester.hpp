#pragma once

#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/TestDataFactory.hpp>

namespace rhbm_gem::rhbm_tester
{

struct BetaReplicaResidual
{
    Eigen::VectorXd ols_residual;
    Eigen::VectorXd mdpde_residual;
};

struct MuReplicaResidual
{
    Eigen::VectorXd median_residual;
    Eigen::VectorXd mdpde_residual;
};

struct NeighborhoodReplicaResidual
{
    Eigen::VectorXd no_cut_ols_residual;
    Eigen::VectorXd no_cut_mdpde_residual;
    Eigen::VectorXd cut_mdpde_residual;
    double trained_alpha_r{ 0.0 };
};

struct ResidualStatistics
{
    Eigen::VectorXd mean;
    Eigen::VectorXd sigma;
};

struct ResidualStatisticsSeries
{
    std::vector<ResidualStatistics> requested_alpha;
    ResidualStatistics trained_alpha;
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

bool RunSingleBetaMDPDETest(
    BetaReplicaResidual & result,
    const RHBMMemberDataset & dataset,
    const GaussianParameterVector & gaussian_truth,
    double alpha_r,
    int thread_size = 1
);

bool RunBetaMDPDETest(
    BetaMDPDETestResidual & result,
    const std::vector<double> & alpha_r_list,
    const test_data_factory::RHBMBetaTestInput & test_input,
    int thread_size = 1
);

bool RunMuMDPDETest(
    MuMDPDETestResidual & result,
    const std::vector<double> & alpha_g_list,
    const test_data_factory::RHBMMuTestInput & test_input,
    int thread_size = 1
);

bool RunBetaMDPDEWithNeighborhoodTest(
    NeighborhoodMDPDETestResidual & result,
    const test_data_factory::RHBMNeighborhoodTestInput & test_input,
    int thread_size = 1,
    double angle = 0.0
);

} // namespace rhbm_gem::rhbm_tester
