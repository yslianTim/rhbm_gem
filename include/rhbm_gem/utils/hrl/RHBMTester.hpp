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

bool RunSingleBetaMDPDETest(
    int gaussian_parameter_size,
    BetaReplicaResidual & result,
    const RHBMMemberDataset & dataset,
    const Eigen::VectorXd & gaussian_truth,
    double alpha_r,
    int thread_size = 1
);

bool RunBetaMDPDETest(
    int gaussian_parameter_size,
    const std::vector<double> & alpha_r_list,
    std::vector<Eigen::VectorXd> & residual_mean_ols_list,
    std::vector<Eigen::VectorXd> & residual_mean_mdpde_list,
    std::vector<Eigen::VectorXd> & residual_sigma_ols_list,
    std::vector<Eigen::VectorXd> & residual_sigma_mdpde_list,
    const rhbm_gem::test_data_factory::RHBMBetaTestInput & test_input,
    int thread_size = 1
);

bool RunMuMDPDETest(
    int gaussian_parameter_size,
    const std::vector<double> & alpha_g_list,
    std::vector<Eigen::VectorXd> & residual_mean_median_list,
    std::vector<Eigen::VectorXd> & residual_mean_mdpde_list,
    std::vector<Eigen::VectorXd> & residual_sigma_median_list,
    std::vector<Eigen::VectorXd> & residual_sigma_mdpde_list,
    const rhbm_gem::test_data_factory::RHBMMuTestInput & test_input,
    int thread_size = 1
);

bool RunBetaMDPDEWithNeighborhoodTest(
    int gaussian_parameter_size,
    std::vector<Eigen::VectorXd> & residual_mean_list,
    std::vector<Eigen::VectorXd> & residual_sigma_list,
    const rhbm_gem::test_data_factory::RHBMNeighborhoodTestInput & test_input,
    double & training_alpha_r_average,
    int thread_size = 1,
    double angle = 0.0
);

} // namespace rhbm_gem::rhbm_tester
