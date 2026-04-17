#pragma once

#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/HRLModelTestDataFactory.hpp>

class HRLModelTester
{
    int m_gaus_par_size;

public:
    HRLModelTester() = delete;
    explicit HRLModelTester(int gaus_par_size);
    ~HRLModelTester() = default;

    bool RunBetaMDPDETest(
        const std::vector<double> & alpha_r_list,
        std::vector<Eigen::VectorXd> & residual_mean_ols_list,
        std::vector<Eigen::VectorXd> & residual_mean_mdpde_list,
        std::vector<Eigen::VectorXd> & residual_sigma_ols_list,
        std::vector<Eigen::VectorXd> & residual_sigma_mdpde_list,
        const HRLBetaTestInput & test_input,
        int thread_size = 1
    );
    bool RunMuMDPDETest(
        const std::vector<double> & alpha_g_list,
        std::vector<Eigen::VectorXd> & residual_mean_median_list,
        std::vector<Eigen::VectorXd> & residual_mean_mdpde_list,
        std::vector<Eigen::VectorXd> & residual_sigma_median_list,
        std::vector<Eigen::VectorXd> & residual_sigma_mdpde_list,
        const HRLMuTestInput & test_input,
        int thread_size = 1
    );
    bool RunBetaMDPDEWithNeighborhoodTest(
        std::vector<Eigen::VectorXd> & residual_mean_list,
        std::vector<Eigen::VectorXd> & residual_sigma_list,
        const HRLNeighborhoodTestInput & test_input,
        double & training_alpha_r_average,
        int thread_size = 1,
        double angle = 0.0
    );

private:
    bool CheckGausParametersDimension(const Eigen::VectorXd & gaus_par);
};
