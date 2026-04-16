#pragma once

#include <vector>
#include <string>
#include <random>
#include <Eigen/Dense>

#include <rhbm_gem/utils/math/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/SimulationDataGenerator.hpp>

class SimulationDataGenerator;

class HRLModelTester
{
    int m_gaus_par_size;
    int m_linear_basis_size;
    int m_replica_size;
    SimulationDataGenerator m_data_generator;

public:
    HRLModelTester() = delete;
    explicit HRLModelTester(int gaus_par_size, int linear_basis_size, int replica_size);
    ~HRLModelTester() = default;

    void SetFittingRange(double x_min, double x_max);
    LocalPotentialSampleList RunDataEntryWithNeighborhoodTest(
        const Eigen::VectorXd & gaus_true,
        int sampling_entry_size,
        double radius_min = 0.0,
        double radius_max = 1.0,
        double neighbor_distance = 2.0,
        size_t neighbor_count = 1,
        double angle = 0.0
    );
    bool RunBetaMDPDETest(
        const std::vector<double> & alpha_r_list,
        std::vector<Eigen::VectorXd> & residual_mean_ols_list,
        std::vector<Eigen::VectorXd> & residual_mean_mdpde_list,
        std::vector<Eigen::VectorXd> & residual_sigma_ols_list,
        std::vector<Eigen::VectorXd> & residual_sigma_mdpde_list,
        const Eigen::VectorXd & gaus_true,
        int sampling_entry_size,
        double data_error_sigma = 1.0,
        double outlier_ratio = 0.0,
        int thread_size = 1
    );
    bool RunMuMDPDETest(
        const std::vector<double> & alpha_g_list,
        std::vector<Eigen::VectorXd> & residual_mean_median_list,
        std::vector<Eigen::VectorXd> & residual_mean_mdpde_list,
        std::vector<Eigen::VectorXd> & residual_sigma_median_list,
        std::vector<Eigen::VectorXd> & residual_sigma_mdpde_list,
        int member_size,
        const Eigen::VectorXd & gaus_prior,
        const Eigen::VectorXd & gaus_sigma,
        const Eigen::VectorXd & outlier_prior,
        const Eigen::VectorXd & outlier_sigma,
        double outlier_ratio = 0.0,
        int thread_size = 1
    );
    bool RunBetaMDPDEWithNeighborhoodTest(
        std::vector<Eigen::VectorXd> & residual_mean_list,
        std::vector<Eigen::VectorXd> & residual_sigma_list,
        const Eigen::VectorXd & gaus_true,
        double & training_alpha_r_average,
        int sampling_entry_size,
        double data_error_sigma = 1.0,
        double neighbor_distance = 2.0,
        size_t neighbor_count = 1,
        int thread_size = 1,
        double angle = 0.0
    );

private:
    bool CheckGausParametersDimension(const Eigen::VectorXd & gaus_par);
    Eigen::MatrixXd BuildRandomGausParameters(
        int member_size,
        const Eigen::VectorXd & gaus_prior,
        const Eigen::VectorXd & gaus_sigma,
        const Eigen::VectorXd & outlier_prior,
        const Eigen::VectorXd & outlier_sigma,
        double outlier_ratio,
        std::mt19937 & generator
    );
    Eigen::MatrixXd BuildBetaMatrix(const Eigen::MatrixXd & gaus_array);
    Eigen::VectorXd CalculateNormalizedResidual(
        const Eigen::VectorXd & estimate,
        const Eigen::VectorXd & true_value
    );

};
