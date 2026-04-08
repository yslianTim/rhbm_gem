#pragma once

#include <vector>
#include <string>
#include <tuple>
#include <random>
#include <Eigen/Dense>

class HRLModelTester
{
    std::random_device m_random_device;
    std::mt19937 m_generator{ m_random_device() };
    int m_gaus_par_size;
    int m_linear_basis_size;
    int m_replica_size;
    double m_x_min, m_x_max;

public:
    HRLModelTester() = delete;
    explicit HRLModelTester(int gaus_par_size, int linear_basis_size, int replica_size);
    ~HRLModelTester() = default;

    void SetFittingRange(double x_min, double x_max) { m_x_min = x_min; m_x_max = x_max; }
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
        const std::vector<double> & alpha_r_list,
        std::vector<Eigen::VectorXd> & residual_mean_ols_list,
        std::vector<Eigen::VectorXd> & residual_mean_mdpde_list,
        std::vector<Eigen::VectorXd> & residual_sigma_ols_list,
        std::vector<Eigen::VectorXd> & residual_sigma_mdpde_list,
        const Eigen::VectorXd & gaus_true,
        int sampling_entry_size,
        double data_error_sigma = 1.0,
        double neighbor_distance = 2.0,
        int thread_size = 1
    );

private:
    bool CheckGausParametersDimension(const Eigen::VectorXd & gaus_par);
    Eigen::MatrixXd BuildRandomGausParameters(
        int member_size,
        const Eigen::VectorXd & gaus_prior,
        const Eigen::VectorXd & gaus_sigma,
        const Eigen::VectorXd & outlier_prior,
        const Eigen::VectorXd & outlier_sigma,
        double outlier_ratio = 0.0
    );
    Eigen::MatrixXd BuildBetaMatrix(const Eigen::MatrixXd & gaus_array);
    std::vector<std::tuple<float, float>> BuildRandomGausSamplingEntry(
        size_t sampling_entry_size,
        const Eigen::VectorXd & gaus_par,
        double outlier_ratio = 0.0
    );
    std::vector<std::tuple<float, float>> BuildRandomGausSamplingEntryWithNeighborhood(
        size_t sampling_entry_size,
        const Eigen::VectorXd & gaus_par,
        double neighbor_distance = 2.0
    );
    std::vector<Eigen::VectorXd> BuildRandomLinearDataEntry(
        size_t sampling_entry_size,
        const Eigen::VectorXd & gaus_par,
        double error_sigma,
        double outlier_ratio = 0.0
    );
    std::vector<Eigen::VectorXd> BuildRandomLinearDataEntryWithNeighborhood(
        size_t sampling_entry_size,
        const Eigen::VectorXd & gaus_par,
        double error_sigma,
        double neighbor_distance = 2.0
    );
    Eigen::VectorXd CalculateNormalizedResidual(
        const Eigen::VectorXd & estimate,
        const Eigen::VectorXd & true_value
    );

};
