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
    int m_member_size;
    int m_replica_size;
    size_t m_sampling_entry_size;
    double m_x_min, m_x_max;
    double m_data_error_sigma;
    double m_outlier_ratio;
    Eigen::VectorXd m_gaus_par_prior;
    Eigen::VectorXd m_gaus_par_sigma;

public:
    HRLModelTester(void) = delete;
    explicit HRLModelTester(int gaus_par_size, int linear_basis_size, int replica_size);
    ~HRLModelTester() = default;

    void SetSamplingEntrySize(size_t sampling_entry_size) { m_sampling_entry_size = sampling_entry_size; }
    void SetDataErrorSigma(double data_error_sigma) { m_data_error_sigma = data_error_sigma; }
    void SetOutlierRatio(double outlier_ratio) { m_outlier_ratio = outlier_ratio; }
    void SetFittingRange(double x_min, double x_max) { m_x_min = x_min; m_x_max = x_max; }
    void SetGausParametersPrior(const Eigen::VectorXd & gaus_par_prior);
    void SetGausParametersSigma(const Eigen::VectorXd & gaus_par_sigma);
    bool RunBetaMDPDETest(
        double alpha_r,
        Eigen::VectorXd & residual_mean_ols,
        Eigen::VectorXd & residual_mean_mdpde,
        Eigen::VectorXd & residual_sigma_ols,
        Eigen::VectorXd & residual_sigma_mdpde,
        int thread_size = 1
    );

private:
    Eigen::MatrixXd BuildRandomGausParameters(int replica_size);
    std::vector<std::tuple<float, float>> BuildRandomGausSamplingEntry(
        size_t sampling_entry_size,
        const Eigen::VectorXd & gaus_par,
        double outlier_ratio = 0.0
    );
    std::vector<Eigen::VectorXd> BuildRandomLinearDataEntry(
        size_t sampling_entry_size,
        const Eigen::VectorXd & gaus_par,
        double error_sigma,
        double outlier_ratio = 0.0
    );
    bool CheckGausParametersDimension(const Eigen::VectorXd & gaus_par);
    Eigen::VectorXd CalculateResidual(const Eigen::VectorXd & estimate, const Eigen::VectorXd & true_value);

};
