#pragma once

#include <vector>
#include <string>
#include <tuple>
#include <Eigen/Dense>

class HRLModelTester
{
    int m_model_par_size;
    int m_linear_basis_size;
    int m_member_size;
    size_t m_sampling_entry_size;
    double m_x_min, m_x_max;
    double m_data_error_sigma;
    Eigen::VectorXd m_model_par_prior;
    Eigen::VectorXd m_model_par_sigma;
    std::vector<Eigen::VectorXd> m_model_par_local_list;
    std::vector<Eigen::VectorXd> m_model_par_0_list;
    std::vector<std::vector<std::tuple<double, double>>> m_sampling_entries_list;
    std::vector<std::tuple<std::vector<Eigen::VectorXd>, std::string>> m_data_array;
    std::vector<Eigen::VectorXd> m_local_estimate_residual_ols_list;
    std::vector<Eigen::VectorXd> m_local_estimate_residual_mdpde_list;
    std::vector<Eigen::VectorXd> m_local_estimate_residual_posterior_list;
    std::vector<Eigen::VectorXd> m_local_estimate_deviation_ols_list;
    std::vector<Eigen::VectorXd> m_local_estimate_deviation_mdpde_list;
    std::vector<Eigen::VectorXd> m_local_estimate_deviation_posterior_list;
    Eigen::VectorXd m_local_estimate_mean_ols;
    Eigen::VectorXd m_local_estimate_mean_mdpde;
    Eigen::VectorXd m_local_estimate_mean_posterior;
    Eigen::VectorXd m_group_estimate_residual_mean;
    Eigen::VectorXd m_group_estimate_residual_mdpde;
    Eigen::VectorXd m_group_estimate_residual_hrl;

public:
    HRLModelTester(void) = delete;
    explicit HRLModelTester(int model_par_size, int linear_basis_size, int member_size);
    ~HRLModelTester() = default;

    void SetSamplingEntrySize(size_t sampling_entry_size) { m_sampling_entry_size = sampling_entry_size; }
    void SetDataErrorSigma(double data_error_sigma) { m_data_error_sigma = data_error_sigma; }
    void SetFittingRange(double x_min, double x_max) { m_x_min = x_min; m_x_max = x_max; }
    void SetModelParametersPrior(const Eigen::VectorXd & model_par_prior);
    void SetModelParametersSigma(const Eigen::VectorXd & model_par_sigma);
    bool RunTest(double alpha_r, double alpha_g);
    Eigen::VectorXd GetGroupEstimateResidualMean(void) const { return m_group_estimate_residual_mean; }
    Eigen::VectorXd GetGroupEstimateResidualMDPDE(void) const { return m_group_estimate_residual_mdpde; }
    Eigen::VectorXd GetGroupEstimateResidualHRL(void) const { return m_group_estimate_residual_hrl; }
    Eigen::VectorXd GetLocalEstimateMeanOLS(void) const { return m_local_estimate_mean_ols; }
    Eigen::VectorXd GetLocalEstimateMeanMDPDE(void) const { return m_local_estimate_mean_mdpde; }
    Eigen::VectorXd GetLocalEstimateMeanPosterior(void) const { return m_local_estimate_mean_posterior; }
    Eigen::VectorXd GetLocalEstimateResidualOLS(size_t index) const { return m_local_estimate_residual_ols_list.at(index); }
    Eigen::VectorXd GetLocalEstimateResidualMDPDE(size_t index) const { return m_local_estimate_residual_mdpde_list.at(index); }
    Eigen::VectorXd GetLocalEstimateResidualPosterior(size_t index) const { return m_local_estimate_residual_posterior_list.at(index); }
    Eigen::VectorXd GetLocalEstimateDeviationOLS(size_t index) const { return m_local_estimate_deviation_ols_list.at(index); }
    Eigen::VectorXd GetLocalEstimateDeviationMDPDE(size_t index) const { return m_local_estimate_deviation_mdpde_list.at(index); }
    Eigen::VectorXd GetLocalEstimateDeviationPosterior(size_t index) const { return m_local_estimate_deviation_posterior_list.at(index); }
    const std::vector<Eigen::VectorXd> & GetLocalEstimateResidualOLSList(void) const { return m_local_estimate_residual_ols_list; }
    const std::vector<Eigen::VectorXd> & GetLocalEstimateResidualMDPDEList(void) const { return m_local_estimate_residual_mdpde_list; }
    const std::vector<Eigen::VectorXd> & GetLocalEstimateResidualPosteriorList(void) const { return m_local_estimate_residual_posterior_list; }
    const std::vector<Eigen::VectorXd> & GetLocalEstimateDeviationOLSList(void) const { return m_local_estimate_deviation_ols_list; }
    const std::vector<Eigen::VectorXd> & GetLocalEstimateDeviationMDPDEList(void) const { return m_local_estimate_deviation_mdpde_list; }
    const std::vector<Eigen::VectorXd> & GetLocalEstimateDeviationPosteriorList(void) const { return m_local_estimate_deviation_posterior_list; }

private:
    void BuildDataArray(void);
    std::vector<std::tuple<double, double>> BuildRandomSamplingEntry(
        size_t sampling_entry_size, const Eigen::VectorXd & model_par);
    bool CheckModelParametersDimension(const Eigen::VectorXd & model_par);
    Eigen::VectorXd CalculateResidual(const Eigen::VectorXd & estimate, const Eigen::VectorXd & true_value);
    Eigen::VectorXd CalculateMoment(const std::vector<std::tuple<double, double>> & sampling_entries) const;

};
