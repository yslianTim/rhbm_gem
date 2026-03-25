#pragma once

#include <vector>
#include <tuple>
#include <string>
#include <unordered_map>
#include <Eigen/Dense>

namespace rhbm_gem {

class DataObjectBase;

class LocalPotentialEntry
{
    double m_alpha_r;
    std::vector<std::tuple<float, float>> m_distance_and_map_value_list;
    std::tuple<double, double> m_gaus_estimate_ols;
    std::tuple<double, double> m_gaus_estimate_mdpde;
    std::unordered_map<std::string, std::tuple<double, double>> m_gaus_estimate_posterior_map;
    std::unordered_map<std::string, std::tuple<double, double>> m_gaus_variance_posterior_map;
    std::unordered_map<std::string, bool> m_outlier_tag_map;
    std::unordered_map<std::string, double> m_statistical_distance_map;
    Eigen::VectorXd m_beta_ols_tmp;
    Eigen::VectorXd m_beta_mdpde_tmp;
    double m_sigma_square_tmp;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> m_data_weight_tmp;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> m_data_covariance_tmp;
    std::vector<Eigen::VectorXd> m_basis_and_response_entry_list_tmp;

public:
    LocalPotentialEntry();
    ~LocalPotentialEntry();

    void SetAlphaR(double value) { m_alpha_r = value; }
    void SetBetaEstimateOLS(const Eigen::VectorXd & value) { m_beta_ols_tmp = value; }
    void SetBetaEstimateMDPDE(const Eigen::VectorXd & value) { m_beta_mdpde_tmp = value; }
    void SetSigmaSquare(double value) { m_sigma_square_tmp = value; }
    void SetDataWeight(const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & value) { m_data_weight_tmp = value; }
    void SetDataCovariance(const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & value) { m_data_covariance_tmp = value; }
    void AddDistanceAndMapValueList(std::vector<std::tuple<float, float>> && list);
    void AddBasisAndResponseEntryList(std::vector<Eigen::VectorXd> && list);
    void ClearDistanceAndMapValueList();
    void AddGausEstimateOLS(double v0, double v1);
    void AddGausEstimateMDPDE(double v0, double v1);
    void AddGausEstimatePosterior(const std::string & key, double v0, double v1);
    void AddGausVariancePosterior(const std::string & key, double v0, double v1);
    void AddOutlierTag(const std::string & key, bool value);
    void AddStatisticalDistance(const std::string & key, double value);

    int GetDistanceAndMapValueListSize() const;
    int GetBasisAndResponseEntryListSize() const;
    double GetAlphaR() const { return m_alpha_r; }
    double GetSigmaSquare() const { return m_sigma_square_tmp; }
    const Eigen::VectorXd & GetBetaEstimateOLS() const { return m_beta_ols_tmp; }
    const Eigen::VectorXd & GetBetaEstimateMDPDE() const { return m_beta_mdpde_tmp; }
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & GetDataWeight() const { return m_data_weight_tmp; }
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & GetDataCovariance() const { return m_data_covariance_tmp; }
    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList() const;
    const std::vector<Eigen::VectorXd> & GetBasisAndResponseEntryList() const;
    double GetMapValueNearCenter() const;
    double GetMomentZeroEstimate() const;
    double GetMomentTwoEstimate() const;
    double GetBetaEstimateOLS(int par_id) const;
    double GetBetaEstimateMDPDE(int par_id) const;
    double GetGausEstimateOLS(int par_id) const;
    double GetAmplitudeEstimateOLS() const;
    double GetWidthEstimateOLS() const;
    double GetIntensityEstimateOLS() const;
    double GetGausEstimateMDPDE(int par_id) const;
    double GetAmplitudeEstimateMDPDE() const;
    double GetWidthEstimateMDPDE() const;
    double GetIntensityEstimateMDPDE() const;
    double GetGausEstimatePosterior(const std::string & key, int par_id) const;
    double GetAmplitudeEstimatePosterior(const std::string & key) const;
    double GetWidthEstimatePosterior(const std::string & key) const;
    double GetIntensityEstimatePosterior(const std::string & key) const;
    double GetGausVariancePosterior(const std::string & key, int par_id) const;
    double GetAmplitudeVariancePosterior(const std::string & key) const;
    double GetWidthVariancePosterior(const std::string & key) const;
    double GetIntensityVariancePosterior(const std::string & key) const;
    bool GetOutlierTag(const std::string & key) const;
    double GetStatisticalDistance(const std::string & key) const;
    double CalculateQScore(int par_choice) const;

private:
    double CalculateIntensityEstimate(double amplitude, double width) const;
    double CalculateIntensityVariance(double amplitude, double sigma_amplitude, double width, double sigma_width) const;

};

} // namespace rhbm_gem
