#pragma once

#include <vector>
#include <tuple>
#include <string>
#include <unordered_map>
#include <Eigen/Dense>

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
    LocalPotentialEntry(void);
    ~LocalPotentialEntry();

    void SetAlphaR(double value) { m_alpha_r = value; }
    void SetBetaEstimateOLS(const Eigen::VectorXd & value) { m_beta_ols_tmp = value; }
    void SetBetaEstimateMDPDE(const Eigen::VectorXd & value) { m_beta_mdpde_tmp = value; }
    void SetSigmaSquare(double value) { m_sigma_square_tmp = value; }
    void SetDataWeight(const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & value) { m_data_weight_tmp = value; }
    void SetDataCovariance(const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & value) { m_data_covariance_tmp = value; }
    void AddDistanceAndMapValueList(std::vector<std::tuple<float, float>> && list);
    void AddBasisAndResponseEntryList(std::vector<Eigen::VectorXd> && list);
    void ClearDistanceAndMapValueList(void);
    void AddGausEstimateOLS(double v0, double v1);
    void AddGausEstimateMDPDE(double v0, double v1);
    void AddGausEstimatePosterior(const std::string & key, double v0, double v1);
    void AddGausVariancePosterior(const std::string & key, double v0, double v1);
    void AddOutlierTag(const std::string & key, bool value);
    void AddStatisticalDistance(const std::string & key, double value);

    int GetDistanceAndMapValueListSize(void) const;
    int GetBasisAndResponseEntryListSize(void) const;
    double GetAlphaR(void) const { return m_alpha_r; }
    double GetSigmaSquare(void) const { return m_sigma_square_tmp; }
    const Eigen::VectorXd & GetBetaEstimateOLS(void) const { return m_beta_ols_tmp; }
    const Eigen::VectorXd & GetBetaEstimateMDPDE(void) const { return m_beta_mdpde_tmp; }
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & GetDataWeight(void) const { return m_data_weight_tmp; }
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & GetDataCovariance(void) const { return m_data_covariance_tmp; }
    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList(void) const;
    const std::vector<Eigen::VectorXd> & GetBasisAndResponseEntryList(void) const;
    double GetMomentZeroEstimate(void) const;
    double GetMomentTwoEstimate(void) const;
    double GetBetaEstimateOLS(int par_id) const;
    double GetBetaEstimateMDPDE(int par_id) const;
    double GetGausEstimateOLS(int par_id) const;
    double GetAmplitudeEstimateOLS(void) const;
    double GetWidthEstimateOLS(void) const;
    double GetIntensityEstimateOLS(void) const;
    double GetGausEstimateMDPDE(int par_id) const;
    double GetAmplitudeEstimateMDPDE(void) const;
    double GetWidthEstimateMDPDE(void) const;
    double GetIntensityEstimateMDPDE(void) const;
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

private:
    double CalculateIntensityEstimate(double amplitude, double width) const;
    double CalculateIntensityVariance(double amplitude, double sigma_amplitude, double width, double sigma_width) const;

};
