#pragma once

#include <vector>
#include <tuple>
#include <string>
#include <unordered_map>
#include <Eigen/Dense>

#include <rhbm_gem/data/object/GaussianStatistics.hpp>

namespace rhbm_gem {

class LocalPotentialEntry
{
    double m_alpha_r;
    std::vector<std::tuple<float, float>> m_distance_and_map_value_list;
    GaussianEstimate m_gaus_estimate_ols;
    GaussianEstimate m_gaus_estimate_mdpde;
    std::unordered_map<std::string, GaussianPosterior> m_gaus_posterior_map;
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
    void SetEstimateOLS(const GaussianEstimate & estimate) { m_gaus_estimate_ols = estimate; }
    void SetEstimateMDPDE(const GaussianEstimate & estimate) { m_gaus_estimate_mdpde = estimate; }
    void SetPosterior(const std::string & key, const GaussianPosterior & posterior);
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
    const GaussianEstimate & GetEstimateOLS() const { return m_gaus_estimate_ols; }
    const GaussianEstimate & GetEstimateMDPDE() const { return m_gaus_estimate_mdpde; }
    const GaussianPosterior & GetPosterior(const std::string & key) const;
    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList() const;
    const std::vector<Eigen::VectorXd> & GetBasisAndResponseEntryList() const;
    double GetMapValueNearCenter() const;
    double GetMomentZeroEstimate() const;
    double GetMomentTwoEstimate() const;
    double GetBetaEstimateOLS(int par_id) const;
    double GetBetaEstimateMDPDE(int par_id) const;
    bool GetOutlierTag(const std::string & key) const;
    double GetStatisticalDistance(const std::string & key) const;
    double CalculateQScore(int par_choice) const;
};

} // namespace rhbm_gem
