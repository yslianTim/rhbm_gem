#pragma once

#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem {

class LocalPotentialFitState
{
    Eigen::VectorXd m_beta_ols;
    Eigen::VectorXd m_beta_mdpde;
    double m_sigma_square;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> m_data_weight;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> m_data_covariance;
    std::vector<Eigen::VectorXd> m_basis_and_response_entry_list;

public:
    LocalPotentialFitState();
    ~LocalPotentialFitState();

    void SetBetaEstimateOLS(const Eigen::VectorXd & value) { m_beta_ols = value; }
    void SetBetaEstimateMDPDE(const Eigen::VectorXd & value) { m_beta_mdpde = value; }
    void SetSigmaSquare(double value) { m_sigma_square = value; }
    void SetDataWeight(const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & value)
    {
        m_data_weight = value;
    }
    void SetDataCovariance(const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & value)
    {
        m_data_covariance = value;
    }
    void SetBasisAndResponseEntryList(std::vector<Eigen::VectorXd> value);

    int GetBasisAndResponseEntryListSize() const;
    double GetSigmaSquare() const { return m_sigma_square; }
    const Eigen::VectorXd & GetBetaEstimateOLS() const { return m_beta_ols; }
    const Eigen::VectorXd & GetBetaEstimateMDPDE() const { return m_beta_mdpde; }
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & GetDataWeight() const
    {
        return m_data_weight;
    }
    const Eigen::DiagonalMatrix<double, Eigen::Dynamic> & GetDataCovariance() const
    {
        return m_data_covariance;
    }
    const std::vector<Eigen::VectorXd> & GetBasisAndResponseEntryList() const;
};

} // namespace rhbm_gem
