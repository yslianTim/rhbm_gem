#pragma once

#include <optional>
#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem {

class LocalPotentialFitState
{
public:
    struct Dataset
    {
        std::vector<Eigen::VectorXd> basis_and_response_entry_list;
    };

    struct FitResult
    {
        Eigen::VectorXd beta_ols;
        Eigen::VectorXd beta_mdpde;
        double sigma_square{ 0.0 };
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> data_weight;
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> data_covariance;
    };

private:
    std::optional<Dataset> m_dataset;
    std::optional<FitResult> m_fit_result;

public:
    LocalPotentialFitState();
    ~LocalPotentialFitState();

    void SetDataset(std::vector<Eigen::VectorXd> basis_and_response_entry_list);
    void SetFitResult(FitResult value);

    bool HasDataset() const { return m_dataset.has_value(); }
    bool HasFitResult() const { return m_fit_result.has_value(); }
    const Dataset & GetDataset() const;
    const FitResult & GetFitResult() const;
};

} // namespace rhbm_gem
