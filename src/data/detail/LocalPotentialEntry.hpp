#pragma once

#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/data/object/GaussianStatistics.hpp>

namespace rhbm_gem {

struct LocalPotentialAnnotation
{
    GaussianPosterior posterior{};
    bool is_outlier{ false };
    double statistical_distance{ 0.0 };
};

class LocalPotentialEntry
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
    double m_alpha_r;
    std::vector<std::tuple<float, float>> m_distance_and_map_value_list;
    GaussianEstimate m_gaus_estimate_ols;
    GaussianEstimate m_gaus_estimate_mdpde;
    std::unordered_map<std::string, LocalPotentialAnnotation> m_annotation_map;
    std::optional<Dataset> m_dataset;
    std::optional<FitResult> m_fit_result;

public:
    LocalPotentialEntry();
    ~LocalPotentialEntry();

    void SetAlphaR(double value) { m_alpha_r = value; }
    void SetDistanceAndMapValueList(std::vector<std::tuple<float, float>> value);
    void SetEstimateOLS(const GaussianEstimate & estimate) { m_gaus_estimate_ols = estimate; }
    void SetEstimateMDPDE(const GaussianEstimate & estimate) { m_gaus_estimate_mdpde = estimate; }
    void SetDataset(std::vector<Eigen::VectorXd> basis_and_response_entry_list);
    void SetFitResult(FitResult value);
    void SetAnnotation(const std::string & key, LocalPotentialAnnotation annotation);
    void ClearTransientFitState();

    int GetDistanceAndMapValueListSize() const;
    double GetAlphaR() const { return m_alpha_r; }
    bool HasDataset() const { return m_dataset.has_value(); }
    bool HasFitResult() const { return m_fit_result.has_value(); }
    const GaussianEstimate & GetEstimateOLS() const { return m_gaus_estimate_ols; }
    const GaussianEstimate & GetEstimateMDPDE() const { return m_gaus_estimate_mdpde; }
    const Dataset & GetDataset() const;
    const FitResult & GetFitResult() const;
    LocalPotentialAnnotation * FindAnnotation(const std::string & key);
    const LocalPotentialAnnotation * FindAnnotation(const std::string & key) const;
    const std::unordered_map<std::string, LocalPotentialAnnotation> & Annotations() const;
    const std::vector<std::tuple<float, float>> & GetDistanceAndMapValueList() const;
    double GetMapValueNearCenter() const;
    double GetMomentZeroEstimate() const;
    double GetMomentTwoEstimate() const;
    double CalculateQScore(int par_choice) const;
};

} // namespace rhbm_gem
