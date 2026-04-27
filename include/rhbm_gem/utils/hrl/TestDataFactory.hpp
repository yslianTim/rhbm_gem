#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem::test_data_factory
{

struct RHBMBetaTestInput
{
    Eigen::VectorXd gaus_true;
    std::vector<RHBMMemberDataset> replica_datasets;
};

struct RHBMMuTestInput
{
    Eigen::VectorXd gaus_true;
    std::vector<Eigen::MatrixXd> replica_beta_matrices;
};

struct RHBMNeighborhoodTestInput
{
    Eigen::VectorXd gaus_true;
    std::vector<RHBMMemberDataset> no_cut_datasets;
    std::vector<RHBMMemberDataset> cut_datasets;
    std::vector<LocalPotentialSampleList> sampling_summaries;
};

class TestDataFactory
{
public:
    struct BetaScenario
    {
        Eigen::VectorXd gaus_true;
        int sampling_entry_size{ 1 };
        double data_error_sigma{ 1.0 };
        double outlier_ratio{ 0.0 };
        int replica_size{ 1 };
        std::optional<std::uint32_t> random_seed{};
    };

    struct MuScenario
    {
        int member_size{ 1 };
        Eigen::VectorXd gaus_prior;
        Eigen::VectorXd gaus_sigma;
        Eigen::VectorXd outlier_prior;
        Eigen::VectorXd outlier_sigma;
        double outlier_ratio{ 0.0 };
        int replica_size{ 1 };
        std::optional<std::uint32_t> random_seed{};
    };

    struct NeighborhoodScenario
    {
        Eigen::VectorXd gaus_true;
        int sampling_entry_size{ 1 };
        double data_error_sigma{ 1.0 };
        double radius_min{ 0.0 };
        double radius_max{ 1.0 };
        double neighbor_distance{ 2.0 };
        size_t neighbor_count{ 1 };
        double rejected_angle{ 0.0 };
        bool include_sampling_summary{ false };
        double summary_radius_min{ 0.0 };
        double summary_radius_max{ 1.0 };
        int replica_size{ 1 };
        std::optional<std::uint32_t> random_seed{};
    };

private:
    linearization_service::LinearizationSpec m_linearization_spec;
    double m_fit_range_min;
    double m_fit_range_max;

public:
    TestDataFactory() = delete;
    explicit TestDataFactory(
        linearization_service::LinearizationSpec linearization_spec);
    ~TestDataFactory() = default;

    void SetFittingRange(double x_min, double x_max);

    RHBMBetaTestInput BuildBetaTestInput(const BetaScenario & scenario) const;
    RHBMMuTestInput BuildMuTestInput(const MuScenario & scenario) const;
    RHBMNeighborhoodTestInput BuildNeighborhoodTestInput(
        const NeighborhoodScenario & scenario) const;
};

} // namespace rhbm_gem::test_data_factory
