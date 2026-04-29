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
    std::vector<double> requested_alpha_r_list;
    bool alpha_training{ true };
};

struct RHBMMuTestInput
{
    Eigen::VectorXd gaus_true;
    std::vector<Eigen::MatrixXd> replica_beta_matrices;
    std::vector<double> requested_alpha_g_list;
    bool alpha_training{ true };
};

struct RHBMNeighborhoodTestInput
{
    RHBMBetaTestInput no_cut_input;
    RHBMBetaTestInput cut_input;
    std::vector<LocalPotentialSampleList> sampling_summaries;
};

struct TestDataBuildOptions
{
    linearization_service::LinearizationSpec linearization_spec{
        linearization_service::LinearizationSpec::AtomDecode()
    };
    linearization_service::LinearizationRange fitting_range{ 0.0, 1.0 };
};

struct BetaScenario
{
    GaussianModel3D gaus_true;
    int sampling_entry_size{ 1 };
    double data_error_sigma{ 1.0 };
    double outlier_ratio{ 0.0 };
    int replica_size{ 1 };
    std::optional<std::uint32_t> random_seed{};
    std::vector<double> requested_alpha_r_list{};
    bool alpha_training{ true };
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
    std::vector<double> requested_alpha_g_list{};
    bool alpha_training{ true };
};

struct NeighborhoodScenario
{
    GaussianModel3D gaus_true;
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

RHBMBetaTestInput BuildBetaTestInput(
    const BetaScenario & scenario,
    const TestDataBuildOptions & options = {});

RHBMMuTestInput BuildMuTestInput(
    const MuScenario & scenario,
    const TestDataBuildOptions & options = {});

RHBMNeighborhoodTestInput BuildNeighborhoodTestInput(
    const NeighborhoodScenario & scenario,
    const TestDataBuildOptions & options = {});

} // namespace rhbm_gem::test_data_factory
