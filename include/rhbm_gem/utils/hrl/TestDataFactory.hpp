#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem::test_data_factory
{

enum class AtomNeighborType : int
{
    None = 0, O = 1, N = 2, C = 3, CA = 4
};

struct RHBMBetaTestInput
{
    GaussianModel3D gaus_true;
    std::vector<LocalPotentialSampleList> replica_sampling_entries;
    std::vector<RHBMMemberDataset> replica_datasets;
    std::vector<double> requested_alpha_r_list;
    bool alpha_training{ true };
    LocalGaussianFitModel local_fit_model{ LocalGaussianFitModel::LogQuadratic };
};

struct RHBMMuTestInput
{
    GaussianModel3D gaus_true;
    std::vector<std::vector<LocalPotentialSampleList>> replica_member_sampling_entries;
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
    double fit_range_min{ 0.0 };
    double fit_range_max{ 1.0 };
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
    LocalGaussianFitModel local_fit_model{ LocalGaussianFitModel::LogQuadratic };
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

struct AtomNeighborhoodScenario
{
    AtomNeighborType neighbor_type{ AtomNeighborType::None };
    GaussianModel3D gaus_true;
    int sampling_entry_size{ 1 };
    double data_error_sigma{ 1.0 };
    double radius_min{ 0.0 };
    double radius_max{ 1.0 };
    double rejected_angle{ 0.0 };
    bool include_sampling_summary{ false };
    double summary_radius_min{ 0.0 };
    double summary_radius_max{ 1.0 };
    int replica_size{ 1 };
    std::optional<std::uint32_t> random_seed{};
    std::vector<double> requested_alpha_r_list{};
    bool alpha_training{ true };
    LocalGaussianFitModel local_fit_model{ LocalGaussianFitModel::LogQuadratic };
};

RHBMBetaTestInput BuildBetaTestInput(
    const BetaScenario & scenario,
    const TestDataBuildOptions & options = {});

RHBMMuTestInput BuildMuTestInput(const MuScenario & scenario);

RHBMNeighborhoodTestInput BuildAtomNeighborhoodTestInput(
    const AtomNeighborhoodScenario & scenario,
    const TestDataBuildOptions & options = {});

} // namespace rhbm_gem::test_data_factory
