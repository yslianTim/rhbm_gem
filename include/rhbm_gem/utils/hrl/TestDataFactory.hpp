#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem::test_data_factory
{
struct LocalTestInput
{
    GaussianModel3D gaus_true;
    std::vector<LocalPotentialSampleList> replica_sampling_entries;
    std::vector<double> requested_alpha_r_list;
    bool alpha_training{ true };
    LocalGaussianFitModel local_fit_model{ LocalGaussianFitModel::LogQuadratic };
};

struct GroupTestInput
{
    GaussianModel3D gaus_true;
    std::vector<std::vector<LocalPotentialSampleList>> replica_member_sampling_entries;
    std::vector<double> requested_alpha_g_list;
    bool alpha_training{ true };
};

struct AtomModelTestInput
{
    LocalTestInput no_cut_input;
    LocalTestInput cut_input;
    std::vector<LocalPotentialSampleList> sampling_summaries;
};

struct LocalScenario
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

struct GroupScenario
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

struct AtomModelScenario
{
    Spot spot{ Spot::UNK };
    GaussianModel3D gaus_true;
    double data_error_sigma{ 1.0 };
    double rejected_angle{ 0.0 };
    bool include_sampling_summary{ false };
    int replica_size{ 1 };
    std::optional<std::uint32_t> random_seed{};
    std::vector<double> requested_alpha_r_list{};
    bool alpha_training{ true };
    LocalGaussianFitModel local_fit_model{ LocalGaussianFitModel::LogQuadratic };
};

LocalTestInput BuildLocalTestInput(const LocalScenario & scenario);
GroupTestInput BuildGroupTestInput(const GroupScenario & scenario);
AtomModelTestInput BuildAtomModelTestInput(const AtomModelScenario & scenario);

} // namespace rhbm_gem::test_data_factory
