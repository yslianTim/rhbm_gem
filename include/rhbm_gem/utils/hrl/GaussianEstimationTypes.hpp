#pragma once

#include <optional>
#include <string>
#include <rhbm_gem/utils/domain/JointEstimationTypes.hpp>
#include <vector>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem {

enum class FittingStage
{
    First,
    Second
};

enum class EstimateMethod { Unspecified, LocalMDPDE, Peeling, JointComponents };
enum class FittingRole { NotRecorded, Target, Halo };
enum class EvidenceStatus { NotRun, Available, Unavailable, Ineligible };

struct EstimateSource
{
    EstimateMethod method{ EstimateMethod::Unspecified };
    std::string atom_id, component_id, run_id;
    FittingRole role{ FittingRole::NotRecorded };
};

struct StageUncertainty
{
    EvidenceStatus status{ EvidenceStatus::NotRun };
    std::string method, reason;
    // Coordinates are (A, C, log B); includes nuisance-parameter coupling.
    std::optional<Eigen::Matrix3d> covariance;
    std::optional<double> residual_variance;
    std::size_t rank{}, degrees_of_freedom{};
    double rank_threshold{};
};

struct LocalStageEstimate
{
    std::optional<GaussianModel3D> point;
    EstimateSource source;
    std::string reason{ "not-fitted" };
    JointCheckStatus convergence{ JointCheckStatus::NotRun };
    StageUncertainty uncertainty;
};

struct PeelingSampleEstimate
{
    std::optional<double> response;
    std::string reason;
};

struct PostFitPeelingResult
{
    EstimateSource source;
    std::string mode{ "grid-consistent" };
    std::size_t neighbor_count{};
    // One entry per raw sample, including samples without coverage.
    std::vector<PeelingSampleEstimate> samples;
};

struct LocalGaussianResult
{
    double alpha_r{ 0.0 };
    GaussianModel3DWithUncertainty ols{
        GaussianModel3D{ 0.0, 0.0 },
        GaussianModel3DUncertainty{}
    };
    GaussianModel3DWithUncertainty mdpde{
        GaussianModel3D{ 0.0, 0.0 },
        GaussianModel3DUncertainty{}
    };
    std::optional<RHBMBetaEstimateResult> fit_result{};
};

struct GroupGaussianMemberResult
{
    GaussianModel3DWithUncertainty posterior{};
    bool is_outlier{ false };
    double statistical_distance{ 0.0 };
};

struct GroupGaussianMemberInput
{
    LocalPotentialSampleList sample_entries{};
    double alpha_r{ 0.0 };
    GaussianModel3D local_model{};
};

struct GroupGaussianResult
{
    double alpha_g{ 0.0 };
    GaussianModel3D mean{ 0.0, 0.0 };
    GaussianModel3D mdpde{ 0.0, 0.0 };
    GaussianModel3DWithUncertainty prior{
        GaussianModel3D{ 0.0, 0.0 },
        GaussianModel3DUncertainty{}
    };
    std::vector<GroupGaussianMemberResult> member_results{};
};

} // namespace rhbm_gem
