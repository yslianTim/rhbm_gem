#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include <rhbm_gem/utils/domain/JointEstimationTypes.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem {

enum class EstimateMethod { Unspecified = 0, LocalMDPDE = 1, JointComponents = 3 };
enum class FittingRole { NotRecorded, Target, Halo };

struct EstimateSource
{
    EstimateMethod method{ EstimateMethod::Unspecified };
    std::string atom_id, component_id, run_id;
    FittingRole role{ FittingRole::NotRecorded };
    bool fixed_from_initial_stage{};
};

struct StageUncertainty
{
    EvidenceStatus status{ EvidenceStatus::NotRun };
    std::string method, reason;
    // Coordinates are (A, C, log tau); includes nuisance-parameter coupling.
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

// An Original-grid voxel observation, independent of InitialStage interpolation.
struct PostJointGridObservation
{
    std::size_t original_row_index{};
    std::size_t global_grid_index{};
    double squared_distance{};
    double map_value{}; // Current stored voxel value; JSON key remains "raw_response".
    std::array<double,3> position{};
};

struct PeelingSampleEstimate
{
    std::optional<double> response;
    std::string reason;
    PostJointGridObservation observation;
};

struct PostFitPeelingResult
{
    EstimateSource source;
    std::string mode{ "original-grid-v1" };
    std::size_t neighbor_count{};
    // Original-grid observations within 2 A, sorted by Original row identity.
    std::vector<PeelingSampleEstimate> samples;
};

struct PeelingRatioAssessment
{
    std::optional<double> ratio;
    std::size_t total_samples{};
    std::size_t valid_samples{};
    std::size_t peeling_samples{};
    double raw_response_sum{};
    double peeled_response_sum{};
    double valid_raw_response_sum{};
    std::string reason;
};

struct LocalFitDiagnostics
{
    RHBMEstimationStatus status{ RHBMEstimationStatus::SUCCESS };
    double sigma_square{};
    RHBMBetaDiagnostics iterations;
    std::optional<RHBMEndpointRefinementDiagnostics> refinement;
    std::optional<RHBMFailedOnlyRecoveryDiagnostics> failed_only_recovery;
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
    std::optional<LocalFitDiagnostics> diagnostics{};
    bool uncertainty_recorded{ true };
};

struct JointFixedChargeProfile
{
    double coefficient{};
    double width{};
};

struct PostJointLocalSampleFailure
{
    std::size_t original_row_index{};
    double distance{};
    std::string reason;
    std::size_t global_grid_index{};
};

// Persisted Post-Joint Local estimation contract, in angstroms.
inline constexpr double PostJointLocalFittingRadius{1.0};

struct PostJointLocalResult
{
    EstimateSource source;
    EvidenceStatus status{ EvidenceStatus::NotRun };
    std::string reason{ "not-fitted" };
    LocalGaussianResult fit;
    JointFixedChargeProfile fixed_charge;
    std::optional<RHBMInformation> information;
    std::size_t sample_count{};
    std::size_t effective_sample_count{};
    std::size_t corrected_sample_count{}, positive_sample_count{};
    double fitting_radius{PostJointLocalFittingRadius};
    std::vector<PostJointLocalSampleFailure> sample_failures;
};

struct GroupGaussianMemberResult
{
    GaussianModel3DWithUncertainty posterior{};
    bool is_outlier{ false };
    double statistical_distance{ 0.0 };
    std::string evidence_source_id;
    bool charge_inferred{ true };
    std::optional<Eigen::Matrix2d> parameter_covariance;
    std::optional<JointFixedChargeProfile> fixed_charge;
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

struct GroupParameterSummary
{
    EvidenceStatus status{ EvidenceStatus::NotRun };
    std::string reason, source_id;
    std::string correlation_policy{ "block-diagonal-by-atom" };
    std::size_t point_count{}, eligible_count{}, excluded_count{};
    std::optional<GaussianModel3D> descriptive_mean;
    std::optional<double> descriptive_charge_median;
    std::optional<GroupGaussianResult> inference;
    std::vector<int> member_ids;
};

} // namespace rhbm_gem
