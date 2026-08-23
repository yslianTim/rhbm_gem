#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/utils/algorithm/Convergence.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

namespace rhbm_gem {
class AtomObject;
}

namespace rhbm_gem::core::detail {

constexpr std::size_t kLogPeakHeightChangeIndex{ 0 };
constexpr std::size_t kLogWidthChangeIndex{ 1 };
constexpr std::size_t kOffsetToPeakRatioChangeIndex{ 2 };
constexpr std::size_t kTransformedChangeSize{ 3 };

std::optional<Eigen::Vector3d> EncodeTransformedCoordinates(const GaussianModel3D & model);

std::optional<GaussianModel3D> DecodeTransformedCoordinates(const Eigen::Vector3d & coordinates);

bool IsValidSecondStageGaussianModel(const GaussianModel3D & model);

std::optional<GaussianModel3D> BuildGaussianParameterMedian(
    const std::vector<GaussianModel3D> & model_list);

std::vector<GaussianModel3D> BuildGroupMedianModelList(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & model_list);

std::vector<double> BuildGroupMedianOffsetList(
    const std::vector<std::size_t> & group_id_by_atom_position,
    const std::vector<GaussianModel3D> & model_list);

bool TryBuildSharedOffsetDampedModelList(
    const std::vector<GaussianModel3D> & previous_model_list,
    const std::vector<GaussianModel3D> & raw_model_list,
    const std::vector<double> & previous_shared_offset_list,
    const std::vector<double> & raw_shared_offset_list,
    double damping,
    std::vector<GaussianModel3D> & candidate_model_list);

struct SharedOffsetResponse
{
    double response{ 0.0 };
    Eigen::Vector2d shape_jacobian{ Eigen::Vector2d::Zero() };
    double offset_jacobian{ 0.0 };
};

std::optional<SharedOffsetResponse> EvaluateSharedOffsetResponse(
    const GaussianModel3D & model,
    double distance);

struct TransformedModelInvariants
{
    GaussianModel3D model{};
    double peak_height{ 0.0 };
};

std::optional<TransformedModelInvariants> BuildTransformedModelInvariants(const GaussianModel3D & model);

std::optional<SharedOffsetResponse> EvaluateSharedOffsetResponse(
    const TransformedModelInvariants & invariants,
    double distance);

std::optional<Eigen::Vector3d> EvaluateTransformedJacobian(
    const TransformedModelInvariants & invariants,
    double distance);


struct LocalGaussianDesignTemplate
{
    std::size_t source_sample_count{ 0 };
    std::vector<std::size_t> source_sample_index_list{};
    std::vector<double> distance_list{};
    RHBMDesignMatrix design_matrix{};
};

RHBMExecutionOptions MakeExecutionOptions(const FitOptions & options);

LocalPotentialSampleList BuildSamplesForZeroOffsetGaussianFit(
    const LocalPotentialSampleList & sample_entries,
    const GaussianModel3D & model);

LocalGaussianDesignTemplate BuildLocalGaussianDesignTemplate(
    const LocalPotentialSampleList & sample_entries,
    double range_min,
    double range_max);

RHBMMemberDataset BuildLocalGaussianPreparedDataset(
    const LocalGaussianDesignTemplate & design_template,
    const std::vector<double> & sample_response_list,
    const GaussianModel3D & offset_model);

LocalGaussianResult EstimateLocalGaussianPrepared(
    const LocalGaussianDesignTemplate & design_template,
    const std::vector<double> & sample_response_list,
    double alpha_r,
    const FitOptions & options,
    const GaussianModel3D & offset_model);


using ClusterKey = std::vector<std::size_t>;
using ResidueKey = std::pair<std::string, int>;

struct SampleRef
{
    std::size_t atom_index{ 0 };
    std::size_t sample_index{ 0 };

    friend auto operator<=>(const SampleRef &, const SampleRef &) = default;
};

struct NeighborAtomSample
{
    bool is_selected{ true };
    std::size_t atom_index{ 0 };
    double distance{ 0.0 };
};

struct UnselectedAtomContributor
{
    int atom_serial_id{ 0 };
    std::optional<std::size_t> selected_group_id{};
    GaussianModel3DWithUncertainty initial_seed{};
};

struct AtomContext
{
    const AtomObject * atom{ nullptr };
    std::size_t group_id{ 0 };
    LocalPotentialSampleList raw_sampling_entries{};
    LocalGaussianResult initial_result{};
    std::optional<GaussianModel3DWithUncertainty> group_prior{};
    std::vector<NeighborAtomSample> neighbor_atom_sample_list{};
    std::vector<std::size_t> neighbor_atom_sample_offset_list{};
    LocalGaussianDesignTemplate refit_design_template{};
    double alpha_r{ 0.0 };
    int neighbor_count_for_peeling{ 0 };

    std::span<const NeighborAtomSample> Neighbors(std::size_t sample_index) const
    {
        const auto begin{ neighbor_atom_sample_offset_list.at(sample_index) };
        const auto end{ neighbor_atom_sample_offset_list.at(sample_index + 1) };
        return std::span<const NeighborAtomSample>{ neighbor_atom_sample_list }
            .subspan(begin, end - begin);
    }
};

struct SecondStageContext
{
    std::vector<AtomContext> selected_atom_list{};
    std::vector<UnselectedAtomContributor> unselected_atom_list{};
    std::vector<std::vector<std::size_t>> selected_atom_index_list_by_group{};

    std::size_t size() const { return selected_atom_list.size(); }
    AtomContext & at(std::size_t index) { return selected_atom_list.at(index); }
    const AtomContext & at(std::size_t index) const { return selected_atom_list.at(index); }
    auto begin() { return selected_atom_list.begin(); }
    auto end() { return selected_atom_list.end(); }
    auto begin() const { return selected_atom_list.begin(); }
    auto end() const { return selected_atom_list.end(); }
};


using FitState = std::vector<LocalGaussianResult>;

struct FitStatePatch
{
    ClusterKey atom_index_list{};
    std::vector<GaussianModel3DWithUncertainty> mdpde_list{};

    static FitStatePatch FromState(
        const FitState & state,
        ClusterKey atom_index_list);

    const GaussianModel3DWithUncertainty * Find(
        std::size_t atom_index) const;

    void ApplyTo(FitState & state) const;
};

struct FitStateProposal
{
    FitStatePatch patch{};
    double effective_damping{ 0.0 };
    double step_norm{ 0.0 };
};

class FitStateView
{
    const FitState & m_base_state;
    const FitStatePatch & m_patch;

public:
    FitStateView(const FitState & base_state, const FitStatePatch & patch)
        : m_base_state{ base_state }, m_patch{ patch }
    {
    }

    const GaussianModel3DWithUncertainty & GetMdpde(std::size_t atom_index) const
    {
        const auto * value{ FindOverride(atom_index) };
        if (value != nullptr) return *value;
        return m_base_state.at(atom_index).mdpde;
    }

    const GaussianModel3D & GetModel(std::size_t atom_index) const
    {
        return GetMdpde(atom_index).GetModel();
    }

    const GaussianModel3D & GetBaseModel(std::size_t atom_index) const
    {
        return m_base_state.at(atom_index).mdpde.GetModel();
    }

    const GaussianModel3DWithUncertainty * FindOverride(std::size_t atom_index) const
    {
        return m_patch.Find(atom_index);
    }

    const ClusterKey & GetOverrideAtomIndexList() const
    {
        return m_patch.atom_index_list;
    }

    FitState Materialize() const;

    std::size_t size() const { return m_base_state.size(); }
};

const GaussianModel3D & GetFitModel(const FitState & state, std::size_t atom_index);

const GaussianModel3D & GetFitModel(const FitStateView & state, std::size_t atom_index);

using PolishProvenance = std::vector<char>;

bool UsesPolish(const PolishProvenance & provenance);


using SecondStageAdjustedResponseCache = std::vector<std::vector<double>>;
using FittedGaussianSnapshot = std::vector<GaussianModel3D>;

const GaussianModel3D & GetFitModel(
    const FittedGaussianSnapshot & state,
    std::size_t atom_index);

FittedGaussianSnapshot BuildFittedGaussianSnapshot(const FitState & state);

FittedGaussianSnapshot BuildFittedGaussianSnapshot(const FitStateView & state);

struct SecondStageModelSnapshot
{
    FittedGaussianSnapshot selected{};
    FittedGaussianSnapshot unselected{};
};

const GaussianModel3D & ResolveNeighborAtomModel(
    const NeighborAtomSample & neighbor_atom_sample,
    const SecondStageModelSnapshot & model_snapshot);

SecondStageModelSnapshot BuildSecondStageModelSnapshot(
    const SecondStageContext & context,
    FittedGaussianSnapshot selected_snapshot);

SecondStageModelSnapshot BuildSecondStageModelSnapshot(
    const SecondStageContext & context,
    const FitState & state);

struct ResidualSample
{
    double adjusted_response{ 0.0 };
    double residual{ 0.0 };
};

struct ResidualBaseline
{
    SecondStageModelSnapshot model_snapshot{};
    std::vector<std::vector<std::optional<ResidualSample>>> sample_list{};

    std::optional<ResidualSample> operator()(const SampleRef & sample_ref) const
    {
        return sample_list.at(sample_ref.atom_index).at(sample_ref.sample_index);
    }

    const FittedGaussianSnapshot & GetState() const
    {
        return model_snapshot.selected;
    }
};

SecondStageAdjustedResponseCache BuildSecondStageAdjustedResponseCache(
    const SecondStageContext & context,
    const SecondStageModelSnapshot & model_snapshot);

LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const AtomContext & atom_context,
    const std::vector<double> & adjusted_response_list);

LocalPotentialSampleList BuildSecondStageAdjustedSamples(
    const AtomContext & atom_context,
    const SecondStageModelSnapshot & model_snapshot);

std::optional<ResidualSample> EvaluateResidualSample(
    const SecondStageContext & context,
    const FitStateView & state,
    const SampleRef & sample_ref,
    const SecondStageModelSnapshot & model_snapshot);

std::optional<ResidualSample> EvaluateResidualSample(
    const SecondStageContext & context,
    const FittedGaussianSnapshot & state,
    const SampleRef & sample_ref,
    const SecondStageModelSnapshot & model_snapshot);

struct SnapshotResidualEvaluator
{
    const SecondStageContext & context;
    const SecondStageModelSnapshot & model_snapshot;

    std::optional<ResidualSample> operator()(
        const SampleRef & sample_ref) const;

    const FittedGaussianSnapshot & GetState() const { return model_snapshot.selected; }
};

ResidualBaseline BuildResidualBaseline(const SecondStageContext & context, const FitState & state);


constexpr double kTransformedChangeTolerance{ 1.0e-4 };

struct TransformedChangeSummary
{
    algorithm::ParameterChangeStats percentile_stats{};
    std::vector<double> maximum_list{};
};

algorithm::ParameterChange CalculateTransformedChange(
    const GaussianModel3D & current,
    const GaussianModel3D & previous);

double GetMaximumTransformedChange(const std::vector<double> & value_list);

bool IsTransformedChangeMaterial(
    const algorithm::ParameterChange & change,
    double minimum_change);

std::vector<double> SummarizeMaximumTransformedChanges(
    const std::vector<algorithm::ParameterChange> & change_list,
    const std::vector<std::size_t> & index_list);

bool IsTransformedChangeConverged(
    const algorithm::ParameterChangeStats & percentile_stats,
    const std::vector<double> & maximum_list);

TransformedChangeSummary SummarizeTransformedChanges(
    const FitState & current_state,
    const FitState & previous_state,
    const std::vector<std::size_t> & index_list);

TransformedChangeSummary SummarizeTransformedChanges(
    const FitStateView & current_state,
    const FittedGaussianSnapshot & previous_state,
    const std::vector<std::size_t> & index_list);

double GetMaximumTransformedChange(const TransformedChangeSummary & summary);

bool IsTransformedChangeConverged(const TransformedChangeSummary & summary);

bool IsTransformedPercentileConverged(const TransformedChangeSummary & summary);


bool IsTrustRegionStepWithinRadius(double step_norm, double radius);

bool IsTrustRegionStepAtGrowthBoundary(double step_norm, double radius);

struct TrustRegionDamping
{
    double effective_damping{ 1.0 };
    double step_norm{ 0.0 };
};

TrustRegionDamping LimitTrustRegionDamping(
    const std::vector<Eigen::Vector3d> & previous_estimation_list,
    const std::vector<Eigen::Vector3d> & candidate_estimation_list,
    double requested_damping,
    double radius);

std::optional<double> CalculateModelTrustRegionStepNorm(
    const std::vector<GaussianModel3D> & previous_model_list,
    const std::vector<GaussianModel3D> & candidate_model_list);

} // namespace rhbm_gem::core::detail
