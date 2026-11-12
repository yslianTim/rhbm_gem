#pragma once

#include <rhbm_gem/data/object/PostJointQScoreResult.hpp>

#include <cmath>
#include <cstddef>
#include <optional>
#include <utility>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {

class LocalPotentialEntry
{
    struct InitialStageRecord
    {
        LocalStageEstimate estimate;
        double alpha_r{};
        bool uncertainty_recorded{ false };
        GaussianModel3DWithUncertainty ols{GaussianModel3D{0, 0}, {}};
        GaussianModel3DUncertainty mdpde_uncertainty;
        // Seeds/invalid native outcomes are not published final points.
        std::optional<GaussianModel3D> unfitted_model;
        std::optional<LocalFitDiagnostics> diagnostics;
        std::optional<RHBMBetaEstimateResult> fit_result;
    };
    LocalPotentialSampleList m_sampling_entries;
    InitialStageRecord m_initial_stage;
    LocalStageEstimate m_joint_stage_estimate;
    std::optional<PostFitPeelingResult> m_post_fit_peeling;
    std::optional<PostJointLocalResult> m_post_joint_local_result;
    std::optional<PostJointQScoreResult> m_post_joint_qscore;
    bool m_sample_geometry_available{ true };
    std::optional<GroupGaussianMemberResult> m_group_member_result{};

public:
    void SetInitialStageAlphaR(double value) { m_initial_stage.alpha_r = value; }
    double InitialStageAlphaR() const { return m_initial_stage.alpha_r; }
    // Workflow replacements have geometry; persistence supplies its recorded availability.
    void SetSamplingEntries(LocalPotentialSampleList value, bool geometry_available = true)
    {
        m_sampling_entries = std::move(value);
        m_sample_geometry_available = geometry_available;
    }
    void SetInitialStageGaussianResult(LocalGaussianResult value)
    {
        auto & record = m_initial_stage;
        record = InitialStageRecord{};
        const auto & model = value.mdpde.GetModel();
        if (std::isfinite(model.GetAmplitude()) && std::isfinite(model.GetOffset()) &&
            std::isfinite(model.GetWidth()) && model.GetWidth() > 0)
        {
            record.estimate.point = model;
            record.estimate.source.method = EstimateMethod::LocalMDPDE;
            record.estimate.reason.clear();
        }
        else record.unfitted_model = model;
        record.alpha_r = value.alpha_r;
        record.uncertainty_recorded = value.uncertainty_recorded;
        record.ols = std::move(value.ols);
        record.mdpde_uncertainty = value.mdpde.GetStandardDeviationModel();
        record.diagnostics = std::move(value.diagnostics);
        if (value.fit_result)
            record.diagnostics = LocalFitDiagnostics{value.fit_result->status, value.fit_result->sigma_square,
                value.fit_result->diagnostics, value.fit_result->refinement, value.fit_result->failed_only_recovery};
        record.fit_result = std::move(value.fit_result);
    }
    const LocalStageEstimate & InitialStageEstimate() const { return m_initial_stage.estimate; }
    const LocalStageEstimate & JointStageEstimate() const { return m_joint_stage_estimate; }
    void SetJointStageEstimate(LocalStageEstimate value) { m_joint_stage_estimate = std::move(value); }
    void SetInitialStageEstimate(LocalStageEstimate value)
    {
        auto & record = m_initial_stage;
        if (!value.point && record.estimate.point) record.unfitted_model = record.estimate.point;
        if (value.point && (record.estimate.point || record.unfitted_model) &&
            InitialStageMDPDEModel().ToVector() != value.point->ToVector())
        {
            record.diagnostics.reset(); record.fit_result.reset(); record.mdpde_uncertainty = {};
        }
        if (value.point) record.unfitted_model.reset();
        record.estimate = std::move(value);
    }
    LocalGaussianResult InitialStageGaussianResult(bool include_transient = true) const
    {
        const auto & r = m_initial_stage;
        return {r.alpha_r, r.ols, {InitialStageMDPDEModel(), r.mdpde_uncertainty}, include_transient ? r.fit_result : std::nullopt, r.diagnostics, r.uncertainty_recorded};
    }
    const GaussianModel3D & InitialStageMDPDEModel() const
    {
        const auto & r = m_initial_stage;
        if (r.estimate.point) return *r.estimate.point;
        if (r.unfitted_model) return *r.unfitted_model;
        static const GaussianModel3D missing{0, 0};
        return missing;
    }
    const GaussianModel3D & InitialStageOLSModel() const { return m_initial_stage.ols.GetModel(); }
    const std::optional<PostFitPeelingResult> & PostFitPeeling() const { return m_post_fit_peeling; }
    void SetPostFitPeeling(PostFitPeelingResult value)
    {
        m_post_joint_local_result.reset();
        m_post_fit_peeling = std::move(value);
    }
    const std::optional<PostJointLocalResult> & PostJointLocalResult() const { return m_post_joint_local_result; }
    void SetPostJointLocalResult(rhbm_gem::PostJointLocalResult value) { m_post_joint_local_result = std::move(value); }
    void ClearPostJointLocalResult() { m_post_joint_local_result.reset(); }
    bool SampleGeometryAvailable() const { return m_sample_geometry_available; }
    void SetGroupMemberResult(GroupGaussianMemberResult value) { m_group_member_result = std::move(value); }
    const std::optional<PostJointQScoreResult> & PostJointQScore() const { return m_post_joint_qscore; }
    void SetPostJointQScore(PostJointQScoreResult value) { m_post_joint_qscore = std::move(value); }
    void ClearPostJointQScore() { m_post_joint_qscore.reset(); }
    void ClearPeeling()
    {
        m_post_joint_qscore.reset();
        m_post_fit_peeling.reset();
        m_post_joint_local_result.reset();
    }
    void ClearGroupMemberResult() { m_group_member_result.reset(); }
    void ClearInitialStageTransientFitState() { m_initial_stage.fit_result.reset(); }
    void ClearPostJointTransientFitState()
    {
        if (m_post_joint_local_result) m_post_joint_local_result->fit.fit_result.reset();
    }
    const std::optional<GroupGaussianMemberResult> & GroupMemberResult() const { return m_group_member_result; }
    // All interpolated sampling records, including unselected points, in stored order.
    const LocalPotentialSampleList & SamplingEntries() const { return m_sampling_entries; }
    // Projects available Original-grid peeled observations, without interpolation.
    // The shared container's is_selected=true is a projection convention only;
    // PostFitPeeling() retains full voxel identity, unavailable reasons and provenance.
    LocalPotentialSampleList PeeledObservationEntries() const
    {
        LocalPotentialSampleList samples;
        if (!m_post_fit_peeling) return samples;
        for (const auto & estimate : m_post_fit_peeling->samples)
        {
            if (!estimate.response) continue;
            const auto & observation = estimate.observation;
            LocalPotentialSample sample{*estimate.response,
                {std::sqrt(observation.squared_distance),observation.position,true}};
            samples.push_back(std::move(sample));
        }
        return samples;
    }
};

} // namespace rhbm_gem
