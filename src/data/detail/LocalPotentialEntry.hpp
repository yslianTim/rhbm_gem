#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>

#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rhbm_gem {

class LocalPotentialEntry
{
    struct StageRecord
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
    LocalPotentialSampleList m_raw_sampling_entries;
    LocalPotentialSampleList m_peeling_sampling_entries;
    std::array<StageRecord, 2> m_stages{};
    std::optional<PostFitPeelingResult> m_post_fit_peeling;
    bool m_sample_geometry_available{ true };
    std::optional<GroupGaussianMemberResult> m_group_member_result{};
    std::optional<GroupParameterEvidence> m_group_evidence;
    int m_neighbor_count_for_peeling{ 0 };

public:
    void SetAlphaR(FittingStage stage, double value) { Record(stage).alpha_r = value; }
    double AlphaR(FittingStage stage) const { return Record(stage).alpha_r; }
    void SetRawSamplingEntries(LocalPotentialSampleList value) { m_raw_sampling_entries = std::move(value); }
    void SetPeelingSamplingEntries(LocalPotentialSampleList value)
    {
        if (StageEstimate(FittingStage::Second).source.method == EstimateMethod::JointComponents)
            throw std::invalid_argument("Joint peeling requires paired raw-sample responses.");
        m_peeling_sampling_entries = std::move(value);
    }
    void SetNeighborCountForPeeling(int value)
    {
        if (StageEstimate(FittingStage::Second).source.method == EstimateMethod::JointComponents)
            throw std::invalid_argument("Joint neighbor count belongs to paired peeling.");
        m_neighbor_count_for_peeling = value;
    }
    void SetGaussianResult(FittingStage stage, LocalGaussianResult value)
    {
        auto & record = Record(stage);
        record = StageRecord{};
        const auto & model = value.mdpde.GetModel();
        if (std::isfinite(model.GetAmplitude()) && std::isfinite(model.GetOffset()) &&
            std::isfinite(model.GetWidth()) && model.GetWidth() > 0)
        {
            record.estimate.point = model;
            record.estimate.source.method = stage == FittingStage::First ? EstimateMethod::LocalMDPDE : EstimateMethod::Peeling;
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
                value.fit_result->diagnostics, value.fit_result->refinement};
        record.fit_result = std::move(value.fit_result);
    }
    const LocalStageEstimate & StageEstimate(FittingStage stage) const { return Record(stage).estimate; }
    void SetStageEstimate(FittingStage stage, LocalStageEstimate value)
    {
        auto & record = Record(stage);
        if (value.source.method == EstimateMethod::JointComponents) record = StageRecord{};
        else if (!value.point && record.estimate.point) record.unfitted_model = record.estimate.point;
        if (record.estimate.point && value.point && record.estimate.point->ToVector() != value.point->ToVector())
        {
            record.diagnostics.reset(); record.fit_result.reset(); record.mdpde_uncertainty = {};
        }
        if (value.point) record.unfitted_model.reset();
        record.estimate = std::move(value);
    }
    LocalGaussianResult GaussianResult(FittingStage stage, bool include_transient = true) const
    {
        const auto & r = Record(stage);
        return {r.alpha_r, r.ols, {MDPDEModel(stage), r.mdpde_uncertainty}, include_transient ? r.fit_result : std::nullopt, r.diagnostics, r.uncertainty_recorded};
    }
    const GaussianModel3D & MDPDEModel(FittingStage stage) const
    {
        const auto & r = Record(stage);
        if (r.estimate.point) return *r.estimate.point;
        if (r.unfitted_model) return *r.unfitted_model;
        static const GaussianModel3D missing{0, 0};
        return missing;
    }
    const GaussianModel3D & OLSModel(FittingStage stage) const { return Record(stage).ols.GetModel(); }
    const std::optional<PostFitPeelingResult> & PostFitPeeling() const { return m_post_fit_peeling; }
    void SetPostFitPeeling(PostFitPeelingResult value)
    {
        m_peeling_sampling_entries.clear();
        m_neighbor_count_for_peeling = 0;
        m_post_fit_peeling = std::move(value);
    }
    void ClearPostFitPeeling() { m_post_fit_peeling.reset(); }
    bool SampleGeometryAvailable() const { return m_sample_geometry_available; }
    void SetSampleGeometryAvailable(bool value) { m_sample_geometry_available = value; }
    void SetGroupMemberResult(GroupGaussianMemberResult value) { m_group_member_result = std::move(value); }
    const std::optional<GroupParameterEvidence> & GroupEvidence() const { return m_group_evidence; }
    void SetGroupEvidence(GroupParameterEvidence value) { m_group_evidence = std::move(value); }
    void ClearGroupEvidence() { m_group_evidence.reset(); }
    void ClearPeeling()
    {
        m_post_fit_peeling.reset();
        m_peeling_sampling_entries.clear();
        m_neighbor_count_for_peeling = 0;
    }
    void ClearGroupMemberResult() { m_group_member_result.reset(); }
    void ClearTransientFitState(FittingStage stage) { Record(stage).fit_result.reset(); }
    int RawSamplingEntryCount() const { return static_cast<int>(m_raw_sampling_entries.size()); }
    int PeelingSamplingEntryCount() const
    {
        if (m_post_fit_peeling) return static_cast<int>(std::count_if(m_post_fit_peeling->samples.begin(),
            m_post_fit_peeling->samples.end(), [](const auto & s) { return s.response.has_value(); }));
        return static_cast<int>(m_peeling_sampling_entries.size());
    }
    int NeighborCountForPeeling() const
    { return m_post_fit_peeling ? static_cast<int>(m_post_fit_peeling->neighbor_count) : m_neighbor_count_for_peeling; }
    const std::optional<GroupGaussianMemberResult> & GroupMemberResult() const { return m_group_member_result; }
    const LocalPotentialSampleList & RawSamplingEntries() const { return m_raw_sampling_entries; }
    LocalPotentialSampleList PeelingSamplingEntries() const
    {
        if (!m_post_fit_peeling) return m_peeling_sampling_entries;
        LocalPotentialSampleList samples;
        for (std::size_t i = 0; i < m_post_fit_peeling->samples.size(); ++i)
        {
            if (!m_post_fit_peeling->samples[i].response) continue;
            auto sample = m_raw_sampling_entries.at(i);
            sample.response = *m_post_fit_peeling->samples[i].response;
            samples.push_back(std::move(sample));
        }
        return samples;
    }

private:
    StageRecord & Record(FittingStage stage) { return m_stages.at(StageIndex(stage)); }
    const StageRecord & Record(FittingStage stage) const { return m_stages.at(StageIndex(stage)); }
    static std::size_t StageIndex(FittingStage stage)
    {
        const auto index = static_cast<std::size_t>(stage);
        if (index >= 2) throw std::invalid_argument("Unknown local fitting stage.");
        return index;
    }
};

} // namespace rhbm_gem
