#pragma once

#include <optional>

#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/GaussianModel3D.hpp>

namespace rhbm_gem {

class AtomObject;
class LocalPotentialEntry;

class AtomLocalPotentialView
{
    const AtomObject * m_atom_object{ nullptr };

public:
    static AtomLocalPotentialView For(const AtomObject & atom_object);
    bool IsAvailable() const;
    const LocalStageEstimate & GetInitialStageEstimate() const;
    bool HasInitialStageModel() const;
    const GaussianModel3D & GetInitialStageModel() const;
    const LocalStageEstimate & GetJointStageEstimate() const;
    bool HasJointStageModel() const;
    const GaussianModel3D & GetJointStageModel() const;
    const std::optional<PostFitPeelingResult> & GetPostFitPeeling() const;
    const std::optional<PostJointLocalResult> & GetPostJointLocalResult() const;
    bool HasSampleGeometry() const;
    LocalGaussianResult GetInitialStageGaussianResult() const;
    const std::optional<GroupGaussianMemberResult> & GetGroupMemberResult() const;
    const GaussianModel3D & GetInitialStageOLS() const;
    const GaussianModel3D & GetInitialStageMDPDE() const;
    LocalPotentialSampleList GetRawSamplingEntries(bool apply_selection = true) const;
    const LocalPotentialSampleList & GetRawSamplingEntriesView() const;
    LocalPotentialSampleList GetPeelingSamplingEntries(bool apply_selection = false) const;
    PeelingRatioAssessment AssessLocalFittingPeelingRatio(
        double distance_min, double distance_max) const;
    std::optional<double> GetLocalFittingPeelingRatio(
        double distance_min,
        double distance_max) const;
    int GetNeighborCountForPeeling() const;
    double GetInitialStageAlphaR() const;

private:
    explicit AtomLocalPotentialView(const AtomObject * atom_object);
    const LocalPotentialEntry * FindEntry() const;
    const LocalPotentialEntry & RequireEntry(const char * context) const;
    const LocalPotentialEntry & RequireInitialStageEntry(const char * diagnostic_error) const;
};

} // namespace rhbm_gem
