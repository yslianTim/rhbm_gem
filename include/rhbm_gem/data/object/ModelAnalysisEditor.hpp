#pragma once
#include <rhbm_gem/data/object/JointAnalysisResult.hpp>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

namespace rhbm_gem {

class AtomObject;
class ModelObject;

class ModelAnalysisEditor
{
    ModelObject & m_model_object;

public:
    // Attach diagnostics to already-published stages; saving validates consistency.
    void SetJointResult(JointAnalysisResult result);
    // Publish snapshot endpoints atomically, then retain the diagnostic snapshot.
    void ApplyJointResult(JointAnalysisResult result, const std::string & run_id);
    void ClearJointResult();
    void UpdateJointMetadata(JointAnalysisMetadata metadata);
    void ApplySecondStageEstimates(std::map<int, LocalStageEstimate> estimates);
    void SetAtomStageEstimate(FittingStage stage, const AtomObject & atom, LocalStageEstimate value);
    void SetAtomPostFitPeeling(const AtomObject & atom, PostFitPeelingResult value);
    void SetAtomPostJointLocalResult(const AtomObject & atom, PostJointLocalResult value);
    void SetAtomGroupEvidence(const AtomObject & atom, GroupParameterEvidence value);
    void ApplyJointPostFitEvidenceBatch(
        std::vector<std::pair<int, PostFitPeelingResult>> peelings,
        std::vector<std::pair<int, GroupParameterEvidence>> evidence);
    void ApplyAtomGroupParameterSummary(GroupKey key, GroupParameterSummary value);
    void Clear();
    void ClearTransientFitStates();
    void InitializeFromSelection();
    void InitializeLocalFittingSeedModels();
    void EnsureSelectedAtomLocalPotentials();
    void EnsureAtomGroupLocalPotentials(GroupKey group_key);
    void SetAtomLocalRawSamplingEntries(
        const AtomObject & atom_object,
        LocalPotentialSampleList value);
    void SetAtomLocalPeelingSamplingEntries(
        const AtomObject & atom_object,
        LocalPotentialSampleList value);
    void SetAtomLocalGaussianResult(
        FittingStage stage,
        const AtomObject & atom_object,
        LocalGaussianResult result);
    void SetAtomLocalAlphaR(
        FittingStage stage,
        const AtomObject & atom_object,
        double alpha_r);
    void RebuildAtomGroupsFromSelection();
    void InitializeLocalAlpha(FittingStage stage, double alpha_r);
    void InitializeGroupAlpha(double alpha_g);
    void CopyLocalFittingStageResult(FittingStage source_stage, FittingStage destination_stage);
    void ApplyAtomGroupGaussianResult(
        GroupKey group_key,
        const GroupGaussianResult & group_result);
    void SetAtomLocalNeighborCountForPeeling(
        const AtomObject & atom_object,
        int neighbor_count);
    void ApplyAtomLocalSecondStageResult(
        const AtomObject & atom_object,
        LocalGaussianResult result,
        LocalPotentialSampleList peeling_sampling_entries);
    void SetAtomGroupAlphaG(GroupKey group_key, double alpha_g);

private:
    explicit ModelAnalysisEditor(ModelObject & model_object);
    friend class ModelObject;
};

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
struct GroupInvalidationStats
{
    std::size_t membership_traversals{};
    std::size_t groups_invalidated{};
};
GroupInvalidationStats GetGroupInvalidationStatsForTesting();
void ResetGroupInvalidationStatsForTesting();
#endif

} // namespace rhbm_gem
