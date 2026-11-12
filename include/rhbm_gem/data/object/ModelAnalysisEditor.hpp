#pragma once
#include <rhbm_gem/data/object/JointAnalysisResult.hpp>
#include <rhbm_gem/data/object/PostJointQScoreResult.hpp>
#include <map>
#include <utility>
#include <vector>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

namespace rhbm_gem {

class AtomObject;
class ModelObject;
class MapObject;

class ModelAnalysisEditor
{
    ModelObject & m_model_object;

public:
    void BindJointMapObservations(const MapObject &);
    // Replacing diagnostics invalidates downstream grid evidence.
    void SetJointResult(JointAnalysisResult result);
    // Publish snapshot endpoints atomically, then retain the diagnostic snapshot.
    void ApplyJointResult(JointAnalysisResult result, const std::string & run_id);
    void ClearJointResult();
    void ClearPostJointQScores();
    void ApplyPostJointQScoresBatch(std::vector<std::pair<int, PostJointQScoreResult>>);
    void UpdateJointMetadata(JointAnalysisMetadata metadata);
    void ApplyJointStageEstimates(std::map<int, LocalStageEstimate> estimates);
    void SetAtomInitialStageEstimate(const AtomObject & atom, LocalStageEstimate value);
    void SetAtomJointStageEstimate(const AtomObject & atom, LocalStageEstimate value);
    void SetAtomPostFitPeeling(const AtomObject & atom, PostFitPeelingResult value);
    void ApplyPostFitPeelingsBatch(std::vector<std::pair<int, PostFitPeelingResult>> peelings);
    void SetAtomPostJointLocalResult(const AtomObject & atom, PostJointLocalResult value);
    void ApplyPostJointLocalResultsBatch(std::vector<std::pair<int, PostJointLocalResult>> results);
    void ApplyAtomGroupParameterSummary(GroupKey key, GroupParameterSummary value);
    void Clear();
    void ClearTransientFitStates();
    void InitializeFromSelection();
    void EnsureSelectedAtomLocalPotentials();
    // Replaces interpolated sampling entries with available coordinates.
    void SetAtomLocalSamplingEntries(
        const AtomObject & atom_object,
        LocalPotentialSampleList value);
    void SetAtomInitialStageGaussianResult(
        const AtomObject & atom_object,
        LocalGaussianResult result);
    void SetAtomInitialStageAlphaR(
        const AtomObject & atom_object,
        double alpha_r);
    void RebuildAtomGroupsFromSelection();
    void SetSelectedAtomsInitialStageAlphaR(double alpha_r);
    void InitializeGroupAlpha(double alpha_g);
    void ApplyAtomGroupGaussianResult(
        GroupKey group_key,
        const GroupGaussianResult & group_result);
    void SetAtomGroupAlphaG(GroupKey group_key, double alpha_g);

private:
    explicit ModelAnalysisEditor(ModelObject & model_object);
    friend class ModelObject;
};

} // namespace rhbm_gem
