#pragma once

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

namespace rhbm_gem {

class AtomObject;
class ModelObject;

class ModelAnalysisEditor
{
    ModelObject & m_model_object;

public:
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
    void SetAtomGroupAlphaR(FittingStage stage, GroupKey group_key, double alpha_r);
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

} // namespace rhbm_gem
