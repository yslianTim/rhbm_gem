#pragma once

#include <rhbm_gem/data/object/AtomLocalPotentialEditor.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>

namespace rhbm_gem {

class AtomObject;
class ModelObject;

class ModelAnalysisEditor
{
    ModelObject & m_model_object;

public:
    explicit ModelAnalysisEditor(ModelObject & model_object);
    void Clear();
    void ClearTransientFitStates();
    AtomLocalPotentialEditor EnsureAtomLocalPotential(const AtomObject & atom_object);
    void RebuildAtomGroupsFromSelection();
    void InitializeLocalAlpha(FittingStage stage, double alpha_r);
    void InitializeGroupAlpha(FittingStage stage, double alpha_g);
    void CopyLocalFittingStageResult(FittingStage source_stage, FittingStage destination_stage);
    void CopyGroupFittingStageResult(FittingStage source_stage, FittingStage destination_stage);
    void CopyFittingStageState(FittingStage source_stage, FittingStage destination_stage);
    void ApplyAtomGroupGaussianResult(
        FittingStage stage,
        GroupKey group_key,
        const GroupGaussianResult & group_result);
    void SetAtomGroupAlphaG(FittingStage stage, GroupKey group_key, double alpha_g);
    
};

} // namespace rhbm_gem
