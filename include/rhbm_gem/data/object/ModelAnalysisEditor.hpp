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
    void InitializeLocalAlpha(LocalFittingStage stage, double alpha_r);
    void InitializeGroupAlpha(LocalFittingStage stage, double alpha_g);
    void CopyAtomGroupPotentialStage(
        LocalFittingStage source_stage,
        LocalFittingStage destination_stage);
    void ApplyAtomGroupGaussianResult(
        LocalFittingStage stage,
        GroupKey group_key,
        const GroupGaussianResult & group_result);
    void SetAtomGroupAlphaG(
        LocalFittingStage stage,
        GroupKey group_key,
        double alpha_g);
    
};

} // namespace rhbm_gem
