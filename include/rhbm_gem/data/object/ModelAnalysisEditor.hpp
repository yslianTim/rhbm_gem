#pragma once

#include <string>

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
    void ApplyAtomGroupGaussianResult(
        GroupKey group_key,
        const std::string & class_key,
        const GroupGaussianResult & group_result);
    void SetAtomGroupAlphaG(GroupKey group_key, const std::string & class_key, double alpha_g);
    
};

} // namespace rhbm_gem
