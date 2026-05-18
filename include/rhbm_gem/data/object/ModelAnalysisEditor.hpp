#pragma once

#include <string>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {

class AtomObject;
class LocalPotentialEntry;
class ModelObject;

class LocalPotentialEditor
{
    AtomObject * m_atom_object{ nullptr };
    LocalPotentialEntry * m_entry{ nullptr };
    
public:
    LocalPotentialEditor() = default;
    void SetSamplingEntries(LocalPotentialSampleList value);
    void SetGaussianResult(LocalGaussianResult value);
    void SetAlphaR(double value);

private:
    explicit LocalPotentialEditor(AtomObject * atom_object);
    ModelObject * GetOwner() const;
    LocalPotentialEntry & EnsureEntry() const;
    friend class ModelAnalysisEditor;

};

class ModelAnalysisEditor
{
    ModelObject & m_model_object;

public:
    explicit ModelAnalysisEditor(ModelObject & model_object);
    void Clear();
    void ClearTransientFitStates();
    LocalPotentialEditor EnsureAtomLocalPotential(const AtomObject & atom_object);
    void RebuildAtomGroupsFromSelection();
    void ApplyAtomGroupGaussianResult(
        GroupKey group_key,
        const std::string & class_key,
        const GroupGaussianResult & group_result);
    void SetAtomGroupAlphaG(GroupKey group_key, const std::string & class_key, double alpha_g);
    
};

} // namespace rhbm_gem
