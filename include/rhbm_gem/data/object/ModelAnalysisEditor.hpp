#pragma once

#include <optional>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/GaussianEstimationTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class ModelObject;

struct LocalPotentialAnnotationData
{
    GaussianModel3DWithUncertainty gaussian{
        GaussianModel3D{ 0.0, 0.0 },
        GaussianModel3DUncertainty{}
    };
    bool is_outlier{ false };
    double statistical_distance{ 0.0 };
};

class MutableLocalPotentialView
{
    AtomObject * m_atom_object{ nullptr };
    BondObject * m_bond_object{ nullptr };
    void * m_entry_ptr{ nullptr };
    
public:
    MutableLocalPotentialView() = default;
    void SetSamplingEntries(LocalPotentialSampleList value);
    void SetGaussianResult(LocalGaussianResult value);
    void SetAlphaR(double value);
    void SetAnnotation(const std::string & key, const LocalPotentialAnnotationData & value);
    double GetAlphaR() const;
    const LocalGaussianResult & GetGaussianResult() const;
    const LocalPotentialSampleList & GetSamplingEntries() const;
    int GetSamplingEntryCount() const;
    const AtomObject * GetAtomObjectPtr() const { return m_atom_object; }
    const BondObject * GetBondObjectPtr() const { return m_bond_object; }
    const void * GetEntryHandle() const { return m_entry_ptr; }

private:
    explicit MutableLocalPotentialView(AtomObject * atom_object);
    explicit MutableLocalPotentialView(BondObject * bond_object);
    friend class ModelAnalysisEditor;

};

class ModelAnalysisEditor
{
    ModelObject & m_model_object;

public:
    explicit ModelAnalysisEditor(ModelObject & model_object);
    static ModelAnalysisEditor Of(ModelObject & model_object);
    void Clear();
    void ClearTransientFitStates();
    MutableLocalPotentialView EnsureAtomLocalPotential(const AtomObject & atom_object);
    MutableLocalPotentialView EnsureBondLocalPotential(const BondObject & bond_object);
    void RebuildAtomGroupsFromSelection();
    void RebuildBondGroupsFromSelection();
    void ApplyAtomGroupGaussianResult(
        GroupKey group_key,
        const std::string & class_key,
        const GroupGaussianResult & group_result);
    void SetAtomGroupAlphaG(GroupKey group_key, const std::string & class_key, double alpha_g);
    void SetBondGroupAlphaG(GroupKey group_key, const std::string & class_key, double alpha_g);
    
};

} // namespace rhbm_gem
