#pragma once

#include <optional>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
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
    mutable RHBMMemberDataset m_dataset_cache{};
    mutable RHBMBetaEstimateResult m_fit_result_cache{};
    
public:
    MutableLocalPotentialView() = default;
    void SetSamplingEntries(LocalPotentialSampleList value);
    void SetDataset(RHBMMemberDataset value);
    void SetFitResult(RHBMBetaEstimateResult value);
    void SetAlphaR(double value);
    void SetAnnotation(const std::string & key, const LocalPotentialAnnotationData & value);
    bool HasDataset() const;
    bool HasFitResult() const;
    const RHBMMemberDataset & GetDataset() const;
    const RHBMBetaEstimateResult & GetFitResult() const;
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
    void ApplyAtomGroupEstimateResult(
        GroupKey group_key,
        const std::string & class_key,
        const RHBMGroupEstimationResult & result,
        double alpha_g);
    void ApplyBondGroupEstimateResult(
        GroupKey group_key,
        const std::string & class_key,
        const RHBMGroupEstimationResult & result,
        double alpha_g);
    void SetAtomGroupAlphaG(GroupKey group_key, const std::string & class_key, double alpha_g);
    void SetBondGroupAlphaG(GroupKey group_key, const std::string & class_key, double alpha_g);
    
};

} // namespace rhbm_gem
