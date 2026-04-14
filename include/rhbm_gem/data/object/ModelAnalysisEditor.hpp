#pragma once

#include <optional>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <rhbm_gem/data/object/GaussianStatistics.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/hrl/HRLModelTypes.hpp>
#include <rhbm_gem/utils/math/SamplingTypes.hpp>

namespace rhbm_gem {

class AtomObject;
class BondObject;
class ModelObject;

struct LocalPotentialFitResult
{
    Eigen::VectorXd beta_ols;
    Eigen::VectorXd beta_mdpde;
    double sigma_square{ 0.0 };
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> data_weight;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> data_covariance;
};

struct LocalPotentialEstimates
{
    GaussianEstimate ols{};
    GaussianEstimate mdpde{};
};

struct LocalPotentialAnnotationData
{
    GaussianPosterior posterior{};
    bool is_outlier{ false };
    double statistical_distance{ 0.0 };
};

struct GroupPotentialStatistics
{
    GaussianEstimate mean{};
    GaussianEstimate mdpde{};
    GaussianEstimate prior{};
    GaussianEstimate prior_variance{};
    double alpha_g{ 0.0 };
};

class MutableLocalPotentialView
{
public:
    MutableLocalPotentialView() = default;

    void SetSamplingEntries(LocalPotentialSampleList value);
    void SetDataset(HRLMemberDataset value);
    void SetFitResult(LocalPotentialFitResult value);
    void SetEstimates(const LocalPotentialEstimates & value);
    void SetAlphaR(double value);
    void SetAnnotation(const std::string & key, const LocalPotentialAnnotationData & value);

    bool HasDataset() const;
    bool HasFitResult() const;
    const HRLMemberDataset & GetDataset() const;
    const LocalPotentialFitResult & GetFitResult() const;
    const AtomObject * GetAtomObjectPtr() const { return m_atom_object; }
    const BondObject * GetBondObjectPtr() const { return m_bond_object; }
    const void * GetEntryHandle() const { return m_entry_ptr; }

private:
    AtomObject * m_atom_object{ nullptr };
    BondObject * m_bond_object{ nullptr };
    void * m_entry_ptr{ nullptr };
    mutable HRLMemberDataset m_dataset_cache{};
    mutable LocalPotentialFitResult m_fit_result_cache{};

    explicit MutableLocalPotentialView(AtomObject * atom_object);
    explicit MutableLocalPotentialView(BondObject * bond_object);

    friend class ModelAnalysisEditor;
};

class ModelAnalysisEditor
{
public:
    explicit ModelAnalysisEditor(ModelObject & model_object);

    static ModelAnalysisEditor Of(ModelObject & model_object);

    void Clear();
    void ClearTransientFitStates();

    MutableLocalPotentialView EnsureAtomLocalPotential(const AtomObject & atom_object);
    MutableLocalPotentialView EnsureBondLocalPotential(const BondObject & bond_object);

    void RebuildAtomGroupsFromSelection();
    void RebuildBondGroupsFromSelection();

    void SetAtomGroupStatistics(
        GroupKey group_key,
        const std::string & class_key,
        const GroupPotentialStatistics & statistics);
    void SetBondGroupStatistics(
        GroupKey group_key,
        const std::string & class_key,
        const GroupPotentialStatistics & statistics);

    void SetAtomGroupAlphaG(GroupKey group_key, const std::string & class_key, double alpha_g);
    void SetBondGroupAlphaG(GroupKey group_key, const std::string & class_key, double alpha_g);

private:
    ModelObject & m_model_object;
};

} // namespace rhbm_gem
