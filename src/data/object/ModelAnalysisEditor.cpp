#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

#include <stdexcept>
#include <string>

namespace rhbm_gem {

namespace {

ModelObject * OwnerOf(const MutableLocalPotentialView & view)
{
    if (view.GetAtomObjectPtr() != nullptr)
    {
        return ModelAnalysisData::OwnerOf(*view.GetAtomObjectPtr());
    }
    if (view.GetBondObjectPtr() != nullptr)
    {
        return ModelAnalysisData::OwnerOf(*view.GetBondObjectPtr());
    }
    return nullptr;
}

LocalPotentialEntry * FindResolvedLocalEntry(const MutableLocalPotentialView & view)
{
    if (view.GetEntryHandle() != nullptr)
    {
        return static_cast<LocalPotentialEntry *>(const_cast<void *>(view.GetEntryHandle()));
    }
    auto * owner{ OwnerOf(view) };
    if (owner == nullptr)
    {
        return nullptr;
    }
    auto & analysis_data{ ModelAnalysisData::Of(*owner) };
    if (view.GetAtomObjectPtr() != nullptr)
    {
        return analysis_data.FindAtomLocalEntry(*view.GetAtomObjectPtr());
    }
    if (view.GetBondObjectPtr() != nullptr)
    {
        return analysis_data.FindBondLocalEntry(*view.GetBondObjectPtr());
    }
    return nullptr;
}

LocalPotentialEntry & EnsureResolvedLocalEntry(const MutableLocalPotentialView & view)
{
    if (view.GetEntryHandle() != nullptr)
    {
        return *static_cast<LocalPotentialEntry *>(const_cast<void *>(view.GetEntryHandle()));
    }
    auto * owner{ OwnerOf(view) };
    if (owner == nullptr)
    {
        throw std::runtime_error("Local analysis owner model is not available.");
    }
    auto & analysis_data{ ModelAnalysisData::Of(*owner) };
    if (view.GetAtomObjectPtr() != nullptr)
    {
        return analysis_data.EnsureAtomLocalEntry(*view.GetAtomObjectPtr());
    }
    if (view.GetBondObjectPtr() != nullptr)
    {
        return analysis_data.EnsureBondLocalEntry(*view.GetBondObjectPtr());
    }
    throw std::runtime_error("Mutable local potential target is not available.");
}

const LocalPotentialEntry & RequireResolvedLocalEntry(
    const MutableLocalPotentialView & view,
    const char * context)
{
    const auto * entry{ FindResolvedLocalEntry(view) };
    if (entry == nullptr)
    {
        throw std::runtime_error(std::string(context) + " is not available.");
    }
    return *entry;
}

LocalPotentialEntry::FitResult ToDetailFitResult(LocalPotentialFitResult value)
{
    return LocalPotentialEntry::FitResult{
        std::move(value.beta_ols),
        std::move(value.beta_mdpde),
        value.sigma_square,
        std::move(value.data_weight),
        std::move(value.data_covariance)
    };
}

LocalPotentialAnnotation ToDetailAnnotation(const LocalPotentialAnnotationData & value)
{
    return LocalPotentialAnnotation{
        value.posterior,
        value.is_outlier,
        value.statistical_distance
    };
}

} // namespace

MutableLocalPotentialView::MutableLocalPotentialView(AtomObject * atom_object) :
    m_atom_object{ atom_object }
{
}

MutableLocalPotentialView::MutableLocalPotentialView(BondObject * bond_object) :
    m_bond_object{ bond_object }
{
}

void MutableLocalPotentialView::SetSamplingEntries(LocalPotentialSampleList value)
{
    EnsureResolvedLocalEntry(*this).SetSamplingEntries(std::move(value));
}

void MutableLocalPotentialView::SetDataset(LocalPotentialDataset value)
{
    EnsureResolvedLocalEntry(*this).SetDataset(std::move(value.basis_and_response_entry_list));
}

void MutableLocalPotentialView::SetFitResult(LocalPotentialFitResult value)
{
    EnsureResolvedLocalEntry(*this).SetFitResult(ToDetailFitResult(std::move(value)));
}

void MutableLocalPotentialView::SetEstimates(const LocalPotentialEstimates & value)
{
    auto & entry{ EnsureResolvedLocalEntry(*this) };
    entry.SetEstimateOLS(value.ols);
    entry.SetEstimateMDPDE(value.mdpde);
}

void MutableLocalPotentialView::SetAlphaR(double value)
{
    EnsureResolvedLocalEntry(*this).SetAlphaR(value);
}

void MutableLocalPotentialView::SetAnnotation(
    const std::string & key,
    const LocalPotentialAnnotationData & value)
{
    EnsureResolvedLocalEntry(*this).SetAnnotation(key, ToDetailAnnotation(value));
}

bool MutableLocalPotentialView::HasDataset() const
{
    const auto * entry{ FindResolvedLocalEntry(*this) };
    return entry != nullptr && entry->HasDataset();
}

bool MutableLocalPotentialView::HasFitResult() const
{
    const auto * entry{ FindResolvedLocalEntry(*this) };
    return entry != nullptr && entry->HasFitResult();
}

const LocalPotentialDataset & MutableLocalPotentialView::GetDataset() const
{
    const auto & dataset{ RequireResolvedLocalEntry(*this, "Local dataset").GetDataset() };
    m_dataset_cache.basis_and_response_entry_list = dataset.basis_and_response_entry_list;
    return m_dataset_cache;
}

const LocalPotentialFitResult & MutableLocalPotentialView::GetFitResult() const
{
    const auto & fit_result{ RequireResolvedLocalEntry(*this, "Local fit result").GetFitResult() };
    m_fit_result_cache.beta_ols = fit_result.beta_ols;
    m_fit_result_cache.beta_mdpde = fit_result.beta_mdpde;
    m_fit_result_cache.sigma_square = fit_result.sigma_square;
    m_fit_result_cache.data_weight = fit_result.data_weight;
    m_fit_result_cache.data_covariance = fit_result.data_covariance;
    return m_fit_result_cache;
}

ModelAnalysisEditor::ModelAnalysisEditor(ModelObject & model_object) :
    m_model_object{ model_object }
{
}

ModelAnalysisEditor ModelAnalysisEditor::Of(ModelObject & model_object)
{
    return ModelAnalysisEditor(model_object);
}

void ModelAnalysisEditor::Clear()
{
    ModelAnalysisData::Of(m_model_object).Clear();
}

void ModelAnalysisEditor::ClearTransientFitStates()
{
    ModelAnalysisData::Of(m_model_object).ClearTransientFitStates();
}

MutableLocalPotentialView ModelAnalysisEditor::EnsureAtomLocalPotential(const AtomObject & atom_object)
{
    auto & mutable_atom{ const_cast<AtomObject &>(atom_object) };
    auto & entry{ ModelAnalysisData::Of(m_model_object).EnsureAtomLocalEntry(mutable_atom) };
    auto view{ MutableLocalPotentialView(&mutable_atom) };
    view.m_entry_ptr = &entry;
    return view;
}

MutableLocalPotentialView ModelAnalysisEditor::EnsureBondLocalPotential(const BondObject & bond_object)
{
    auto & mutable_bond{ const_cast<BondObject &>(bond_object) };
    auto & entry{ ModelAnalysisData::Of(m_model_object).EnsureBondLocalEntry(mutable_bond) };
    auto view{ MutableLocalPotentialView(&mutable_bond) };
    view.m_entry_ptr = &entry;
    return view;
}

void ModelAnalysisEditor::RebuildAtomGroupsFromSelection()
{
    ModelAnalysisData::Of(m_model_object).RebuildAtomGroupEntriesFromSelection(m_model_object);
}

void ModelAnalysisEditor::RebuildBondGroupsFromSelection()
{
    ModelAnalysisData::Of(m_model_object).RebuildBondGroupEntriesFromSelection(m_model_object);
}

void ModelAnalysisEditor::SetAtomGroupStatistics(
    GroupKey group_key,
    const std::string & class_key,
    const GroupPotentialStatistics & statistics)
{
    ModelAnalysisData::Of(m_model_object).EnsureAtomGroupEntry(class_key).SetGroupStatistics(
        group_key,
        statistics.mean,
        statistics.mdpde,
        statistics.prior,
        statistics.prior_variance,
        statistics.alpha_g);
}

void ModelAnalysisEditor::SetBondGroupStatistics(
    GroupKey group_key,
    const std::string & class_key,
    const GroupPotentialStatistics & statistics)
{
    ModelAnalysisData::Of(m_model_object).EnsureBondGroupEntry(class_key).SetGroupStatistics(
        group_key,
        statistics.mean,
        statistics.mdpde,
        statistics.prior,
        statistics.prior_variance,
        statistics.alpha_g);
}

void ModelAnalysisEditor::SetAtomGroupAlphaG(
    GroupKey group_key,
    const std::string & class_key,
    double alpha_g)
{
    ModelAnalysisData::Of(m_model_object).EnsureAtomGroupEntry(class_key).SetAlphaG(group_key, alpha_g);
}

void ModelAnalysisEditor::SetBondGroupAlphaG(
    GroupKey group_key,
    const std::string & class_key,
    double alpha_g)
{
    ModelAnalysisData::Of(m_model_object).EnsureBondGroupEntry(class_key).SetAlphaG(group_key, alpha_g);
}

} // namespace rhbm_gem
