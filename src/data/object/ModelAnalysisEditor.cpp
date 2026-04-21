#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/hrl/GaussianLinearizationService.hpp>

#include <stdexcept>
#include <string>

namespace rhbm_gem {

namespace {

struct LocalPotentialEstimates
{
    GaussianEstimate ols{};
    GaussianEstimate mdpde{};
};

GaussianLinearizationContext BuildLocalDecodeContext(const LocalPotentialEntry & entry)
{
    Eigen::VectorXd model_par_init{ Eigen::VectorXd::Zero(3) };
    model_par_init(0) = entry.GetMomentZeroEstimate();
    model_par_init(1) = entry.GetMomentTwoEstimate();
    return GaussianLinearizationContext::FromModelParameters(model_par_init);
}

const GaussianLinearizationService & LocalDecodeService()
{
    static const GaussianLinearizationService service{
        GaussianLinearizationSpec::AtomLocalDecode()
    };
    return service;
}

const GaussianLinearizationService & AtomGroupDecodeService()
{
    static const GaussianLinearizationService service{
        GaussianLinearizationSpec::AtomGroupDecode()
    };
    return service;
}

const GaussianLinearizationService & BondGroupDecodeService()
{
    static const GaussianLinearizationService service{
        GaussianLinearizationSpec::BondGroupDecode()
    };
    return service;
}

template <typename EntryT>
const EntryT & RequireGroupEntry(const EntryT * entry, const char * context)
{
    if (entry == nullptr)
    {
        throw std::runtime_error(std::string(context) + " is not available.");
    }
    return *entry;
}

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

LocalPotentialAnnotation ToDetailAnnotation(const LocalPotentialAnnotationData & value)
{
    return LocalPotentialAnnotation{
        value.posterior,
        value.is_outlier,
        value.statistical_distance
    };
}

LocalPotentialEstimates BuildLocalPotentialEstimates(
    const LocalPotentialEntry & entry,
    const HRLBetaEstimateResult & value)
{
    const auto context{ BuildLocalDecodeContext(entry) };
    const auto gaus_ols{ LocalDecodeService().DecodeLocalEstimate(value.beta_ols, context) };
    const auto gaus_mdpde{ LocalDecodeService().DecodeLocalEstimate(value.beta_mdpde, context) };

    return LocalPotentialEstimates{
        gaus_ols,
        gaus_mdpde
    };
}

void ValidateGroupEstimateResult(
    const HRLGroupEstimationResult & result,
    std::size_t member_count,
    const char * context)
{
    const auto expected_member_count{ static_cast<Eigen::Index>(member_count) };
    if (result.beta_posterior_matrix.cols() != expected_member_count)
    {
        throw std::invalid_argument(
            std::string(context) + " beta_posterior_matrix member count is inconsistent.");
    }
    if (result.capital_sigma_posterior_list.size() != member_count)
    {
        throw std::invalid_argument(
            std::string(context) + " capital_sigma_posterior_list member count is inconsistent.");
    }
    if (result.outlier_flag_array.rows() != expected_member_count)
    {
        throw std::invalid_argument(
            std::string(context) + " outlier_flag_array member count is inconsistent.");
    }
    if (result.statistical_distance_array.rows() != expected_member_count)
    {
        throw std::invalid_argument(
            std::string(context) + " statistical_distance_array member count is inconsistent.");
    }
}

void ApplyAtomGroupStatistics(
    ModelAnalysisData & analysis_data,
    GroupKey group_key,
    const std::string & class_key,
    const HRLGroupEstimationResult & result,
    double alpha_g)
{
    const auto gaus_group_mean{ AtomGroupDecodeService().DecodeGroupEstimate(result.mu_mean) };
    const auto gaus_group_mdpde{ AtomGroupDecodeService().DecodeGroupEstimate(result.mu_mdpde) };
    const auto gaus_prior{
        AtomGroupDecodeService().DecodePosteriorEstimate(
            result.mu_prior,
            result.capital_lambda)
    };

    analysis_data.EnsureAtomGroupEntry(class_key).SetGroupStatistics(
        group_key,
        gaus_group_mean,
        gaus_group_mdpde,
        gaus_prior.estimate,
        gaus_prior.variance,
        alpha_g);
}

void ApplyBondGroupStatistics(
    ModelAnalysisData & analysis_data,
    GroupKey group_key,
    const std::string & class_key,
    const HRLGroupEstimationResult & result,
    double alpha_g)
{
    const auto gaus_group_mean{ BondGroupDecodeService().DecodeGroupEstimate(result.mu_mean) };
    const auto gaus_group_mdpde{ BondGroupDecodeService().DecodeGroupEstimate(result.mu_mdpde) };
    const auto gaus_prior{
        BondGroupDecodeService().DecodePosteriorEstimate(
            result.mu_prior,
            result.capital_lambda)
    };

    analysis_data.EnsureBondGroupEntry(class_key).SetGroupStatistics(
        group_key,
        gaus_group_mean,
        gaus_group_mdpde,
        gaus_prior.estimate,
        gaus_prior.variance,
        alpha_g);
}

LocalPotentialAnnotationData BuildAtomAnnotationData(
    const HRLGroupEstimationResult & result,
    Eigen::Index member_index)
{
    const auto beta_vector_posterior{ result.beta_posterior_matrix.col(member_index) };
    const auto & sigma_matrix_posterior{
        result.capital_sigma_posterior_list.at(static_cast<std::size_t>(member_index))
    };
    const auto gaus_posterior{
        AtomGroupDecodeService().DecodePosteriorEstimate(
            beta_vector_posterior,
            sigma_matrix_posterior)
    };
    return LocalPotentialAnnotationData{
        gaus_posterior,
        static_cast<bool>(result.outlier_flag_array(member_index)),
        result.statistical_distance_array(member_index)
    };
}

LocalPotentialAnnotationData BuildBondAnnotationData(
    const HRLGroupEstimationResult & result,
    Eigen::Index member_index)
{
    const auto beta_vector_posterior{ result.beta_posterior_matrix.col(member_index) };
    const auto & sigma_matrix_posterior{
        result.capital_sigma_posterior_list.at(static_cast<std::size_t>(member_index))
    };
    const auto gaus_posterior{
        BondGroupDecodeService().DecodePosteriorEstimate(
            beta_vector_posterior,
            sigma_matrix_posterior)
    };
    return LocalPotentialAnnotationData{
        gaus_posterior,
        static_cast<bool>(result.outlier_flag_array(member_index)),
        result.statistical_distance_array(member_index)
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

void MutableLocalPotentialView::SetDataset(HRLMemberDataset value)
{
    EnsureResolvedLocalEntry(*this).SetDataset(std::move(value));
}

void MutableLocalPotentialView::SetFitResult(HRLBetaEstimateResult value)
{
    auto & entry{ EnsureResolvedLocalEntry(*this) };
    const auto estimates{ BuildLocalPotentialEstimates(entry, value) };

    entry.SetFitResult(std::move(value));
    entry.SetEstimateOLS(estimates.ols);
    entry.SetEstimateMDPDE(estimates.mdpde);
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

const HRLMemberDataset & MutableLocalPotentialView::GetDataset() const
{
    const auto & dataset{ RequireResolvedLocalEntry(*this, "Local dataset").GetDataset() };
    m_dataset_cache = dataset;
    return m_dataset_cache;
}

const HRLBetaEstimateResult & MutableLocalPotentialView::GetFitResult() const
{
    const auto & fit_result{ RequireResolvedLocalEntry(*this, "Local fit result").GetFitResult() };
    m_fit_result_cache = fit_result;
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

void ModelAnalysisEditor::ApplyAtomGroupEstimateResult(
    GroupKey group_key,
    const std::string & class_key,
    const HRLGroupEstimationResult & result,
    double alpha_g)
{
    auto & analysis_data{ ModelAnalysisData::Of(m_model_object) };
    const auto & group_entry{
        RequireGroupEntry(analysis_data.FindAtomGroupEntry(class_key), "Atom group entry")
    };
    const auto & atom_list{ group_entry.GetMembers(group_key) };
    ValidateGroupEstimateResult(result, atom_list.size(), "Atom group result");

    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(atom_list.size()); i++)
    {
        auto & atom_entry{ analysis_data.EnsureAtomLocalEntry(*atom_list[static_cast<std::size_t>(i)]) };
        atom_entry.SetAnnotation(class_key, ToDetailAnnotation(BuildAtomAnnotationData(result, i)));
    }

    ApplyAtomGroupStatistics(analysis_data, group_key, class_key, result, alpha_g);
}

void ModelAnalysisEditor::ApplyBondGroupEstimateResult(
    GroupKey group_key,
    const std::string & class_key,
    const HRLGroupEstimationResult & result,
    double alpha_g)
{
    auto & analysis_data{ ModelAnalysisData::Of(m_model_object) };
    const auto & group_entry{
        RequireGroupEntry(analysis_data.FindBondGroupEntry(class_key), "Bond group entry")
    };
    const auto & bond_list{ group_entry.GetMembers(group_key) };
    ValidateGroupEstimateResult(result, bond_list.size(), "Bond group result");

    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(bond_list.size()); i++)
    {
        auto & bond_entry{ analysis_data.EnsureBondLocalEntry(*bond_list[static_cast<std::size_t>(i)]) };
        bond_entry.SetAnnotation(class_key, ToDetailAnnotation(BuildBondAnnotationData(result, i)));
    }

    ApplyBondGroupStatistics(analysis_data, group_key, class_key, result, alpha_g);
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
