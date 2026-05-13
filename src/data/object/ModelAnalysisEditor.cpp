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
        value.gaussian,
        value.is_outlier,
        value.statistical_distance
    };
}

LocalPotentialAnnotation ToDetailAnnotation(const LocalGaussianResult & value)
{
    return LocalPotentialAnnotation{
        value.mdpde,
        value.is_outlier,
        value.statistical_distance
    };
}

void ValidateGroupGaussianResult(
    const std::vector<LocalGaussianResult> & member_results,
    std::size_t member_count,
    const char * context)
{
    if (member_results.size() != member_count)
    {
        throw std::invalid_argument(
            std::string(context) + " member result count is inconsistent.");
    }
}

void ApplyAtomGroupGaussianResultToEntry(
    ModelAnalysisData & analysis_data,
    GroupKey group_key,
    const std::string & class_key,
    const GroupGaussianResult & result)
{
    analysis_data.EnsureAtomGroupEntry(class_key).SetGroupStatistics(
        group_key,
        result.mean,
        result.mdpde,
        result.prior.GetModel(),
        result.prior.GetStandardDeviationModel(),
        result.alpha_g);
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

void MutableLocalPotentialView::SetGaussianResult(LocalGaussianResult value)
{
    EnsureResolvedLocalEntry(*this).SetGaussianResult(std::move(value));
}

void MutableLocalPotentialView::SetAlphaR(double value)
{
    EnsureResolvedLocalEntry(*this).SetAlphaR(value);
}

double MutableLocalPotentialView::GetAlphaR() const
{
    return RequireResolvedLocalEntry(*this, "Local alpha-r").GetAlphaR();
}

void MutableLocalPotentialView::SetAnnotation(
    const std::string & key,
    const LocalPotentialAnnotationData & value)
{
    EnsureResolvedLocalEntry(*this).SetAnnotation(key, ToDetailAnnotation(value));
}

const LocalGaussianResult & MutableLocalPotentialView::GetGaussianResult() const
{
    return RequireResolvedLocalEntry(*this, "Local Gaussian result").GetGaussianResult();
}

const LocalPotentialSampleList & MutableLocalPotentialView::GetSamplingEntries() const
{
    return RequireResolvedLocalEntry(*this, "Local sampling entries").GetSamplingEntries();
}

int MutableLocalPotentialView::GetSamplingEntryCount() const
{
    return RequireResolvedLocalEntry(*this, "Local sampling entry count").GetSamplingEntryCount();
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

void ModelAnalysisEditor::ApplyAtomGroupGaussianResult(
    GroupKey group_key,
    const std::string & class_key,
    const GroupGaussianResult & group_result)
{
    auto & analysis_data{ ModelAnalysisData::Of(m_model_object) };
    const auto & group_entry{
        RequireGroupEntry(analysis_data.FindAtomGroupEntry(class_key), "Atom group entry")
    };
    const auto & atom_list{ group_entry.GetMembers(group_key) };
    ValidateGroupGaussianResult(group_result.member_results, atom_list.size(), "Atom group result");

    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        auto & atom_entry{ analysis_data.EnsureAtomLocalEntry(*atom_list.at(i)) };
        atom_entry.SetAnnotation(class_key, ToDetailAnnotation(group_result.member_results.at(i)));
    }

    ApplyAtomGroupGaussianResultToEntry(analysis_data, group_key, class_key, group_result);
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
