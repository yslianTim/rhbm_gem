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

LocalPotentialEditor::LocalPotentialEditor(AtomObject * atom_object) :
    m_atom_object{ atom_object }
{
}

LocalPotentialEditor::LocalPotentialEditor(BondObject * bond_object) :
    m_bond_object{ bond_object }
{
}

ModelObject * LocalPotentialEditor::GetOwner() const
{
    if (m_atom_object != nullptr)
    {
        return ModelAnalysisData::OwnerOf(*m_atom_object);
    }
    if (m_bond_object != nullptr)
    {
        return ModelAnalysisData::OwnerOf(*m_bond_object);
    }
    return nullptr;
}

LocalPotentialEntry & LocalPotentialEditor::EnsureEntry() const
{
    if (m_entry != nullptr)
    {
        return *m_entry;
    }
    auto * owner{ GetOwner() };
    if (owner == nullptr)
    {
        throw std::runtime_error("Local analysis owner model is not available.");
    }
    auto & analysis_data{ ModelAnalysisData::Of(*owner) };
    if (m_atom_object != nullptr)
    {
        return analysis_data.EnsureAtomLocalEntry(*m_atom_object);
    }
    if (m_bond_object != nullptr)
    {
        return analysis_data.EnsureBondLocalEntry(*m_bond_object);
    }
    throw std::runtime_error("Local potential edit target is not available.");
}

void LocalPotentialEditor::SetSamplingEntries(LocalPotentialSampleList value)
{
    EnsureEntry().SetSamplingEntries(std::move(value));
}

void LocalPotentialEditor::SetGaussianResult(LocalGaussianResult value)
{
    EnsureEntry().SetGaussianResult(std::move(value));
}

void LocalPotentialEditor::SetAlphaR(double value)
{
    EnsureEntry().SetAlphaR(value);
}

void LocalPotentialEditor::SetAnnotation(
    const std::string & key,
    const LocalPotentialAnnotation & value)
{
    EnsureEntry().SetAnnotation(key, value);
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

LocalPotentialEditor ModelAnalysisEditor::EnsureAtomLocalPotential(const AtomObject & atom_object)
{
    auto & mutable_atom{ const_cast<AtomObject &>(atom_object) };
    auto & entry{ ModelAnalysisData::Of(m_model_object).EnsureAtomLocalEntry(mutable_atom) };
    auto editor{ LocalPotentialEditor(&mutable_atom) };
    editor.m_entry = &entry;
    return editor;
}

LocalPotentialEditor ModelAnalysisEditor::EnsureBondLocalPotential(const BondObject & bond_object)
{
    auto & mutable_bond{ const_cast<BondObject &>(bond_object) };
    auto & entry{ ModelAnalysisData::Of(m_model_object).EnsureBondLocalEntry(mutable_bond) };
    auto editor{ LocalPotentialEditor(&mutable_bond) };
    editor.m_entry = &entry;
    return editor;
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
