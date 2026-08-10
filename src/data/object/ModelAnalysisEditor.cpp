#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/AtomClassifier.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <stdexcept>

namespace rhbm_gem {

ModelAnalysisEditor::ModelAnalysisEditor(ModelObject & model_object) :
    m_model_object{ model_object }
{
}

void ModelAnalysisEditor::Clear()
{
    ModelAnalysisData::Of(m_model_object).Clear();
}

void ModelAnalysisEditor::ClearTransientFitStates()
{
    for (auto & [serial_id, entry] : ModelAnalysisData::Of(m_model_object).AtomLocalEntries())
    {
        (void)serial_id;
        if (entry != nullptr)
        {
            entry->ClearTransientFitState(LocalFittingStage::First);
            entry->ClearTransientFitState(LocalFittingStage::Second);
            entry->ClearTransientFitState(LocalFittingStage::Third);
        }
    }
}

AtomLocalPotentialEditor ModelAnalysisEditor::EnsureAtomLocalPotential(const AtomObject & atom_object)
{
    auto & entry{ ModelAnalysisData::Of(m_model_object).EnsureAtomLocalEntry(atom_object) };
    return AtomLocalPotentialEditor(entry);
}

void ModelAnalysisEditor::RebuildAtomGroupsFromSelection()
{
    auto & analysis_data{ ModelAnalysisData::Of(m_model_object) };
    auto & group_entry{ analysis_data.AtomGroupEntry() };
    group_entry = AtomGroupPotentialEntry{};
    for (auto * atom : m_model_object.GetSelectedAtoms())
    {
        const auto group_key{ data_internal::GetGroupKey(atom) };
        group_entry.AddMember(group_key, *atom);
    }
}

void ModelAnalysisEditor::InitializeLocalAlpha(
    LocalFittingStage stage,
    double alpha_r)
{
    for (auto * atom : m_model_object.GetSelectedAtoms())
    {
        EnsureAtomLocalPotential(*atom).SetAlphaR(stage, alpha_r);
    }
}

void ModelAnalysisEditor::InitializeGroupAlpha(double alpha_g)
{
    auto & group_entry{ ModelAnalysisData::Of(m_model_object).AtomGroupEntry() };
    for (const auto group_key : group_entry.CollectGroupKeys())
    {
        group_entry.SetAlphaG(group_key, alpha_g);
    }
}

void ModelAnalysisEditor::ApplyAtomGroupGaussianResult(
    LocalFittingStage stage,
    GroupKey group_key,
    const GroupGaussianResult & group_result)
{
    auto & analysis_data{ ModelAnalysisData::Of(m_model_object) };
    auto & group_entry{ analysis_data.AtomGroupEntry() };
    if (!group_entry.HasGroup(group_key))
    {
        throw std::runtime_error("Atom group entry is not available.");
    }
    const auto & atom_list{ group_entry.GetMembers(group_key) };
    if (group_result.member_results.size() != atom_list.size())
    {
        throw std::invalid_argument("Atom group result member result count is inconsistent.");
    }

    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        const auto & member_result{ group_result.member_results.at(i) };
        auto & atom_entry{ analysis_data.EnsureAtomLocalEntry(*atom_list.at(i)) };
        atom_entry.SetPosteriorResult(
            stage,
            member_result.mdpde,
            member_result.is_outlier,
            member_result.statistical_distance);
    }
    group_entry.SetGaussianResult(group_key, group_result);
}

void ModelAnalysisEditor::SetAtomGroupAlphaG(
    GroupKey group_key,
    double alpha_g)
{
    ModelAnalysisData::Of(m_model_object).AtomGroupEntry().SetAlphaG(group_key, alpha_g);
}

} // namespace rhbm_gem
