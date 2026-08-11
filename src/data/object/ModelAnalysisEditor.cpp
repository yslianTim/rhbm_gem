#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/AtomClassifier.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <stdexcept>
#include <utility>

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
            entry->ClearTransientFitState(FittingStage::First);
            entry->ClearTransientFitState(FittingStage::Second);
            entry->ClearTransientFitState(FittingStage::Third);
        }
    }
}

void ModelAnalysisEditor::InitializeLocalFittingSeedModels()
{
    const auto seed_model{ GaussianModel3D{ 0.0, 1.0, 0.0 } };
    for (auto * atom : m_model_object.GetSelectedAtoms())
    {
        auto local_editor{ EnsureAtomLocalPotential(*atom) };
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        auto result{ local_view.GetGaussianResult(FittingStage::First) };
        result.ols = GaussianModel3DWithUncertainty{
            seed_model,
            GaussianModel3DUncertainty{}
        };
        result.mdpde = GaussianModel3DWithUncertainty{
            seed_model,
            GaussianModel3DUncertainty{}
        };
        result.posterior.reset();
        result.is_outlier = false;
        result.statistical_distance = 0.0;
        result.fit_result.reset();
        local_editor.SetGaussianResult(FittingStage::First, result);
        local_editor.SetGaussianResult(FittingStage::Second, result);
        local_editor.SetGaussianResult(FittingStage::Third, result);
    }
}

AtomLocalPotentialEditor ModelAnalysisEditor::EnsureAtomLocalPotential(const AtomObject & atom_object)
{
    auto & entry{ ModelAnalysisData::Of(m_model_object).EnsureAtomLocalEntry(atom_object) };
    return AtomLocalPotentialEditor(entry);
}

void ModelAnalysisEditor::EnsureSelectedAtomLocalPotentials()
{
    for (auto * atom : m_model_object.GetSelectedAtoms())
    {
        EnsureAtomLocalPotential(*atom);
    }
}

void ModelAnalysisEditor::EnsureAtomGroupLocalPotentials(FittingStage stage, GroupKey group_key)
{
    const auto & atom_list{
        ModelAnalysisData::Of(m_model_object).AtomGroupEntry().GetMembers(stage, group_key)
    };
    for (auto * atom : atom_list)
    {
        EnsureAtomLocalPotential(*atom);
    }
}

AtomLocalPotentialEditor ModelAnalysisEditor::GetAtomLocalPotentialEditor(
    const AtomObject & atom_object) const
{
    auto * entry{ ModelAnalysisData::Of(m_model_object).FindAtomLocalEntry(atom_object) };
    if (entry == nullptr)
    {
        throw std::runtime_error("Atom local potential entry is not available.");
    }
    return AtomLocalPotentialEditor(*entry);
}

void ModelAnalysisEditor::RebuildAtomGroupsFromSelection()
{
    auto & analysis_data{ ModelAnalysisData::Of(m_model_object) };
    auto & group_entry{ analysis_data.AtomGroupEntry() };
    group_entry = AtomGroupPotentialEntry{};
    for (auto * atom : m_model_object.GetSelectedAtoms())
    {
        const auto group_key{ data_internal::GetGroupKey(atom) };
        group_entry.AddMember(FittingStage::First, group_key, *atom);
        group_entry.AddMember(FittingStage::Second, group_key, *atom);
        group_entry.AddMember(FittingStage::Third, group_key, *atom);
    }
}

void ModelAnalysisEditor::InitializeLocalAlpha(FittingStage stage, double alpha_r)
{
    for (auto * atom : m_model_object.GetSelectedAtoms())
    {
        EnsureAtomLocalPotential(*atom).SetAlphaR(stage, alpha_r);
    }
}

void ModelAnalysisEditor::SetAtomGroupAlphaR(
    FittingStage stage,
    GroupKey group_key,
    double alpha_r)
{
    const auto & atom_list{
        ModelAnalysisData::Of(m_model_object).AtomGroupEntry().GetMembers(stage, group_key)
    };
    for (auto * atom : atom_list)
    {
        EnsureAtomLocalPotential(*atom).SetAlphaR(stage, alpha_r);
    }
}

void ModelAnalysisEditor::InitializeGroupAlpha(FittingStage stage, double alpha_g)
{
    auto & group_entry{ ModelAnalysisData::Of(m_model_object).AtomGroupEntry() };
    for (const auto group_key : group_entry.CollectGroupKeys(stage))
    {
        group_entry.SetAlphaG(stage, group_key, alpha_g);
    }
}

void ModelAnalysisEditor::CopyLocalFittingStageResult(
    FittingStage source_stage,
    FittingStage destination_stage)
{
    for (auto * atom : m_model_object.GetSelectedAtoms())
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        EnsureAtomLocalPotential(*atom).SetGaussianResult(
            destination_stage,
            local_view.GetGaussianResult(source_stage));
    }
}

void ModelAnalysisEditor::CopyGroupFittingStageResult(
    FittingStage source_stage,
    FittingStage destination_stage)
{
    ModelAnalysisData::Of(m_model_object).AtomGroupEntry().CopyStage(
        source_stage,
        destination_stage);
}

void ModelAnalysisEditor::CopyFittingStageState(
    FittingStage source_stage,
    FittingStage destination_stage)
{
    CopyLocalFittingStageResult(source_stage, destination_stage);
    CopyGroupFittingStageResult(source_stage, destination_stage);
}

void ModelAnalysisEditor::ApplyAtomGroupGaussianResult(
    FittingStage stage,
    GroupKey group_key,
    const GroupGaussianResult & group_result)
{
    auto & analysis_data{ ModelAnalysisData::Of(m_model_object) };
    auto & group_entry{ analysis_data.AtomGroupEntry() };
    if (!group_entry.HasGroup(stage, group_key))
    {
        throw std::runtime_error("Atom group entry is not available.");
    }
    const auto & atom_list{ group_entry.GetMembers(stage, group_key) };
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
    group_entry.SetGaussianResult(stage, group_key, group_result);
}

void ModelAnalysisEditor::ApplyAtomLocalGaussianResult(
    FittingStage stage,
    const AtomObject & atom_object,
    LocalGaussianResult result)
{
    GetAtomLocalPotentialEditor(atom_object).SetGaussianResult(stage, std::move(result));
}

void ModelAnalysisEditor::SetAtomLocalNeighborCountForPeeling(
    const AtomObject & atom_object,
    int neighbor_count)
{
    EnsureAtomLocalPotential(atom_object).SetNeighborCountForPeeling(neighbor_count);
}

void ModelAnalysisEditor::ApplyAtomLocalSecondStageResult(
    const AtomObject & atom_object,
    LocalGaussianResult result,
    LocalPotentialSampleList peeling_sampling_entries)
{
    auto local_editor{ EnsureAtomLocalPotential(atom_object) };
    local_editor.SetGaussianResult(FittingStage::Second, std::move(result));
    local_editor.SetPeelingSamplingEntries(std::move(peeling_sampling_entries));
}

void ModelAnalysisEditor::SetAtomGroupAlphaG(
    FittingStage stage,
    GroupKey group_key,
    double alpha_g)
{
    ModelAnalysisData::Of(m_model_object).AtomGroupEntry().SetAlphaG(stage, group_key, alpha_g);
}

} // namespace rhbm_gem
