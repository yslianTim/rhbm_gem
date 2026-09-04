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

namespace {

constexpr double kInitialLocalAlpha{ 0.0 };
constexpr double kInitialGroupAlpha{ 0.0 };

LocalPotentialEntry & EnsureAtomLocalPotential(
    ModelObject & model_object,
    const AtomObject & atom_object)
{
    return ModelAnalysisData::Of(model_object).EnsureAtomLocalEntry(atom_object);
}

} // namespace

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

void ModelAnalysisEditor::InitializeFromSelection()
{
    Clear();
    RebuildAtomGroupsFromSelection();
    for (const auto stage : {
            FittingStage::First,
            FittingStage::Second,
            FittingStage::Third })
    {
        InitializeLocalAlpha(stage, kInitialLocalAlpha);
    }
    InitializeGroupAlpha(kInitialGroupAlpha);
}

void ModelAnalysisEditor::InitializeLocalFittingSeedModels()
{
    const auto seed_model{ GaussianModel3D{ 0.0, 1.0, 0.0 } };
    for (auto * atom : m_model_object.GetSelectedAtoms())
    {
        EnsureAtomLocalPotential(m_model_object, *atom).ClearGroupMemberResult();
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        auto result{ local_view.GetGaussianResult(FittingStage::First) };
        result.ols = GaussianModel3DWithUncertainty{
            seed_model,
            GaussianModel3DUncertainty{}
        };
        result.mdpde = GaussianModel3DWithUncertainty{
            seed_model,
            GaussianModel3DUncertainty{}
        };
        result.fit_result.reset();
        SetAtomLocalGaussianResult(FittingStage::First, *atom, result);
        SetAtomLocalGaussianResult(FittingStage::Second, *atom, result);
        SetAtomLocalGaussianResult(FittingStage::Third, *atom, std::move(result));
    }
}

void ModelAnalysisEditor::EnsureSelectedAtomLocalPotentials()
{
    for (auto * atom : m_model_object.GetSelectedAtoms())
    {
        EnsureAtomLocalPotential(m_model_object, *atom);
    }
}

void ModelAnalysisEditor::EnsureAtomGroupLocalPotentials(GroupKey group_key)
{
    const auto & atom_list{
        ModelAnalysisData::Of(m_model_object).AtomGroupEntry().GetMembers(group_key)
    };
    for (auto * atom : atom_list)
    {
        EnsureAtomLocalPotential(m_model_object, *atom);
    }
}

void ModelAnalysisEditor::SetAtomLocalRawSamplingEntries(
    const AtomObject & atom_object,
    LocalPotentialSampleList value)
{
    EnsureAtomLocalPotential(m_model_object, atom_object)
        .SetRawSamplingEntries(std::move(value));
}

void ModelAnalysisEditor::SetAtomLocalPeelingSamplingEntries(
    const AtomObject & atom_object,
    LocalPotentialSampleList value)
{
    EnsureAtomLocalPotential(m_model_object, atom_object)
        .SetPeelingSamplingEntries(std::move(value));
}

void ModelAnalysisEditor::SetAtomLocalGaussianResult(
    FittingStage stage,
    const AtomObject & atom_object,
    LocalGaussianResult result)
{
    EnsureAtomLocalPotential(m_model_object, atom_object)
        .SetGaussianResult(stage, std::move(result));
}

void ModelAnalysisEditor::SetAtomLocalAlphaR(
    FittingStage stage,
    const AtomObject & atom_object,
    double alpha_r)
{
    EnsureAtomLocalPotential(m_model_object, atom_object).SetAlphaR(stage, alpha_r);
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

void ModelAnalysisEditor::InitializeLocalAlpha(FittingStage stage, double alpha_r)
{
    for (auto * atom : m_model_object.GetSelectedAtoms())
    {
        SetAtomLocalAlphaR(stage, *atom, alpha_r);
    }
}

void ModelAnalysisEditor::SetAtomGroupAlphaR(
    FittingStage stage,
    GroupKey group_key,
    double alpha_r)
{
    const auto & atom_list{
        ModelAnalysisData::Of(m_model_object).AtomGroupEntry().GetMembers(group_key)
    };
    for (auto * atom : atom_list)
    {
        SetAtomLocalAlphaR(stage, *atom, alpha_r);
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

void ModelAnalysisEditor::CopyLocalFittingStageResult(
    FittingStage source_stage,
    FittingStage destination_stage)
{
    for (auto * atom : m_model_object.GetSelectedAtoms())
    {
        const auto local_view{ AtomLocalPotentialView::For(*atom) };
        SetAtomLocalGaussianResult(
            destination_stage, *atom, local_view.GetGaussianResult(source_stage));
    }
}

void ModelAnalysisEditor::ApplyAtomGroupGaussianResult(
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
        atom_entry.SetGroupMemberResult(member_result);
    }
    group_entry.SetGaussianResult(group_key, group_result);
}

void ModelAnalysisEditor::SetAtomLocalNeighborCountForPeeling(
    const AtomObject & atom_object,
    int neighbor_count)
{
    EnsureAtomLocalPotential(m_model_object, atom_object)
        .SetNeighborCountForPeeling(neighbor_count);
}

void ModelAnalysisEditor::ApplyAtomLocalSecondStageResult(
    const AtomObject & atom_object,
    LocalGaussianResult result,
    LocalPotentialSampleList peeling_sampling_entries)
{
    SetAtomLocalGaussianResult(FittingStage::Second, atom_object, std::move(result));
    SetAtomLocalPeelingSamplingEntries(atom_object, std::move(peeling_sampling_entries));
}

void ModelAnalysisEditor::SetAtomGroupAlphaG(
    GroupKey group_key,
    double alpha_g)
{
    ModelAnalysisData::Of(m_model_object).AtomGroupEntry().SetAlphaG(group_key, alpha_g);
}

} // namespace rhbm_gem
