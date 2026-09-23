#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/AtomClassifier.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <stdexcept>
#include <algorithm>
#include <set>
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

bool SameSource(const EstimateSource & a, const EstimateSource & b)
{
    return a.method == b.method && a.atom_id == b.atom_id && a.component_id == b.component_id &&
        a.run_id == b.run_id && a.role == b.role;
}

bool SameEndpoint(const LocalStageEstimate & a, const LocalStageEstimate & b)
{
    return SameSource(a.source, b.source) && a.point.has_value() == b.point.has_value() &&
        (!a.point || a.point->ToVector() == b.point->ToVector()) &&
        a.convergence == b.convergence && a.reason == b.reason;
}

bool SameUncertainty(const StageUncertainty & a, const StageUncertainty & b)
{
    return a.status == b.status && a.method == b.method && a.reason == b.reason &&
        a.residual_variance == b.residual_variance && a.rank == b.rank &&
        a.degrees_of_freedom == b.degrees_of_freedom && a.rank_threshold == b.rank_threshold &&
        a.covariance.has_value() == b.covariance.has_value() &&
        (!a.covariance || *a.covariance == *b.covariance);
}

void InvalidateGroup(ModelAnalysisData & data, GroupKey key)
{
    auto & groups = data.AtomGroupEntry();
    for (const auto * atom : groups.GetMembers(key))
        if (auto * entry = data.FindAtomLocalEntry(*atom)) entry->ClearGroupMemberResult();
    groups.ClearResult(key);
}

void InvalidateGroups(ModelAnalysisData & data, const std::set<int> & ids)
{
    auto & groups = data.AtomGroupEntry();
    for (const auto key : groups.CollectGroupKeys())
        if (std::any_of(groups.GetMembers(key).begin(), groups.GetMembers(key).end(),
            [&](const auto * atom) { return ids.contains(atom->GetSerialID()); }))
            InvalidateGroup(data, key);
}

void InvalidateJointRuns(ModelAnalysisData & data, const std::set<std::string> & runs)
{
    if (runs.empty()) return;
    std::set<int> affected;
    for (auto & [id, entry] : data.AtomLocalEntries())
    {
        auto stage = entry->StageEstimate(FittingStage::Second);
        if (stage.source.method != EstimateMethod::JointComponents || !runs.contains(stage.source.run_id)) continue;
        stage.uncertainty = {};
        entry->SetStageEstimate(FittingStage::Second, std::move(stage));
        entry->ClearGroupEvidence();
        entry->ClearGroupMemberResult();
        entry->ClearPeeling();
        affected.insert(id);
    }
    InvalidateGroups(data, affected);
    data.joint_result.reset();
}

void InvalidateStageChange(ModelAnalysisData & data, int id,
    const LocalStageEstimate & before, LocalStageEstimate & after)
{
    if (!SameEndpoint(before, after))
    {
        std::set<std::string> runs;
        for (const auto * stage : {&before, static_cast<const LocalStageEstimate *>(&after)})
            if (stage->source.method == EstimateMethod::JointComponents) runs.insert(stage->source.run_id);
        InvalidateJointRuns(data, runs);
        // An uncertainty copied along with a changed endpoint is no longer evidence for it.
        if (before.source.method == EstimateMethod::JointComponents &&
            SameUncertainty(before.uncertainty, after.uncertainty)) after.uncertainty = {};
        auto & entry = *data.AtomLocalEntries().at(id);
        entry.ClearGroupEvidence();
        if (!runs.empty()) entry.ClearPeeling();
        entry.ClearGroupMemberResult();
        InvalidateGroups(data, {id});
    }
    else if (!SameUncertainty(before.uncertainty, after.uncertainty))
    {
        data.AtomLocalEntries().at(id)->ClearGroupEvidence();
        data.AtomLocalEntries().at(id)->ClearGroupMemberResult();
        InvalidateGroups(data, {id});
    }
}

} // namespace

ModelAnalysisEditor::ModelAnalysisEditor(ModelObject & model_object) :
    m_model_object{ model_object }
{
}

void ModelAnalysisEditor::SetJointResult(JointAnalysisResult result)
{
    ModelAnalysisData::Of(m_model_object).joint_result=std::move(result);
}

void ModelAnalysisEditor::SetAtomStageEstimate(FittingStage stage, const AtomObject & atom, LocalStageEstimate value)
{
    value.source.atom_id = std::to_string(atom.GetSerialID());
    auto & entry = EnsureAtomLocalPotential(m_model_object, atom);
    const auto before = entry.StageEstimate(stage);
    if (stage == FittingStage::Second)
        InvalidateStageChange(ModelAnalysisData::Of(m_model_object), atom.GetSerialID(), before, value);
    entry.SetStageEstimate(stage, std::move(value));
}

void ModelAnalysisEditor::ApplySecondStageEstimates(std::map<int, LocalStageEstimate> estimates)
{
    // Validate the complete mapping before touching any published result.
    for (auto & [id, estimate] : estimates)
    {
        const AtomObject * atom = nullptr;
        try { atom = m_model_object.FindAtomPtr(id); }
        catch (const std::out_of_range &) { throw std::invalid_argument("Joint contributor does not belong to model."); }
        if (!atom || atom->GetElement() == Element::HYDROGEN)
            throw std::invalid_argument("Joint contributor does not belong to model.");
        if (estimate.point) GaussianModel3D::RequireFinitePositiveWidthModel(*estimate.point);
        estimate.source.atom_id = std::to_string(id);
    }
    auto & data = ModelAnalysisData::Of(m_model_object);
    std::set<std::string> runs;
    std::set<int> ids;
    for (const auto & [id, estimate] : estimates)
    {
        auto & entry = EnsureAtomLocalPotential(m_model_object, *m_model_object.FindAtomPtr(id));
        for (const auto * source : {&entry.StageEstimate(FittingStage::Second).source, &estimate.source})
            if (source->method == EstimateMethod::JointComponents) runs.insert(source->run_id);
        ids.insert(id);
    }
    InvalidateJointRuns(data, runs);
    InvalidateGroups(data, ids);
    for (auto & [id, estimate] : estimates)
    {
        auto & entry = *data.AtomLocalEntries().at(id);
        entry.ClearPeeling(); entry.ClearGroupEvidence(); entry.ClearGroupMemberResult();
        entry.SetStageEstimate(FittingStage::Second, std::move(estimate));
    }
}

void ModelAnalysisEditor::SetAtomPostFitPeeling(const AtomObject & atom, PostFitPeelingResult value)
{
    auto & entry = EnsureAtomLocalPotential(m_model_object, atom);
    const auto & raw = entry.RawSamplingEntries();
    if (raw.size() != value.samples.size()) throw std::invalid_argument("Peeling/raw sample count mismatch.");
    if (!SameSource(value.source, entry.StageEstimate(FittingStage::Second).source))
        throw std::invalid_argument("Peeling source mismatch.");
    entry.SetPostFitPeeling(std::move(value));
}

void ModelAnalysisEditor::SetAtomGroupEvidence(const AtomObject & atom, GroupParameterEvidence value)
{
    auto & entry = EnsureAtomLocalPotential(m_model_object, atom);
    entry.ClearGroupMemberResult();
    InvalidateGroups(ModelAnalysisData::Of(m_model_object), {atom.GetSerialID()});
    entry.SetGroupEvidence(std::move(value));
}

void ModelAnalysisEditor::ApplyAtomGroupParameterSummary(GroupKey key, GroupParameterSummary value)
{
    auto & groups = ModelAnalysisData::Of(m_model_object).AtomGroupEntry();
    const auto & members = groups.GetMembers(key);
    if (value.inference && value.inference->member_results.size() != value.member_ids.size())
        throw std::invalid_argument("Parameter posterior identity count mismatch.");
    if (std::set<int>(value.member_ids.begin(), value.member_ids.end()).size() != value.member_ids.size())
        throw std::invalid_argument("Duplicate parameter posterior identity.");
    for (const auto id : value.member_ids)
        if (std::none_of(members.begin(), members.end(), [id](const auto * atom) { return atom->GetSerialID() == id; }))
            throw std::invalid_argument("Parameter posterior identity outside group.");
    for (const auto * atom : members) EnsureAtomLocalPotential(m_model_object, *atom).ClearGroupMemberResult();
    if (value.inference)
        for (std::size_t i = 0; i < value.member_ids.size(); ++i)
            EnsureAtomLocalPotential(m_model_object, *m_model_object.FindAtomPtr(value.member_ids[i]))
                .SetGroupMemberResult(value.inference->member_results[i]);
    groups.SetParameterSummary(key, std::move(value));
}

void ModelAnalysisEditor::UpdateJointMetadata(JointAnalysisMetadata metadata)
{
    auto & result = ModelAnalysisData::Of(m_model_object).joint_result;
    if (!result) throw std::runtime_error("Joint diagnostic snapshot unavailable.");
    result->metadata = std::move(metadata);
}

void ModelAnalysisEditor::ClearJointResult()
{
    ModelAnalysisData::Of(m_model_object).joint_result.reset();
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
        }
    }
}

void ModelAnalysisEditor::InitializeFromSelection()
{
    Clear();
    RebuildAtomGroupsFromSelection();
    for (const auto stage : {
            FittingStage::First,
            FittingStage::Second })
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
        SetAtomLocalGaussianResult(FittingStage::Second, *atom, std::move(result));
        SetAtomStageEstimate(FittingStage::First, *atom, LocalStageEstimate{});
        SetAtomStageEstimate(FittingStage::Second, *atom, LocalStageEstimate{});
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
    auto & entry = EnsureAtomLocalPotential(m_model_object, atom_object);
    entry.ClearPeeling();
    entry.SetRawSamplingEntries(std::move(value));
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
    auto & entry = EnsureAtomLocalPotential(m_model_object, atom_object);
    const auto before = entry.StageEstimate(stage);
    LocalPotentialEntry prepared;
    prepared.SetGaussianResult(stage, result);
    auto estimate = prepared.StageEstimate(stage);
    estimate.source.atom_id = std::to_string(atom_object.GetSerialID());
    if (stage == FittingStage::Second)
        InvalidateStageChange(ModelAnalysisData::Of(m_model_object), atom_object.GetSerialID(), before, estimate);
    entry.SetGaussianResult(stage, std::move(result));
    entry.SetStageEstimate(stage, std::move(estimate));
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
    for (const auto key : group_entry.CollectGroupKeys()) InvalidateGroup(analysis_data, key);
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

void ModelAnalysisEditor::InitializeGroupAlpha(double alpha_g)
{
    auto & group_entry{ ModelAnalysisData::Of(m_model_object).AtomGroupEntry() };
    for (const auto group_key : group_entry.CollectGroupKeys())
    {
        SetAtomGroupAlphaG(group_key, alpha_g);
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
    auto & data = ModelAnalysisData::Of(m_model_object);
    if (data.AtomGroupEntry().HasGroup(group_key) && data.AtomGroupEntry().GetAlphaG(group_key) != alpha_g)
        InvalidateGroup(data, group_key);
    data.AtomGroupEntry().SetAlphaG(group_key, alpha_g);
}

} // namespace rhbm_gem
