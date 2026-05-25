#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

#include "data/detail/GroupPotentialEntry.hpp"
#include "data/detail/AtomClassifier.hpp"
#include "data/detail/LocalPotentialEntry.hpp"
#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>

#include <stdexcept>
#include <string>

namespace rhbm_gem {

AtomLocalPotentialEditor::AtomLocalPotentialEditor(LocalPotentialEntry & entry) :
    m_entry{ &entry }
{
}

void AtomLocalPotentialEditor::SetSamplingEntries(LocalPotentialSampleList value)
{
    m_entry->SetSamplingEntries(std::move(value));
}

void AtomLocalPotentialEditor::SetGaussianResult(LocalGaussianResult value)
{
    m_entry->SetGaussianResult(std::move(value));
}

void AtomLocalPotentialEditor::SetAlphaR(double value)
{
    m_entry->SetAlphaR(value);
}

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
            entry->ClearTransientFitState();
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
    analysis_data.AtomGroupEntries().clear();
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        auto & group_entry{ analysis_data.EnsureAtomGroupEntry(class_key) };
        for (auto * atom : m_model_object.GetSelectedAtoms())
        {
            const auto group_key{ AtomClassifier::GetGroupKeyInClass(atom, class_key) };
            group_entry.AddMember(group_key, *atom);
        }
    }
}

void ModelAnalysisEditor::ApplyAtomGroupGaussianResult(
    GroupKey group_key,
    const std::string & class_key,
    const GroupGaussianResult & group_result)
{
    auto & analysis_data{ ModelAnalysisData::Of(m_model_object) };
    auto group_entry{ analysis_data.FindAtomGroupEntry(class_key) };
    if (group_entry == nullptr)
    {
        throw std::runtime_error("Atom group entry is not available.");
    }
    const auto & atom_list{ group_entry->GetMembers(group_key) };
    if (group_result.member_results.size() != atom_list.size())
    {
        throw std::invalid_argument("Atom group result member result count is inconsistent.");
    }

    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        auto & atom_entry{ analysis_data.EnsureAtomLocalEntry(*atom_list.at(i)) };
        atom_entry.SetAnnotation(class_key,
            LocalPotentialAnnotation{
                group_result.member_results.at(i).mdpde,
                group_result.member_results.at(i).is_outlier,
                group_result.member_results.at(i).statistical_distance
            }
        );
    }
    analysis_data.EnsureAtomGroupEntry(class_key).SetGaussianResult(group_key, group_result);
}

void ModelAnalysisEditor::SetAtomGroupAlphaG(
    GroupKey group_key,
    const std::string & class_key,
    double alpha_g)
{
    ModelAnalysisData::Of(m_model_object).EnsureAtomGroupEntry(class_key).SetAlphaG(group_key, alpha_g);
}

} // namespace rhbm_gem
