#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <stdexcept>

namespace rhbm_gem {

ModelAnalysisView::ModelAnalysisView(const ModelObject & model_object) :
    m_model_object(model_object)
{
}

bool ModelAnalysisView::HasGroupedAnalysisData(FittingStage stage) const
{
    return !ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().CollectGroupKeys(stage).empty();
}

bool ModelAnalysisView::HasAtomGroup(
    FittingStage stage,
    GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().HasGroup(stage, group_key);
}

const GaussianModel3D & ModelAnalysisView::GetAtomGroupMean(
    FittingStage stage,
    GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().GetMean(stage, group_key);
}

const GaussianModel3D & ModelAnalysisView::GetAtomGroupMDPDE(
    FittingStage stage,
    GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().GetMDPDE(stage, group_key);
}

const GaussianModel3D & ModelAnalysisView::GetAtomGroupPrior(
    FittingStage stage,
    GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().GetPrior(stage, group_key);
}

GaussianModel3DWithUncertainty ModelAnalysisView::GetAtomGroupPriorWithUncertainty(
    FittingStage stage,
    GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().GetPriorWithUncertainty(stage, group_key);
}

const std::vector<AtomObject *> & ModelAnalysisView::GetAtomObjectList(
    FittingStage stage,
    GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().GetMembers(stage, group_key);
}

double ModelAnalysisView::GetAtomAlphaR(
    FittingStage stage,
    GroupKey group_key) const
{
    const auto & atom_list{ GetAtomObjectList(stage, group_key) };
    if (atom_list.empty())
    {
        throw std::runtime_error("Atom group has no members.");
    }
    return AtomLocalPotentialView::For(*atom_list.front()).GetAlphaR(stage);
}

double ModelAnalysisView::GetAtomAlphaG(
    FittingStage stage,
    GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().GetAlphaG(stage, group_key);
}

std::vector<GroupKey> ModelAnalysisView::CollectAtomGroupKeys(
    FittingStage stage) const
{
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().CollectGroupKeys(stage);
}

} // namespace rhbm_gem
