#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>

namespace rhbm_gem {

ModelAnalysisView::ModelAnalysisView(const ModelObject & model_object) :
    m_model_object(model_object)
{
}

const std::optional<JointAnalysisResult> & ModelAnalysisView::GetJointResult() const
{
    return ModelAnalysisData::Of(m_model_object).joint_result;
}

bool ModelAnalysisView::HasGroupedAnalysisData() const
{
    return !ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().CollectGroupKeys().empty();
}

bool ModelAnalysisView::HasAtomGroup(GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().HasGroup(group_key);
}

const std::optional<GroupParameterSummary> & ModelAnalysisView::GetGroupParameterSummary(GroupKey key) const
{
    return ModelAnalysisData::Of(m_model_object).AtomGroupEntry().GetParameterSummary(key);
}

bool ModelAnalysisView::HasAtomGroupPrior(GroupKey key) const
{
    if (!HasAtomGroup(key)) return false;
    const auto & summary = GetGroupParameterSummary(key);
    return summary ? summary->inference.has_value() :
        ModelAnalysisData::Of(m_model_object).AtomGroupEntry().GetPrior(key).GetWidth() > 0;
}

const GaussianModel3D & ModelAnalysisView::GetAtomGroupMean(GroupKey group_key) const
{
    const auto & summary=GetGroupParameterSummary(group_key);
    if (summary)
    {
        if (!summary->descriptive_mean) throw std::runtime_error("Group descriptive mean unavailable.");
        return *summary->descriptive_mean;
    }
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().GetMean(group_key);
}

const GaussianModel3D & ModelAnalysisView::GetAtomGroupMDPDE(GroupKey group_key) const
{
    const auto & summary=GetGroupParameterSummary(group_key);
    if (summary && !summary->inference) throw std::runtime_error("Group MDPDE unavailable.");
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().GetMDPDE(group_key);
}

const GaussianModel3D & ModelAnalysisView::GetAtomGroupPrior(GroupKey group_key) const
{
    const auto & summary=GetGroupParameterSummary(group_key);
    if (summary && !summary->inference) throw std::runtime_error("Group prior unavailable.");
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().GetPrior(group_key);
}

GaussianModel3DWithUncertainty ModelAnalysisView::GetAtomGroupPriorWithUncertainty(GroupKey group_key) const
{
    const auto & summary=GetGroupParameterSummary(group_key);
    if (summary && !summary->inference) throw std::runtime_error("Group prior uncertainty unavailable.");
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().GetPriorWithUncertainty(group_key);
}

const std::vector<AtomObject *> & ModelAnalysisView::GetAtomObjectList(GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().GetMembers(group_key);
}

double ModelAnalysisView::GetAtomAlphaG(GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().GetAlphaG(group_key);
}

std::vector<GroupKey> ModelAnalysisView::CollectAtomGroupKeys() const
{
    return ModelAnalysisData::Of(m_model_object)
        .AtomGroupEntry().CollectGroupKeys();
}

} // namespace rhbm_gem
