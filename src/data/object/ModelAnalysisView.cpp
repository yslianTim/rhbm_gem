#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/ModelObject.hpp>

#include <stdexcept>

namespace rhbm_gem {

namespace {

template <typename EntryT>
std::vector<GroupKey> CollectGroupKeys(const EntryT * entry)
{
    return entry == nullptr ? std::vector<GroupKey>{} : entry->CollectGroupKeys();
}

template <typename EntryT>
const EntryT * RequireGroupEntry(const EntryT * entry, const char * context)
{
    if (entry == nullptr)
    {
        throw std::runtime_error(std::string(context) + " is not available.");
    }
    return entry;
}

} // namespace

ModelAnalysisView::ModelAnalysisView(const ModelObject & model_object) :
    m_model_object(model_object)
{
}

bool ModelAnalysisView::HasGroupedAnalysisData() const
{
    const auto & analysis_data{ ModelAnalysisData::Of(m_model_object) };
    return !analysis_data.AtomGroupEntries().empty();
}

bool ModelAnalysisView::HasAtomGroup(GroupKey group_key, const std::string & class_key) const
{
    const auto * entry{ ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key) };
    if (entry == nullptr)
    {
        return false;
    }
    return entry->HasGroup(group_key);
}

const GaussianModel3D & ModelAnalysisView::GetAtomGroupMean(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key),
            "Atom group entry")
    };
    return entry->GetMean(group_key);
}

const GaussianModel3D & ModelAnalysisView::GetAtomGroupMDPDE(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key),
            "Atom group entry")
    };
    return entry->GetMDPDE(group_key);
}

const GaussianModel3D & ModelAnalysisView::GetAtomGroupPrior(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key),
            "Atom group entry")
    };
    return entry->GetPrior(group_key);
}

GaussianModel3DWithUncertainty ModelAnalysisView::GetAtomGroupPriorWithUncertainty(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key),
            "Atom group entry")
    };
    return entry->GetPriorWithUncertainty(group_key);
}

const std::vector<AtomObject *> & ModelAnalysisView::GetAtomObjectList(
    GroupKey group_key,
    const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key),
            "Atom group entry")
    };
    return entry->GetMembers(group_key);
}

double ModelAnalysisView::GetAtomAlphaR(GroupKey group_key, const std::string & class_key) const
{
    const auto & atom_list{ GetAtomObjectList(group_key, class_key) };
    if (atom_list.empty())
    {
        throw std::runtime_error("Atom group has no members.");
    }
    return AtomLocalPotentialView::RequireFor(*atom_list.front()).GetAlphaR();
}

double ModelAnalysisView::GetAtomAlphaG(GroupKey group_key, const std::string & class_key) const
{
    const auto * entry{
        RequireGroupEntry(
            ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key),
            "Atom group entry")
    };
    return entry->GetAlphaG(group_key);
}

std::vector<GroupKey> ModelAnalysisView::CollectAtomGroupKeys(const std::string & class_key) const
{
    return CollectGroupKeys(ModelAnalysisData::Of(m_model_object).FindAtomGroupEntry(class_key));
}

} // namespace rhbm_gem
