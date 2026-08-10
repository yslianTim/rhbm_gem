#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>

#include <map>
#include <stdexcept>
#include <string>

namespace rhbm_gem {

namespace {

} // namespace

ModelAnalysisView::ModelAnalysisView(const ModelObject & model_object) :
    m_model_object(model_object)
{
}

bool ModelAnalysisView::HasGroupedAnalysisData() const
{
    return !ModelAnalysisData::Of(m_model_object).AtomGroupEntry().CollectGroupKeys().empty();
}

bool ModelAnalysisView::HasAtomGroup(GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object).AtomGroupEntry().HasGroup(group_key);
}

const GaussianModel3D & ModelAnalysisView::GetAtomGroupMean(GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object).AtomGroupEntry().GetMean(group_key);
}

const GaussianModel3D & ModelAnalysisView::GetAtomGroupMDPDE(GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object).AtomGroupEntry().GetMDPDE(group_key);
}

const GaussianModel3D & ModelAnalysisView::GetAtomGroupPrior(GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object).AtomGroupEntry().GetPrior(group_key);
}

GaussianModel3DWithUncertainty ModelAnalysisView::GetAtomGroupPriorWithUncertainty(GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object).AtomGroupEntry().GetPriorWithUncertainty(group_key);
}

const std::vector<AtomObject *> & ModelAnalysisView::GetAtomObjectList(GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object).AtomGroupEntry().GetMembers(group_key);
}

double ModelAnalysisView::GetAtomAlphaR(
    LocalFittingStage stage,
    GroupKey group_key) const
{
    const auto & atom_list{ GetAtomObjectList(group_key) };
    if (atom_list.empty())
    {
        throw std::runtime_error("Atom group has no members.");
    }
    return AtomLocalPotentialView::RequireFor(*atom_list.front()).GetAlphaR(stage);
}

double ModelAnalysisView::GetAtomAlphaG(GroupKey group_key) const
{
    return ModelAnalysisData::Of(m_model_object).AtomGroupEntry().GetAlphaG(group_key);
}

std::vector<GroupKey> ModelAnalysisView::CollectAtomGroupKeys() const
{
    return ModelAnalysisData::Of(m_model_object).AtomGroupEntry().CollectGroupKeys();
}

std::string ModelAnalysisView::GetAtomCountingSummary() const
{
    std::map<Element, std::size_t> element_counts;
    for (const auto * atom : m_model_object.GetSelectedAtoms())
    {
        element_counts[atom->GetElement()]++;
    }

    std::string description{
        "Number of selected atom = " + std::to_string(m_model_object.GetSelectedAtomCount())
    };
    for (const auto & [element, count] : element_counts)
    {
        description +=
            "\n - Element type: " + ChemicalDataHelper::GetLabel(element) + " include "
            + std::to_string(count) + " atoms.";
    }
    return description;
}

std::string ModelAnalysisView::GetAtomGroupingSummary() const
{
    std::string description{ "Atomic model includes " };
    description += std::to_string(CollectAtomGroupKeys().size()) + " atom groups.";
    return description;
}

} // namespace rhbm_gem
