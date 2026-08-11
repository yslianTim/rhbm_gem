#include <rhbm_gem/data/object/ModelAnalysisView.hpp>

#include "data/detail/AtomClassifier.hpp"
#include "data/detail/ModelAnalysisData.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include <algorithm>
#include <array>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace rhbm_gem {

namespace {

constexpr std::array<Spot, 5> kGroupPriorSummarySpotList{
    Spot::C, Spot::CA, Spot::CB, Spot::N, Spot::O
};

struct GaussianModelParameterSamples
{
    std::vector<double> amplitude_list{};
    std::vector<double> width_list{};
    std::vector<double> offset_list{};
};

} // namespace

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

std::optional<GaussianModel3DWithUncertainty>
ModelAnalysisView::FindAtomGroupPriorWithUncertainty(
    FittingStage stage,
    const AtomObject & atom_object) const
{
    const auto group_key{ data_internal::GetGroupKey(&atom_object) };
    if (!HasAtomGroup(stage, group_key)) return std::nullopt;
    return GetAtomGroupPriorWithUncertainty(stage, group_key);
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
    return AtomLocalPotentialView::RequireFor(*atom_list.front()).GetAlphaR(stage);
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

std::string ModelAnalysisView::GetAtomGroupingSummary(
    FittingStage stage) const
{
    std::string description{ "Atomic model includes " };
    description += std::to_string(CollectAtomGroupKeys(stage).size()) + " atom groups.";
    return description;
}

std::string ModelAnalysisView::GetGroupPriorSpotSummary(FittingStage stage) const
{
    std::map<Spot, GaussianModelParameterSamples> spot_sample_map;
    for (const auto group_key : CollectAtomGroupKeys(stage))
    {
        const auto & atom_list{ GetAtomObjectList(stage, group_key) };
        if (atom_list.empty()) continue;

        const auto spot{ atom_list.front()->GetSpot() };
        if (std::find(
                kGroupPriorSummarySpotList.begin(), kGroupPriorSummarySpotList.end(),
                spot) == kGroupPriorSummarySpotList.end()) continue;
        const auto & prior{ GetAtomGroupPrior(stage, group_key) };
        auto & sample_list{ spot_sample_map[spot] };
        sample_list.amplitude_list.emplace_back(prior.GetAmplitude());
        sample_list.width_list.emplace_back(prior.GetWidth());
        sample_list.offset_list.emplace_back(prior.GetOffset());
    }

    if (spot_sample_map.empty())
    {
        return "Group fitting prior summary by Spot: no atom groups available.";
    }

    std::ostringstream summary;
    summary << "Group fitting prior summary by Spot:\n"
        << "|---Spot---|------Amplitude------|--------Width--------|-------Offset--------|\n"
        << "|          |   mean   |   s.d.   |   mean   |   s.d.   |   mean   |   s.d.   |";
    for (const auto spot : kGroupPriorSummarySpotList)
    {
        const auto sample_iter{ spot_sample_map.find(spot) };
        if (sample_iter == spot_sample_map.end()) continue;
        const auto & sample_list{ sample_iter->second };

        const auto amplitude_mean{
            array_helper::ComputeMean(
                sample_list.amplitude_list.data(), sample_list.amplitude_list.size())
        };
        const auto width_mean{
            array_helper::ComputeMean(
                sample_list.width_list.data(), sample_list.width_list.size())
        };
        const auto offset_mean{
            array_helper::ComputeMean(
                sample_list.offset_list.data(), sample_list.offset_list.size())
        };

        summary << "\n| " << std::left << std::setw(8)
            << ChemicalDataHelper::GetLabel(spot)
            << " | " << std::right << std::fixed << std::setprecision(2)
            << std::setw(8) << amplitude_mean
            << " | " << std::setw(8)
            << array_helper::ComputeStandardDeviation(
                sample_list.amplitude_list.data(),
                sample_list.amplitude_list.size(),
                amplitude_mean)
            << " | " << std::setw(8) << width_mean
            << " | " << std::setw(8)
            << array_helper::ComputeStandardDeviation(
                sample_list.width_list.data(),
                sample_list.width_list.size(),
                width_mean)
            << " | " << std::setw(8) << offset_mean
            << " | " << std::setw(8)
            << array_helper::ComputeStandardDeviation(
                sample_list.offset_list.data(),
                sample_list.offset_list.size(),
                offset_mean)
            << " |";
    }
    return summary.str();
}

std::string ModelAnalysisView::GetLocalFittingResultCsv(bool peeling_applied) const
{
    auto atom_list{ m_model_object.GetSelectedAtoms() };
    std::sort(
        atom_list.begin(),
        atom_list.end(),
        [](const AtomObject * lhs, const AtomObject * rhs)
        {
            return lhs->GetSerialID() < rhs->GetSerialID();
        });

    std::ostringstream table;
    table << std::fixed << std::setprecision(2);
    table
        << "serial id,residue,spot,neighbor count,peeling ratio,"
        << "amplitude 1st,amplitude 2nd,amplitude 3rd,"
        << "width 1st,width 2nd,width 3rd,"
        << "offset 1st,offset 2nd,offset 3rd";
    for (const auto * atom : atom_list)
    {
        const auto local_view{ AtomLocalPotentialView::RequireFor(*atom) };
        const auto & first_model{ local_view.GetEstimateMDPDE(FittingStage::First) };
        const auto & second_model{ local_view.GetEstimateMDPDE(FittingStage::Second) };
        const auto & third_model{ local_view.GetEstimateMDPDE(FittingStage::Third) };

        table << '\n'
            << atom->GetSerialID() << ','
            << ChemicalDataHelper::GetLabel(atom->GetResidue()) << ','
            << atom->GetAtomID() << ','
            << local_view.GetNeighborCountForPeeling() << ',';
        const auto peeling_ratio{ local_view.GetLocalFittingPeelingRatio(peeling_applied) };
        if (!peeling_ratio.has_value())
        {
            table << "nan";
        }
        else
        {
            table << *peeling_ratio;
        }
        table << ','
            << first_model.GetAmplitude() << ','
            << second_model.GetAmplitude() << ','
            << third_model.GetAmplitude() << ','
            << first_model.GetWidth() << ','
            << second_model.GetWidth() << ','
            << third_model.GetWidth() << ','
            << first_model.GetOffset() << ','
            << second_model.GetOffset() << ','
            << third_model.GetOffset();
    }
    table << '\n';
    return table.str();
}

} // namespace rhbm_gem
