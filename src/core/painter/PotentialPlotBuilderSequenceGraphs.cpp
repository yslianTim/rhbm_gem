#include <rhbm_gem/core/painter/PotentialPlotBuilder.hpp>

#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondClassifier.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TGraphErrors.h>
#endif

namespace rhbm_gem {

#ifdef HAVE_ROOT

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateNormalizedAtomGausEstimateScatterGraph(
    Element element, double reference_amplitude, bool reverse)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto model_object{ m_query.GetModelObject() };
    auto graph{ ROOTHelper::CreateGraphErrors() };
    std::unordered_map<int, double> amplitude_diff_to_carbonyl_oxygen_map;
    for (auto & atom : model_object->GetSelectedAtomList())
    {
        if (atom->GetSpot() != Spot::O) continue;
        if (atom->GetSpecialAtomFlag() == false)
        {
            auto entry{ atom->GetLocalPotentialEntry() };
            auto sequence_id{ atom->GetSequenceID() };
            auto amplitude_estimate{ entry->GetAmplitudeEstimateMDPDE() };
            amplitude_diff_to_carbonyl_oxygen_map[sequence_id] = amplitude_estimate - reference_amplitude;
        }
    }
    auto count{ 0 };
    for (auto & atom : model_object->GetSelectedAtomList())
    {
        if (atom->GetElement() != element) continue;
        auto sequence_id{ atom->GetSequenceID() };
        auto entry{ atom->GetLocalPotentialEntry() };
        auto normalized_amplitude{ entry->GetAmplitudeEstimateMDPDE() };
        if (amplitude_diff_to_carbonyl_oxygen_map.find(sequence_id) != amplitude_diff_to_carbonyl_oxygen_map.end())
        {
            normalized_amplitude -= amplitude_diff_to_carbonyl_oxygen_map.at(sequence_id);
        }
        if (reverse == false)
        {
            graph->SetPoint(count, normalized_amplitude, entry->GetWidthEstimateMDPDE());
        }
        else
        {
            graph->SetPoint(count, entry->GetWidthEstimateMDPDE(), normalized_amplitude);
        }
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateNormalizedBondGausEstimateScatterGraph(
    Element element, double reference_amplitude, bool reverse)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto model_object{ m_query.GetModelObject() };
    auto graph{ ROOTHelper::CreateGraphErrors() };
    std::unordered_map<int, double> amplitude_diff_to_carbonyl_oxygen_map;
    for (auto & bond : model_object->GetSelectedBondList())
    {
        if (bond->GetSpecialBondFlag() == false)
        {
            auto entry{ bond->GetLocalPotentialEntry() };
            auto sequence_id{ bond->GetAtomObject1()->GetSequenceID() };
            auto amplitude_estimate{ entry->GetAmplitudeEstimateMDPDE() };
            amplitude_diff_to_carbonyl_oxygen_map[sequence_id] = amplitude_estimate - reference_amplitude;
        }
    }
    auto count{ 0 };
    for (auto & bond : model_object->GetSelectedBondList())
    {
        if (bond->GetAtomObject1()->GetElement() != element) continue;
        auto sequence_id{ bond->GetAtomObject1()->GetSequenceID() };
        auto entry{ bond->GetLocalPotentialEntry() };
        auto normalized_amplitude{ entry->GetAmplitudeEstimateMDPDE() };
        if (amplitude_diff_to_carbonyl_oxygen_map.find(sequence_id) != amplitude_diff_to_carbonyl_oxygen_map.end())
        {
            normalized_amplitude -= amplitude_diff_to_carbonyl_oxygen_map.at(sequence_id);
        }
        if (reverse == false)
        {
            graph->SetPoint(count, normalized_amplitude, entry->GetWidthEstimateMDPDE());
        }
        else
        {
            graph->SetPoint(count, entry->GetWidthEstimateMDPDE(), normalized_amplitude);
        }
        count++;
    }
    return graph;
}

std::unordered_map<std::string, std::unique_ptr<TGraphErrors>>
PotentialPlotBuilder::CreateAtomMapValueToSequenceIDGraphMap(
    size_t main_chain_element_id, Residue residue)
{
    if (IsModelObjectAvailable() == false)
    {
        return {};
    }
    auto model_object{ m_query.GetModelObject() };

    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> graph_map;
    std::unordered_map<std::string, int> count_map;

    for (auto & atom : model_object->GetSelectedAtomList())
    {
        if (atom->GetElement() != AtomClassifier::GetMainChainElement(main_chain_element_id)) continue;
        if (atom->GetSpot() != AtomClassifier::GetMainChainSpot(main_chain_element_id)) continue;
        if (residue != Residue::UNK && atom->GetResidue() != residue) continue;
        auto entry{ atom->GetLocalPotentialEntry() };
        auto sequence_id{ atom->GetSequenceID() };
        auto chain_id{ atom->GetChainID() };
        if (sequence_id < 0) continue;
        if (graph_map.find(chain_id) == graph_map.end())
        {
            graph_map[chain_id] = ROOTHelper::CreateGraphErrors();
            count_map[chain_id] = 0;
        }
        auto x_value{ static_cast<double>(sequence_id) };
        graph_map[chain_id]->SetPoint(count_map[chain_id], x_value, entry->GetMapValueNearCenter());
        count_map[chain_id]++;
    }
    return graph_map;
}

std::unordered_map<std::string, std::unique_ptr<TGraphErrors>>
PotentialPlotBuilder::CreateAtomQScoreToSequenceIDGraphMap(
    size_t main_chain_element_id, const int par_choice)
{
    if (IsModelObjectAvailable() == false)
    {
        return {};
    }
    auto model_object{ m_query.GetModelObject() };

    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> graph_map;
    std::unordered_map<std::string, int> count_map;

    for (auto & atom : model_object->GetSelectedAtomList())
    {
        if (atom->GetElement() != AtomClassifier::GetMainChainElement(main_chain_element_id)) continue;
        if (atom->GetSpot() != AtomClassifier::GetMainChainSpot(main_chain_element_id)) continue;
        auto entry{ atom->GetLocalPotentialEntry() };
        auto sequence_id{ atom->GetSequenceID() };
        auto chain_id{ atom->GetChainID() };
        if (sequence_id < 0) continue;
        if (graph_map.find(chain_id) == graph_map.end())
        {
            graph_map[chain_id] = ROOTHelper::CreateGraphErrors();
            count_map[chain_id] = 0;
        }
        auto x_value{ static_cast<double>(sequence_id) };
        auto q_score{ entry->CalculateQScore(par_choice) };
        // tmp
        if (sequence_id == 20 || sequence_id == 40 || sequence_id == 60 || sequence_id == 80 || sequence_id == 100)
        {
            if (main_chain_element_id == 0)
            {   
                Logger::Log(LogLevel::Info, Form("Sequence ID %d, Q-Score: %.4f, %d", sequence_id, q_score, par_choice));
            }
        }
        // end tmp
        graph_map[chain_id]->SetPoint(count_map[chain_id], x_value, q_score);
        count_map[chain_id]++;
    }
    return graph_map;
}

std::unordered_map<std::string, std::unique_ptr<TGraphErrors>>
PotentialPlotBuilder::CreateAtomGausEstimateToSequenceIDGraphMap(
    size_t main_chain_element_id, const int par_id, Residue residue)
{
    if (IsModelObjectAvailable() == false)
    {
        return {};
    }
    auto model_object{ m_query.GetModelObject() };

    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> graph_map;
    std::unordered_map<std::string, int> count_map;

    for (auto & atom : model_object->GetSelectedAtomList())
    {
        if (atom->GetElement() != AtomClassifier::GetMainChainElement(main_chain_element_id)) continue;
        if (atom->GetSpot() != AtomClassifier::GetMainChainSpot(main_chain_element_id)) continue;
        if (residue != Residue::UNK && atom->GetResidue() != residue) continue;
        auto entry{ atom->GetLocalPotentialEntry() };
        auto sequence_id{ atom->GetSequenceID() };
        auto chain_id{ atom->GetChainID() };
        if (sequence_id < 0) continue;
        if (graph_map.find(chain_id) == graph_map.end())
        {
            graph_map[chain_id] = ROOTHelper::CreateGraphErrors();
            count_map[chain_id] = 0;
        }
        auto x_value{ static_cast<double>(sequence_id) };
        graph_map[chain_id]->SetPoint(count_map[chain_id], x_value, entry->GetGausEstimateMDPDE(par_id));
        count_map[chain_id]++;
    }
    return graph_map;
}

std::unordered_map<std::string, std::unique_ptr<TGraphErrors>>
PotentialPlotBuilder::CreateBondGausEstimateToSequenceIDGraphMap(
    size_t main_chain_element_id, const int par_id, Residue residue)
{
    if (IsModelObjectAvailable() == false)
    {
        return {};
    }

    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> graph_map;
    std::unordered_map<std::string, int> count_map;
    auto class_key{ ChemicalDataHelper::GetSimpleBondClassKey() };
    auto group_key{ BondClassifier::GetMainChainSimpleBondClassGroupKey(main_chain_element_id) };
    for (auto & bond : m_query.GetBondObjectList(group_key, class_key))
    {
        if (residue != Residue::UNK && bond->GetAtomObject1()->GetResidue() != residue) continue;
        auto entry{ bond->GetLocalPotentialEntry() };
        auto sequence_id{ bond->GetAtomObject1()->GetSequenceID() };
        auto chain_id{ bond->GetAtomObject1()->GetChainID() };
        if (sequence_id < 0) continue;
        if (graph_map.find(chain_id) == graph_map.end())
        {
            graph_map[chain_id] = ROOTHelper::CreateGraphErrors();
            count_map[chain_id] = 0;
        }
        auto x_value{ static_cast<double>(sequence_id) };
        auto y_value{ entry->GetGausEstimateMDPDE(par_id) };
        graph_map[chain_id]->SetPoint(count_map[chain_id], x_value, y_value);
        count_map[chain_id]++;
    }

    for (auto & [chain, graph] : graph_map)
    {
        graph->Sort(); // sort points along x axis
    }
    return graph_map;
}

std::unordered_map<std::string, std::unique_ptr<TGraphErrors>>
PotentialPlotBuilder::CreateAtomGausEstimatePosteriorToSequenceIDGraphMap(
    size_t main_chain_element_id, const std::string & class_key, const int par_id, Residue residue)
{
    if (IsModelObjectAvailable() == false)
    {
        return {};
    }
    auto model_object{ m_query.GetModelObject() };

    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> graph_map;
    std::unordered_map<std::string, int> count_map;

    for (auto & atom : model_object->GetSelectedAtomList())
    {
        if (atom->GetElement() != AtomClassifier::GetMainChainElement(main_chain_element_id)) continue;
        if (atom->GetSpot() != AtomClassifier::GetMainChainSpot(main_chain_element_id)) continue;
        if (residue != Residue::UNK && atom->GetResidue() != residue) continue;
        auto entry{ atom->GetLocalPotentialEntry() };
        auto sequence_id{ atom->GetSequenceID() };
        auto chain_id{ atom->GetChainID() };
        if (sequence_id < 0) continue;
        if (graph_map.find(chain_id) == graph_map.end())
        {
            graph_map[chain_id] = ROOTHelper::CreateGraphErrors();
            count_map[chain_id] = 0;
        }
        auto x_value{ static_cast<double>(sequence_id) };
        graph_map[chain_id]->SetPoint(
            count_map[chain_id], x_value, entry->GetGausEstimatePosterior(class_key, par_id)
        );
        graph_map[chain_id]->SetPointError(
            count_map[chain_id], 0.0, entry->GetGausVariancePosterior(class_key, par_id)
        );
        count_map[chain_id]++;
    }
    return graph_map;
}

#endif

} // namespace rhbm_gem
