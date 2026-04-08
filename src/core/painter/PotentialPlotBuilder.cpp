#include "PotentialPlotBuilder.hpp"

#include "data/detail/ModelAnalysisData.hpp"
#include "data/detail/ModelDerivedState.hpp"
#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include "data/detail/BondClassifier.hpp"
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <detail/PotentialSeriesOps.hpp>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#endif

#include <array>
#include <cmath>
#include <stdexcept>

namespace rhbm_gem {

PotentialPlotBuilder::PotentialPlotBuilder(ModelObject * model_object) :
    m_model_object{ model_object }
{
}

PotentialPlotBuilder::PotentialPlotBuilder(AtomObject * atom_object) :
    m_atom_object{ atom_object }
{
}

PotentialPlotBuilder::PotentialPlotBuilder(BondObject * bond_object) :
    m_bond_object{ bond_object }
{
}

ModelPotentialView PotentialPlotBuilder::GetModelView() const
{
    if (m_model_object == nullptr)
    {
        throw std::runtime_error("Model object is not available.");
    }
    return ModelPotentialView(*m_model_object);
}

const LocalPotentialEntry & PotentialPlotBuilder::GetLocalEntry() const
{
    if (m_atom_object != nullptr)
    {
        return ModelAnalysisData::RequireLocalEntry(*m_atom_object);
    }
    if (m_bond_object != nullptr)
    {
        return ModelAnalysisData::RequireLocalEntry(*m_bond_object);
    }
    throw std::runtime_error("Local entry is not available.");
}

bool PotentialPlotBuilder::IsModelObjectAvailable() const
{
    if (m_model_object == nullptr)
    {
        Logger::Log(LogLevel::Error, "Model object is not available.");
        return false;
    }
    return true;
}

bool PotentialPlotBuilder::IsAtomLocalEntryAvailable() const
{
    return m_atom_object != nullptr && ModelAnalysisData::FindLocalEntry(*m_atom_object) != nullptr;
}

size_t PotentialPlotBuilder::GetAtomResidueCount(
    const std::string & class_key, Residue residue, Structure structure) const
{
    GroupKey group_key{ 0 };
        if (class_key == ChemicalDataHelper::GetComponentAtomClassKey())
    {
        group_key = AtomClassifier::GetMainChainComponentAtomClassGroupKey(0, residue);
    }
    else if (class_key == ChemicalDataHelper::GetStructureAtomClassKey())
    {
        group_key = AtomClassifier::GetMainChainStructureAtomClassGroupKey(0, structure, residue);
    }
    if (IsAvailableAtomGroupKey(group_key, class_key) == false)
    {
        return 0;
    }
    return GetModelView().GetAtomObjectList(group_key, class_key).size();
}

size_t PotentialPlotBuilder::GetBondResidueCount(
    const std::string & class_key, Residue residue) const
{
    GroupKey group_key{ 0 };
        if (class_key == ChemicalDataHelper::GetComponentBondClassKey())
    {
        group_key = BondClassifier::GetMainChainComponentBondClassGroupKey(0, residue);
    }
    if (IsAvailableBondGroupKey(group_key, class_key) == false)
    {
        return 0;
    }
    return GetModelView().GetBondObjectList(group_key, class_key).size();
}

bool PotentialPlotBuilder::IsAvailableAtomGroupKey(
    GroupKey group_key, const std::string & class_key, bool varbose) const
{
    return GetModelView().HasAtomGroup(group_key, class_key, varbose);
}

bool PotentialPlotBuilder::IsAvailableBondGroupKey(
    GroupKey group_key, const std::string & class_key, bool varbose) const
{
    return GetModelView().HasBondGroup(group_key, class_key, varbose);
}

#ifdef HAVE_ROOT
std::unique_ptr<TH1D> PotentialPlotBuilder::CreateComponentCountHistogram(
    std::vector<GroupKey> & group_key_list, const std::string & class_key) const
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto component_size{ static_cast<int>(group_key_list.size()) };
    auto hist{
        ROOTHelper::CreateHist1D(
            "hist_" + class_key, "Component Count", component_size, -0.5, component_size - 0.5)
    };
    for (size_t i = 0; i < group_key_list.size(); i++)
    {
        auto count{ GetModelView().GetAtomObjectList(group_key_list.at(i), class_key).size() };
        hist->SetBinContent(static_cast<int>(i + 1), static_cast<double>(count));
    }
    return hist;
}

std::unique_ptr<TH2D> PotentialPlotBuilder::CreateAtomResidueCountHistogram2D(
    const std::string & class_key)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto hist{ ROOTHelper::CreateHist2D("hist", "Residue Count", 20, -0.5, 19.5, 3, -0.5, 2.5) };
    std::vector<Structure> structure_list{ Structure::FREE, Structure::HELX_P, Structure::SHEET };
    for (size_t i = 0; i < structure_list.size(); i++)
    {
        auto structure{ structure_list.at(i) };
        for (auto & residue : ChemicalDataHelper::GetStandardAminoAcidList())
        {
            auto count{ GetAtomResidueCount(class_key, residue, structure) };
            hist->SetBinContent(
                static_cast<int>(residue), static_cast<int>(i + 1), static_cast<double>(count));
        }
    }
    return hist;
}

std::unique_ptr<TH1D> PotentialPlotBuilder::CreateAtomResidueCountHistogram(
    const std::string & class_key, Structure structure)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto hist{ ROOTHelper::CreateHist1D("hist", "Residue Count", 20, -0.5, 19.5) };
    for (auto & residue : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        auto count{ GetAtomResidueCount(class_key, residue, structure) };
        hist->SetBinContent(static_cast<int>(residue), static_cast<double>(count));
    }
    return hist;
}

std::unique_ptr<TH1D> PotentialPlotBuilder::CreateBondResidueCountHistogram(
    const std::string & class_key)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto hist{ ROOTHelper::CreateHist1D("hist", "Residue Count", 20, -0.5, 19.5) };
    for (auto & residue : ChemicalDataHelper::GetStandardAminoAcidList())
    {
        auto count{ GetBondResidueCount(class_key, residue) };
        hist->SetBinContent(static_cast<int>(residue), static_cast<double>(count));
    }
    return hist;
}

std::unique_ptr<TH1D> PotentialPlotBuilder::CreateAtomGausEstimateHistogram(
    GroupKey group_key, const std::string & class_key, int par_id) const
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    if (IsAvailableAtomGroupKey(group_key, class_key) == false)
    {
        Logger::Log(LogLevel::Error, "Group key is not available.");
        return nullptr;
    }

    const auto & atom_list{ GetModelView().GetAtomObjectList(group_key, class_key) };
    std::vector<double> gaus_estimate_list;
    gaus_estimate_list.reserve(atom_list.size());
    for (auto atom : atom_list)
    {
        auto * local_entry{ ModelAnalysisData::FindLocalEntry(*atom) };
        gaus_estimate_list.emplace_back(local_entry->GetEstimateMDPDE().GetParameter(par_id));
    }

    double x_min{ 0.0 };
    double x_max{ 1.0 };

    auto estimate_range{
        ArrayStats<double>::ComputeScalingPercentileRangeTuple(gaus_estimate_list, 0.1, 0.05, 0.95)
    };

    if (gaus_estimate_list.size() > 1)
    {
        x_min = std::get<0>(estimate_range);
        x_max = std::get<1>(estimate_range);
    }
    else
    {
        x_max = 2.0 * std::ceil(gaus_estimate_list.at(0));
    }
    if (x_max == 0.0)
    {
        x_max = 1.0;
    }

    auto estimate_average{ 0.5 * (x_max + x_min) };
    if (par_id == 0 && (x_max - x_min) < 5.0)
    {
        x_max = estimate_average + 2.5;
        x_min = estimate_average - 2.5;
    }
    else if (par_id == 1 && (x_max - x_min) < 0.1)
    {
        x_max = estimate_average + 0.05;
        x_min = estimate_average - 0.05;
    }

    auto hist_name{ std::to_string(group_key) + "_" + class_key + "_par" + std::to_string(par_id) };
    auto hist{ ROOTHelper::CreateHist1D(hist_name.data(), "", 25, x_min, x_max) };
    for (auto & value : gaus_estimate_list)
    {
        hist->Fill(value);
    }
    return hist;
}

std::unique_ptr<TH1D> PotentialPlotBuilder::CreateLinearModelDataHistogram(int dimension_id) const
{
    auto data_array{ series_ops::BuildLinearModelDistanceAndMapValueList(GetLocalEntry()) };
    std::vector<float> data_list;
    data_list.reserve(data_array.size());
    for (auto & [distance, map_value] : data_array)
    {
        switch (dimension_id)
        {
            case 0:
                data_list.emplace_back(distance);
                break;
            case 1:
                data_list.emplace_back(map_value);
                break;
            default:
                throw std::runtime_error("Dimension id is invalid.");
        }
    }

    auto data_range{ ArrayStats<float>::ComputeScalingRangeTuple(data_list, 0.1f) };
    double x_min{ std::get<0>(data_range) };
    double x_max{ std::get<1>(data_range) };

    auto hist_name{ "data_dim" + std::to_string(dimension_id) };
    auto hist{ ROOTHelper::CreateHist1D(hist_name.data(), "", 25, x_min, x_max) };
    for (auto & value : data_list)
    {
        hist->Fill(value);
    }
    return hist;
}

std::unique_ptr<TH2D> PotentialPlotBuilder::CreateDistanceToMapValueHistogram(
    int x_bin_size, int y_bin_size) const
{
    auto distance_range{ series_ops::ComputeDistanceRange(GetLocalEntry(), 0.0) };
    auto map_value_range{ series_ops::ComputeMapValueRange(GetLocalEntry(), 0.1) };
    auto hist{
        ROOTHelper::CreateHist2D(
            "hist_distance_mapvalue", "Distance vs Map Value",
            x_bin_size, std::get<0>(distance_range), std::get<1>(distance_range),
            y_bin_size, std::get<0>(map_value_range), std::get<1>(map_value_range))
    };
    for (auto & [distance, map_value] : GetLocalEntry().GetDistanceAndMapValueList())
    {
        hist->Fill(distance, map_value);
    }
    return hist;
}

std::vector<std::unique_ptr<TH1D>> PotentialPlotBuilder::CreateMainChainAtomGausRankHistogram(
    int par_id, int & chain_size, Residue residue,
    size_t extra_id, std::vector<Residue> veto_residues_list)
{
    if (IsModelObjectAvailable() == false)
    {
        return {};
    }

    auto model_object{ m_model_object };
    std::unordered_map<std::string, std::unordered_map<int, std::array<double, 4>>> values_map;
    for (auto & atom : model_object->GetSelectedAtoms())
    {
        size_t id;
        if (AtomClassifier::IsMainChainMember(atom->GetSpot(), id) == false)
        {
            continue;
        }
        if (residue != Residue::UNK && atom->GetResidue() != residue)
        {
            continue;
        }
        if (atom->GetSpecialAtomFlag() == true)
        {
            continue;
        }
        bool is_veto_residue{ false };
        for (auto & veto_residue : veto_residues_list)
        {
            if (atom->GetResidue() == veto_residue)
            {
                is_veto_residue = true;
            }
        }
        if (is_veto_residue == true)
        {
            continue;
        }
        auto sequence_id{ atom->GetSequenceID() };
        auto chain_id{ atom->GetChainID() };
        auto * entry{ ModelAnalysisData::FindLocalEntry(*atom) };
        auto gaus_value{ entry->GetEstimateMDPDE().GetParameter(par_id) };
        values_map[chain_id][sequence_id].at(id) = gaus_value;
    }
    chain_size = static_cast<int>(values_map.size());

    std::vector<std::unique_ptr<TH1D>> hist_list;
    for (size_t i = 0; i < 4; i++)
    {
        auto name{
            Form("h%d_%d_%d_%d", static_cast<int>(extra_id), static_cast<int>(i),
                static_cast<int>(residue), par_id)
        };
        auto hist{ ROOTHelper::CreateHist1D(name, "", 4, 0.5, 4.5) };
        for (auto & [chain_id, values_map_tmp] : values_map)
        {
            for (auto & [sequence_id, values] : values_map_tmp)
            {
                hist->Fill(ArrayStats<double>::ComputeRank(values, i));
            }
        }
        hist_list.emplace_back(std::move(hist));
    }
    return hist_list;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateAtomGausEstimateToResidueGraph(
    std::vector<GroupKey> & group_key_list, const std::string & class_key, const int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    const auto model_view{ GetModelView() };

    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (IsAvailableAtomGroupKey(group_key, class_key) == false)
        {
            continue;
        }
        auto x_value{ static_cast<int>(model_view.GetResidueFromAtomGroupKey(group_key, class_key)) - 1 };
        auto y_value{ model_view.GetAtomGroupPrior(group_key, class_key).GetParameter(par_id) };
        auto y_error{ model_view.GetAtomGroupPriorPosterior(group_key, class_key).GetVariance(par_id) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateBondGausEstimateToResidueGraph(
    std::vector<GroupKey> & group_key_list, const std::string & class_key, const int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    const auto model_view{ GetModelView() };

    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (IsAvailableBondGroupKey(group_key, class_key) == false)
        {
            continue;
        }
        auto x_value{ static_cast<int>(model_view.GetResidueFromBondGroupKey(group_key, class_key)) - 1 };
        auto y_value{ model_view.GetBondGroupPrior(group_key, class_key).GetParameter(par_id) };
        auto y_error{ model_view.GetBondGroupPriorPosterior(group_key, class_key).GetVariance(par_id) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateAtomGausEstimateToSpotGraph(
    std::vector<GroupKey> & group_key_list, const std::string & class_key, const int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    const auto model_view{ GetModelView() };

    auto count{ 0 };
    for (size_t i = 0; i < group_key_list.size(); i++)
    {
        auto & group_key{ group_key_list[i] };
        if (IsAvailableAtomGroupKey(group_key, class_key) == false)
        {
            continue;
        }
        auto x_value{ static_cast<double>(i) };
        auto y_value{ model_view.GetAtomGroupPrior(group_key, class_key).GetParameter(par_id) };
        auto y_error{ model_view.GetAtomGroupPriorPosterior(group_key, class_key).GetVariance(par_id) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateAtomGausEstimateToAtomIdGraph(
    const std::map<std::string, GroupKey> & group_key_map,
    const std::vector<std::string> & atom_id_list,
    const std::string & class_key, const int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };
    const auto model_view{ GetModelView() };
    auto count{ 0 };
    for (size_t i = 0; i < atom_id_list.size(); i++)
    {
        auto & atom_id{ atom_id_list[i] };
        if (group_key_map.find(atom_id) == group_key_map.end())
        {
            continue;
        }
        auto & group_key{ group_key_map.at(atom_id) };
        if (IsAvailableAtomGroupKey(group_key, class_key) == false)
        {
            continue;
        }
        auto x_value{ static_cast<double>(i) };
        auto y_value{ model_view.GetAtomGroupPrior(group_key, class_key).GetParameter(par_id) };
        auto y_error{ model_view.GetAtomGroupPriorPosterior(group_key, class_key).GetVariance(par_id) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateAtomGausEstimateScatterGraph(
    GroupKey group_key, const std::string & class_key, int par1_id, int par2_id, bool select_outliers)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto atom_list{ GetModelView().GetAtomObjectList(group_key, class_key) };
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto count{ 0 };
    for (auto atom : atom_list)
    {
        auto * entry{ ModelAnalysisData::FindLocalEntry(*atom) };
        const auto * annotation{ entry->FindAnnotation(class_key) };
        auto is_outlier{ annotation != nullptr && annotation->is_outlier };
        if (select_outliers == true && is_outlier == false)
        {
            continue;
        }
        graph->SetPoint(
            count,
            entry->GetEstimateMDPDE().GetParameter(par1_id),
            entry->GetEstimateMDPDE().GetParameter(par2_id));
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateAtomGausEstimateScatterGraph(
    std::vector<GroupKey> & group_key_list, const std::string & class_key,
    int par1_id, int par2_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };

    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (IsAvailableAtomGroupKey(group_key, class_key) == false)
        {
            continue;
        }
        graph->SetPoint(
            count,
            GetModelView().GetAtomGausEstimatePrior(group_key, class_key, par1_id),
            GetModelView().GetAtomGausEstimatePrior(group_key, class_key, par2_id));
        graph->SetPointError(
            count,
            GetModelView().GetAtomGausVariancePrior(group_key, class_key, par1_id),
            GetModelView().GetAtomGausVariancePrior(group_key, class_key, par2_id));
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateBondGausEstimateScatterGraph(
    std::vector<GroupKey> & group_key_list, const std::string & class_key,
    int par1_id, int par2_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };

    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (IsAvailableBondGroupKey(group_key, class_key) == false)
        {
            continue;
        }
        graph->SetPoint(
            count,
            GetModelView().GetBondGausEstimatePrior(group_key, class_key, par1_id),
            GetModelView().GetBondGausEstimatePrior(group_key, class_key, par2_id));
        graph->SetPointError(
            count,
            GetModelView().GetBondGausVariancePrior(group_key, class_key, par1_id),
            GetModelView().GetBondGausVariancePrior(group_key, class_key, par2_id));
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateAtomGausEstimateScatterGraph(
    Element element, bool reverse)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }

    auto model_object{ m_model_object };
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto count{ 0 };
    for (auto & atom : model_object->GetSelectedAtoms())
    {
        if (atom->GetElement() != element)
        {
            continue;
        }
        auto * entry{ ModelAnalysisData::FindLocalEntry(*atom) };
        if (reverse == false)
        {
            graph->SetPoint(count, entry->GetEstimateMDPDE().amplitude, entry->GetEstimateMDPDE().width);
        }
        else
        {
            graph->SetPoint(count, entry->GetEstimateMDPDE().width, entry->GetEstimateMDPDE().amplitude);
        }
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateBondGausEstimateScatterGraph(
    Element element, bool reverse)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }

    auto model_object{ m_model_object };
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto count{ 0 };
    for (auto & bond : model_object->GetSelectedBonds())
    {
        if (bond->GetAtomObject1()->GetElement() != element)
        {
            continue;
        }
        auto * entry{ ModelAnalysisData::FindLocalEntry(*bond) };
        if (reverse == false)
        {
            graph->SetPoint(count, entry->GetEstimateMDPDE().amplitude, entry->GetEstimateMDPDE().width);
        }
        else
        {
            graph->SetPoint(count, entry->GetEstimateMDPDE().width, entry->GetEstimateMDPDE().amplitude);
        }
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateGausEstimateScatterGraph(
    GroupKey group_key1, GroupKey group_key2, const std::string & class_key, const int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    if (IsAvailableAtomGroupKey(group_key1, class_key) == false ||
        IsAvailableAtomGroupKey(group_key2, class_key) == false)
    {
        Logger::Log(LogLevel::Error, "Group key is not available.");
        return nullptr;
    }
    auto graph{ ROOTHelper::CreateGraphErrors() };

    const auto & atom_list1{ GetModelView().GetAtomObjectList(group_key1, class_key) };
    const auto & atom_list2{ GetModelView().GetAtomObjectList(group_key2, class_key) };

    auto count{ 0 };
    for (auto atom1 : atom_list1)
    {
        auto * entry1{ ModelAnalysisData::FindLocalEntry(*atom1) };
        for (auto atom2 : atom_list2)
        {
            auto * entry2{ ModelAnalysisData::FindLocalEntry(*atom2) };
            if (atom1->GetSequenceID() == atom2->GetSequenceID() &&
                atom1->GetChainID() == atom2->GetChainID())
            {
                graph->SetPoint(
                    count,
                    entry1->GetEstimateMDPDE().GetParameter(par_id),
                    entry2->GetEstimateMDPDE().GetParameter(par_id));
                count++;
                break;
            }
        }
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateDistanceToMapValueGraph()
{
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto count{ 0 };
    for (auto & [distance, map_value] : GetLocalEntry().GetDistanceAndMapValueList())
    {
        graph->SetPoint(count, distance, map_value);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateLinearModelDistanceToMapValueGraph()
{
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto count{ 0 };
    for (auto & [x, y] : series_ops::BuildLinearModelDistanceAndMapValueList(GetLocalEntry()))
    {
        graph->SetPoint(count, x, y);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateBinnedDistanceToMapValueGraph(
    int bin_size, double x_min, double x_max)
{
    auto data_array{ series_ops::BuildBinnedDistanceAndMapValueList(GetLocalEntry(), bin_size, x_min, x_max) };
    auto graph{ ROOTHelper::CreateGraphErrors(bin_size) };
    auto count{ 0 };
    for (auto & [distance, map_value] : data_array)
    {
        graph->SetPoint(count, distance, map_value);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateInRangeAtomsToGausEstimateGraph(
    GroupKey group_key, const std::string & class_key, double range, int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto model_object{ m_model_object };
    auto graph{ ROOTHelper::CreateGraphErrors() };

    auto count{ 0 };
    for (auto & atom : GetModelView().GetAtomObjectList(group_key, class_key))
    {
        auto in_range_atom_list{
            ModelDerivedState::Of(*model_object).FindAtomsInRange(*model_object, *atom, range) };
        auto * atom_entry{ ModelAnalysisData::FindLocalEntry(*atom) };
        graph->SetPoint(
            count,
            static_cast<double>(in_range_atom_list.size()),
            atom_entry->GetEstimateMDPDE().GetParameter(par_id));
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateCOMDistanceToGausEstimateGraph(
    GroupKey group_key, const std::string & class_key, int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto model_object{ m_model_object };
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto center_of_mass_pos{ model_object->GetCenterOfMassPosition() };

    auto count{ 0 };
    for (auto & atom : GetModelView().GetAtomObjectList(group_key, class_key))
    {
        const auto & atom_pos{ atom->GetPositionRef() };
        auto distance{ ArrayStats<float>::ComputeNorm(atom_pos, center_of_mass_pos) };
        auto * atom_entry{ ModelAnalysisData::FindLocalEntry(*atom) };
        graph->SetPoint(
            count,
            static_cast<double>(distance),
            atom_entry->GetEstimateMDPDE().GetParameter(par_id));
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateAtomXYPositionTomographyGraph(
    double normalized_z_pos, double z_ratio_window, bool com_center)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto model_object{ m_model_object };
    auto com_pos{
        com_center ? model_object->GetCenterOfMassPosition() : std::array<float, 3>{ 0.0, 0.0, 0.0 }
    };
    auto graph{ ROOTHelper::CreateGraphErrors() };
    auto z_pos{ model_object->GetModelPosition(2, normalized_z_pos) };
    auto window_width{ 0.5 * model_object->GetModelLength(2) * z_ratio_window };
    auto z_window_min{ z_pos - window_width };
    auto z_window_max{ z_pos + window_width };

    auto count{ 0 };
    for (auto & atom : model_object->GetAtomList())
    {
        auto position{ atom->GetPosition() };
        if (position.at(2) < z_window_min || position.at(2) >= z_window_max)
        {
            continue;
        }
        graph->SetPoint(count, position.at(0) - com_pos.at(0), position.at(1) - com_pos.at(1));
        count++;
    }
    return graph;
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomLocalLinearModelFunctionOLS() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    auto * atom_local_entry{ ModelAnalysisData::FindLocalEntry(*m_atom_object) };
    const auto beta{ atom_local_entry->GetEstimateOLS().ToBeta() };
    auto beta_0{ beta(0) };
    auto beta_1{ beta(1) };
    return ROOTHelper::CreateLinearModelFunction("linear", beta_0, beta_1);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomLocalLinearModelFunctionMDPDE() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    auto * atom_local_entry{ ModelAnalysisData::FindLocalEntry(*m_atom_object) };
    const auto beta{ atom_local_entry->GetEstimateMDPDE().ToBeta() };
    auto beta_0{ beta(0) };
    auto beta_1{ beta(1) };
    return ROOTHelper::CreateLinearModelFunction("linear", beta_0, beta_1);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomLocalGausFunctionOLS() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    auto * atom_local_entry{ ModelAnalysisData::FindLocalEntry(*m_atom_object) };
    auto amplitude{ atom_local_entry->GetEstimateOLS().amplitude };
    auto width{ atom_local_entry->GetEstimateOLS().width };
    return ROOTHelper::CreateGaus3DFunctionIn1D("gaus", amplitude, width);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomLocalGausFunctionMDPDE() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    auto * atom_local_entry{ ModelAnalysisData::FindLocalEntry(*m_atom_object) };
    auto amplitude{ atom_local_entry->GetEstimateMDPDE().amplitude };
    auto width{ atom_local_entry->GetEstimateMDPDE().width };
    return ROOTHelper::CreateGaus3DFunctionIn1D("gaus", amplitude, width);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomGroupLinearModelFunctionMean(
    GroupKey group_key, const std::string & class_key, double x_min, double x_max) const
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    const auto beta{ GetModelView().GetAtomGroupMean(group_key, class_key).ToBeta() };
    auto mu_0{ beta(0) };
    auto mu_1{ beta(1) };
    return ROOTHelper::CreateLinearModelFunction("linear_mean", mu_0, mu_1, x_min, x_max);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomGroupLinearModelFunctionPrior(
    GroupKey group_key, const std::string & class_key, double x_min, double x_max) const
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    const auto beta{ GetModelView().GetAtomGroupPrior(group_key, class_key).ToBeta() };
    auto mu_0{ beta(0) };
    auto mu_1{ beta(1) };
    return ROOTHelper::CreateLinearModelFunction("linear_prior", mu_0, mu_1, x_min, x_max);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomGroupGausFunctionMean(
    GroupKey group_key, const std::string & class_key) const
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    const auto & mean{ GetModelView().GetAtomGroupMean(group_key, class_key) };
    auto amplitude{ mean.amplitude };
    auto width{ mean.width };
    return ROOTHelper::CreateGaus3DFunctionIn1D("group_gaus_mean", amplitude, width);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomGroupGausFunctionPrior(
    GroupKey group_key, const std::string & class_key) const
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    const auto & prior{ GetModelView().GetAtomGroupPrior(group_key, class_key) };
    auto amplitude{ prior.amplitude };
    auto width{ prior.width };
    return ROOTHelper::CreateGaus3DFunctionIn1D("group_gaus_prior", amplitude, width);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateBondGroupGausFunctionPrior(
    GroupKey group_key, const std::string & class_key) const
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    const auto & prior{ GetModelView().GetBondGroupPrior(group_key, class_key) };
    auto amplitude{ prior.amplitude };
    auto width{ prior.width };
    return ROOTHelper::CreateGaus2DFunctionIn1D("group_gaus_prior", amplitude, width);
}

#endif

} // namespace rhbm_gem
#include "PotentialPlotBuilder.hpp"

#include "data/detail/AtomClassifier.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include "data/detail/BondClassifier.hpp"
#include <rhbm_gem/data/object/BondObject.hpp>
#include "data/detail/ModelAnalysisData.hpp"
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
    auto model_object{ m_model_object };
    auto graph{ ROOTHelper::CreateGraphErrors() };
    std::unordered_map<int, double> amplitude_diff_to_carbonyl_oxygen_map;
    for (auto & atom : model_object->GetSelectedAtoms())
    {
        if (atom->GetSpot() != Spot::O) continue;
        if (atom->GetSpecialAtomFlag() == false)
        {
            auto * entry{ ModelAnalysisData::FindLocalEntry(*atom) };
            auto sequence_id{ atom->GetSequenceID() };
            auto amplitude_estimate{ entry->GetEstimateMDPDE().amplitude };
            amplitude_diff_to_carbonyl_oxygen_map[sequence_id] = amplitude_estimate - reference_amplitude;
        }
    }
    auto count{ 0 };
    for (auto & atom : model_object->GetSelectedAtoms())
    {
        if (atom->GetElement() != element) continue;
        auto sequence_id{ atom->GetSequenceID() };
        auto * entry{ ModelAnalysisData::FindLocalEntry(*atom) };
        auto normalized_amplitude{ entry->GetEstimateMDPDE().amplitude };
        if (amplitude_diff_to_carbonyl_oxygen_map.find(sequence_id) != amplitude_diff_to_carbonyl_oxygen_map.end())
        {
            normalized_amplitude -= amplitude_diff_to_carbonyl_oxygen_map.at(sequence_id);
        }
        if (reverse == false)
        {
            graph->SetPoint(count, normalized_amplitude, entry->GetEstimateMDPDE().width);
        }
        else
        {
            graph->SetPoint(count, entry->GetEstimateMDPDE().width, normalized_amplitude);
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
    auto model_object{ m_model_object };
    auto graph{ ROOTHelper::CreateGraphErrors() };
    std::unordered_map<int, double> amplitude_diff_to_carbonyl_oxygen_map;
    for (auto & bond : model_object->GetSelectedBonds())
    {
        if (bond->GetSpecialBondFlag() == false)
        {
            auto * entry{ ModelAnalysisData::FindLocalEntry(*bond) };
            auto sequence_id{ bond->GetAtomObject1()->GetSequenceID() };
            auto amplitude_estimate{ entry->GetEstimateMDPDE().amplitude };
            amplitude_diff_to_carbonyl_oxygen_map[sequence_id] = amplitude_estimate - reference_amplitude;
        }
    }
    auto count{ 0 };
    for (auto & bond : model_object->GetSelectedBonds())
    {
        if (bond->GetAtomObject1()->GetElement() != element) continue;
        auto sequence_id{ bond->GetAtomObject1()->GetSequenceID() };
        auto * entry{ ModelAnalysisData::FindLocalEntry(*bond) };
        auto normalized_amplitude{ entry->GetEstimateMDPDE().amplitude };
        if (amplitude_diff_to_carbonyl_oxygen_map.find(sequence_id) != amplitude_diff_to_carbonyl_oxygen_map.end())
        {
            normalized_amplitude -= amplitude_diff_to_carbonyl_oxygen_map.at(sequence_id);
        }
        if (reverse == false)
        {
            graph->SetPoint(count, normalized_amplitude, entry->GetEstimateMDPDE().width);
        }
        else
        {
            graph->SetPoint(count, entry->GetEstimateMDPDE().width, normalized_amplitude);
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
    auto model_object{ m_model_object };

    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> graph_map;
    std::unordered_map<std::string, int> count_map;

    for (auto & atom : model_object->GetSelectedAtoms())
    {
        if (atom->GetElement() != AtomClassifier::GetMainChainElement(main_chain_element_id)) continue;
        if (atom->GetSpot() != AtomClassifier::GetMainChainSpot(main_chain_element_id)) continue;
        if (residue != Residue::UNK && atom->GetResidue() != residue) continue;
        auto * entry{ ModelAnalysisData::FindLocalEntry(*atom) };
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
    auto model_object{ m_model_object };

    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> graph_map;
    std::unordered_map<std::string, int> count_map;

    for (auto & atom : model_object->GetSelectedAtoms())
    {
        if (atom->GetElement() != AtomClassifier::GetMainChainElement(main_chain_element_id)) continue;
        if (atom->GetSpot() != AtomClassifier::GetMainChainSpot(main_chain_element_id)) continue;
        auto * entry{ ModelAnalysisData::FindLocalEntry(*atom) };
        if (entry == nullptr) continue;
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
    auto model_object{ m_model_object };

    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> graph_map;
    std::unordered_map<std::string, int> count_map;

    for (auto & atom : model_object->GetSelectedAtoms())
    {
        if (atom->GetElement() != AtomClassifier::GetMainChainElement(main_chain_element_id)) continue;
        if (atom->GetSpot() != AtomClassifier::GetMainChainSpot(main_chain_element_id)) continue;
        if (residue != Residue::UNK && atom->GetResidue() != residue) continue;
        auto * entry{ ModelAnalysisData::FindLocalEntry(*atom) };
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
            count_map[chain_id], x_value, entry->GetEstimateMDPDE().GetParameter(par_id));
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
    for (auto & bond : GetModelView().GetBondObjectList(group_key, class_key))
    {
        if (residue != Residue::UNK && bond->GetAtomObject1()->GetResidue() != residue) continue;
        auto * entry{ ModelAnalysisData::FindLocalEntry(*bond) };
        auto sequence_id{ bond->GetAtomObject1()->GetSequenceID() };
        auto chain_id{ bond->GetAtomObject1()->GetChainID() };
        if (sequence_id < 0) continue;
        if (graph_map.find(chain_id) == graph_map.end())
        {
            graph_map[chain_id] = ROOTHelper::CreateGraphErrors();
            count_map[chain_id] = 0;
        }
        auto x_value{ static_cast<double>(sequence_id) };
        auto y_value{ entry->GetEstimateMDPDE().GetParameter(par_id) };
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
    auto model_object{ m_model_object };

    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> graph_map;
    std::unordered_map<std::string, int> count_map;

    for (auto & atom : model_object->GetSelectedAtoms())
    {
        if (atom->GetElement() != AtomClassifier::GetMainChainElement(main_chain_element_id)) continue;
        if (atom->GetSpot() != AtomClassifier::GetMainChainSpot(main_chain_element_id)) continue;
        if (residue != Residue::UNK && atom->GetResidue() != residue) continue;
        auto * entry{ ModelAnalysisData::FindLocalEntry(*atom) };
        auto sequence_id{ atom->GetSequenceID() };
        auto chain_id{ atom->GetChainID() };
        if (sequence_id < 0) continue;
        if (graph_map.find(chain_id) == graph_map.end())
        {
            graph_map[chain_id] = ROOTHelper::CreateGraphErrors();
            count_map[chain_id] = 0;
        }
        auto x_value{ static_cast<double>(sequence_id) };
        const auto * annotation{ entry->FindAnnotation(class_key) };
        if (annotation == nullptr)
        {
            continue;
        }
        graph_map[chain_id]->SetPoint(
            count_map[chain_id], x_value, annotation->posterior.GetEstimate(par_id)
        );
        graph_map[chain_id]->SetPointError(
            count_map[chain_id], 0.0, annotation->posterior.GetVariance(par_id)
        );
        count_map[chain_id]++;
    }
    return graph_map;
}

#endif

} // namespace rhbm_gem
