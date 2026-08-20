#include "painter/detail/PotentialPlotBuilder.hpp"

#include "data/detail/AtomClassifier.hpp"
#include "data/detail/ModelDerivedState.hpp"

#include <rhbm_gem/core/QScoreHelper.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/KeyPacker.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/LocalPotentialSeries.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#endif

#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <unordered_map>

namespace rhbm_gem {

namespace {

LocalPotentialSampleList LoadSamplingEntries(
    const AtomLocalPotentialView & view,
    bool apply_selection,
    bool use_peeling_sampling_entries)
{
    return use_peeling_sampling_entries
        ? view.GetPeelingSamplingEntries(apply_selection)
        : view.GetRawSamplingEntries(apply_selection);
}

SeriesPointList BuildLocalDatasetSeries(
    const AtomLocalPotentialView & view,
    bool apply_selection,
    bool use_peeling_sampling_entries)
{
    auto model_prior{ view.GetEstimateMDPDE(FittingStage::Third) };
    auto offset{ model_prior.GetOffset() };
    double range_max{ 0.0 };
    auto sampling_entries{
        LoadSamplingEntries(view, apply_selection, use_peeling_sampling_entries)
    };
    for (auto & sample : sampling_entries)
    {
        auto distance{ static_cast<double>(sample.point.distance) };
        sample.response -= static_cast<float>(offset * model_prior.OffsetBasisAtDistance(distance));
        if (std::isfinite(distance) && distance >= 0.0 && distance > range_max)
        {
            range_max = distance;
        }
    }
    return linearization_service::BuildDatasetSeries(sampling_entries, 0.0, range_max);
}

std::vector<GroupKey> CollectComponentAtomGroupKeys(
    const ModelAnalysisView & model_view,
    AtomKey atom_key)
{
    std::vector<GroupKey> result;
    for (const auto group_key :
        model_view.CollectAtomGroupKeys(FittingStage::Third))
    {
        const auto unpacked_key{ KeyPackerComponentAtomClass::Unpack(group_key) };
        if (std::get<1>(unpacked_key) == atom_key)
        {
            result.emplace_back(group_key);
        }
    }
    return result;
}

} // namespace

PotentialPlotBuilder::PotentialPlotBuilder(ModelObject * model_object) :
    m_model_object{ model_object }
{
}

PotentialPlotBuilder::PotentialPlotBuilder(AtomObject * atom_object) :
    m_atom_object{ atom_object }
{
}

ModelAnalysisView PotentialPlotBuilder::GetModelView() const
{
    if (m_model_object == nullptr)
    {
        throw std::runtime_error("Model object is not available.");
    }
    return ModelAnalysisView(*m_model_object);
}

AtomLocalPotentialView PotentialPlotBuilder::GetLocalEntry() const
{
    if (m_atom_object != nullptr)
    {
        return AtomLocalPotentialView::RequireFor(*m_atom_object);
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
    return m_atom_object != nullptr && AtomLocalPotentialView::For(*m_atom_object).IsAvailable();
}

bool PotentialPlotBuilder::IsAvailableAtomGroupKey(GroupKey group_key) const
{
    return GetModelView().HasAtomGroup(FittingStage::Third, group_key);
}

#ifdef HAVE_ROOT
std::unique_ptr<TH1D> PotentialPlotBuilder::CreateComponentCountHistogram(
    const std::vector<GroupKey> & group_key_list) const
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto component_size{ static_cast<int>(group_key_list.size()) };
    auto hist{
        root_helper::CreateHist1D(
            "hist_component", "Component Count", component_size, -0.5, component_size - 0.5)
    };
    for (size_t i = 0; i < group_key_list.size(); i++)
    {
        auto count{ GetModelView().GetAtomObjectList(
            FittingStage::Third, group_key_list.at(i)).size() };
        hist->SetBinContent(static_cast<int>(i + 1), static_cast<double>(count));
    }
    return hist;
}

std::unique_ptr<TH1D> PotentialPlotBuilder::CreateAtomGausEstimateHistogram(
    GroupKey group_key, int par_id) const
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    if (IsAvailableAtomGroupKey(group_key) == false)
    {
        Logger::Log(LogLevel::Error, "Group key is not available.");
        return nullptr;
    }

    const auto & atom_list{
        GetModelView().GetAtomObjectList(FittingStage::Third, group_key)
    };
    std::vector<double> gaus_estimate_list;
    gaus_estimate_list.reserve(atom_list.size());
    for (auto atom : atom_list)
    {
        const auto local_entry{ AtomLocalPotentialView::RequireFor(*atom) };
        gaus_estimate_list.emplace_back(
            local_entry.GetEstimateMDPDE(FittingStage::Third)
                .GetDisplayParameter(par_id));
    }

    double x_min{ 0.0 };
    double x_max{ 1.0 };

    auto estimate_range{
        array_helper::ComputeScalingPercentileRangeTuple(gaus_estimate_list, 0.1, 0.05, 0.95)
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

    auto hist_name{ std::to_string(group_key) + "_par" + std::to_string(par_id) };
    auto hist{ root_helper::CreateHist1D(hist_name.data(), "", 25, x_min, x_max) };
    for (auto & value : gaus_estimate_list)
    {
        hist->Fill(value);
    }
    return hist;
}

std::unique_ptr<TH1D> PotentialPlotBuilder::CreateLinearModelDataHistogram(
    int dimension_id,
    bool apply_selection,
    bool use_peeling_sampling_entries) const
{
    auto data_array{
        BuildLocalDatasetSeries(GetLocalEntry(), apply_selection, use_peeling_sampling_entries)
    };
    std::vector<float> data_list;
    data_list.reserve(data_array.size());
    for (const auto & point : data_array)
    {
        switch (dimension_id)
        {
            case 0:
                data_list.emplace_back(static_cast<float>(point.GetBasisValue(1)));
                break;
            case 1:
                data_list.emplace_back(static_cast<float>(point.response));
                break;
            default:
                throw std::runtime_error("Dimension id is invalid.");
        }
    }

    auto data_range{ array_helper::ComputeScalingRangeTuple(data_list, 0.1f) };
    double x_min{ std::get<0>(data_range) };
    double x_max{ std::get<1>(data_range) };

    auto hist_name{ "data_dim" + std::to_string(dimension_id) };
    auto hist{ root_helper::CreateHist1D(hist_name.data(), "", 25, x_min, x_max) };
    for (auto & value : data_list)
    {
        hist->Fill(value);
    }
    return hist;
}

std::unique_ptr<TH2D> PotentialPlotBuilder::CreateDistanceToMapValueHistogram(
    int x_bin_size,
    int y_bin_size,
    bool apply_selection,
    bool use_peeling_sampling_entries) const
{
    const auto local_entry{ GetLocalEntry() };
    const auto sampling_entries{
        LoadSamplingEntries(local_entry, apply_selection, use_peeling_sampling_entries)
    };
    auto map_value_range{
        local_potential_series::ComputeResponseRange(sampling_entries, 0.1)
    };
    auto hist{
        root_helper::CreateHist2D(
            "hist_distance_mapvalue", "Distance vs Map Value",
            x_bin_size, 0.0, 2.0,
            y_bin_size, std::get<0>(map_value_range), std::get<1>(map_value_range))
    };
    for (const auto & sample : sampling_entries)
    {
        hist->Fill(sample.point.distance, sample.response);
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
        if (data_internal::IsMainChainMember(atom->GetSpot(), id) == false)
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
        const auto entry{ AtomLocalPotentialView::RequireFor(*atom) };
        auto gaus_value{
            entry.GetEstimateMDPDE(FittingStage::Third)
                .GetDisplayParameter(par_id)
        };
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
        auto hist{ root_helper::CreateHist1D(name, "", 4, 0.5, 4.5) };
        for (auto & [chain_id, values_map_tmp] : values_map)
        {
            for (auto & [sequence_id, values] : values_map_tmp)
            {
                hist->Fill(array_helper::ComputeRank(values, i));
            }
        }
        hist_list.emplace_back(std::move(hist));
    }
    return hist_list;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateAtomGausEstimateToResidueGraph(
    const std::vector<GroupKey> & group_key_list, const int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ root_helper::CreateGraphErrors() };
    const auto model_view{ GetModelView() };

    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (IsAvailableAtomGroupKey(group_key) == false)
        {
            continue;
        }
        auto x_value{ static_cast<int>(data_internal::GetResidueFromGroupKey(group_key)) - 1 };
        auto y_value{ model_view.GetAtomGroupPrior(
            FittingStage::Third, group_key).GetDisplayParameter(par_id) };
        auto y_error{ model_view.GetAtomGroupPriorWithUncertainty(
            FittingStage::Third, group_key).GetDisplayStandardDeviation(par_id) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateAtomGausEstimateToAtomIdGraph(
    const std::map<std::string, GroupKey> & group_key_map,
    const std::vector<std::string> & atom_id_list,
    const int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ root_helper::CreateGraphErrors() };
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
        if (IsAvailableAtomGroupKey(group_key) == false)
        {
            continue;
        }
        auto x_value{ static_cast<double>(i) };
        auto y_value{ model_view.GetAtomGroupPrior(
            FittingStage::Third, group_key).GetDisplayParameter(par_id) };
        auto y_error{ model_view.GetAtomGroupPriorWithUncertainty(
            FittingStage::Third, group_key).GetDisplayStandardDeviation(par_id) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, 0.0, y_error);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateAtomGausEstimateScatterGraph(
    GroupKey group_key, int par1_id, int par2_id, bool select_outliers)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto atom_list{ GetModelView().GetAtomObjectList(
        FittingStage::Third, group_key) };
    auto graph{ root_helper::CreateGraphErrors() };
    auto count{ 0 };
    for (auto atom : atom_list)
    {
        const auto entry{ AtomLocalPotentialView::RequireFor(*atom) };
        const auto & result{
            entry.GetGaussianResult(FittingStage::Third)
        };
        auto is_outlier{ result.posterior.has_value() && result.is_outlier };
        if (select_outliers == true && is_outlier == false)
        {
            continue;
        }
        graph->SetPoint(
            count,
            entry.GetEstimateMDPDE(FittingStage::Third)
                .GetDisplayParameter(par1_id),
            entry.GetEstimateMDPDE(FittingStage::Third)
                .GetDisplayParameter(par2_id));
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateAtomGausEstimateScatterGraph(
    const std::vector<GroupKey> & group_key_list,
    int par1_id, int par2_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto graph{ root_helper::CreateGraphErrors() };

    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (IsAvailableAtomGroupKey(group_key) == false)
        {
            continue;
        }
        graph->SetPoint(
            count,
            GetModelView().GetAtomGroupPrior(
                FittingStage::Third, group_key).GetDisplayParameter(par1_id),
            GetModelView().GetAtomGroupPrior(
                FittingStage::Third, group_key).GetDisplayParameter(par2_id));
        graph->SetPointError(
            count,
            GetModelView().GetAtomGroupPriorWithUncertainty(
                FittingStage::Third, group_key)
                .GetDisplayStandardDeviation(par1_id),
            GetModelView().GetAtomGroupPriorWithUncertainty(
                FittingStage::Third, group_key)
                .GetDisplayStandardDeviation(par2_id));
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
    auto graph{ root_helper::CreateGraphErrors() };
    auto count{ 0 };
    for (auto & atom : model_object->GetSelectedAtoms())
    {
        if (atom->GetElement() != element)
        {
            continue;
        }
        const auto entry{ AtomLocalPotentialView::RequireFor(*atom) };
        if (reverse == false)
        {
            graph->SetPoint(
                count,
                entry.GetEstimateMDPDE(FittingStage::Third).GetAmplitude(),
                entry.GetEstimateMDPDE(FittingStage::Third).GetWidth());
        }
        else
        {
            graph->SetPoint(
                count,
                entry.GetEstimateMDPDE(FittingStage::Third).GetWidth(),
                entry.GetEstimateMDPDE(FittingStage::Third).GetAmplitude());
        }
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateDistanceToMapValueGraph(
    bool apply_selection,
    bool use_peeling_sampling_entries)
{
    auto graph{ root_helper::CreateGraphErrors() };
    auto count{ 0 };
    for (const auto & sample : LoadSamplingEntries(GetLocalEntry(), apply_selection, use_peeling_sampling_entries))
    {
        graph->SetPoint(count, sample.point.distance, sample.response);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateLinearModelDistanceToMapValueGraph(
    bool apply_selection,
    bool use_peeling_sampling_entries)
{
    auto graph{ root_helper::CreateGraphErrors() };
    auto count{ 0 };
    for (const auto & point : BuildLocalDatasetSeries(GetLocalEntry(), apply_selection, use_peeling_sampling_entries))
    {
        graph->SetPoint(count, point.GetBasisValue(1), point.response);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateBinnedDistanceToRawMapValueGraph(
    int bin_size, double x_min, double x_max)
{
    auto data_array{
        local_potential_series::BuildBinnedDistanceResponseSeries(
            GetLocalEntry().GetRawSamplingEntries(true), bin_size, x_min, x_max)
    };
    auto graph{ root_helper::CreateGraphErrors(bin_size) };
    auto count{ 0 };
    for (const auto & point : data_array)
    {
        graph->SetPoint(count, point.GetBasisValue(0), point.response);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateBinnedDistanceToPeelingMapValueGraph(
    int bin_size, double x_min, double x_max)
{
    auto data_array{
        local_potential_series::BuildBinnedDistanceResponseSeries(
            GetLocalEntry().GetPeelingSamplingEntries(false), bin_size, x_min, x_max)
    };
    auto graph{ root_helper::CreateGraphErrors(bin_size) };
    auto count{ 0 };
    for (const auto & point : data_array)
    {
        graph->SetPoint(count, point.GetBasisValue(0), point.response);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateBinnedDistanceToNeighborMapValueGraph(
    int bin_size, double x_min, double x_max)
{
    auto raw_sampling_entries{ GetLocalEntry().GetRawSamplingEntries(false) };
    auto peeling_sampling_entries{ GetLocalEntry().GetPeelingSamplingEntries(false) };
    auto neighbor_sampling_entries{ raw_sampling_entries };
    for (size_t i = 0; i < raw_sampling_entries.size(); i++)
    {
        neighbor_sampling_entries.at(i).response = raw_sampling_entries.at(i).response - peeling_sampling_entries.at(i).response;
    }
    auto data_array{
        local_potential_series::BuildBinnedDistanceResponseSeries(neighbor_sampling_entries, bin_size, x_min, x_max)
    };
    auto graph{ root_helper::CreateGraphErrors(bin_size) };
    auto count{ 0 };
    for (const auto & point : data_array)
    {
        graph->SetPoint(count, point.GetBasisValue(0), point.response);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateInRangeAtomsToGausEstimateGraph(
    GroupKey group_key, double range, int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto model_object{ m_model_object };
    auto graph{ root_helper::CreateGraphErrors() };

    auto count{ 0 };
    for (auto & atom : GetModelView().GetAtomObjectList(
        FittingStage::Third, group_key))
    {
        auto in_range_atom_list{
            ModelDerivedState::Of(*model_object).FindAtomsInRange(*model_object, *atom, range) };
        const auto atom_entry{ AtomLocalPotentialView::RequireFor(*atom) };
        graph->SetPoint(
            count,
            static_cast<double>(in_range_atom_list.size()),
            atom_entry.GetEstimateMDPDE(FittingStage::Third)
                .GetDisplayParameter(par_id));
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateCOMDistanceToGausEstimateGraph(
    GroupKey group_key, int par_id)
{
    if (IsModelObjectAvailable() == false)
    {
        return nullptr;
    }
    auto model_object{ m_model_object };
    auto graph{ root_helper::CreateGraphErrors() };
    auto center_of_mass_pos{ model_object->GetCenterOfMassPosition() };

    auto count{ 0 };
    for (auto & atom : GetModelView().GetAtomObjectList(
        FittingStage::Third, group_key))
    {
        const auto & atom_pos{ atom->GetPositionRef() };
        auto distance{ array_helper::ComputeNorm(atom_pos, center_of_mass_pos) };
        const auto atom_entry{ AtomLocalPotentialView::RequireFor(*atom) };
        graph->SetPoint(
            count,
            static_cast<double>(distance),
            atom_entry.GetEstimateMDPDE(FittingStage::Third)
                .GetDisplayParameter(par_id));
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
    auto graph{ root_helper::CreateGraphErrors() };
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

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateMapValueScatterGraph(
    AtomKey atom_key,
    ModelObject * model1,
    ModelObject * model2,
    int bin_size,
    double x_min,
    double x_max)
{
    auto graph{ root_helper::CreateGraphErrors() };
    const ModelAnalysisView entry1_view{ *model1 };
    const ModelAnalysisView entry2_view{ *model2 };
    const auto group1{ CollectComponentAtomMembers(entry1_view, atom_key) };
    const auto group2{ CollectComponentAtomMembers(entry2_view, atom_key) };
    if (group1.empty() || group2.empty())
    {
        return graph;
    }

    std::unordered_map<int, AtomObject *> model1_atom_map;
    model1_atom_map.reserve(group1.size());
    for (auto * atom_object : group1)
    {
        model1_atom_map[atom_object->GetSerialID()] = atom_object;
    }

    std::unordered_map<int, AtomObject *> model2_atom_map;
    model2_atom_map.reserve(group2.size());
    for (auto * atom_object : group2)
    {
        model2_atom_map[atom_object->GetSerialID()] = atom_object;
    }

    auto count{ 0 };
    for (auto & [atom_id, atom_object1] : model1_atom_map)
    {
        if (model2_atom_map.find(atom_id) == model2_atom_map.end())
        {
            continue;
        }
        auto * atom_object2{ model2_atom_map.at(atom_id) };
        auto data1_array{
            local_potential_series::BuildBinnedDistanceResponseSeries(
                AtomLocalPotentialView::RequireFor(*atom_object1).GetRawSamplingEntries(),
                bin_size,
                x_min,
                x_max)
        };
        auto data2_array{
            local_potential_series::BuildBinnedDistanceResponseSeries(
                AtomLocalPotentialView::RequireFor(*atom_object2).GetRawSamplingEntries(),
                bin_size,
                x_min,
                x_max)
        };
        for (size_t i = 0; i < static_cast<size_t>(bin_size); i++)
        {
            graph->SetPoint(count, data1_array.at(i).response, data2_array.at(i).response);
            count++;
        }
    }
    return graph;
}

std::vector<AtomObject *> PotentialPlotBuilder::CollectComponentAtomMembers(
    const ModelAnalysisView & model_view,
    AtomKey atom_key)
{
    std::vector<AtomObject *> result;
    for (const auto group_key : CollectComponentAtomGroupKeys(model_view, atom_key))
    {
        const auto & members{ model_view.GetAtomObjectList(
            FittingStage::Third, group_key) };
        result.insert(result.end(), members.begin(), members.end());
    }
    return result;
}

std::optional<GaussianModel3DWithUncertainty> PotentialPlotBuilder::ComputeComponentAtomAveragePrior(
    const ModelAnalysisView & model_view,
    AtomKey atom_key)
{
    double amplitude{ 0.0 };
    double width{ 0.0 };
    double offset{ 0.0 };
    double amplitude_sd{ 0.0 };
    double width_sd{ 0.0 };
    double offset_sd{ 0.0 };
    std::size_t count{ 0 };
    for (const auto group_key : CollectComponentAtomGroupKeys(model_view, atom_key))
    {
        const auto prior{ model_view.GetAtomGroupPriorWithUncertainty(
            FittingStage::Third, group_key) };
        const auto & model{ prior.GetModel() };
        const auto & uncertainty{ prior.GetStandardDeviationModel() };
        amplitude += model.GetAmplitude();
        width += model.GetWidth();
        offset += model.GetOffset();
        amplitude_sd += uncertainty.GetAmplitude();
        width_sd += uncertainty.GetWidth();
        offset_sd += uncertainty.GetOffset();
        count++;
    }
    if (count == 0)
    {
        return std::nullopt;
    }
    const auto scale{ 1.0 / static_cast<double>(count) };
    return GaussianModel3DWithUncertainty{
        GaussianModel3D{ amplitude * scale, width * scale, offset * scale },
        GaussianModel3DUncertainty{ amplitude_sd * scale, width_sd * scale, offset_sd * scale }
    };
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomLocalLinearModelFunctionOLS() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    const auto atom_local_entry{ AtomLocalPotentialView::RequireFor(*m_atom_object) };
    const auto beta{
        linearization_service::EncodeGaussianToParameterVector(
            atom_local_entry.GetEstimateOLS(FittingStage::Third))
    };
    auto beta_0{ beta(0) };
    auto beta_1{ beta(1) };
    return root_helper::CreateLinearModelFunction("linear", beta_0, beta_1);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomLocalLinearModelFunctionMDPDE() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    const auto atom_local_entry{ AtomLocalPotentialView::RequireFor(*m_atom_object) };
    const auto beta{
        linearization_service::EncodeGaussianToParameterVector(
            atom_local_entry.GetEstimateMDPDE(FittingStage::Third))
    };
    auto beta_0{ beta(0) };
    auto beta_1{ beta(1) };
    return root_helper::CreateLinearModelFunction("linear", beta_0, beta_1);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomLocalGausFunctionOLS() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    const auto atom_local_entry{ AtomLocalPotentialView::RequireFor(*m_atom_object) };
    const auto & model{
        atom_local_entry.GetEstimateOLS(FittingStage::Third)
    };
    auto amplitude{ model.GetAmplitude() };
    auto width{ model.GetWidth() };
    auto offset{ model.GetOffset() };
    return root_helper::CreateGaus3DFunctionIn1D("gaus", amplitude, width, offset);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomLocalGausFunctionMDPDE() const
{
    if (IsAtomLocalEntryAvailable() == false) return nullptr;
    const auto atom_local_entry{ AtomLocalPotentialView::RequireFor(*m_atom_object) };
    const auto & model{
        atom_local_entry.GetEstimateMDPDE(FittingStage::Third)
    };
    auto amplitude{ model.GetAmplitude() };
    auto width{ model.GetWidth() };
    auto offset{ model.GetOffset() };
    return root_helper::CreateGaus3DFunctionIn1D("gaus", amplitude, width, offset);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomGroupGausFunctionMean(
    FittingStage stage, GroupKey group_key) const
{
    if (IsModelObjectAvailable() == false) return nullptr;
    const auto & mean{ GetModelView().GetAtomGroupMean(stage, group_key) };
    auto amplitude{ mean.GetAmplitude() };
    auto width{ mean.GetWidth() };
    auto offset{ mean.GetOffset() };
    return root_helper::CreateGaus3DFunctionIn1D("group_gaus_mean", amplitude, width, offset);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomGroupGausFunctionPrior(
    FittingStage stage, GroupKey group_key) const
{
    if (IsModelObjectAvailable() == false) return nullptr;
    const auto & prior{ GetModelView().GetAtomGroupPrior(stage, group_key) };
    auto amplitude{ prior.GetAmplitude() };
    auto width{ prior.GetWidth() };
    auto offset{ prior.GetOffset() };
    return root_helper::CreateGaus3DFunctionIn1D("group_gaus_prior", amplitude, width, offset);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateComponentAtomAverageGausFunctionPrior(AtomKey atom_key) const
{
    if (IsModelObjectAvailable() == false) return nullptr;
    const auto prior{ ComputeComponentAtomAveragePrior(GetModelView(), atom_key) };
    if (!prior.has_value()) return nullptr;
    const auto & model{ prior->GetModel() };
    return root_helper::CreateGaus3DFunctionIn1D(
        "component_atom_average_gaus_prior",
        model.GetAmplitude(),
        model.GetWidth(),
        model.GetOffset());
}

#endif

} // namespace rhbm_gem

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
    auto graph{ root_helper::CreateGraphErrors() };
    std::unordered_map<int, double> amplitude_diff_to_carbonyl_oxygen_map;
    for (auto & atom : model_object->GetSelectedAtoms())
    {
        if (atom->GetSpot() != Spot::O) continue;
        if (atom->GetSpecialAtomFlag() == false)
        {
            const auto entry{ AtomLocalPotentialView::RequireFor(*atom) };
            auto sequence_id{ atom->GetSequenceID() };
            auto amplitude_estimate{
                entry.GetEstimateMDPDE(FittingStage::Third).GetAmplitude()
            };
            amplitude_diff_to_carbonyl_oxygen_map[sequence_id] = amplitude_estimate - reference_amplitude;
        }
    }
    auto count{ 0 };
    for (auto & atom : model_object->GetSelectedAtoms())
    {
        if (atom->GetElement() != element) continue;
        auto sequence_id{ atom->GetSequenceID() };
        const auto entry{ AtomLocalPotentialView::RequireFor(*atom) };
        auto normalized_amplitude{
            entry.GetEstimateMDPDE(FittingStage::Third).GetAmplitude()
        };
        if (amplitude_diff_to_carbonyl_oxygen_map.find(sequence_id) != amplitude_diff_to_carbonyl_oxygen_map.end())
        {
            normalized_amplitude -= amplitude_diff_to_carbonyl_oxygen_map.at(sequence_id);
        }
        if (reverse == false)
        {
            graph->SetPoint(
                count,
                normalized_amplitude,
                entry.GetEstimateMDPDE(FittingStage::Third).GetWidth());
        }
        else
        {
            graph->SetPoint(
                count,
                entry.GetEstimateMDPDE(FittingStage::Third).GetWidth(),
                normalized_amplitude);
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
        if (atom->GetElement() != data_internal::GetMainChainElement(main_chain_element_id)) continue;
        if (atom->GetSpot() != data_internal::GetMainChainSpot(main_chain_element_id)) continue;
        if (residue != Residue::UNK && atom->GetResidue() != residue) continue;
        const auto entry{ AtomLocalPotentialView::RequireFor(*atom) };
        auto sequence_id{ atom->GetSequenceID() };
        auto chain_id{ atom->GetChainID() };
        if (sequence_id < 0) continue;
        if (graph_map.find(chain_id) == graph_map.end())
        {
            graph_map[chain_id] = root_helper::CreateGraphErrors();
            count_map[chain_id] = 0;
        }
        auto x_value{ static_cast<double>(sequence_id) };
        graph_map[chain_id]->SetPoint(
            count_map[chain_id],
            x_value,
            local_potential_series::ComputeMapValueNearCenter(entry.GetRawSamplingEntries()));
        count_map[chain_id]++;
    }
    return graph_map;
}

std::unordered_map<std::string, std::unique_ptr<TGraphErrors>>
PotentialPlotBuilder::CreateAverageQScoreToSequenceIDGraphMap(
    bool use_standard,
    bool use_peeling_sampling_entries)
{
    if (IsModelObjectAvailable() == false) return {};
    auto apply_selection{ false };
    auto model_object{ m_model_object };
    std::vector<std::string> chain_id_list;
    for (auto & [entity_id, chain_ids] : model_object->GetChainIDListMap())
    {
        (void) entity_id;
        if (chain_ids.empty()) continue;
        chain_id_list.insert(chain_id_list.end(), chain_ids.begin(), chain_ids.end());
    }

    auto reference_height{ model_object->GetReferenceHeight() };
    auto reference_offset{ model_object->GetReferenceOffset() };
    auto reference_width{ 0.6 };
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> graph_map;
    for (auto & chain_id : chain_id_list)
    {
        std::unordered_map<int, std::vector<double>> q_scores_map;
        for (auto & atom : model_object->GetSelectedAtoms())
        {
            if (atom->GetChainID() != chain_id) continue;
            const auto entry{ AtomLocalPotentialView::For(*atom) };
            if (!entry.IsAvailable()) continue;
            auto sequence_id{ atom->GetSequenceID() };
            if (sequence_id < 0) continue;
            if (use_peeling_sampling_entries)
            {
                reference_height = entry.GetGaussianResult(FittingStage::Third)
                    .mdpde.GetModel().GetHeight();
                reference_offset = entry.GetGaussianResult(FittingStage::Third)
                    .mdpde.GetModel().GetOffset();
                reference_width = entry.GetGaussianResult(FittingStage::Third)
                    .mdpde.GetModel().GetWidth();
            }
            if (!use_peeling_sampling_entries) apply_selection = true;
            auto q_score{ use_standard ?
                atom->GetStandardQScore() :
                core::CalculateQScoreForAtom(
                    LoadSamplingEntries(
                        entry,
                        apply_selection,
                        use_peeling_sampling_entries),
                    reference_height,
                    reference_offset,
                    reference_width
                )
            };
            q_scores_map[sequence_id].emplace_back(q_score);
        }
        if (q_scores_map.empty()) continue;

        graph_map[chain_id] = root_helper::CreateGraphErrors();
        int count{ 0 };
        for (const auto & [seq_id, q_score_vec] : q_scores_map)
        {
            auto q_score_average{ array_helper::ComputeMean(q_score_vec.data(), q_score_vec.size()) };
            graph_map[chain_id]->SetPoint(count, static_cast<double>(seq_id), q_score_average);
            count++;
        }
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
        if (atom->GetElement() != data_internal::GetMainChainElement(main_chain_element_id)) continue;
        if (atom->GetSpot() != data_internal::GetMainChainSpot(main_chain_element_id)) continue;
        if (residue != Residue::UNK && atom->GetResidue() != residue) continue;
        const auto entry{ AtomLocalPotentialView::RequireFor(*atom) };
        auto sequence_id{ atom->GetSequenceID() };
        auto chain_id{ atom->GetChainID() };
        if (sequence_id < 0) continue;
        if (graph_map.find(chain_id) == graph_map.end())
        {
            graph_map[chain_id] = root_helper::CreateGraphErrors();
            count_map[chain_id] = 0;
        }
        auto x_value{ static_cast<double>(sequence_id) };
        graph_map[chain_id]->SetPoint(
            count_map[chain_id],
            x_value,
            entry.GetEstimateMDPDE(FittingStage::Third)
                .GetDisplayParameter(par_id));
        count_map[chain_id]++;
    }
    return graph_map;
}

#endif

} // namespace rhbm_gem
