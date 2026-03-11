#include <rhbm_gem/core/painter/PotentialPlotBuilder.hpp>

#include "internal/plot/IPlotRenderBackend.hpp"

#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/GroupPotentialEntry.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/math/KDTreeAlgorithm.hpp>

#ifdef HAVE_ROOT
#include <TGraphErrors.h>
#include <TF1.h>
#endif

#include <array>

namespace rhbm_gem {

#ifdef HAVE_ROOT

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateDistanceToMapValueGraph()
{
    auto graph{ m_plot_render_backend->CreateGraphErrors() };
    auto count{ 0 };
    for (auto & [distance, map_value] : GetDistanceAndMapValueList())
    {
        graph->SetPoint(count, distance, map_value);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateLinearModelDistanceToMapValueGraph()
{
    auto graph{ m_plot_render_backend->CreateGraphErrors() };
    auto count{ 0 };
    for (auto & [x, y] : GetLinearModelDistanceAndMapValueList())
    {
        graph->SetPoint(count, x, y);
        count++;
    }
    return graph;
}

std::unique_ptr<TGraphErrors> PotentialPlotBuilder::CreateBinnedDistanceToMapValueGraph(
    int bin_size, double x_min, double x_max)
{
    auto data_array{ GetBinnedDistanceAndMapValueList(bin_size, x_min, x_max) };
    auto graph{ m_plot_render_backend->CreateGraphErrors(bin_size) };
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
    auto graph{ m_plot_render_backend->CreateGraphErrors() };
    auto kd_tree_root{ m_model_object->GetKDTreeRoot() };
    if (kd_tree_root == nullptr)
    {
        Logger::Log(LogLevel::Error, "KDTree is not available for the model object.");
        return nullptr;
    }

    auto count{ 0 };
    for (auto & atom : GetAtomObjectList(group_key, class_key))
    {
        auto in_range_atom_list{ KDTreeAlgorithm<AtomObject>::RangeSearch(kd_tree_root, atom, range) };
        auto atom_entry{ atom->GetLocalPotentialEntry() };
        graph->SetPoint(count,
            static_cast<double>(in_range_atom_list.size()),
            atom_entry->GetGausEstimateMDPDE(par_id));
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
    auto graph{ m_plot_render_backend->CreateGraphErrors() };
    auto center_of_mass_pos{ m_model_object->GetCenterOfMassPosition() };

    auto count{ 0 };
    for (auto & atom : GetAtomObjectList(group_key, class_key))
    {
        const auto & atom_pos{ atom->GetPositionRef() };
        auto distance{ ArrayStats<float>::ComputeNorm(atom_pos, center_of_mass_pos) };
        auto atom_entry{ atom->GetLocalPotentialEntry() };
        graph->SetPoint(count,
            static_cast<double>(distance), atom_entry->GetGausEstimateMDPDE(par_id));
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
    auto com_pos{ com_center ? m_model_object->GetCenterOfMassPosition() : std::array<float, 3>{0.0, 0.0, 0.0} };
    auto graph{ m_plot_render_backend->CreateGraphErrors() };
    auto z_pos{ m_model_object->GetModelPosition(2, normalized_z_pos) };
    auto window_width{ 0.5 * m_model_object->GetModelLength(2) * z_ratio_window };
    auto z_window_min{ z_pos - window_width };
    auto z_window_max{ z_pos + window_width };

    auto count{ 0 };
    for (auto & atom : m_model_object->GetAtomList())
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
    auto beta_0{ m_atom_local_entry->GetBetaEstimateOLS(0) };
    auto beta_1{ m_atom_local_entry->GetBetaEstimateOLS(1) };
    return m_plot_render_backend->CreateLinearModelFunction("linear", beta_0, beta_1);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomLocalLinearModelFunctionMDPDE() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    auto beta_0{ m_atom_local_entry->GetBetaEstimateMDPDE(0) };
    auto beta_1{ m_atom_local_entry->GetBetaEstimateMDPDE(1) };
    return m_plot_render_backend->CreateLinearModelFunction("linear", beta_0, beta_1);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomLocalGausFunctionOLS() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    auto amplitude{ m_atom_local_entry->GetGausEstimateOLS(0) };
    auto width{ m_atom_local_entry->GetGausEstimateOLS(1) };
    return m_plot_render_backend->CreateGaus3DFunctionIn1D("gaus", amplitude, width);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomLocalGausFunctionMDPDE() const
{
    if (IsAtomLocalEntryAvailable() == false)
    {
        return nullptr;
    }
    auto amplitude{ m_atom_local_entry->GetGausEstimateMDPDE(0) };
    auto width{ m_atom_local_entry->GetGausEstimateMDPDE(1) };
    return m_plot_render_backend->CreateGaus3DFunctionIn1D("gaus", amplitude, width);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomGroupLinearModelFunctionMean(
    GroupKey group_key, const std::string & class_key, double x_min, double x_max) const
{
    auto mu_0
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetMuEstimateMean(group_key, 0)
    };
    auto mu_1
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetMuEstimateMean(group_key, 1)
    };
    return m_plot_render_backend->CreateLinearModelFunction("linear_mean", mu_0, mu_1, x_min, x_max);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomGroupLinearModelFunctionPrior(
    GroupKey group_key, const std::string & class_key, double x_min, double x_max) const
{
    auto mu_0
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetMuEstimatePrior(group_key, 0)
    };
    auto mu_1
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetMuEstimatePrior(group_key, 1)
    };
    return m_plot_render_backend->CreateLinearModelFunction("linear_prior", mu_0, mu_1, x_min, x_max);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomGroupGausFunctionMean(
    GroupKey group_key, const std::string & class_key) const
{
    auto amplitude
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimateMean(group_key, 0)
    };
    auto width
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimateMean(group_key, 1)
    };
    return m_plot_render_backend->CreateGaus3DFunctionIn1D("group_gaus_mean", amplitude, width);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateAtomGroupGausFunctionPrior(
    GroupKey group_key, const std::string & class_key) const
{
    auto amplitude
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, 0)
    };
    auto width
    {
        m_model_object->GetAtomGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, 1)
    };
    return m_plot_render_backend->CreateGaus3DFunctionIn1D("group_gaus_prior", amplitude, width);
}

std::unique_ptr<TF1> PotentialPlotBuilder::CreateBondGroupGausFunctionPrior(
    GroupKey group_key, const std::string & class_key) const
{
    auto amplitude
    {
        m_model_object->GetBondGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, 0)
    };
    auto width
    {
        m_model_object->GetBondGroupPotentialEntry(class_key)->GetGausEstimatePrior(group_key, 1)
    };
    return m_plot_render_backend->CreateGaus2DFunctionIn1D("group_gaus_prior", amplitude, width);
}

#endif

} // namespace rhbm_gem
