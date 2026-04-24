#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>

#ifdef HAVE_ROOT
#include <TAxis.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TPad.h>
#include <TPaveText.h>
#endif

namespace rhbm_gem::painter_internal {

inline std::string BuildPainterOutputLabel(const ModelObject & model_object)
{
    auto is_simulation{ model_object.GetEmdID().find("Simulation") != std::string::npos };
    auto label{ model_object.GetPdbID() };
    if (is_simulation)
    {
        label += "_" + model_object.GetEmdID()
                 + "_bw" + string_helper::ToStringWithPrecision(model_object.GetResolution(), 2);
    }
    label += ".pdf";
    return label;
}

#ifdef HAVE_ROOT
inline void PrintGausTitlePad(
    ::TPad * pad,
    ::TPaveText * text,
    const std::string & title,
    float text_size)
{
    pad->cd();
    auto left_margin{ 0.001 + pad->GetLeftMargin() * pad->GetAbsWNDC() };
    auto right_margin{ 0.001 + pad->GetRightMargin() * pad->GetAbsWNDC() };
    auto bottom_margin{ 0.005 + (1.0 - pad->GetTopMargin()) * pad->GetAbsHNDC() };
    root_helper::SetPaveTextMarginInCanvas(pad, text, left_margin, right_margin, bottom_margin, 0.005);
    root_helper::SetPaveTextDefaultStyle(text);
    root_helper::SetPaveAttribute(text, 0, 0.2);
    root_helper::SetFillAttribute(text, 1001, kAzure - 7);
    root_helper::SetTextAttribute(text, text_size, 23, 22, 0.0, kYellow - 10);
    root_helper::SetLineAttribute(text, 1, 0);
    text->AddText(title.data());
    pad->Update();
}

inline void PrintGausResultGlobalPad(
    ::TPad * pad,
    ::TH2 * hist,
    double left_margin,
    double right_margin,
    double bottom_margin,
    double top_margin,
    bool is_right_side_pad)
{
    pad->cd();
    root_helper::SetPadMarginInCanvas(pad, left_margin, right_margin, bottom_margin, top_margin);
    root_helper::SetPadLayout(pad, 1, 1, 0, 0);
    root_helper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0f);
    root_helper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0f);
    root_helper::SetAxisLabelAttribute(hist->GetXaxis(), 55.0f, 0.11f, 103, kCyan + 3);
    root_helper::SetAxisLabelAttribute(hist->GetYaxis(), 50.0f, 0.01f);

    auto x_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(pad, 0.008, 1) };
    root_helper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 21);
    root_helper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
    hist->GetXaxis()->SetLimits(-1.0, 20.0);
    hist->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    hist->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);
    for (size_t i = 0; i < ChemicalDataHelper::GetStandardAminoAcidCount(); i++)
    {
        auto residue{ ChemicalDataHelper::GetStandardAminoAcidList().at(i) };
        auto label{ ChemicalDataHelper::GetLabel(residue) };
        auto label_index{ static_cast<int>(i) + 2 };
        hist->GetXaxis()->ChangeLabel(label_index, 90.0, -1, 12, -1, -1, label.data());
    }
    hist->SetStats(0);
    hist->Draw((is_right_side_pad) ? "Y+" : "");
}

inline void BuildMapValueScatterGraph(
    GroupKey group_key,
    ::TGraphErrors * graph,
    ModelObject * model1,
    ModelObject * model2,
    int bin_size = 15,
    double x_min = 0.0,
    double x_max = 1.5)
{
    const ModelAnalysisView entry1_view{ *model1 };
    const ModelAnalysisView entry2_view{ *model2 };
    const auto & class_key{ ChemicalDataHelper::GetSimpleAtomClassKey() };
    if (!entry1_view.HasAtomGroup(group_key, class_key) ||
        !entry2_view.HasAtomGroup(group_key, class_key))
    {
        return;
    }
    const auto & group1{ entry1_view.GetAtomObjectList(group_key, class_key) };
    const auto & group2{ entry2_view.GetAtomObjectList(group_key, class_key) };

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
            LocalPotentialView::RequireFor(*atom_object1)
                .GetBinnedDistanceResponseSeries(bin_size, x_min, x_max)
        };
        auto data2_array{
            LocalPotentialView::RequireFor(*atom_object2)
                .GetBinnedDistanceResponseSeries(bin_size, x_min, x_max)
        };
        for (size_t i = 0; i < static_cast<size_t>(bin_size); i++)
        {
            graph->SetPoint(count, data1_array.at(i).response, data2_array.at(i).response);
            count++;
        }
    }
}
#endif

} // namespace rhbm_gem::painter_internal
