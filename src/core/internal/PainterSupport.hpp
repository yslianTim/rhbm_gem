#pragma once

#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/PotentialEntryQuery.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>

#include "PainterTypeCheck.hpp"

#ifdef HAVE_ROOT
#include <TAxis.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TPad.h>
#include <TPaveText.h>
#endif

namespace rhbm_gem::painter_internal {

template <typename ExpectedType>
void AppendPainterObject(
    DataObjectBase * data_object,
    std::string_view painter_name,
    std::string_view method_name,
    std::vector<ExpectedType *> & object_list)
{
    auto & typed_data_object{
        RequirePainterObject<ExpectedType>(data_object, painter_name, method_name) };
    object_list.push_back(&typed_data_object);
}

template <typename ExpectedType>
void AppendPainterReferenceObject(
    DataObjectBase * data_object,
    const std::string & label,
    std::string_view painter_name,
    std::string_view method_name,
    std::unordered_map<std::string, std::vector<ExpectedType *>> & object_map)
{
    auto & typed_data_object{
        RequirePainterObject<ExpectedType>(data_object, painter_name, method_name) };
    object_map[label].push_back(&typed_data_object);
}

inline std::string BuildPainterOutputLabel(const ModelObject & model_object)
{
    auto is_simulation{ model_object.GetEmdID().find("Simulation") != std::string::npos };
    auto label{ model_object.GetPdbID() };
    if (is_simulation)
    {
        label += "_" + model_object.GetEmdID()
                 + "_bw" + StringHelper::ToStringWithPrecision(model_object.GetResolution(), 2);
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
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, left_margin, right_margin, bottom_margin, 0.005);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.2);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure - 7);
    ROOTHelper::SetTextAttribute(text, text_size, 23, 22, 0.0, kYellow - 10);
    ROOTHelper::SetLineAttribute(text, 1, 0);
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
    ROOTHelper::SetPadMarginInCanvas(pad, left_margin, right_margin, bottom_margin, top_margin);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 55.0f, 0.11f, 103, kCyan + 3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 50.0f, 0.01f);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.008, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 21);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
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
    PotentialEntryQuery entry1_iter{ model1 };
    PotentialEntryQuery entry2_iter{ model2 };
    if (!entry1_iter.IsAvailableAtomGroupKey(group_key, ChemicalDataHelper::GetSimpleAtomClassKey()))
    {
        return;
    }
    if (!entry2_iter.IsAvailableAtomGroupKey(group_key, ChemicalDataHelper::GetSimpleAtomClassKey()))
    {
        return;
    }

    auto model1_atom_map{
        entry1_iter.GetAtomObjectMap(group_key, ChemicalDataHelper::GetSimpleAtomClassKey()) };
    auto model2_atom_map{
        entry2_iter.GetAtomObjectMap(group_key, ChemicalDataHelper::GetSimpleAtomClassKey()) };
    auto count{ 0 };
    for (auto & [atom_id, atom_object1] : model1_atom_map)
    {
        if (model2_atom_map.find(atom_id) == model2_atom_map.end())
        {
            continue;
        }
        auto * atom_object2{ model2_atom_map.at(atom_id) };
        PotentialEntryQuery atom1_iter{ atom_object1 };
        PotentialEntryQuery atom2_iter{ atom_object2 };
        auto data1_array{ atom1_iter.GetBinnedDistanceAndMapValueList(bin_size, x_min, x_max) };
        auto data2_array{ atom2_iter.GetBinnedDistanceAndMapValueList(bin_size, x_min, x_max) };
        for (size_t i = 0; i < static_cast<size_t>(bin_size); i++)
        {
            graph->SetPoint(count, std::get<1>(data1_array.at(i)), std::get<1>(data2_array.at(i)));
            count++;
        }
    }
}
#endif

} // namespace rhbm_gem::painter_internal
