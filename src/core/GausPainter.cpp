#include "GausPainter.hpp"
#include "ModelObject.hpp"
#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "DataObjectBase.hpp"
#include "PotentialEntryIterator.hpp"
#include "ChemicalComponentEntry.hpp"
#include "FilePathHelper.hpp"
#include "ChemicalDataHelper.hpp"
#include "ComponentHelper.hpp"
#include "AtomClassifier.hpp"
#include "BondClassifier.hpp"
#include "ArrayStats.hpp"
#include "GlobalEnumClass.hpp"
#include "AtomKeySystem.hpp"
#include "StringHelper.hpp"
#include "Logger.hpp"

#ifdef HAVE_ROOT
#include "ROOTHelper.hpp"
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TColor.h>
#include <TMarker.h>
#include <TAxis.h>
#include <TH2.h>
#include <TF1.h>
#include <TLine.h>
#include <TH3F.h>
#include <TH1.h>
#endif

GausPainter::GausPainter(void) :
    m_folder_path{ "./" },
    m_atom_classifier{ std::make_unique<AtomClassifier>() },
    m_bond_classifier{ std::make_unique<BondClassifier>() }
{

}

GausPainter::~GausPainter()
{

}

void GausPainter::SetFolder(const std::string & folder_path)
{
    m_folder_path = FilePathHelper::EnsureTrailingSlash(folder_path);
}

void GausPainter::AddDataObject(DataObjectBase * data_object)
{
    auto model_object{ dynamic_cast<ModelObject *>(data_object) };
    if (model_object == nullptr)
    {
        throw std::runtime_error(
            "GausPainter::AddDataObject(): invalid data_object type");
    }
    m_model_object_list.push_back(model_object);
}

void GausPainter::AddReferenceDataObject(DataObjectBase * data_object, const std::string & label)
{
    auto model_object{ dynamic_cast<ModelObject *>(data_object) };
    if (model_object == nullptr)
    {
        throw std::runtime_error(
            "GausPainter::AddReferenceDataObject(): invalid data_object type");
    }
    m_ref_model_object_list_map[label].push_back(model_object);
}

void GausPainter::Painting(void)
{
    Logger::Log(LogLevel::Info, "GausPainter::Painting() called.");
    Logger::Log(LogLevel::Info, "Folder path: " + m_folder_path);
    Logger::Log(LogLevel::Info, "Number of atom objects to be painted: "
                + std::to_string(m_model_object_list.size()));

    for (auto model_object : m_model_object_list)
    {
        auto is_simulation{ model_object->GetEmdID().find("Simulation") != std::string::npos };
        auto label{ model_object->GetPdbID() };
        if (is_simulation == true)
        {
            label += "_"+ model_object->GetEmdID() +
                     "_bw"+ StringHelper::ToStringWithPrecision(model_object->GetResolution(), 2);
        }
        label += ".pdf";
        
        PaintAtomGausInGroupSummary(model_object, "atom_gaus_in_group_"+ label);
        PaintAtomGroupGausSummary(model_object, "atom_group_gaus_summary_"+ label);
    }
}

void GausPainter::PaintAtomGausInGroupSummary(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " GausPainter::PaintAtomGausInGroup");

    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    const auto & chemical_component_map{ model_object->GetChemicalComponentEntryMap() };
    auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.0);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 750) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size{ 6 };

    std::unique_ptr<TPad> pad[pad_size];
    pad[0] = ROOTHelper::CreatePad("pad0","", 0.00, 0.80, 1.00, 1.00); // The title pad
    pad[1] = ROOTHelper::CreatePad("pad1","", 0.00, 0.00, 0.55, 0.80); // The left pad
    pad[2] = ROOTHelper::CreatePad("pad2","", 0.55, 0.00, 0.90, 0.60); // The bottom-middle pad
    pad[3] = ROOTHelper::CreatePad("pad3","", 0.55, 0.60, 0.90, 0.80); // The top-middle pad
    pad[4] = ROOTHelper::CreatePad("pad4","", 0.90, 0.00, 1.00, 0.60); // The bottom-right pad
    pad[5] = ROOTHelper::CreatePad("pad5","", 0.90, 0.60, 1.00, 0.80); // The top-right pad

    for (auto & [component_key, component_entry] : chemical_component_map)
    {
        const auto & component_atom_map{ component_entry->GetComponentAtomEntryMap() };
        auto component_id{ component_entry->GetComponentId() };
        for (auto & [atom_key, atom_entry] : component_atom_map)
        {
            if (atom_entry.element_type == Element::HYDROGEN) continue;
            if (atom_entry.atom_id == "OXT") continue;
            auto group_key{ AtomClassifier::GetGroupKeyInClass(component_key, atom_key) };
            if (entry_iter->IsAvailableAtomGroupKey(group_key, class_key) == false) continue;
            if (entry_iter->GetAtomObjectMap(group_key, class_key).size() < 1) continue;

            auto amplitude_hist{ entry_iter->CreateAtomGausEstimateHistogram(group_key, class_key, 0) };
            auto width_hist{ entry_iter->CreateAtomGausEstimateHistogram(group_key, class_key, 1) };

            canvas->cd();
            for (int i = 0; i < pad_size; i++)
            {
                ROOTHelper::SetPadDefaultStyle(pad[i].get());
                pad[i]->Draw();
            }

            pad[0]->cd();
            ROOTHelper::SetPadMarginInCanvas(gPad, 0.01, 0.01, 0.01, 0.01);
            auto resolution_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
            ROOTHelper::SetPaveTextMarginInCanvas(gPad, resolution_text.get(), 0.80, 0.01, 0.02, 0.02);
            ROOTHelper::SetPaveTextDefaultStyle(resolution_text.get());
            ROOTHelper::SetPaveAttribute(resolution_text.get(), 0, 0.1);
            ROOTHelper::SetFillAttribute(resolution_text.get(), 1001, kAzure-7);
            ROOTHelper::SetLineAttribute(resolution_text.get(), 1, 0);
            ROOTHelper::SetTextAttribute(resolution_text.get(), 75.0f, 133, 22, 0.0, kYellow-10);
            if (model_object->GetEmdID().find("Simulation") != std::string::npos)
            {
                resolution_text->AddText(Form("#sigma_{B}=%.2f", model_object->GetResolution()));
            }
            else
            {
                resolution_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
            }
            resolution_text->Draw();

            auto info_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
            ROOTHelper::SetPaveTextMarginInCanvas(gPad, info_text.get(), 0.55, 0.21, 0.02, 0.02);
            ROOTHelper::SetPaveTextDefaultStyle(info_text.get());
            ROOTHelper::SetPaveAttribute(info_text.get(), 0, 0.1);
            ROOTHelper::SetFillAttribute(info_text.get(), 1001, kAzure-7, 0.5);
            ROOTHelper::SetTextAttribute(info_text.get(), 55, 133, 12);
            info_text->AddText(("#font[102]{PDB-" + model_object->GetPdbID() +"}").data());
            info_text->AddText(("#font[102]{"+ model_object->GetEmdID() +"}").data());
            info_text->Draw();

            auto component_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
            ROOTHelper::SetPaveTextMarginInCanvas(gPad, component_text.get(), 0.01, 0.75, 0.02, 0.02);
            ROOTHelper::SetPaveTextDefaultStyle(component_text.get());
            ROOTHelper::SetPaveAttribute(component_text.get(), 0, 0.1);
            ROOTHelper::SetFillAttribute(component_text.get(), 1001, kAzure-7);
            ROOTHelper::SetTextAttribute(component_text.get(), 75.0f, 103, 22, 0.0, kYellow-10);
            component_text->AddText((component_id + "/" + atom_entry.atom_id).data());
            component_text->Draw();

            auto result_text{ ROOTHelper::CreatePaveText(0.26, 0.10, 0.55, 0.90, "nbNDC", false) };
            ROOTHelper::SetPaveTextDefaultStyle(result_text.get());
            ROOTHelper::SetTextAttribute(result_text.get(), 50.0f, 133, 12, 0.0, kRed);
            ROOTHelper::SetFillAttribute(result_text.get(), 4000);
            auto amplitude_prior{ entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 0) };
            auto amplitude_error{ entry_iter->GetAtomGausVariancePrior(group_key, class_key, 0) };
            auto width_prior{ entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 1) };
            auto width_error{ entry_iter->GetAtomGausVariancePrior(group_key, class_key, 1) };
            result_text->AddText(Form("#font[2]{A} = %.2f #pm %.2f", amplitude_prior, amplitude_error));
            result_text->AddText(Form("#tau = %.2f #pm %.2f", width_prior, width_error));
            result_text->Draw();

            pad[2]->cd();
            ROOTHelper::SetPadMarginInCanvas(gPad, 0.09, 0.01, 0.12, 0.02);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto frame_scatter{ ROOTHelper::CreateHist2D("hist_scatter","", 100, 0.0, 1.0, 100, 0.0, 1.0) };
            ROOTHelper::SetAxisTitleAttribute(frame_scatter->GetXaxis(), 40.0f, 1.1f);
            ROOTHelper::SetAxisTitleAttribute(frame_scatter->GetYaxis(), 40.0f, 1.2f);
            ROOTHelper::SetAxisLabelAttribute(frame_scatter->GetXaxis(), 40.0f, 0.01f);
            ROOTHelper::SetAxisLabelAttribute(frame_scatter->GetYaxis(), 40.0f, 0.01f);
            auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.03, 0) };
            auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.015, 1) };
            ROOTHelper::SetAxisTickAttribute(frame_scatter->GetXaxis(), static_cast<float>(x_tick_length), 505);
            ROOTHelper::SetAxisTickAttribute(frame_scatter->GetYaxis(), static_cast<float>(y_tick_length), 505);
            frame_scatter->GetXaxis()->SetTitle("Local Amplitude");
            frame_scatter->GetYaxis()->SetTitle("Local Width");
            frame_scatter->GetXaxis()->SetLimits(
                amplitude_hist->GetXaxis()->GetXmin(),
                amplitude_hist->GetXaxis()->GetXmax()
            );
            frame_scatter->GetYaxis()->SetLimits(
                width_hist->GetXaxis()->GetXmin(),
                width_hist->GetXaxis()->GetXmax()
            );
            frame_scatter->GetXaxis()->CenterTitle();
            frame_scatter->GetYaxis()->CenterTitle();
            frame_scatter->SetStats(0);
            frame_scatter->Draw();
            auto scatter_graph{ entry_iter->CreateAtomGausEstimateScatterGraph(group_key, class_key, 0, 1) };
            ROOTHelper::SetMarkerAttribute(scatter_graph.get(), 24, 1.5, kAzure-7, 1.0f);
            scatter_graph->Draw("P");

            pad[3]->cd();
            ROOTHelper::SetPadMarginInCanvas(gPad, 0.09, 0.01, 0.02, 0.02);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            ROOTHelper::SetAxisTitleAttribute(amplitude_hist->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisTitleAttribute(amplitude_hist->GetYaxis(), 40.0, 1.2f);
            ROOTHelper::SetAxisLabelAttribute(amplitude_hist->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(amplitude_hist->GetYaxis(), 40.0, 0.01f);
            ROOTHelper::SetAxisTickAttribute(amplitude_hist->GetXaxis(), 0.0f, 505);
            ROOTHelper::SetAxisTickAttribute(amplitude_hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
            ROOTHelper::SetFillAttribute(amplitude_hist.get(), 1001, kAzure-7, 0.5f);
            amplitude_hist->GetYaxis()->SetTitle("Counts");
            amplitude_hist->GetYaxis()->CenterTitle();
            amplitude_hist->SetStats(0);
            amplitude_hist->Draw("BAR");

            pad[4]->cd();
            ROOTHelper::SetPadMarginInCanvas(gPad, 0.01, 0.01, 0.12, 0.02);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            ROOTHelper::SetAxisTitleAttribute(width_hist->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisTitleAttribute(width_hist->GetYaxis(), 40.0, 1.0f);
            ROOTHelper::SetAxisLabelAttribute(width_hist->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(width_hist->GetYaxis(), 40.0, 0.01f);
            ROOTHelper::SetAxisTickAttribute(width_hist->GetXaxis(), 0.0f, 505);
            ROOTHelper::SetAxisTickAttribute(width_hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
            ROOTHelper::SetFillAttribute(width_hist.get(), 1001, kAzure-7, 0.5f);
            width_hist->GetYaxis()->SetTitle("Counts");
            width_hist->GetYaxis()->CenterTitle();
            width_hist->SetStats(0);
            width_hist->Draw("HBAR");

            pad[1]->cd();
            ROOTHelper::SetPadMarginInCanvas(gPad, 0.10, 0.00, 0.12, 0.10);
            auto member_size{ entry_iter->GetAtomObjectList(group_key, class_key).size() };
            std::vector<std::unique_ptr<TGraphErrors>> map_value_graph_list;
            std::vector<double> y_array;
            map_value_graph_list.reserve(member_size);
            y_array.reserve(member_size * 2);
            for (auto atom : entry_iter->GetAtomObjectList(group_key, class_key))
            {
                auto atom_iter{ std::make_unique<PotentialEntryIterator>(atom) };
                auto graph{ atom_iter->CreateBinnedDistanceToMapValueGraph() };
                ROOTHelper::SetLineAttribute(graph.get(), 1, 3, static_cast<short>(kAzure-7), 0.3f);
                map_value_graph_list.emplace_back(std::move(graph));
                auto map_value_range{ atom_iter->GetMapValueRange(0.0) };
                y_array.emplace_back(std::get<0>(map_value_range));
                y_array.emplace_back(std::get<1>(map_value_range));
            }
            auto gaus_function{ entry_iter->CreateAtomGroupGausFunctionPrior(group_key, class_key) };

            auto frame{ ROOTHelper::CreateHist2D("hist_0","", 100, 0.0, 1.0, 100, 0.0, 1.0) };
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 50.0f, 0.75f);
            ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), 50.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame->GetXaxis(), 0.06f, 505);
            ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 60.0f, 1.25f);
            ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 60.0f, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame->GetYaxis(), 0.03f, 506);
            ROOTHelper::SetLineAttribute(frame.get(), 1, 0);
            frame->SetStats(0);
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetTitle("Radial Distance #[]{#AA}");
            frame->GetYaxis()->SetTitle("Map Value");
            auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.15) };
            auto x_min{ 0.01 };
            auto x_max{ 1.49 };
            auto y_min{ std::get<0>(y_range) };
            auto y_max{ std::get<1>(y_range) };
            frame->GetXaxis()->SetLimits(x_min, x_max);
            frame->GetYaxis()->SetLimits(y_min, y_max);
            frame->SetStats(0);
            frame->Draw();
            for (auto & graph : map_value_graph_list)
            {
                graph->Draw("L X0");
            }
            auto line_ref{ ROOTHelper::CreateLine(0.01, 0.0, 1.49, 0.0) };
            ROOTHelper::SetLineAttribute(line_ref.get(), 2, 2, kBlack);
            line_ref->Draw();
            ROOTHelper::SetLineAttribute(gaus_function.get(), 2, 4, kRed+1);
            gaus_function->Draw("SAME");

            auto legend{ ROOTHelper::CreateLegend(0.20, 0.90, 1.00, 1.00, false) };
            ROOTHelper::SetLegendDefaultStyle(legend.get());
            ROOTHelper::SetFillAttribute(legend.get(), 4000);
            ROOTHelper::SetTextAttribute(legend.get(), 40.0f, 133, 12, 0.0);
            legend->SetMargin(0.15f);
            legend->SetNColumns(2);
            legend->AddEntry(gaus_function.get(),
                "Gaussian Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})}", "l");
            legend->AddEntry(map_value_graph_list.at(0).get(),
                "Map Value", "l");
            legend->Draw();

            auto count_text{ ROOTHelper::CreatePaveText(0.50, 0.80, 1.00, 0.90, "nbNDC", false) };
            ROOTHelper::SetPaveTextDefaultStyle(count_text.get());
            ROOTHelper::SetTextAttribute(count_text.get(), 40.0f, 133, 12, 0.0, kGray+1);
            ROOTHelper::SetFillAttribute(count_text.get(), 4000);
            count_text->AddText(Form("Number of Members: %ld ", member_size));
            count_text->Draw();
            
            ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
        }
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void GausPainter::PaintAtomGroupGausSummary(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " GausPainter::PaintAtomGroupGausSummary");

    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    const auto & chemical_component_map{ model_object->GetChemicalComponentEntryMap() };
    
    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.0);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 700) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size{ 7 };
    
    std::unique_ptr<TPad> pad[pad_size];
    std::unique_ptr<TH2> frame[pad_size];
    pad[0] = ROOTHelper::CreatePad("pad0","", 0.00, 0.00, 0.52, 0.45); // The bottom-left pad
    pad[1] = ROOTHelper::CreatePad("pad1","", 0.00, 0.45, 0.52, 0.80); // The top-left pad
    pad[2] = ROOTHelper::CreatePad("pad2","", 0.52, 0.00, 0.65, 0.45); // The bottom-middle pad
    pad[3] = ROOTHelper::CreatePad("pad3","", 0.52, 0.45, 0.65, 0.80); // The top-middle pad
    pad[4] = ROOTHelper::CreatePad("pad4","", 0.65, 0.00, 1.00, 0.80); // The bottom-right pad
    pad[5] = ROOTHelper::CreatePad("pad5","", 0.52, 0.80, 1.00, 1.00); // The title pad
    pad[6] = ROOTHelper::CreatePad("pad6","", 0.00, 0.80, 0.52, 1.00); // The histogram pad
    frame[1] = ROOTHelper::CreateHist2D("hist_1","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[0] = ROOTHelper::CreateHist2D("hist_0","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[3] = ROOTHelper::CreateHist2D("hist_3","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[2] = ROOTHelper::CreateHist2D("hist_2","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[4] = ROOTHelper::CreateHist2D("hist_4","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[5] = ROOTHelper::CreateHist2D("hist_5","", 100, 0.0, 1.0, 100, 0.0, 1.0);

    auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    for (auto & [component_key, component_entry] : chemical_component_map)
    {
        if (component_entry->GetComponentType() == "non-polymer") continue;
        const auto & component_atom_map{ component_entry->GetComponentAtomEntryMap() };
        auto component_id{ component_entry->GetComponentId() };
        std::map<Element, std::map<std::string, GroupKey>> group_key_list_map;
        std::vector<std::string> atom_id_list;
        atom_id_list.reserve(component_atom_map.size());
        GroupKey group_key_tmp{ 0 };
        for (auto & [atom_key, atom_entry] : component_atom_map)
        {
            if (atom_entry.element_type == Element::HYDROGEN) continue;
            if (atom_entry.atom_id == "OXT") continue;
            auto group_key{ AtomClassifier::GetGroupKeyInClass(component_key, atom_key) };
            group_key_list_map[atom_entry.element_type].emplace(atom_entry.atom_id, group_key);
            if (group_key_tmp == 0) group_key_tmp = group_key;
            atom_id_list.emplace_back(atom_entry.atom_id);
        }
        if (entry_iter->IsAvailableAtomGroupKey(group_key_tmp, class_key) == false) continue;

        std::map<Element, std::unique_ptr<TGraphErrors>> amplitude_graph_map;
        std::map<Element, std::unique_ptr<TGraphErrors>> width_graph_map;
        std::map<Element, std::unique_ptr<TGraphErrors>> correlation_graph_map;
        std::map<Element, std::unique_ptr<TH2D>> amplitude_hist_map;
        std::map<Element, std::unique_ptr<TH2D>> width_hist_map;
        std::vector<Element> element_list;
        std::vector<double> amplitude_array, width_array;
        amplitude_array.reserve(component_atom_map.size());
        width_array.reserve(component_atom_map.size());
        
        auto component_size{ entry_iter->GetAtomObjectList(group_key_tmp, class_key).size() };
        for (auto & [element, group_key_list] : group_key_list_map)
        {
            auto amplitude_graph{
                entry_iter->CreateAtomGausEstimateToAtomIdGraph(
                    group_key_list, atom_id_list, class_key, 0)
            };
            auto width_graph{
                entry_iter->CreateAtomGausEstimateToAtomIdGraph(
                    group_key_list, atom_id_list, class_key, 1)
            };
            std::vector<GroupKey> group_key_list_vector;
            group_key_list_vector.reserve(group_key_list.size());
            for (auto & [atom_id, group_key] : group_key_list)
            {
                group_key_list_vector.emplace_back(group_key);
            }
            auto correlation_graph{
                entry_iter->CreateAtomGausEstimateScatterGraph(group_key_list_vector, class_key, 0, 1)
            };
            for (int p = 0; p < amplitude_graph->GetN(); p++)
            {
                amplitude_array.push_back(amplitude_graph->GetPointY(p));
                width_array.push_back(width_graph->GetPointY(p));
            }
            auto component_color{ ChemicalDataHelper::GetDisplayColor(element) };
            auto component_marker{ ChemicalDataHelper::GetDisplayMarker(element) };
            ROOTHelper::SetMarkerAttribute(amplitude_graph.get(), component_marker, 1.3f, component_color);
            ROOTHelper::SetMarkerAttribute(width_graph.get(), component_marker, 1.3f, component_color);
            ROOTHelper::SetMarkerAttribute(correlation_graph.get(), component_marker, 1.3f, component_color);
            ROOTHelper::SetLineAttribute(amplitude_graph.get(), 1, 1, component_color);
            ROOTHelper::SetLineAttribute(width_graph.get(), 1, 1, component_color);
            ROOTHelper::SetLineAttribute(correlation_graph.get(), 1, 1, component_color);

            amplitude_graph_map[element] = std::move(amplitude_graph);
            width_graph_map[element] = std::move(width_graph);
            correlation_graph_map[element] = std::move(correlation_graph);
            element_list.emplace_back(element);
        }

        auto scaling{ 0.3 };
        auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, scaling) };
        auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, scaling) };
        auto element_count{ element_list.size() };
        std::vector<std::string> element_label_list;
        element_label_list.reserve(element_count);
        for (size_t i = 0; i < element_count; i++)
        {
            auto element{ element_list.at(i) };
            auto element_label{ ChemicalDataHelper::GetLabel(element) };
            auto x_value{ static_cast<double>(i) };
            element_label_list.emplace_back(element_label);
            std::string name_amplitude{ "amplitude_hist_"+ element_label };
            std::string name_width{ "width_hist_"+ element_label };
            auto amplitude_hist{
                ROOTHelper::CreateHist2D(
                    name_amplitude.data(),"",
                    static_cast<int>(element_count), -0.5, static_cast<int>(element_count)-0.5,
                    100, std::get<0>(amplitude_range), std::get<1>(amplitude_range))
            };
            auto width_hist{
                ROOTHelper::CreateHist2D(
                    name_width.data(),"",
                    static_cast<int>(element_count), -0.5, static_cast<int>(element_count)-0.5,
                    100, std::get<0>(width_range), std::get<1>(width_range))
            };
            for (int p = 0; p < amplitude_graph_map.at(element)->GetN(); p++)
            {
                amplitude_hist->Fill(x_value, amplitude_graph_map.at(element)->GetPointY(p));
            }
            for (int p = 0; p < width_graph_map.at(element)->GetN(); p++)
            {
                width_hist->Fill(x_value, width_graph_map.at(element)->GetPointY(p));
            }
            auto element_color{ ChemicalDataHelper::GetDisplayColor(element) };
            ROOTHelper::SetLineAttribute(amplitude_hist.get(), 1, 1, element_color);
            ROOTHelper::SetLineAttribute(width_hist.get(), 1, 1, element_color);
            ROOTHelper::SetFillAttribute(amplitude_hist.get(), 1001, element_color, 0.3f);
            ROOTHelper::SetFillAttribute(width_hist.get(), 1001, element_color, 0.3f);
            amplitude_hist->SetBarWidth(0.5f);
            width_hist->SetBarWidth(0.5f);
            amplitude_hist_map[element] = std::move(amplitude_hist);
            width_hist_map[element] = std::move(width_hist);
        }

        canvas->cd();
        for (int i = 0; i < pad_size; i++)
        {
            ROOTHelper::SetPadDefaultStyle(pad[i].get());
            pad[i]->Draw();
        }

        frame[1]->GetYaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        frame[0]->GetYaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));
        frame[3]->GetYaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        frame[2]->GetYaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));
        frame[4]->GetXaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        frame[4]->GetYaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));

        pad[1]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.07, 0.005, 0.01, 0.01);
        RemodelFrameInPad(frame[1].get(), pad[1].get(), 0.0, 0.015);
        RemodelAxisLabels(frame[1].get(), atom_id_list, 90.0, 12);
        ROOTHelper::SetAxisTitleAttribute(frame[1]->GetXaxis(), 0.0);
        ROOTHelper::SetAxisTitleAttribute(frame[1]->GetYaxis(), 35.0, 1.4f);
        ROOTHelper::SetAxisLabelAttribute(frame[1]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(frame[1]->GetYaxis(), 30.0, 0.01f);
        frame[1]->GetYaxis()->SetTitle("Amplitude");
        frame[1]->Draw();
        for (auto & [element, graph] : amplitude_graph_map) graph->Draw("P");

        pad[0]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.07, 0.005, 0.11, 0.01);
        RemodelFrameInPad(frame[0].get(), pad[0].get(), 0.0, 0.015);
        RemodelAxisLabels(frame[0].get(), atom_id_list, 90.0, 12);
        ROOTHelper::SetAxisTitleAttribute(frame[0]->GetXaxis(), 0.0);
        ROOTHelper::SetAxisTitleAttribute(frame[0]->GetYaxis(), 35.0, 1.4f);
        ROOTHelper::SetAxisLabelAttribute(frame[0]->GetXaxis(), 32.0, 0.11f, 103, kCyan+3);
        ROOTHelper::SetAxisLabelAttribute(frame[0]->GetYaxis(), 30.0, 0.01f);
        frame[0]->GetYaxis()->SetTitle("Width");
        frame[0]->Draw();
        for (auto & [element, graph] : width_graph_map) graph->Draw("P");

        pad[3]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.005, 0.005, 0.01, 0.01);
        RemodelFrameInPad(frame[3].get(), pad[3].get(), 0.0, 0.015);
        RemodelAxisLabels(frame[3].get(), element_label_list, 0.0, -1);
        ROOTHelper::SetAxisTitleAttribute(frame[3]->GetXaxis(), 35.0f, 1.0f);
        ROOTHelper::SetAxisTitleAttribute(frame[3]->GetYaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(frame[3]->GetXaxis(), 35.0f, 0.005f, 103);
        ROOTHelper::SetAxisLabelAttribute(frame[3]->GetYaxis(), 0.0f);
        frame[3]->GetXaxis()->SetTitle("Element");
        frame[3]->Draw();
        for (auto & [element, hist] : amplitude_hist_map) hist->Draw("CANDLE3 SAME");

        pad[2]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.005, 0.005, 0.11, 0.01);
        RemodelFrameInPad(frame[2].get(), pad[2].get(), 0.0, 0.015);
        RemodelAxisLabels(frame[2].get(), element_label_list, 0.0, -1);
        ROOTHelper::SetAxisTitleAttribute(frame[2]->GetXaxis(), 35.0f, 1.0f);
        ROOTHelper::SetAxisTitleAttribute(frame[2]->GetYaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(frame[2]->GetXaxis(), 35.0f, 0.005f, 103);
        ROOTHelper::SetAxisLabelAttribute(frame[2]->GetYaxis(), 0.0f);
        frame[2]->GetXaxis()->SetTitle("Element");
        frame[2]->Draw();
        for (auto & [element, hist] : width_hist_map) hist->Draw("CANDLE3 SAME");

        pad[4]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.005, 0.07, 0.11, 0.01);
        RemodelFrameInPad(frame[4].get(), pad[4].get(), 0.03, 0.015);
        ROOTHelper::SetAxisTitleAttribute(frame[4]->GetXaxis(), 35.0f, 1.0f);
        ROOTHelper::SetAxisTitleAttribute(frame[4]->GetYaxis(), 35.0f, 1.5f);
        ROOTHelper::SetAxisLabelAttribute(frame[4]->GetXaxis(), 35.0f);
        ROOTHelper::SetAxisLabelAttribute(frame[4]->GetYaxis(), 35.0f, 0.01f);
        frame[4]->GetXaxis()->SetTitle("Amplitude Estimate");
        frame[4]->GetYaxis()->SetTitle("Width Estimate");
        frame[4]->Draw("Y+");
        for (auto & [element, graph] : correlation_graph_map) graph->Draw("P X0");

        pad[5]->cd();
        auto info_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, info_text.get(), 0.005, 0.24, 0.01, 0.02);
        ROOTHelper::SetPaveTextDefaultStyle(info_text.get());
        ROOTHelper::SetPaveAttribute(info_text.get(), 0, 0.1);
        ROOTHelper::SetFillAttribute(info_text.get(), 1001, kAzure-7, 0.5);
        ROOTHelper::SetTextAttribute(info_text.get(), 55, 133, 12);
        info_text->AddText(("#font[102]{PDB-" + model_object->GetPdbID() +"}").data());
        info_text->AddText(("#font[102]{"+ model_object->GetEmdID() +"}").data());
        info_text->Draw();

        auto resolution_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, resolution_text.get(), 0.25, 0.07, 0.01, 0.02);
        ROOTHelper::SetPaveTextDefaultStyle(resolution_text.get());
        ROOTHelper::SetPaveAttribute(resolution_text.get(), 0, 0.1);
        ROOTHelper::SetFillAttribute(resolution_text.get(), 1001, kAzure-7);
        ROOTHelper::SetLineAttribute(resolution_text.get(), 1, 0);
        ROOTHelper::SetTextAttribute(resolution_text.get(), 75.0f, 133, 22, 0.0, kYellow-10);
        if (model_object->GetEmdID().find("Simulation") != std::string::npos)
        {
            resolution_text->AddText(Form("#sigma_{B}=%.2f", model_object->GetResolution()));
        }
        else
        {
            resolution_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
        }
        resolution_text->Draw();

        pad[6]->cd();
        auto component_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, component_text.get(), 0.01, 0.33, 0.01, 0.02);
        ROOTHelper::SetPaveTextDefaultStyle(component_text.get());
        ROOTHelper::SetPaveAttribute(component_text.get(), 0, 0.1);
        ROOTHelper::SetFillAttribute(component_text.get(), 1001, kAzure-7);
        ROOTHelper::SetTextAttribute(component_text.get(), 75.0f, 103, 22, 0.0, kYellow-10);
        component_text->AddText(component_entry->GetComponentId().data());
        component_text->Draw();

        auto component_info_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC", false) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, component_info_text.get(), 0.19, 0.01, 0.01, 0.02);
        ROOTHelper::SetPaveTextDefaultStyle(component_info_text.get());
        ROOTHelper::SetPaveAttribute(component_info_text.get(), 0);
        ROOTHelper::SetFillAttribute(component_info_text.get(), 4000);
        ROOTHelper::SetTextAttribute(component_info_text.get(), 30.0f, 133);
        auto component_name{ component_entry->GetComponentName() };
        StringHelper::EraseCharFromString(component_name, '\"');
        component_info_text->AddText(component_name.data());
        component_info_text->AddText(("Formula: " + component_entry->GetComponentFormula()).data());
        component_info_text->AddText(Form("Number of members: %zu", component_size));
        component_info_text->Draw();

        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

#ifdef HAVE_ROOT
void GausPainter::RemodelFrameInPad(
    TH2 * frame, TPad * pad, double x_tick, double y_tick)
{
    pad->cd();
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, x_tick, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, y_tick, 1) };
    ROOTHelper::SetAxisTickAttribute(frame->GetXaxis(), static_cast<float>(x_tick_length), 505);
    ROOTHelper::SetAxisTickAttribute(frame->GetYaxis(), static_cast<float>(y_tick_length), 505);
    frame->SetStats(0);
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
}

void GausPainter::RemodelAxisLabels(
    TH2 * frame, const std::vector<std::string> & label_list, double angle, int align)
{
    auto label_size{ static_cast<int>(label_list.size()) };
    ROOTHelper::SetAxisTickAttribute(frame->GetXaxis(), 0.0, label_size + 1);
    frame->GetXaxis()->SetLimits(-1.0, static_cast<float>(label_size));
    frame->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    frame->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);
    for (size_t i = 0; i < label_list.size(); i++)
    {
        auto label_index{ static_cast<int>(i) + 2 };
        frame->GetXaxis()->ChangeLabel(label_index, angle, -1, align, -1, -1, label_list.at(i).data());
    }
}
#endif
