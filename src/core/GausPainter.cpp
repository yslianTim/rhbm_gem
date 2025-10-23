#include "GausPainter.hpp"
#include "ModelObject.hpp"
#include "AtomObject.hpp"
#include "BondObject.hpp"
#include "DataObjectBase.hpp"
#include "PotentialEntryIterator.hpp"
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

#include <vector>
#include <tuple>

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
        
        PaintAtomGroupGausSummary(model_object, "atom_group_gaus_summary_"+ label);
    }
}

void GausPainter::PaintAtomGroupGausSummary(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " GausPainter::PaintAtomGroupGausSummary");

    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    
    #ifdef HAVE_ROOT

    const int component_size{ 4 };
    Residue component_id[component_size]{ Residue::A, Residue::C, Residue::G, Residue::U };

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

    std::vector<GroupKey> group_key_list[component_size];
    for (size_t i = 0; i < component_size; i++)
    {
        group_key_list[i] =
            m_atom_classifier->GetNucleotideMainChainComponentAtomClassGroupKeyList(component_id[i]);
    }

    std::unique_ptr<TGraphErrors> amplitude_graph[component_size];
    std::unique_ptr<TGraphErrors> width_graph[component_size];
    std::unique_ptr<TGraphErrors> correlation_graph[component_size];
    std::unique_ptr<TH2D> amplitude_hist[component_size];
    std::unique_ptr<TH2D> width_hist[component_size];
    std::vector<double> amplitude_array, width_array;
    amplitude_array.reserve(80);
    width_array.reserve(80);
    auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    for (size_t i = 0; i < component_size; i++)
    {
        amplitude_graph[i] = entry_iter->CreateAtomGausEstimateToSpotGraph(group_key_list[i], class_key, 0);
        width_graph[i] = entry_iter->CreateAtomGausEstimateToSpotGraph(group_key_list[i], class_key, 1);
        correlation_graph[i] = entry_iter->CreateAtomGausEstimateScatterGraph(group_key_list[i], class_key, 0, 1);
        for (int p = 0; p < amplitude_graph[i]->GetN(); p++)
        {
            amplitude_array.push_back(amplitude_graph[i]->GetPointY(p));
            width_array.push_back(width_graph[i]->GetPointY(p));
        }
        auto component_color{ AtomClassifier::GetMainChainElementColor(i) }; // TODO: nucleotide color
        auto component_marker{ AtomClassifier::GetMainChainElementSolidMarker(i) }; // TODO: nucleotide marker
        ROOTHelper::SetMarkerAttribute(amplitude_graph[i].get(), component_marker, 1.3f, component_color);
        ROOTHelper::SetMarkerAttribute(width_graph[i].get(), component_marker, 1.3f, component_color);
        ROOTHelper::SetMarkerAttribute(correlation_graph[i].get(), component_marker, 1.3f, component_color);
        ROOTHelper::SetLineAttribute(amplitude_graph[i].get(), 1, 1, component_color);
        ROOTHelper::SetLineAttribute(width_graph[i].get(), 1, 1, component_color);
        ROOTHelper::SetLineAttribute(correlation_graph[i].get(), 1, 1, component_color);
    }

    auto scaling{ 0.3 };
    auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, scaling) };
    auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, scaling) };
    for (int i = 0; i < component_size; i++)
    {
        std::string name_amplitude{ "amplitude_hist_"+ std::to_string(i) };
        std::string name_width{ "width_hist_"+ std::to_string(i) };
        amplitude_hist[i] = ROOTHelper::CreateHist2D(name_amplitude.data(),"", 4, -0.5, 3.5, 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        width_hist[i] = ROOTHelper::CreateHist2D(name_width.data(),"", 4, -0.5, 3.5, 100, std::get<0>(width_range), std::get<1>(width_range));
        for (int p = 0; p < amplitude_graph[i]->GetN(); p++)
        {
            amplitude_hist[i]->Fill(i, amplitude_graph[i]->GetPointY(p));
        }
        for (int p = 0; p < width_graph[i]->GetN(); p++)
        {
            width_hist[i]->Fill(i, width_graph[i]->GetPointY(p));
        }
        auto element_color{ AtomClassifier::GetMainChainElementColor(static_cast<size_t>(i)) };
        ROOTHelper::SetLineAttribute(amplitude_hist[i].get(), 1, 1, element_color);
        ROOTHelper::SetLineAttribute(width_hist[i].get(), 1, 1, element_color);
        ROOTHelper::SetFillAttribute(amplitude_hist[i].get(), 1001, element_color, 0.3f);
        ROOTHelper::SetFillAttribute(width_hist[i].get(), 1001, element_color, 0.3f);
        amplitude_hist[i]->SetBarWidth(0.5f);
        width_hist[i]->SetBarWidth(0.5f);
    }

    //auto count_hist{ entry_iter->CreateAtomResidueCountHistogram(class_key, Structure::FREE) };
    //auto max_count{ static_cast<int>(count_hist->GetMaximum()) };

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
    //frame[5]->GetYaxis()->SetLimits(0.0, max_count*1.1);

    pad[1]->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.07, 0.005, 0.01, 0.01);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(frame[1]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTitleAttribute(frame[1]->GetYaxis(), 35.0f, 1.4f);
    ROOTHelper::SetAxisLabelAttribute(frame[1]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(frame[1]->GetYaxis(), 30.0f, 0.01f);
    auto x_tick_length_1{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.0, 0) };
    auto y_tick_length_1{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.015, 1) };
    ROOTHelper::SetAxisTickAttribute(frame[1]->GetXaxis(), static_cast<float>(x_tick_length_1), 13);
    ROOTHelper::SetAxisTickAttribute(frame[1]->GetYaxis(), static_cast<float>(y_tick_length_1), 506);
    frame[1]->GetXaxis()->SetLimits(-1.0, 12.0);
    frame[1]->SetStats(0);
    frame[1]->GetYaxis()->SetTitle("Amplitude");
    frame[1]->GetYaxis()->CenterTitle();
    frame[1]->Draw();
    for (int i = 0; i < component_size; i++) amplitude_graph[i]->Draw("P X0");

    pad[0]->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.07, 0.005, 0.11, 0.01);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(frame[0]->GetXaxis(), 0.0);
    ROOTHelper::SetAxisTitleAttribute(frame[0]->GetYaxis(), 35.0, 1.4f);
    ROOTHelper::SetAxisLabelAttribute(frame[0]->GetXaxis(), 32.0, 0.11f, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(frame[0]->GetYaxis(), 30.0, 0.01f);
    auto x_tick_length_2{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.0, 0) };
    auto y_tick_length_2{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.015, 1) };
    ROOTHelper::SetAxisTickAttribute(frame[0]->GetXaxis(), static_cast<float>(x_tick_length_2), 13);
    ROOTHelper::SetAxisTickAttribute(frame[0]->GetYaxis(), static_cast<float>(y_tick_length_2), 505);
    frame[0]->GetXaxis()->SetLimits(-1.0, 12.0);
    frame[0]->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    frame[0]->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);
    for (size_t i = 0; i < group_key_list[0].size(); i++)
    {
        auto spot{ AtomClassifier::GetNucleotideMainChainSpot(i) };
        auto label{ ChemicalDataHelper::GetLabel(spot) };
        auto label_index{ static_cast<int>(i) + 2 };
        frame[0]->GetXaxis()->ChangeLabel(label_index, 90.0, -1, 12, -1, -1, label.data());
    }
    frame[0]->SetStats(0);
    frame[0]->GetYaxis()->SetTitle("Width");
    frame[0]->GetYaxis()->CenterTitle();
    frame[0]->Draw();
    for (int i = 0; i < component_size; i++) width_graph[i]->Draw("P X0");

    pad[3]->cd();
    //PrintAmplitudeSummaryPad(pad[3].get(), frame[3].get());
    for (int i = 0; i < component_size; i++) amplitude_hist[i]->Draw("CANDLE3 SAME");

    pad[2]->cd();
    //PrintAtomWidthSummaryPad(pad[2].get(), frame[2].get());
    for (int i = 0; i < component_size; i++) width_hist[i]->Draw("CANDLE3 SAME");

    pad[4]->cd();
    //PrintGausSummaryPad(pad[4].get(), frame[4].get());
    for (int i = 0; i < component_size; i++) correlation_graph[i]->Draw("P X0");

    pad[5]->cd();
    auto info_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    auto resolution_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextMarginInCanvas(gPad, info_text.get(), 0.005, 0.24, 0.01, 0.02);
    ROOTHelper::SetPaveTextDefaultStyle(info_text.get());
    ROOTHelper::SetPaveAttribute(info_text.get(), 0, 0.1);
    ROOTHelper::SetFillAttribute(info_text.get(), 1001, kAzure-7, 0.5);
    ROOTHelper::SetTextAttribute(info_text.get(), 55, 133, 12);
    info_text->AddText(("#font[102]{PDB-" + model_object->GetPdbID() +"}").data());
    info_text->AddText(("#font[102]{"+ model_object->GetEmdID() +"}").data());
    info_text->Draw();

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
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.07, 0.005, 0.01, 0.02);
    ROOTHelper::SetPadLayout(gPad, 1, 0, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(frame[5]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTitleAttribute(frame[5]->GetYaxis(), 32.0f, 1.2f);
    ROOTHelper::SetAxisLabelAttribute(frame[5]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(frame[5]->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisTickAttribute(frame[5]->GetXaxis(), 0.0f, 21);
    ROOTHelper::SetAxisTickAttribute(frame[5]->GetYaxis(), 0.0f);
    frame[5]->GetXaxis()->SetLimits(-1.0, 20.0);
    frame[5]->SetMinimum(0.5);
    frame[5]->SetStats(0);
    frame[5]->GetYaxis()->SetTitle("#splitline{Residue}{Counts}");
    frame[5]->GetYaxis()->CenterTitle();
    frame[5]->Draw();
    //gStyle->SetTextFont(132);
    //ROOTHelper::SetFillAttribute(count_hist.get(), 1001, kAzure-7, 0.5f);
    //ROOTHelper::SetLineAttribute(count_hist.get(), 1, 1, kAzure-7);
    //ROOTHelper::SetMarkerAttribute(count_hist.get(), 20, 7.0f, kAzure);
    //count_hist->Draw("HIST TEXT0 SAME");

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);

    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}
