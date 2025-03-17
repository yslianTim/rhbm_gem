#include "ModelPainter.hpp"
#include "ModelObject.hpp"
#include "GroupPotentialEntry.hpp"
#include "FilePathHelper.hpp"
#include "AtomicInfoHelper.hpp"
#include "ArrayStats.hpp"

#ifdef HAVE_ROOT
#include "ROOTHelper.hpp"
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TColor.h>
#include <TMarker.h>
#include <TAxis.h>
#include <TH2.h>
#endif

#include <vector>
#include <tuple>

using ElementKeyType = GroupKeyMapping<ElementGroupClassifierTag>::type;
using ResidueKeyType = GroupKeyMapping<ResidueGroupClassifierTag>::type;

ModelPainter::ModelPainter(ModelObject * model) :
    m_model_object{ model }, m_folder_path{ "./" }
{

}

ModelPainter::~ModelPainter()
{

}

void ModelPainter::SetFolder(const std::string & folder_path)
{
    m_folder_path = FilePathHelper::EnsureTrailingSlash(folder_path);
}

void ModelPainter::Painting(void)
{
    std::cout <<"- ModelPainter::Painting"<<std::endl;
    std::cout <<"  Folder path: "<< m_folder_path << std::endl;
    PaintResidueClassGroupGaus("test1.pdf");
}

void ModelPainter::PaintResidueClassGroupGaus(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintResidueGroupGaus"<<std::endl;

    auto group_entry{ m_model_object->GetGroupPotentialEntry(AtomicInfoHelper::GetResidueClassKey()) };
    auto residue_class_group_entry{ dynamic_cast<GroupPotentialEntry<ResidueKeyType> *>(group_entry) };
    
    #ifdef HAVE_ROOT

    const int primary_element_size{ 4 };
    const char * element_label[primary_element_size]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };
    int element_index[primary_element_size]   { 6, 6, 7, 8 };
    int remoteness_index[primary_element_size]{ 1, 0, 0, 0 };
    int color_element[primary_element_size]   { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    int marker_element[primary_element_size]  { kFullSquare, kFullCircle, kFullTriangleUp, kFullTriangleDown };

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1200, 600) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size{ 6 };
    std::unique_ptr<TPad> pad[pad_size];
    std::unique_ptr<TH2> frame[pad_size];
    pad[0] = ROOTHelper::CreatePad("pad1","", 0.00, 0.00, 0.52, 0.55);
    pad[1] = ROOTHelper::CreatePad("pad2","", 0.00, 0.55, 0.52, 1.00);
    pad[2] = ROOTHelper::CreatePad("pad3","", 0.52, 0.00, 0.65, 0.55);
    pad[3] = ROOTHelper::CreatePad("pad4","", 0.52, 0.55, 0.65, 1.00);
    pad[4] = ROOTHelper::CreatePad("pad5","", 0.65, 0.00, 1.00, 0.80);
    pad[5] = ROOTHelper::CreatePad("pad6","", 0.65, 0.80, 1.00, 1.00);

    std::unique_ptr<TGraphErrors> amplitude_graph[primary_element_size];
    std::unique_ptr<TGraphErrors> width_graph[primary_element_size];
    std::unique_ptr<TGraphErrors> correlation_graph[primary_element_size];
    std::unique_ptr<TH2D> amplitude_hist[primary_element_size];
    std::unique_ptr<TH2D> width_hist[primary_element_size];
    std::vector<double> amplitude_array, width_array;
    amplitude_array.reserve(primary_element_size * AtomicInfoHelper::GetStandardResidueCount());
    width_array.reserve(primary_element_size * AtomicInfoHelper::GetStandardResidueCount());
    for (int i = 0; i < primary_element_size; i++)
    {
        auto element{ element_index[i] };
        auto remoteness{ remoteness_index[i] };
        amplitude_graph[i] = ROOTHelper::CreateGraphErrors();
        width_graph[i] = ROOTHelper::CreateGraphErrors();
        correlation_graph[i] = ROOTHelper::CreateGraphErrors();
        for (int residue = 0; residue < AtomicInfoHelper::GetStandardResidueCount(); residue++)
        {
            auto group_key{ std::make_tuple(residue, element, remoteness, 0, false) };
            auto gaus_estimate{ residue_class_group_entry->GetGausEstimatePrior(&group_key) };
            auto gaus_variance{ residue_class_group_entry->GetGausVariancePrior(&group_key) };
            amplitude_array.push_back(std::get<0>(gaus_estimate));
            width_array.push_back(std::get<1>(gaus_estimate));
            amplitude_graph[i]->SetPoint(residue, residue, std::get<0>(gaus_estimate));
            amplitude_graph[i]->SetPointError(residue, 0.0, std::get<0>(gaus_variance));
            width_graph[i]->SetPoint(residue, residue, std::get<1>(gaus_estimate));
            width_graph[i]->SetPointError(residue, 0.0, std::get<1>(gaus_variance));
            correlation_graph[i]->SetPoint(residue, std::get<0>(gaus_estimate), std::get<1>(gaus_estimate));
            correlation_graph[i]->SetPointError(residue, std::get<0>(gaus_variance), std::get<1>(gaus_variance));
        }
        ROOTHelper::SetMarkerAttribute(amplitude_graph[i].get(), marker_element[i], 1.3, color_element[i]);
        ROOTHelper::SetMarkerAttribute(width_graph[i].get(), marker_element[i], 1.3, color_element[i]);
        ROOTHelper::SetMarkerAttribute(correlation_graph[i].get(), marker_element[i], 1.3, color_element[i]);
        ROOTHelper::SetLineAttribute(amplitude_graph[i].get(), 1, 1, color_element[i]);
        ROOTHelper::SetLineAttribute(width_graph[i].get(), 1, 1, color_element[i]);
        ROOTHelper::SetLineAttribute(correlation_graph[i].get(), 1, 1, color_element[i]);
    }

    auto scaling{ 0.3 };
    auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, scaling) };
    auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, scaling) };

    for (int i = 0; i < primary_element_size; i++)
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
        ROOTHelper::SetLineAttribute(amplitude_hist[i].get(), 1, 1, color_element[i]);
        ROOTHelper::SetLineAttribute(width_hist[i].get(), 1, 1, color_element[i]);
        ROOTHelper::SetFillAttribute(amplitude_hist[i].get(), 1001, color_element[i], 0.3);
        ROOTHelper::SetFillAttribute(width_hist[i].get(), 1001, color_element[i], 0.3);

        amplitude_hist[i]->SetBarWidth(0.6);
        width_hist[i]->SetBarWidth(0.6);
    }

    canvas->cd();
    for (int i = 0; i < pad_size; i++)
    {
        ROOTHelper::SetPadDefaultStyle(pad[i].get());
        pad[i]->Draw();
    }

    frame[1] = ROOTHelper::CreateHist2D("hist_1","", 100, -1.0, 20.0, 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    frame[0] = ROOTHelper::CreateHist2D("hist_0","", 100, -1.0, 20.0, 100, std::get<0>(width_range), std::get<1>(width_range));
    frame[3] = ROOTHelper::CreateHist2D("hist_3","", 100, 0.0, 1.0, 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    frame[2] = ROOTHelper::CreateHist2D("hist_2","", 100, 0.0, 1.0, 100, std::get<0>(width_range), std::get<1>(width_range));
    frame[4] = ROOTHelper::CreateHist2D("hist_4","", 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range), 100, std::get<0>(width_range), std::get<1>(width_range));
    auto data_info_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };

    PrintAmplitudePad(pad[1].get(), frame[1].get());
    PrintWidthPad(pad[0].get(), frame[0].get());
    PrintAmplitudeSummaryPad(pad[3].get(), frame[3].get());
    PrintWidthSummaryPad(pad[2].get(), frame[2].get());
    PrintGausSummaryPad(pad[4].get(), frame[4].get());
    PrintTitlePad(pad[5].get(), data_info_text.get(), m_model_object->GetPdbID() , m_model_object->GetEmdID());

    pad[1]->cd();
    for (int i = 0; i < primary_element_size; i++) amplitude_graph[i]->Draw("PL X0");

    pad[0]->cd();
    for (int i = 0; i < primary_element_size; i++) width_graph[i]->Draw("PL X0");

    pad[3]->cd();
    for (int i = 0; i < primary_element_size; i++) amplitude_hist[i]->Draw("CANDLE3 SAME");

    pad[2]->cd();
    for (int i = 0; i < primary_element_size; i++) width_hist[i]->Draw("CANDLE3 SAME");

    pad[4]->cd();
    for (int i = 0; i < primary_element_size; i++) correlation_graph[i]->Draw("P X0");

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

#ifdef HAVE_ROOT

void ModelPainter::PrintTitlePad(
    TPad * pad, TPaveText * text, const std::string & pdb_id, const std::string & emd_id)
{
    pad->cd();
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, 0.005, 0.005, 0.01, 0.02);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.1);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7, 0.5);
    ROOTHelper::SetTextAttribute(text, 50, 133, 12);
    text->AddText(("#font[102]{PDB-" + pdb_id +"}").data());
    text->AddText(("#font[102]{"+ emd_id +"}").data());
    pad->Update();
}

void ModelPainter::PrintAmplitudePad(TPad * pad, TH2 * hist)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.08, 0.005, 0.01, 0.02);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0, 1.4);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0, 0.01);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 15, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 15, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, 21);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 506);
    hist->SetStats(0);
    hist->GetYaxis()->SetTitle("Amplitude");
    hist->GetYaxis()->CenterTitle();
    hist->Draw();
}

void ModelPainter::PrintWidthPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.08, 0.005, 0.15, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0, 1.4);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 30.0, -0.05, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0, 0.01);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 15, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 15, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, 21);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 505);
    hist->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    hist->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);
    for (int i = 0; i < AtomicInfoHelper::GetStandardResidueCount(); i++)
    {
        auto label{ AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::ResidueTag>(i) };
        hist->GetXaxis()->ChangeLabel(i+2, 90.0, -1, 32, -1, -1, label.data());
    }

    hist->SetStats(0);
    hist->GetYaxis()->SetTitle("Width");
    hist->GetYaxis()->CenterTitle();
    hist->Draw();
}

void ModelPainter::PrintAmplitudeSummaryPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.005, 0.01, 0.02);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 0.0);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 15, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 15, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, 5);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 506);
    hist->GetXaxis()->SetLimits(-1.0, 4.0);
    hist->SetStats(0);
    hist->Draw();
}

void ModelPainter::PrintWidthSummaryPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.005, 0.15, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 35.0, 1.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 30.0);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 0.0);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 15, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 15, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, 5);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 505);
    hist->GetXaxis()->SetLimits(-1.0, 4.0);
    hist->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    hist->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);

    std::cout <<"N labels = "<< hist->GetXaxis()->GetNlabels() << std::endl;
    const char * label_name[4]{ "C_{#alpha}", "C", "N", "O" };
    for (int i = 0; i < 4; i++)
    {
        hist->GetXaxis()->ChangeLabel(i+2, 0.0, -1, -1, -1, -1, label_name[i]);
    }
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("Element");
    hist->GetXaxis()->CenterTitle();
    hist->Draw();
}

void ModelPainter::PrintGausSummaryPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.07, 0.15, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 35.0, 1.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0, 1.2);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 30.0);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0, 0.01);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 15, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 15, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, 505);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 505);
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("Amplitude");
    hist->GetYaxis()->SetTitle("Width");
    hist->Draw("Y+");
}

#endif