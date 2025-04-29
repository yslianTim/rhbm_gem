#include "DemoPainter.hpp"
#include "ModelObject.hpp"
#include "AtomObject.hpp"
#include "DataObjectBase.hpp"
#include "PotentialEntryIterator.hpp"
#include "FilePathHelper.hpp"
#include "AtomicInfoHelper.hpp"
#include "AtomClassifier.hpp"
#include "ArrayStats.hpp"
#include "KeyPacker.hpp"

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
#endif

#include <iostream>
#include <vector>
#include <tuple>

DemoPainter::DemoPainter(void) :
    m_folder_path{ "./" }, m_atom_classifier{ std::make_unique<AtomClassifier>() }
{

}

DemoPainter::~DemoPainter()
{

}

void DemoPainter::SetFolder(const std::string & folder_path)
{
    m_folder_path = FilePathHelper::EnsureTrailingSlash(folder_path);
}

void DemoPainter::AddDataObject(DataObjectBase * data_object)
{
    m_model_object_list.push_back(dynamic_cast<ModelObject *>(data_object));
}

void DemoPainter::AddReferenceDataObject(DataObjectBase * data_object, const std::string & label)
{
    m_ref_model_object_map[label].push_back(dynamic_cast<ModelObject *>(data_object));
}

void DemoPainter::Painting(void)
{
    std::cout <<"- DemoPainter::Painting"<<std::endl;
    std::cout <<"  Folder path: "<< m_folder_path << std::endl;

    if (m_model_object_list.size() == 4)
    {
        PaintResidueClassGroupGausMainChainSummary("figure_1_c.pdf");
        PaintWidthToBfactorScatterPlotSummary("test_width_bfactor.pdf");
    }

    if (m_model_object_list.size() > 5)
    {
        PaintElementClassGroupGausToFSC("figure_5.pdf");
    }
}

void DemoPainter::PaintResidueClassGroupGausMainChainSummary(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    auto residue_class{ AtomicInfoHelper::GetResidueClassKey() };
    std::cout <<"- DemoPainter::PaintResidueClassGroupGausMainChainSummary"<< std::endl;

    #ifdef HAVE_ROOT

    const int primary_element_size{ 4 };
    short color_element[primary_element_size] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    short marker_element[primary_element_size]{ 25, 24, 26, 32 };
    float marker_size{ 1.5f };

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 3000, 1800) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size_x{ 6 };
    const int pad_size_y{ 4 };
    std::unique_ptr<TPad> pad[pad_size_x][pad_size_y];
    std::unique_ptr<TH2> frame[pad_size_x][pad_size_y];
    
    const double x_pos[pad_size_x + 1]{ 0.00, 0.10, 0.40, 0.47, 0.77, 0.84, 1.00 };
    const double y_pos[pad_size_y + 1]{ 0.00, 0.28, 0.50, 0.72, 1.00 };
    for (int i = 0; i < pad_size_x; i++)
    {
        for (int j = 0; j < pad_size_y; j++)
        {
            auto pad_name{ "pad"+ std::to_string(i) + std::to_string(j) };
            pad[i][j] = ROOTHelper::CreatePad(pad_name.data(),"", x_pos[i], y_pos[j], x_pos[i+1], y_pos[j+1]);
        }
    }

    std::unique_ptr<TPaveText> info_text[pad_size_y];
    std::unique_ptr<TPaveText> icon_text[pad_size_y];
    std::unique_ptr<TGraphErrors> amplitude_graph[pad_size_y][primary_element_size];
    std::unique_ptr<TGraphErrors> width_graph[pad_size_y][primary_element_size];
    std::unique_ptr<TGraphErrors> correlation_graph[pad_size_y][primary_element_size];
    std::unique_ptr<TH2D> amplitude_hist[pad_size_y][primary_element_size];
    std::unique_ptr<TH2D> width_hist[pad_size_y][primary_element_size];
    std::vector<double> width_array_total;
    for (size_t j = 0; j < pad_size_y; j++)
    {
        std::vector<double> amplitude_array, width_array;
        info_text[j] = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
        icon_text[j] = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(m_model_object_list.at(static_cast<size_t>(j))) };
        for (size_t k = 0; k < primary_element_size; k++)
        {
            auto group_key_list{ m_atom_classifier->GetMainChainResidueClassGroupKeyList(k) };
            amplitude_graph[j][k] = entry_iter->CreateGausEstimateToResidueGraph(group_key_list, residue_class, 0);
            width_graph[j][k] = entry_iter->CreateGausEstimateToResidueGraph(group_key_list, residue_class, 1);
            correlation_graph[j][k] = entry_iter->CreateGausEstimateScatterGraph(group_key_list, residue_class, true);
            for (int p = 0; p < amplitude_graph[j][k]->GetN(); p++)
            {
                amplitude_array.push_back(amplitude_graph[j][k]->GetPointY(p));
                width_array.push_back(width_graph[j][k]->GetPointY(p));
                width_array_total.push_back(width_graph[j][k]->GetPointY(p));
            }
            ROOTHelper::SetMarkerAttribute(amplitude_graph[j][k].get(), marker_element[k], marker_size, color_element[k]);
            ROOTHelper::SetMarkerAttribute(width_graph[j][k].get(), marker_element[k], marker_size, color_element[k]);
            ROOTHelper::SetMarkerAttribute(correlation_graph[j][k].get(), marker_element[k], marker_size, color_element[k]);
            ROOTHelper::SetLineAttribute(amplitude_graph[j][k].get(), 1, 1, color_element[k]);
            ROOTHelper::SetLineAttribute(width_graph[j][k].get(), 1, 1, color_element[k]);
            ROOTHelper::SetLineAttribute(correlation_graph[j][k].get(), 1, 1, color_element[k]);
        }
        auto scaling{ 0.3 };
        auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, scaling) };
        auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, scaling) };

        for (int k = 0; k < primary_element_size; k++)
        {
            std::string name_amplitude{ "amplitude_hist_"+ std::to_string(j) + std::to_string(k) };
            std::string name_width{ "width_hist_"+ std::to_string(j) + std::to_string(k) };
            amplitude_hist[j][k] = ROOTHelper::CreateHist2D(name_amplitude.data(),"", 4, -0.5, 3.5, 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
            width_hist[j][k] = ROOTHelper::CreateHist2D(name_width.data(),"", 4, -0.5, 3.5, 100, std::get<0>(width_range), std::get<1>(width_range));
            for (int p = 0; p < amplitude_graph[j][k]->GetN(); p++)
            {
                amplitude_hist[j][k]->Fill(k, amplitude_graph[j][k]->GetPointY(p));
            }
            for (int p = 0; p < width_graph[j][k]->GetN(); p++)
            {
                width_hist[j][k]->Fill(k, width_graph[j][k]->GetPointY(p));
            }
            ROOTHelper::SetLineAttribute(amplitude_hist[j][k].get(), 1, 1, color_element[k]);
            ROOTHelper::SetLineAttribute(width_hist[j][k].get(), 1, 1, color_element[k]);
            ROOTHelper::SetFillAttribute(amplitude_hist[j][k].get(), 1001, color_element[k], 0.3f);
            ROOTHelper::SetFillAttribute(width_hist[j][k].get(), 1001, color_element[k], 0.3f);

            amplitude_hist[j][k]->SetBarWidth(0.6f);
            width_hist[j][k]->SetBarWidth(0.6f);
        }
        frame[0][j] = ROOTHelper::CreateHist2D(("frame0"+std::to_string(j)).data(),"", 100, 0.0, 1.0, 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        frame[1][j] = ROOTHelper::CreateHist2D(("frame1"+std::to_string(j)).data(),"", 100, 0.0, 1.0, 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        frame[2][j] = ROOTHelper::CreateHist2D(("frame2"+std::to_string(j)).data(),"", 100, 0.0, 1.0, 100, std::get<0>(width_range), std::get<1>(width_range));
        frame[3][j] = ROOTHelper::CreateHist2D(("frame3"+std::to_string(j)).data(),"", 100, 0.0, 1.0, 100, std::get<0>(width_range), std::get<1>(width_range));
        frame[4][j] = ROOTHelper::CreateHist2D(("frame4"+std::to_string(j)).data(),"", 100, std::get<0>(width_range), std::get<1>(width_range), 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    }

    auto amplitude_title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    auto width_title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    auto correlation_title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };

    auto width_range_total{ ArrayStats<double>::ComputeScalingRangeTuple(width_array_total, 0.3) };
    for (int j = 0; j < pad_size_y; j++)
    {
        frame[4][j]->GetXaxis()->SetLimits(std::get<0>(width_range_total), std::get<1>(width_range_total));
    }

    canvas->cd();
    for (int i = 0; i < pad_size_x; i++)
    {
        for (int j = 0; j < pad_size_y; j++)
        {
            ROOTHelper::SetPadDefaultStyle(pad[i][j].get());
            pad[i][j]->Draw();
        }
    }

    for (int j = 0; j < pad_size_y; j++)
    {
        auto draw_axis_flag{ (j == 0) ? true : false };
        auto draw_title_flag{ (j == pad_size_y - 1) ? true : false };

        pad[0][j]->cd();
        auto pdb_text{ m_model_object_list.at(static_cast<size_t>(j))->GetPdbID() };
        auto emd_text{ m_model_object_list.at(static_cast<size_t>(j))->GetEmdID() };
        auto resolution{ m_model_object_list.at(static_cast<size_t>(j))->GetResolution() };
        PrintInfoMainChainPad(pad[0][j].get(), info_text[j].get(), pdb_text, emd_text, draw_axis_flag, draw_title_flag);
        PrintIconMainChainPad(pad[0][j].get(), icon_text[j].get(), resolution, draw_axis_flag, draw_title_flag);

        pad[1][j]->cd();
        PrintResultPad(pad[1][j].get(), frame[0][j].get(), draw_axis_flag, draw_title_flag);
        for (int k = 0; k < primary_element_size; k++) amplitude_graph[j][k]->Draw("PL X0");
        if (j == pad_size_y - 1) PrintTitlePad(pad[1][j].get(), amplitude_title_text.get(), "Amplitude");

        pad[2][j]->cd();
        PrintSummaryPad(pad[2][j].get(), frame[1][j].get(), draw_axis_flag, draw_title_flag);
        for (int k = 0; k < primary_element_size; k++) amplitude_hist[j][k]->Draw("CANDLE3 SAME");

        pad[3][j]->cd();
        PrintResultPad(pad[3][j].get(), frame[2][j].get(), draw_axis_flag, draw_title_flag);
        for (int k = 0; k < primary_element_size; k++) width_graph[j][k]->Draw("PL X0");
        if (j == pad_size_y - 1) PrintTitlePad(pad[3][j].get(), width_title_text.get(), "Width");

        pad[4][j]->cd();
        PrintSummaryPad(pad[4][j].get(), frame[3][j].get(), draw_axis_flag, draw_title_flag);
        for (int k = 0; k < primary_element_size; k++) width_hist[j][k]->Draw("CANDLE3 SAME");

        pad[5][j]->cd();
        PrintCorrelationPad(pad[5][j].get(), frame[4][j].get(), draw_axis_flag, draw_title_flag);
        for (int k = 0; k < primary_element_size; k++) correlation_graph[j][k]->Draw("P X0");
        if (j == pad_size_y - 1) PrintTitlePad(pad[5][j].get(), correlation_title_text.get(), "Scatter Plot");
    }

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void DemoPainter::PaintResidueClassGroupGausSideChainSummary(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    auto residue_class{ AtomicInfoHelper::GetResidueClassKey() };
    std::cout <<"- DemoPainter::PaintResidueClassGroupGausSideChainSummary"<< std::endl;

    #ifdef HAVE_ROOT

    float marker_size{ 1.5f };

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 3000, 2000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_extra_size{ 4 };
    std::unique_ptr<TPad> pad_extra[pad_extra_size];
    std::unique_ptr<TH2> frame_extra[pad_extra_size];

    const double x_extra_pos[pad_extra_size + 1]{ 0.00, 0.25, 0.50, 0.75, 1.00 };
    for (int i = 0; i < pad_extra_size; i++)
    {
        auto pad_name{ "pad_extra"+ std::to_string(i) };
        auto frame_name{ "frame_extra"+ std::to_string(i) };
        pad_extra[i] = ROOTHelper::CreatePad(pad_name.data(),"", x_extra_pos[i], 0.00, x_extra_pos[i+1], 0.20);
        frame_extra[i] = ROOTHelper::CreateHist2D(frame_name.data(),"", 100, 0.0, 1.0, 100, 0.0, 1.0);
    }

    std::unique_ptr<TPaveText> side_chain_title_text[pad_extra_size];
    std::vector<std::unique_ptr<TGraphErrors>> amplitude_graph_list[2];
    std::vector<std::unique_ptr<TGraphErrors>> width_graph_list[2];
    std::vector<std::string> label_list[2];
    std::vector<double> amplitude_array, width_array;
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(m_model_object_list.at(3)) };
    Residue residue_list[2]{ Residue::ASP, Residue::GLU };
    for (int i = 0; i < 2; i++)
    {
        auto residue{ residue_list[i] };
        auto count_total{ 0 };
        for (auto element : AtomicInfoHelper::GetStandardElementList())
        {
            auto amplitude_graph_tmp{ ROOTHelper::CreateGraphErrors() };
            auto width_graph_tmp{ ROOTHelper::CreateGraphErrors() };
            auto count_element{ 0 };
            for (auto remoteness : AtomicInfoHelper::GetStandardRemotenessList())
            {
                for (auto branch : AtomicInfoHelper::GetStandardBranchList())
                {
                    auto group_key{ KeyPackerResidueClass::Pack(residue, element, remoteness, branch, false) };
                    if (entry_iter->IsAvailableGroupKey(group_key, residue_class) == false) continue;
                    auto gaus_estimate{ entry_iter->GetGausEstimatePrior(group_key, residue_class) };
                    auto gaus_variance{ entry_iter->GetGausVariancePrior(group_key, residue_class) };
                    amplitude_array.push_back(std::get<0>(gaus_estimate));
                    width_array.push_back(std::get<1>(gaus_estimate));
                    amplitude_graph_tmp->SetPoint(count_element, count_total, std::get<0>(gaus_estimate));
                    amplitude_graph_tmp->SetPointError(count_element, 0.0, std::get<0>(gaus_variance));
                    width_graph_tmp->SetPoint(count_element, count_total, std::get<1>(gaus_estimate));
                    width_graph_tmp->SetPointError(count_element, 0.0, std::get<1>(gaus_variance));
                    auto element_label{ AtomicInfoHelper::GetLabel(element) };
                    auto remoteness_label{ AtomicInfoHelper::GetLabel(remoteness) };
                    auto branch_label{ AtomicInfoHelper::GetLabel(branch) };
                    label_list[i].emplace_back(element_label +"_{"+ remoteness_label + branch_label+"}");
                    count_element++;
                    count_total++;
                }
            }
            auto color{ AtomicInfoHelper::GetDisplayColor(element) };
            ROOTHelper::SetMarkerAttribute(amplitude_graph_tmp.get(), kFullCircle, marker_size, color);
            ROOTHelper::SetMarkerAttribute(width_graph_tmp.get(), kFullCircle, marker_size, color);
            amplitude_graph_list[i].push_back(std::move(amplitude_graph_tmp));
            width_graph_list[i].push_back(std::move(width_graph_tmp));
        }
    }
    auto scaling{ 0.3 };
    auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, scaling) };
    auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, scaling) };
    frame_extra[0]->GetYaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    frame_extra[1]->GetYaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    frame_extra[2]->GetYaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));
    frame_extra[3]->GetYaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));

    for (int i = 0; i < pad_extra_size; i++)
    {
        side_chain_title_text[i] = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
        ROOTHelper::SetPadDefaultStyle(pad_extra[i].get());
        pad_extra[i]->Draw();
    }

    PrintLeftSideChainPad(pad_extra[0].get(), frame_extra[0].get(), residue_list[0],"Amplitude", label_list[0]);
    PrintRightSideChainPad(pad_extra[1].get(), frame_extra[1].get(), residue_list[1], label_list[1]);
    PrintLeftSideChainPad(pad_extra[2].get(), frame_extra[2].get(), residue_list[0],"Width", label_list[0]);
    PrintRightSideChainPad(pad_extra[3].get(), frame_extra[3].get(), residue_list[1], label_list[1]);

    pad_extra[0]->cd();
    for (auto & graph : amplitude_graph_list[0]) graph->Draw("P X0");
    PrintTitleSideChainPad(pad_extra[0].get(), side_chain_title_text[0].get(), AtomicInfoHelper::GetLabel(residue_list[0]));

    pad_extra[1]->cd();
    for (auto & graph : amplitude_graph_list[1]) graph->Draw("P X0");
    PrintTitleSideChainPad(pad_extra[1].get(), side_chain_title_text[1].get(), AtomicInfoHelper::GetLabel(residue_list[1]));

    pad_extra[2]->cd();
    for (auto & graph : width_graph_list[0]) graph->Draw("P X0");
    PrintTitleSideChainPad(pad_extra[2].get(), side_chain_title_text[2].get(), AtomicInfoHelper::GetLabel(residue_list[0]));

    pad_extra[3]->cd();
    for (auto & graph : width_graph_list[1]) graph->Draw("P X0");
    PrintTitleSideChainPad(pad_extra[3].get(), side_chain_title_text[3].get(), AtomicInfoHelper::GetLabel(residue_list[1]));

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void DemoPainter::PaintElementClassGroupGausToFSC(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- DemoPainter::PaintElementClassGroupGausToFSC"<< std::endl;

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int column_size{ 4 };
    const int row_size{ 1 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 400) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), column_size, row_size, 0.08f, 0.01f, 0.17f, 0.12f, 0.01f, 0.005f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    
    const int primary_element_size{ 4 };
    const char * element_label[primary_element_size]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };
    short color_element[primary_element_size] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };

    std::unique_ptr<TGraphErrors> graph[primary_element_size];
    std::unique_ptr<TF1> fit_function[column_size];
    std::vector<double> x_array, y_array;
    double r_square[column_size]{ 0.0 };
    double slope[column_size]{ 0.0 };
    double intercept[column_size]{ 0.0 };
    for (size_t i = 0; i < primary_element_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        graph[i] = ROOTHelper::CreateGraphErrors();
        auto count{ 0 };
        for (auto model : m_model_object_list)
        {
            auto entry_iter{ std::make_unique<PotentialEntryIterator>(model) };
            auto gaus_estimate{ entry_iter->GetGausEstimatePrior(group_key, AtomicInfoHelper::GetElementClassKey()) };
            auto gaus_variance{ entry_iter->GetGausVariancePrior(group_key, AtomicInfoHelper::GetElementClassKey()) };
            graph[i]->SetPoint(count, model->GetResolution(), std::get<1>(gaus_estimate));
            graph[i]->SetPointError(count, 0.0, std::get<1>(gaus_variance));
            count++;
        }
        for (int p = 0; p < graph[i]->GetN(); p++)
        {
            x_array.push_back(graph[i]->GetPointX(p));
            y_array.push_back(graph[i]->GetPointY(p));
        }
        ROOTHelper::SetMarkerAttribute(graph[i].get(), 20, 1.2f, color_element[i]);
        ROOTHelper::SetLineAttribute(graph[i].get(), 1, 2, color_element[i]);

        r_square[i] = ROOTHelper::PerformLinearRegression(graph[i].get(), slope[i], intercept[i]);
        fit_function[i] = ROOTHelper::CreateFunction1D(Form("fit_%d", static_cast<int>(i)), "x*[1]+[0]");
        fit_function[i]->SetParameters(intercept[i], slope[i]);
        ROOTHelper::SetLineAttribute(fit_function[i].get(), 2, 2, kRed);
    }

    auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.2) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.3) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };

    std::unique_ptr<TH2> frame[column_size][row_size];
    std::unique_ptr<TPaveText> title_text[column_size];
    std::unique_ptr<TPaveText> r_square_text[column_size];
    std::unique_ptr<TPaveText> fit_info_text[column_size];
    for (int i = 0; i < column_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min, x_max, 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 28, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 30, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 35, 1.2f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("FSC Resolution #[]{#AA}");
            frame[i][j]->GetYaxis()->SetTitle("Width");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            graph[i]->Draw("P");
            fit_function[i]->SetRange(x_min, x_max);
            fit_function[i]->Draw("SAME");
        }
        title_text[i] = ROOTHelper::CreatePaveText(0.01, 1.01, 0.99, 1.16, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
        ROOTHelper::SetPaveAttribute(title_text[i].get(), 0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[i].get(), 30, 133, 22);
        ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, color_element[i], 0.5f);
        title_text[i]->AddText(element_label[i]);
        title_text[i]->Draw();

        r_square_text[i] = ROOTHelper::CreatePaveText(0.50, 0.05, 0.95, 0.18, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(r_square_text[i].get());
        ROOTHelper::SetPaveAttribute(r_square_text[i].get(), 0, 0.5);
        ROOTHelper::SetLineAttribute(r_square_text[i].get(), 1, 0);
        ROOTHelper::SetTextAttribute(r_square_text[i].get(), 25, 133, 22);
        ROOTHelper::SetFillAttribute(r_square_text[i].get(), 1001, kAzure-7, 0.20f);
        r_square_text[i]->AddText(Form("R^{2} = %.2f", r_square[i]));
        r_square_text[i]->Draw();

        fit_info_text[i] = ROOTHelper::CreatePaveText(0.12, 0.88, 0.70, 0.99, "nbNDC", true);
        ROOTHelper::SetPaveTextDefaultStyle(fit_info_text[i].get());
        ROOTHelper::SetTextAttribute(fit_info_text[i].get(), 25, 133, 22, 0.0, kGray);
        fit_info_text[i]->AddText(Form("#font[1]{y} = %.2f#font[1]{x}%+.2f", slope[i], intercept[i]));
        fit_info_text[i]->Draw();
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void DemoPainter::PaintWidthToBfactorScatterPlotSummary(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- DemoPainter::PaintWidthToBfactorScatterPlotSummary"<< std::endl;

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int column_size{ 4 };
    const int row_size{ 4 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 900) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), column_size, row_size, 0.08f, 0.01f, 0.08f, 0.05f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    
    const int primary_element_size{ 4 };
    const char * element_label[primary_element_size]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };
    short color_element[primary_element_size] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };

    std::unique_ptr<TGraphErrors> graph[primary_element_size][row_size];
    std::unique_ptr<TF1> fit_function[column_size][row_size];
    std::vector<double> x_array[column_size], y_array[row_size];
    double r_square[column_size][row_size]{ 0.0 };
    double slope[column_size][row_size]{ 0.0 };
    double intercept[column_size][row_size]{ 0.0 };
    for (size_t i = 0; i < primary_element_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        for (size_t j = 0; j < row_size; j++)
        {
            auto entry_iter{ std::make_unique<PotentialEntryIterator>(m_model_object_list.at(j)) };
            graph[i][j] = entry_iter->CreateBfactorToWidthScatterGraph(group_key, AtomicInfoHelper::GetElementClassKey());
            for (int p = 0; p < graph[i][j]->GetN(); p++)
            {
                x_array[i].push_back(graph[i][j]->GetPointX(p));
                y_array[j].push_back(graph[i][j]->GetPointY(p));
            }
            r_square[i][j] = ROOTHelper::PerformLinearRegression(graph[i][j].get(), slope[i][j], intercept[i][j]);
            fit_function[i][j] = ROOTHelper::CreateFunction1D(Form("fit_%d", static_cast<int>(i)), "x*[1]+[0]");
            fit_function[i][j]->SetParameters(intercept[i][j], slope[i][j]);
        }
    }

    double x_min[column_size]{ 0.0 };
    double x_max[column_size]{ 0.0 };
    double y_min[row_size]{ 0.0 };
    double y_max[row_size]{ 0.0 };
    for (size_t i = 0; i < column_size; i++)
    {
        auto x_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(x_array[i], 0.3) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);
    }
    for (size_t j = 0; j < row_size; j++)
    {
        auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array[j], 0.3) };
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }

    std::unique_ptr<TH2> frame[column_size][row_size];
    std::unique_ptr<TPaveText> title_text[column_size];
    std::unique_ptr<TPaveText> r_square_text[column_size][row_size];
    std::unique_ptr<TPaveText> fit_info_text[column_size][row_size];
    for (int i = 0; i < column_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min[j], y_max[j]);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 30, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 30, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 35, 1.2f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("B-factor #[]{#AA^{2}}");
            frame[i][j]->GetYaxis()->SetTitle("Width");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();

            ROOTHelper::SetMarkerAttribute(graph[i][j].get(), 24, 1.0f, color_element[i]);
            ROOTHelper::SetLineAttribute(graph[i][j].get(), 1, 2, color_element[i]);
            graph[i][j]->Draw("P");

            ROOTHelper::SetLineAttribute(fit_function[i][j].get(), 2, 2, kRed);
            fit_function[i][j]->SetRange(x_min[i], x_max[i]);
            fit_function[i][j]->Draw("SAME");

            r_square_text[i][j] = ROOTHelper::CreatePaveText(0.50, 0.05, 0.95, 0.21, "nbNDC ARC", true);
            ROOTHelper::SetPaveTextDefaultStyle(r_square_text[i][j].get());
            ROOTHelper::SetPaveAttribute(r_square_text[i][j].get(), 0, 0.5);
            ROOTHelper::SetLineAttribute(r_square_text[i][j].get(), 1, 0);
            ROOTHelper::SetTextAttribute(r_square_text[i][j].get(), 25, 133, 22);
            ROOTHelper::SetFillAttribute(r_square_text[i][j].get(), 1001, kAzure-7, 0.20f);
            r_square_text[i][j]->AddText(Form("R^{2} = %.2f", r_square[i][j]));
            r_square_text[i][j]->Draw();

            fit_info_text[i][j] = ROOTHelper::CreatePaveText(0.12, 0.88, 0.70, 0.99, "nbNDC", true);
            ROOTHelper::SetPaveTextDefaultStyle(fit_info_text[i][j].get());
            ROOTHelper::SetTextAttribute(fit_info_text[i][j].get(), 25, 133, 22, 0.0, kGray);
            fit_info_text[i][j]->AddText(Form("#font[1]{y} = %.2f#font[1]{x}%+.2f", slope[i][j], intercept[i][j]));
            fit_info_text[i][j]->Draw();
        }
        title_text[i] = ROOTHelper::CreatePaveText(0.01, 1.01, 0.99, 1.22, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
        ROOTHelper::SetPaveAttribute(title_text[i].get(), 0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[i].get(), 30, 133, 22);
        ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, color_element[i], 0.5f);
        title_text[i]->AddText(element_label[i]);
        title_text[i]->Draw();
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

#ifdef HAVE_ROOT

void DemoPainter::ModifyAxisLabelSideChain(
    TPad * pad, TH2 * hist, Residue residue, const std::vector<std::string> & label_list)
{
    if (AtomicInfoHelper::IsStandardResidue(residue) == false) return;

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    auto label_size{ static_cast<int>(label_list.size()) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), label_size+1);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 506);

    hist->GetXaxis()->SetLimits(-1.0, static_cast<double>(label_list.size()));
    hist->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    hist->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);

    for (int i = 0; i < label_size; i++)
    {
        hist->GetXaxis()->ChangeLabel(i+2, 0.0, -1, -1, -1, -1, label_list.at(static_cast<size_t>(i)).data());
    }
}

void DemoPainter::PrintIconMainChainPad(
    TPad * pad, TPaveText * text, double resolution, bool is_bottom_pad, bool is_top_pad)
{
    pad->cd();
    auto bottom_margin{ (is_bottom_pad) ? 0.20 : 0.14 };
    auto top_margin{ (is_top_pad) ? 0.06 : 0.01 };
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, 0.005, 0.005, bottom_margin, top_margin);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.2);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(text, 80, 133, 22, 0.0, kYellow-10);
    text->AddText(Form("%.2f #AA", resolution));
    pad->Update();
}

void DemoPainter::PrintInfoMainChainPad(
    TPad * pad, TPaveText * text,
    const std::string & pdb_id, const std::string & emd_id, bool is_bottom_pad, bool is_top_pad)
{
    pad->cd();
    auto bottom_margin{ (is_bottom_pad) ? 0.06 : 0.01 };
    auto top_margin{ (is_top_pad) ? 0.13 : 0.07 };
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, 0.005, 0.005, bottom_margin, top_margin);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.1);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7, 0.5);
    ROOTHelper::SetTextAttribute(text, 50, 133, 13);
    ROOTHelper::SetLineAttribute(text, 1, 0);
    text->AddText(("#font[102]{PDB-" + pdb_id +"}").data());
    text->AddText(("#font[102]{"+ emd_id +"}").data());
    pad->Update();
}

void DemoPainter::PrintTitlePad(TPad * pad, TPaveText * text, const std::string & title)
{
    pad->cd();
    auto left_margin{ (0.01 + pad->GetLeftMargin()) * pad->GetAbsWNDC() };
    auto right_margin{ (0.01 + pad->GetRightMargin()) * pad->GetAbsWNDC() };
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, left_margin, right_margin, 0.225, 0.001);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.2);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(text, 60, 23, 22, 0.0, kYellow-10);
    ROOTHelper::SetLineAttribute(text, 1, 0);
    text->AddText(title.data());
    pad->Update();
}

void DemoPainter::PrintResultPad(TPad * pad, TH2 * hist, bool draw_x_axis, bool draw_title_label)
{
    pad->cd();
    auto bottom_margin{ (draw_x_axis) ? 0.06 : 0.005 };
    auto top_margin{ (draw_title_label) ? 0.06 : 0.005 };
    auto title_size{ (draw_x_axis) ? 40.0f : 0.0f };
    auto label_size{ (draw_x_axis) ? 40.0f : 0.0f };
    ROOTHelper::SetPadMarginInCanvas(pad, 0.025, 0.005, bottom_margin, top_margin);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), title_size);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), label_size, 0.1f, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 40.0f, 0.01f);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.008, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 21);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
    hist->GetXaxis()->SetLimits(-1.0, 20.0);
    hist->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    hist->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);
    for (size_t i = 0; i < AtomicInfoHelper::GetStandardResidueCount(); i++)
    {
        auto residue{ AtomicInfoHelper::GetStandardResidueList().at(i) };
        auto label{ AtomicInfoHelper::GetLabel(residue) };
        auto label_index{ static_cast<int>(i) + 2 };
        hist->GetXaxis()->ChangeLabel(label_index, 90.0, -1, 12, -1, -1, label.data());
    }

    hist->SetStats(0);
    hist->Draw();
}

void DemoPainter::PrintSummaryPad(TPad * pad, TH2 * hist, bool draw_x_axis, bool draw_title_label)
{
    pad->cd();
    auto bottom_margin{ (draw_x_axis) ? 0.06 : 0.005 };
    auto top_margin{ (draw_title_label) ? 0.06 : 0.005 };
    auto title_size{ (draw_x_axis) ? 45.0f : 0.0f };
    auto label_size{ (draw_x_axis) ? 40.0f : 0.0f };
    ROOTHelper::SetPadMarginInCanvas(pad, 0.001, 0.001, bottom_margin, top_margin);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), title_size, 0.9f);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), label_size, 0.005f, 103);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 0.0);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.008, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 5);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
    hist->GetXaxis()->SetLimits(-1.0, 4.0);
    hist->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    hist->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);

    int color[4]   { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    const char * label_name[4]{ "C_{#alpha}", "C", "N", "O" };
    for (int i = 0; i < 4; i++)
    {
        hist->GetXaxis()->ChangeLabel(i+2, 0.0, -1, -1, color[i], -1, label_name[i]);
    }
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("Element");
    hist->GetXaxis()->CenterTitle();
    hist->Draw();
}

void DemoPainter::PrintCorrelationPad(TPad * pad, TH2 * hist, bool draw_x_axis, bool draw_title_label)
{
    pad->cd();
    auto bottom_margin{ (draw_x_axis) ? 0.06 : 0.005 };
    auto top_margin{ (draw_title_label) ? 0.06 : 0.005 };
    auto title_size{ (draw_x_axis) ? 45.0f : 0.0f };
    auto label_size{ (draw_x_axis) ? 40.0f : 0.0f };
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.04, bottom_margin, top_margin);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), title_size, 0.9f);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 45.0f, 1.3f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), label_size);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 40.0f, 0.02f);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.014, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.008, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 505);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("Width");
    hist->GetYaxis()->SetTitle("Amplitude");
    hist->Draw("Y+");
}

void DemoPainter::PrintLeftSideChainPad(
    TPad * pad, TH2 * hist, Residue residue, const std::string & y_title,
    const std::vector<std::string> & label_list)
{
    pad->cd();
    pad->Clear();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.04, 0.005, 0.03, 0.02);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 45.0, 1.2f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 40.0, 0.005f, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 40.0, 0.01f);
    ModifyAxisLabelSideChain(pad, hist, residue, label_list);
    hist->SetStats(0);
    hist->GetYaxis()->SetTitle(y_title.data());
    hist->GetYaxis()->CenterTitle();
    hist->Draw();
}

void DemoPainter::PrintRightSideChainPad(
    TPad * pad, TH2 * hist, Residue residue,
    const std::vector<std::string> & label_list)
{
    pad->cd();
    pad->Clear();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.04, 0.03, 0.02);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 40.0, 0.005f, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 0.0);
    ModifyAxisLabelSideChain(pad, hist, residue, label_list);
    hist->SetStats(0);
    hist->Draw();
}

void DemoPainter::PrintTitleSideChainPad(
    TPad * pad, TPaveText * text, const std::string & residue_name)
{
    pad->cd();
    auto left_margin{ (0.33 + pad->GetLeftMargin()) * pad->GetAbsWNDC() };
    auto right_margin{ (0.33 + pad->GetRightMargin()) * pad->GetAbsWNDC() };
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, left_margin, right_margin, 0.165, 0.01);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.2);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(text, 50, 103, 22, 0.0, kYellow-10);
    text->AddText(residue_name.data());
    pad->Update();
}

#endif