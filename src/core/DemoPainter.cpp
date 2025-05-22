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
#include "GlobalEnumClass.hpp"

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

    for (auto & model : m_model_object_list)
    {
        model->BuildKDTreeRoot();
    }

    if (m_model_object_list.size() == 4)
    {
        PaintAtomMapValueExample("figure_model_example.pdf");
        PaintGroupGausMainChainSummary("figure_1_a.pdf");
        PaintAtomGausMainChainDemo("figure_2_a.pdf", 0);
        PaintAtomGausMainChain("figure_1_b.pdf", 0);
        PaintAtomGausMainChain("figure_4_a.pdf", 1);
        PaintWidthToBfactorScatterPlotSummary("test_width_bfactor.pdf");
        PaintResidueClassWidthScatterPlot("figure_4_b.pdf", 1, true);
        PaintAtomRankMainChain("test_rank_a.pdf", 0);
        PaintAtomRankMainChain("test_rank_b.pdf", 1);
    }

    if (m_model_object_list.size() > 5)
    {
        PaintElementClassGroupGausToFSC("figure_5.pdf");
    }
}

void DemoPainter::PaintAtomMapValueExample(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    auto class_key{ AtomicInfoHelper::GetResidueClassKey() };
    std::cout <<"- DemoPainter::PaintAtomMapValueExample"<< std::endl;

    ModelObject * model_object;
    for (auto model : m_model_object_list)
    {
        if (model->GetPdbID() == "6Z6U")
        {
            model_object = model;
            break;
        }
    };
    if (model_object == nullptr) return;
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::vector<std::unique_ptr<TGraphErrors>> map_value_graph_list;
    std::unique_ptr<TF1> gaus_function;
    double amplitude_prior;
    double width_prior;
    std::vector<double> y_array;
    auto element{ AtomClassifier::GetMainChainElement(0) };
    auto remoteness{ AtomClassifier::GetMainChainRemoteness(0) };
    auto group_key{ KeyPackerResidueClass::Pack(Residue::ALA, element, remoteness, Branch::NONE, false) };
    if (entry_iter->IsAvailableGroupKey(group_key, class_key) == false) return;
    for (auto atom : entry_iter->GetAtomObjectList(group_key, class_key))
    {
        auto atom_iter{ std::make_unique<PotentialEntryIterator>(atom) };
        auto graph{ atom_iter->CreateBinnedDistanceToMapValueGraph() };
        ROOTHelper::SetLineAttribute(graph.get(), 1, 5, static_cast<short>(kAzure-7), 0.3f);
        map_value_graph_list.emplace_back(std::move(graph));
        auto map_value_range{ atom_iter->GetMapValueRange(0.0) };
        y_array.emplace_back(std::get<0>(map_value_range));
        y_array.emplace_back(std::get<1>(map_value_range));
    }
    gaus_function = entry_iter->CreateGroupGausFunctionPrior(group_key, class_key);
    amplitude_prior = entry_iter->GetGausEstimatePrior(group_key, class_key, 0);
    width_prior = entry_iter->GetGausEstimatePrior(group_key, class_key, 1);


    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.15) };
    auto x_min{ 0.01 };
    auto x_max{ 1.49 };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };
    auto line_ref{ ROOTHelper::CreateLine(x_min, 0.0, x_max, 0.0) };

    auto frame{ ROOTHelper::CreateHist2D("frame","", 100, x_min, x_max, 100, y_min, y_max) };
    ROOTHelper::SetPadMarginAttribute(gPad, 0.20f, 0.05f, 0.20f, 0.05f);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 90.0f, 1.0f);
    ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 90.0f, 1.0f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), 80.0f, 0.005f, 133);
    ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 80.0f, 0.01f, 133);
    ROOTHelper::SetAxisTickAttribute(frame->GetXaxis(), 0.05f, 505);
    ROOTHelper::SetAxisTickAttribute(frame->GetYaxis(), 0.05f, 505);
    ROOTHelper::SetLineAttribute(frame.get(), 1, 0);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("Radial Distance #[]{#AA}");
    frame->GetYaxis()->SetTitle("Normalized Map Value");
    frame->Draw();

    for (auto & graph : map_value_graph_list) graph->Draw("L X0");

    ROOTHelper::SetLineAttribute(line_ref.get(), 2, 3, kBlack);
    line_ref->Draw();

    ROOTHelper::SetLineAttribute(gaus_function.get(), 9, 5, kRed+1);
    gaus_function->Draw("SAME");

    auto result_text{ ROOTHelper::CreatePaveText(0.25, 0.20, 0.95, 0.40, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(result_text.get());
    ROOTHelper::SetTextAttribute(result_text.get(), 80.0f, 133, 22, 0.0, kRed+1);
    ROOTHelper::SetFillAttribute(result_text.get(), 4000);
    result_text->AddText(Form("#font[2]{A} = %.2f  ;  #font[2]{#tau} = %.2f", amplitude_prior, width_prior));
    result_text->Draw();

    auto legend{ ROOTHelper::CreateLegend(0.25, 0.82, 0.99, 0.95, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    ROOTHelper::SetTextAttribute(legend.get(), 50.0f, 133, 12, 0.0);
    legend->SetMargin(0.15f);
    legend->AddEntry(gaus_function.get(),
        "Single Gaussian Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})}", "l");
    legend->AddEntry(map_value_graph_list.at(0).get(),
        "Map Value Distribution (Median)", "l");
    legend->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void DemoPainter::PaintGroupGausMainChainSummary(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    auto residue_class{ AtomicInfoHelper::GetResidueClassKey() };
    std::cout <<"- DemoPainter::PaintGroupGausMainChainSummary"<< std::endl;

    #ifdef HAVE_ROOT

    const int main_chain_element_count{ 4 };
    short color_element[main_chain_element_count] { kRed+1, kViolet+1, kGreen+2, kAzure+2 };
    short marker_element[main_chain_element_count]{ 25, 24, 26, 32 };
    float marker_size{ 1.5f };

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 3000, 1200) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size_x{ 4 };
    const int pad_size_y{ 4 };
    std::unique_ptr<TPad> pad[pad_size_x][pad_size_y];
    std::unique_ptr<TH2> frame[pad_size_x][pad_size_y];
    
    const double x_pos[pad_size_x + 1]{ 0.00, 0.16, 0.52, 0.64, 1.00 };
    const double y_pos[pad_size_y + 1]{ 0.00, 0.30, 0.50, 0.70, 1.00 };
    for (int i = 0; i < pad_size_x; i++)
    {
        for (int j = 0; j < pad_size_y; j++)
        {
            auto pad_name{ "pad"+ std::to_string(i) + std::to_string(j) };
            pad[i][j] = ROOTHelper::CreatePad(pad_name.data(),"", x_pos[i], y_pos[j], x_pos[i+1], y_pos[j+1]);
        }
    }

    std::unique_ptr<TPaveText> info_text[pad_size_y];
    std::unique_ptr<TPaveText> resolution_text[pad_size_y];
    std::unique_ptr<TGraphErrors> amplitude_graph[pad_size_y][main_chain_element_count];
    std::unique_ptr<TGraphErrors> width_graph[pad_size_y][main_chain_element_count];
    std::unique_ptr<TGraphErrors> correlation_graph[pad_size_y][main_chain_element_count];
    std::vector<double> width_array_total;
    for (size_t j = 0; j < pad_size_y; j++)
    {
        std::vector<double> amplitude_array, width_array;
        info_text[j] = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
        resolution_text[j] = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(m_model_object_list.at(static_cast<size_t>(j))) };
        for (size_t k = 0; k < main_chain_element_count; k++)
        {
            auto group_key_list{ m_atom_classifier->GetMainChainResidueClassGroupKeyList(k) };
            amplitude_graph[j][k] = entry_iter->CreateGausEstimateToResidueGraph(group_key_list, residue_class, 0);
            width_graph[j][k] = entry_iter->CreateGausEstimateToResidueGraph(group_key_list, residue_class, 1);
            correlation_graph[j][k] = entry_iter->CreateGausEstimateScatterGraph(group_key_list, residue_class, 1, 0);
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
        auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, 0.2) };
        auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, 0.1) };

        frame[0][j] = ROOTHelper::CreateHist2D(("frame0"+std::to_string(j)).data(),"", 100, 0.0, 1.0, 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        frame[1][j] = ROOTHelper::CreateHist2D(("frame1"+std::to_string(j)).data(),"", 100, std::get<0>(width_range), std::get<1>(width_range), 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        frame[2][j] = ROOTHelper::CreateHist2D(("frame2"+std::to_string(j)).data(),"", 100, 0.0, 1.0, 100, std::get<0>(width_range), std::get<1>(width_range));
    }

    auto amplitude_title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    auto width_title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    auto correlation_title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };

    auto width_range_total{ ArrayStats<double>::ComputeScalingRangeTuple(width_array_total, 0.2) };
    for (int j = 0; j < pad_size_y; j++)
    {
        frame[2][j]->GetXaxis()->SetLimits(std::get<0>(width_range_total), std::get<1>(width_range_total));
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

    std::unique_ptr<TLegend> legend;
    for (int j = 0; j < pad_size_y; j++)
    {
        auto draw_axis_flag{ (j == 0) ? true : false };
        auto draw_title_flag{ (j == pad_size_y - 1) ? true : false };

        pad[0][j]->cd();
        auto pdb_id{ m_model_object_list.at(static_cast<size_t>(j))->GetPdbID() };
        auto emd_id{ m_model_object_list.at(static_cast<size_t>(j))->GetEmdID() };
        auto resolution{ m_model_object_list.at(static_cast<size_t>(j))->GetResolution() };
    
        // PDB & EMD Text
        auto info_top_margin{ (draw_title_flag) ? 0.12 : 0.025 };
        auto info_bottom_margin{ (draw_axis_flag) ? 0.12 : 0.025 };
        ROOTHelper::SetPaveTextMarginInCanvas(pad[0][j].get(), info_text[j].get(), 0.003, 0.00, info_bottom_margin, info_top_margin);
        ROOTHelper::SetPaveTextDefaultStyle(info_text[j].get());
        ROOTHelper::SetPaveAttribute(info_text[j].get(), 0, 0.1);
        ROOTHelper::SetFillAttribute(info_text[j].get(), 1001, kAzure-7, 0.3f);
        ROOTHelper::SetTextAttribute(info_text[j].get(), 80.0f, 103, 32);
        ROOTHelper::SetLineAttribute(info_text[j].get(), 1, 0);
        info_text[j]->AddText(pdb_id.data());
        info_text[j]->AddText(emd_id.data());
        info_text[j]->Draw();
        
        // Resolution Text
        auto top_margin{ (draw_title_flag) ? 0.10 : 0.01 };
        auto bottom_margin{ (draw_axis_flag) ? 0.20 : 0.10 };
        ROOTHelper::SetPaveTextMarginInCanvas(pad[0][j].get(), resolution_text[j].get(), 0.001, 0.07, bottom_margin, top_margin);
        ROOTHelper::SetPaveTextDefaultStyle(resolution_text[j].get());
        ROOTHelper::SetPaveAttribute(resolution_text[j].get(), 0, 0.2);
        ROOTHelper::SetFillAttribute(resolution_text[j].get(), 1001, kAzure-7);
        ROOTHelper::SetTextAttribute(resolution_text[j].get(), 80.0f, 133, 22, 0.0, kYellow-10);
        resolution_text[j]->AddText(Form("%.2f #AA", resolution));
        resolution_text[j]->Draw();

        pad[1][j]->cd();
        PrintGausResultPad(pad[1][j].get(), frame[0][j].get(), draw_axis_flag, draw_title_flag, false);
        for (int k = 0; k < main_chain_element_count; k++) amplitude_graph[j][k]->Draw("PL X0");
        if (j == pad_size_y - 1) PrintGausTitlePad(pad[1][j].get(), amplitude_title_text.get(), "Amplitude #font[2]{A}", 80.0f);

        pad[2][j]->cd();
        PrintGausCorrelationPad(pad[2][j].get(), frame[1][j].get(), draw_axis_flag, draw_title_flag);
        for (int k = 0; k < main_chain_element_count; k++) correlation_graph[j][k]->Draw("P X0");
        if (j == pad_size_y - 1) PrintGausTitlePad(pad[2][j].get(), correlation_title_text.get(), "Scatter Plot", 60.0f);

        pad[3][j]->cd();
        PrintGausResultPad(pad[3][j].get(), frame[2][j].get(), draw_axis_flag, draw_title_flag, true);
        for (int k = 0; k < main_chain_element_count; k++) width_graph[j][k]->Draw("PL X0");
        if (j == pad_size_y - 1) PrintGausTitlePad(pad[3][j].get(), width_title_text.get(), "Width #font[2]{#tau}", 80.0f);
        
        if (j == 0)
        {
            pad[0][j]->cd();
            legend = ROOTHelper::CreateLegend(0.15, 0.00, 1.00, 0.40, false);
            ROOTHelper::SetLegendDefaultStyle(legend.get());
            ROOTHelper::SetTextAttribute(legend.get(), 60.0f, 103, 12);
            ROOTHelper::SetFillAttribute(legend.get(), 4000);
            const char * label_name[4]{ "C#alpha", "C", "N", "O" };
            for (size_t k = 0; k < main_chain_element_count; k++)
            {
                legend->AddEntry(correlation_graph[j][k].get(), label_name[k], "p");
            }
            legend->SetNColumns(4);
            legend->Draw();
        }
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
                    auto amplitude_value{ entry_iter->GetGausEstimatePrior(group_key, residue_class, 0) };
                    auto width_value{ entry_iter->GetGausEstimatePrior(group_key, residue_class, 1) };
                    auto amplitude_error{ entry_iter->GetGausVariancePrior(group_key, residue_class, 0) };
                    auto width_error{ entry_iter->GetGausVariancePrior(group_key, residue_class, 1) };
                    amplitude_array.emplace_back(amplitude_value);
                    width_array.emplace_back(width_value);
                    amplitude_graph_tmp->SetPoint(count_element, count_total, amplitude_value);
                    amplitude_graph_tmp->SetPointError(count_element, 0.0, amplitude_error);
                    width_graph_tmp->SetPoint(count_element, count_total, width_value);
                    width_graph_tmp->SetPointError(count_element, 0.0, width_error);
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
    const int col_size{ 4 };
    const int row_size{ 1 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 400) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.08f, 0.01f, 0.17f, 0.12f, 0.01f, 0.005f);
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
    std::unique_ptr<TF1> fit_function[col_size];
    std::vector<double> x_array, y_array;
    double r_square[col_size]{ 0.0 };
    double slope[col_size]{ 0.0 };
    double intercept[col_size]{ 0.0 };
    for (size_t i = 0; i < primary_element_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        graph[i] = ROOTHelper::CreateGraphErrors();
        auto count{ 0 };
        for (auto model : m_model_object_list)
        {
            auto entry_iter{ std::make_unique<PotentialEntryIterator>(model) };
            auto width_value{ entry_iter->GetGausEstimatePrior(group_key, AtomicInfoHelper::GetElementClassKey(), 1) };
            auto width_error{ entry_iter->GetGausVariancePrior(group_key, AtomicInfoHelper::GetElementClassKey(), 1) };
            graph[i]->SetPoint(count, model->GetResolution(), width_value);
            graph[i]->SetPointError(count, 0.0, width_error);
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

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> title_text[col_size];
    std::unique_ptr<TPaveText> r_square_text[col_size];
    std::unique_ptr<TPaveText> fit_info_text[col_size];
    for (int i = 0; i < col_size; i++)
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
    const int col_size{ 4 };
    const int row_size{ 4 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 900) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.08f, 0.01f, 0.08f, 0.05f, 0.01f, 0.01f);
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
    std::unique_ptr<TF1> fit_function[col_size][row_size];
    std::vector<double> x_array[col_size], y_array[row_size];
    double r_square[col_size][row_size]{ {0.0} };
    double slope[col_size][row_size]{ {0.0} };
    double intercept[col_size][row_size]{ {0.0} };
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

    double x_min[col_size]{ 0.0 };
    double x_max[col_size]{ 0.0 };
    double y_min[row_size]{ 0.0 };
    double y_max[row_size]{ 0.0 };
    for (size_t i = 0; i < col_size; i++)
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

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> title_text[col_size];
    std::unique_ptr<TPaveText> r_square_text[col_size][row_size];
    std::unique_ptr<TPaveText> fit_info_text[col_size][row_size];
    for (int i = 0; i < col_size; i++)
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

void DemoPainter::PaintResidueClassWidthScatterPlot(
    const std::string & name, int par_id, bool draw_box_plot)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- DemoPainter::PaintResidueClassWidthScatterPlot"<< std::endl;
    auto residue_class{ AtomicInfoHelper::GetResidueClassKey() };

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 4 };
    const int row_size{ 4 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.08f, 0.05f, 0.07f, 0.07f, 0.01f, 0.01f);
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

    std::vector<std::unique_ptr<TGraphErrors>> graph_list[col_size][row_size];
    std::vector<double> x_array[col_size][row_size];
    std::vector<double> y_array[col_size][row_size];
    std::vector<double> global_x_array[col_size];
    std::vector<double> global_y_array[row_size];
    for (size_t i = 0; i < col_size; i++)
    {
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(m_model_object_list.at(i)) };
        for (size_t j = 0; j < row_size; j++)
        {
            for (auto residue : AtomicInfoHelper::GetStandardResidueList())
            {
                auto group_key{ m_atom_classifier->GetMainChainResidueClassGroupKey(j, residue) };
                if (entry_iter->IsAvailableGroupKey(group_key, residue_class) == false) continue;
                auto graph{ (par_id == 0) ?
                    entry_iter->CreateCOMDistanceToGausEstimateGraph(group_key, residue_class, 1) :
                    entry_iter->CreateInRangeAtomsToGausEstimateGraph(group_key, residue_class, 5.0, 1) };
                for (int p = 0; p < graph->GetN(); p++)
                {
                    x_array[i][j].emplace_back(graph->GetPointX(p));
                    y_array[i][j].emplace_back(graph->GetPointY(p));
                    global_x_array[i].emplace_back(graph->GetPointX(p));
                    global_y_array[j].emplace_back(graph->GetPointY(p));
                }
                graph_list[i][j].emplace_back(std::move(graph));
            }
        }
    }

    double x_min[col_size]{ 0.0 };
    double x_max[col_size]{ 0.0 };
    double y_min[row_size]{ 0.0 };
    double y_max[row_size]{ 0.0 };
    for (size_t i = 0; i < col_size; i++)
    {
        auto x_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(global_x_array[i], 0.05, 0.005, 0.995) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);
    }
    for (size_t j = 0; j < row_size; j++)
    {
        auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(global_y_array[j], 0.2, 0.005, 0.995) };
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }

    std::unique_ptr<TH2D> summary_hist[col_size][row_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto x_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(x_array[i][j], 0.0, 0.005, 0.995) };
            auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array[i][j], 0.0) };
            summary_hist[i][j] = ROOTHelper::CreateHist2D(
                Form("summary_hist_%d_%d", i, j), "",
                5, std::get<0>(x_range), std::get<1>(x_range),
                100, std::get<0>(y_range), std::get<1>(y_range));
            for (auto & graph : graph_list[i][j])
            {
                for (int p = 0; p < graph->GetN(); p++)
                {
                    summary_hist[i][j]->Fill(graph->GetPointX(p), graph->GetPointY(p));
                }
            }
        }
    }

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> title_x_text[col_size];
    std::unique_ptr<TPaveText> title_y_text[row_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min[j], y_max[j]);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 30.0f, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 30.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 35.0f, 1.2f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 30.0f, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            auto title{ (par_id == 0) ? "To C.M. Distance #[]{#AA}" : "In-Range Atoms" };
            frame[i][j]->GetXaxis()->SetTitle(title);
            frame[i][j]->GetYaxis()->SetTitle("Width");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();

            if (draw_box_plot == true)
            {
                ROOTHelper::SetLineAttribute(summary_hist[i][j].get(), 1, 1, color_element[j]);
                ROOTHelper::SetFillAttribute(summary_hist[i][j].get(), 1001, color_element[j], 0.3f);
                summary_hist[i][j]->Draw("CANDLE3 SAME");
            }
            else
            {
                for (auto & graph : graph_list[i][j])
                {
                    ROOTHelper::SetMarkerAttribute(graph.get(), 24, 1.0f, color_element[j]);
                    ROOTHelper::SetLineAttribute(graph.get(), 1, 2, color_element[j]);
                    graph->Draw("P");
                }
            }

            if (i == col_size - 1)
            {
                title_y_text[j] = ROOTHelper::CreatePaveText(1.02, 0.01, 1.22, 0.99, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(title_y_text[j].get());
                ROOTHelper::SetPaveAttribute(title_y_text[j].get(), 0, 0.2);
                ROOTHelper::SetTextAttribute(title_y_text[j].get(), 25.0f, 133, 22);
                ROOTHelper::SetFillAttribute(title_y_text[j].get(), 1001, color_element[j], 0.5f);
                title_y_text[j]->AddText(element_label[j]);
                title_y_text[j]->Draw();
                auto text{ title_y_text[j]->GetLineWith(element_label[j]) };
                text->SetTextAngle(90.0f);
            }
        }

        title_x_text[i] = ROOTHelper::CreatePaveText(0.01, 1.02, 0.99, 1.32, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_x_text[i].get());
        ROOTHelper::SetPaveAttribute(title_x_text[i].get(), 0, 0.2);
        ROOTHelper::SetTextAttribute(title_x_text[i].get(), 30.0f, 103, 22);
        ROOTHelper::SetFillAttribute(title_x_text[i].get(), 1001, kAzure-7, 0.5f);
        title_x_text[i]->AddText(Form("PDB-%s", m_model_object_list.at(static_cast<size_t>(i))->GetPdbID().data()));
        title_x_text[i]->AddText(m_model_object_list.at(static_cast<size_t>(i))->GetEmdID().data());
        title_x_text[i]->Draw();
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void DemoPainter::PaintAtomGausMainChainDemo(const std::string & name, int par_id)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- DemoPainter::PaintAtomGausMainChainDemo"<< std::endl;

    ModelObject * model_object;
    for (auto model : m_model_object_list)
    {
        if (model->GetPdbID() == "6Z6U")
        {
            model_object = model;
            break;
        }
    };
    if (model_object == nullptr) return;
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    
    #ifdef HAVE_ROOT

    const int main_chain_element_count{ 4 };
    short color_element[main_chain_element_count] { kRed+1, kViolet+1, kGreen+2, kAzure+2 };

    gStyle->SetLineScalePS(1.0);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 3000, 500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TGraphErrors> gaus_graph[main_chain_element_count];
    std::unique_ptr<TGraphErrors> gaus_gly_graph;
    std::unique_ptr<TGraphErrors> gaus_pro_graph;
    std::vector<double> x_array, y_array;

    for (size_t i = 0; i < main_chain_element_count; i++)
    {
        gaus_graph[i] = std::move(entry_iter->CreateGausEstimateToResidueIDGraphMap(i, par_id).at("A"));
        for (int p = 0; p < gaus_graph[i]->GetN(); p++)
        {
            x_array.emplace_back(gaus_graph[i]->GetPointX(p));
            y_array.emplace_back(gaus_graph[i]->GetPointY(p));
        }
    }
    auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.01) };
    auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array, 0.2) };
    auto frame{ ROOTHelper::CreateHist2D(
        "frame","",
        100, 0.0, std::get<1>(x_range), 100, std::get<0>(y_range), std::get<1>(y_range)) };

    gaus_gly_graph = std::move(entry_iter->CreateGausEstimateToResidueIDGraphMap(0, par_id, Residue::GLY).at("A"));
    gaus_pro_graph = std::move(entry_iter->CreateGausEstimateToResidueIDGraphMap(0, par_id, Residue::PRO).at("A"));
    for (size_t j = 0; j < main_chain_element_count; j++)
    {
        for (int p = 0; p < gaus_gly_graph->GetN(); p++)
        {
            gaus_gly_graph->SetPointY(p, std::get<1>(y_range));
        }
        for (int p = 0; p < gaus_pro_graph->GetN(); p++)
        {
            gaus_pro_graph->SetPointY(p, std::get<1>(y_range));
        }
    }

    canvas->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.07, 0.08, 0.25, 0.05);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 70.0f, 0.75f, 133);
    ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 80.0f, 1.30f, 133);
    ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), 60.0f, 0.001f, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 70.0f, 0.005f, 133);
    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.055, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.010, 1) };
    ROOTHelper::SetAxisTickAttribute(frame->GetXaxis(), static_cast<float>(x_tick_length), 520);
    ROOTHelper::SetAxisTickAttribute(frame->GetYaxis(), static_cast<float>(y_tick_length), 505);
    frame->SetStats(0);
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->SetTitle(Form("Residue ID (Chain #font[102]{A} in #font[102]{%s})", model_object->GetPdbID().data()));
    frame->GetYaxis()->SetTitle((par_id == 0) ? "Amplitude #font[2]{A}" : "Width #font[2]{#tau}");
    frame->Draw();
    for (size_t i = 0; i < main_chain_element_count; i++)
    {
        ROOTHelper::SetMarkerAttribute(gaus_graph[i].get(), 20, 1.0f, color_element[i]);
        ROOTHelper::SetLineAttribute(gaus_graph[i].get(), 1, 2, color_element[i]);
        gaus_graph[i]->Draw("PL X0");
    }
    ROOTHelper::SetMarkerAttribute(gaus_gly_graph.get(), 23, 2.5f, kOrange-6);
    gaus_gly_graph->Draw("P X0");
    ROOTHelper::SetMarkerAttribute(gaus_pro_graph.get(), 23, 2.5f, kCyan+1);
    gaus_pro_graph->Draw("P X0");

    auto legend{ ROOTHelper::CreateLegend(0.93, 0.20, 1.00, 1.00, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 70.0f, 103, 12);
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    const char * label_name[4]{ "C#alpha", "C", "N", "O" };
    legend->AddEntry(gaus_gly_graph.get(), "GLY", "p");
    legend->AddEntry(gaus_pro_graph.get(), "PRO", "p");
    for (size_t i = 0; i < main_chain_element_count; i++)
    {
        legend->AddEntry(gaus_graph[i].get(), label_name[i], "pl");
    }
    legend->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void DemoPainter::PaintAtomGausMainChain(const std::string & name, int par_id)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- DemoPainter::PaintAtomGausMainChain"<< std::endl;
    
    #ifdef HAVE_ROOT

    const int main_chain_element_count{ 4 };
    short color_element[main_chain_element_count] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };

    gStyle->SetLineScalePS(1.0);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 3000, 1800) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size{ 4 };
    std::unique_ptr<TPad> pad[pad_size + 1];
    std::unique_ptr<TH2> frame[pad_size];
    
    pad[0] = ROOTHelper::CreatePad("pad0","", 0.00, 0.00, 1.00, 0.24);
    pad[1] = ROOTHelper::CreatePad("pad1","", 0.00, 0.24, 1.00, 0.48);
    pad[2] = ROOTHelper::CreatePad("pad2","", 0.00, 0.48, 1.00, 0.72);
    pad[3] = ROOTHelper::CreatePad("pad3","", 0.00, 0.72, 1.00, 0.96);
    pad[4] = ROOTHelper::CreatePad("pad4","", 0.00, 0.96, 1.00, 1.00);

    std::unique_ptr<TPaveText> info_text[pad_size];
    std::unique_ptr<TGraphErrors> gaus_graph[pad_size][main_chain_element_count];
    std::unique_ptr<TGraphErrors> gaus_gly_graph[pad_size];
    std::unique_ptr<TGraphErrors> gaus_pro_graph[pad_size];
    std::vector<double> x_array[pad_size], y_array[pad_size];
    std::string chain_id[pad_size]{ "A", "A", "B", "A" };
    for (size_t i = 0; i < pad_size; i++)
    {
        info_text[i] = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(m_model_object_list.at(static_cast<size_t>(i))) };
        for (size_t j = 0; j < main_chain_element_count; j++)
        {
            gaus_graph[i][j] = std::move(entry_iter->CreateGausEstimateToResidueIDGraphMap(j, par_id).at(chain_id[i]));
            for (int p = 0; p < gaus_graph[i][j]->GetN(); p++)
            {
                x_array[i].emplace_back(gaus_graph[i][j]->GetPointX(p));
                y_array[i].emplace_back(gaus_graph[i][j]->GetPointY(p));
            }
        }
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array[i], 0.01) };
        auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array[i], 0.2) };
        frame[i] = ROOTHelper::CreateHist2D(
            ("frame"+std::to_string(i)).data(),"",
            100, 0.0, std::get<1>(x_range), 100, std::get<0>(y_range), std::get<1>(y_range));

        gaus_gly_graph[i] = std::move(entry_iter->CreateGausEstimateToResidueIDGraphMap(0, par_id, Residue::GLY).at(chain_id[i]));
        gaus_pro_graph[i] = std::move(entry_iter->CreateGausEstimateToResidueIDGraphMap(0, par_id, Residue::PRO).at(chain_id[i]));
        for (size_t j = 0; j < main_chain_element_count; j++)
        {
            for (int p = 0; p < gaus_gly_graph[i]->GetN(); p++)
            {
                gaus_gly_graph[i]->SetPointY(p, std::get<1>(y_range));
            }
            for (int p = 0; p < gaus_pro_graph[i]->GetN(); p++)
            {
                gaus_pro_graph[i]->SetPointY(p, std::get<1>(y_range));
            }
        }
    }

    canvas->cd();
    for (int i = 0; i < pad_size + 1; i++)
    {
        ROOTHelper::SetPadDefaultStyle(pad[i].get());
        pad[i]->Draw();
    }

    for (int i = 0; i < pad_size; i++)
    {
        pad[i]->cd();
        PrintGausResultInResidueIDPad(pad[i].get(), frame[i].get(), par_id);
        for (size_t j = 0; j < main_chain_element_count; j++)
        {
            ROOTHelper::SetMarkerAttribute(gaus_graph[i][j].get(), 20, 0.5f, color_element[j]);
            ROOTHelper::SetLineAttribute(gaus_graph[i][j].get(), 1, 1, color_element[j]);
            gaus_graph[i][j]->Draw("PL X0");
        }
        ROOTHelper::SetMarkerAttribute(gaus_gly_graph[i].get(), 23, 1.0f, kOrange-6);
        gaus_gly_graph[i]->Draw("P X0");
        ROOTHelper::SetMarkerAttribute(gaus_pro_graph[i].get(), 23, 1.0f, kMagenta-2);
        gaus_pro_graph[i]->Draw("P X0");
        PrintInfoInResidueIDPad(
            pad[i].get(), info_text[i].get(),
            m_model_object_list.at(static_cast<size_t>(i)), chain_id[i], gaus_graph[i][0]->GetN());
    }

    pad[4]->cd();
    auto legend{ ROOTHelper::CreateLegend(0.05, 0.00, 0.90, 1.00, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 50.0f, 133, 12);
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    for (size_t j = 0; j < main_chain_element_count; j++)
    {
        legend->AddEntry(gaus_graph[0][j].get(), AtomClassifier::GetMainChainElementLabel(j).data(), "pl");
    }
    legend->AddEntry(gaus_gly_graph[0].get(), "GLY", "p");
    legend->AddEntry(gaus_pro_graph[0].get(), "PRO", "p");
    legend->SetNColumns(6);
    legend->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void DemoPainter::PaintAtomRankMainChain(
    const std::string & name, int par_id)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- DemoPainter::PaintAtomRankMainChain"<< std::endl;

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 4 };
    const int row_size{ 4 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.05f, 0.15f, 0.08f, 0.07f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int main_chain_element_count{ 4 };
    short color_element[main_chain_element_count] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };

    const std::string y_title[3]
    {
        "Amplitude Rank",
        "Width Rank",
        "Intensity Rank"
    };
    
    int chain_size[row_size];
    std::unique_ptr<TH1D> rank_reference[row_size];
    std::vector<std::unique_ptr<TH1D>> rank_hist_list[row_size][3];
    std::vector<Residue> veto_residues_list{ Residue::GLY, Residue::PRO };
    for (int j = 0; j < row_size; j++)
    {
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(m_model_object_list.at(static_cast<size_t>(j))) };
        rank_hist_list[j][0] = entry_iter->CreateMainChainRankHistogram(par_id, chain_size[j], Residue::UNK, j, veto_residues_list);
        rank_hist_list[j][1] = entry_iter->CreateMainChainRankHistogram(par_id, chain_size[j], Residue::GLY, j);
        rank_hist_list[j][2] = entry_iter->CreateMainChainRankHistogram(par_id, chain_size[j], Residue::PRO, j);
        for (auto & hist : rank_hist_list[j][0]) hist->Scale(100.0 / hist->GetEntries());
        for (auto & hist : rank_hist_list[j][1]) hist->Scale(100.0 / hist->GetEntries());
        for (auto & hist : rank_hist_list[j][2]) hist->Scale(100.0 / hist->GetEntries());
        rank_reference[j] = ROOTHelper::CreateHist1D(Form("rank_reference_%d", j),"", 4, 0.5, 4.5);
        for (int p = 1; p <= rank_reference[j]->GetNbinsX(); p++)
        {
            rank_reference[j]->SetBinContent(p, 150.0);
        }
    }


    gStyle->SetTextFont(22);
    gStyle->SetPaintTextFormat(".1f");
    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> element_text[col_size];
    std::unique_ptr<TPaveText> info_text[row_size];
    std::unique_ptr<TLegend> legend;
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 0, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 100, 0.1, 4.9, 100, 0.01, 110.0);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.9f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.005f, 103, kCyan+3);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 5);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 30.0f, 1.3f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 30.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle(y_title[par_id].data());
            frame[i][j]->GetYaxis()->SetTitle("Count Ratio (%)");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();

            
            ROOTHelper::SetFillAttribute(rank_reference[j].get(), 1001, kGray, 0.2f);
            rank_reference[j]->SetBarWidth(0.9f);
            rank_reference[j]->SetBarOffset(0.05f);
            rank_reference[j]->Draw("BAR SAME");

            auto text_size{ 0.2f * 3.0f/static_cast<float>(gPad->GetAbsHNDC()) };
            ROOTHelper::SetFillAttribute(rank_hist_list[j][0].at(static_cast<size_t>(i)).get(), 1001, kAzure-7, 0.5f);
            ROOTHelper::SetLineAttribute(rank_hist_list[j][0].at(static_cast<size_t>(i)).get(), 1, 0, kAzure-7);
            ROOTHelper::SetMarkerAttribute(rank_hist_list[j][0].at(static_cast<size_t>(i)).get(), 1, text_size, kAzure-7);
            rank_hist_list[j][0].at(static_cast<size_t>(i))->SetBarWidth(0.25f);
            rank_hist_list[j][0].at(static_cast<size_t>(i))->SetBarOffset(0.125f);
            rank_hist_list[j][0].at(static_cast<size_t>(i))->Draw("BAR TEXT0 SAME");

            ROOTHelper::SetFillAttribute(rank_hist_list[j][1].at(static_cast<size_t>(i)).get(), 1001, kOrange-6, 0.5f);
            ROOTHelper::SetLineAttribute(rank_hist_list[j][1].at(static_cast<size_t>(i)).get(), 1, 0, kOrange-6);
            ROOTHelper::SetMarkerAttribute(rank_hist_list[j][1].at(static_cast<size_t>(i)).get(), 1, text_size, kOrange-6);
            rank_hist_list[j][1].at(static_cast<size_t>(i))->SetBarWidth(0.25f);
            rank_hist_list[j][1].at(static_cast<size_t>(i))->SetBarOffset(0.375f);
            rank_hist_list[j][1].at(static_cast<size_t>(i))->Draw("BAR TEXT0 SAME");

            ROOTHelper::SetFillAttribute(rank_hist_list[j][2].at(static_cast<size_t>(i)).get(), 1001, kMagenta-2, 0.5f);
            ROOTHelper::SetLineAttribute(rank_hist_list[j][2].at(static_cast<size_t>(i)).get(), 1, 0, kMagenta-2);
            ROOTHelper::SetMarkerAttribute(rank_hist_list[j][2].at(static_cast<size_t>(i)).get(), 1, text_size, kMagenta-2);
            rank_hist_list[j][2].at(static_cast<size_t>(i))->SetBarWidth(0.25f);
            rank_hist_list[j][2].at(static_cast<size_t>(i))->SetBarOffset(0.625f);
            rank_hist_list[j][2].at(static_cast<size_t>(i))->Draw("BAR TEXT0 SAME");

            if (j == row_size - 1)
            {
                element_text[i] = ROOTHelper::CreatePaveText(0.01, 1.02, 0.99, 1.30, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(element_text[i].get());
                ROOTHelper::SetPaveAttribute(element_text[i].get(), 0, 0.2);
                ROOTHelper::SetTextAttribute(element_text[i].get(), 40.0f, 133, 22);
                ROOTHelper::SetFillAttribute(element_text[i].get(), 1001, color_element[i], 0.5f);
                element_text[i]->AddText(AtomClassifier::GetMainChainElementLabel(static_cast<size_t>(i)).data());
                element_text[i]->Draw();
            }

            if (i == col_size - 1)
            {
                auto other_size{ static_cast<int>(rank_hist_list[j][0].at(static_cast<size_t>(i))->GetEntries()) };
                auto gly_size{ static_cast<int>(rank_hist_list[j][1].at(static_cast<size_t>(i))->GetEntries()) };
                auto pro_size{ static_cast<int>(rank_hist_list[j][2].at(static_cast<size_t>(i))->GetEntries()) };
                info_text[j] = ROOTHelper::CreatePaveText(1.01, 0.01, 1.77, 0.99, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(info_text[j].get());
                ROOTHelper::SetPaveAttribute(info_text[j].get(), 0, 0.1);
                ROOTHelper::SetFillAttribute(info_text[j].get(), 1001, kAzure-7, 0.5);
                ROOTHelper::SetTextAttribute(info_text[j].get(), 30.0f, 133, 12);
                ROOTHelper::SetLineAttribute(info_text[j].get(), 1, 0);
                info_text[j]->AddText("---------------------");
                info_text[j]->AddText(("#font[102]{PDB-" + m_model_object_list.at(static_cast<size_t>(j))->GetPdbID() +"}").data());
                info_text[j]->AddText(("#font[102]{"+ m_model_object_list.at(static_cast<size_t>(j))->GetEmdID() +"}").data());
                info_text[j]->AddText(Form("#GLY = %d", gly_size));
                info_text[j]->AddText(Form("#PRO = %d", pro_size));
                info_text[j]->AddText(Form("#Others = %d", other_size));
                info_text[j]->AddText("---------------------");
                info_text[j]->Draw();
            }

            if (i == col_size - 1 && j == row_size - 1)
            {
                legend = ROOTHelper::CreateLegend(1.01, 1.02, 1.77, 1.30, true);
                ROOTHelper::SetLegendDefaultStyle(legend.get());
                ROOTHelper::SetTextAttribute(legend.get(), 20.0f, 133, 12);
                ROOTHelper::SetFillAttribute(legend.get(), 4000);
                legend->AddEntry(rank_hist_list[j][0].at(0).get(), "Other Residues", "f");
                legend->AddEntry(rank_hist_list[j][1].at(0).get(), "GLY", "f");
                legend->AddEntry(rank_hist_list[j][2].at(0).get(), "PRO", "f");
                legend->Draw();
            }
        }
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

void DemoPainter::PrintGausTitlePad(
    TPad * pad, TPaveText * text, const std::string & title, float text_size)
{
    pad->cd();
    auto left_margin{ 0.001 + pad->GetLeftMargin() * pad->GetAbsWNDC() };
    auto right_margin{ 0.001 + pad->GetRightMargin() * pad->GetAbsWNDC() };
    auto bottom_margin{ 0.005 + (1.0 - pad->GetTopMargin()) * pad->GetAbsHNDC() };
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, left_margin, right_margin, bottom_margin, 0.005);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.2);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(text, text_size, 23, 22, 0.0, kYellow-10);
    ROOTHelper::SetLineAttribute(text, 1, 0);
    text->AddText(title.data());
    pad->Update();
}

void DemoPainter::PrintGausResultPad(
    TPad * pad, TH2 * hist, bool draw_x_axis, bool draw_title_label, bool is_right_side_pad)
{
    pad->cd();
    auto left_margin{ (is_right_side_pad) ? 0.005 : 0.030 };
    auto right_margin{ (is_right_side_pad) ? 0.030 : 0.005 };
    auto bottom_margin{ (draw_x_axis) ? 0.105 : 0.010 };
    auto top_margin{ (draw_title_label) ? 0.105 : 0.010 };
    auto label_size{ (draw_x_axis) ? 55.0f : 0.0f };
    ROOTHelper::SetPadMarginInCanvas(pad, left_margin, right_margin, bottom_margin, top_margin);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), label_size, 0.17f, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 50.0f, 0.01f);

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
    hist->Draw((is_right_side_pad) ? "Y+" : "");
}

void DemoPainter::PrintGausCorrelationPad(TPad * pad, TH2 * hist, bool draw_x_axis, bool draw_title_label)
{
    pad->cd();
    auto bottom_margin{ (draw_x_axis) ? 0.105 : 0.010 };
    auto top_margin{ (draw_title_label) ? 0.105 : 0.010 };
    auto title_size{ (draw_x_axis) ? 60.0f : 0.0f };
    auto label_size{ (draw_x_axis) ? 50.0f : 0.0f };
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.005, bottom_margin, top_margin);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), title_size, 0.9f);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), label_size);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 0.0f);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.022, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.008, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 505);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("Width #font[2]{#tau}");
    hist->GetXaxis()->CenterTitle();
    hist->Draw("");
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

void DemoPainter::PrintGausResultInResidueIDPad(TPad * pad, TH2 * hist, int par_id)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.05, 0.10, 0.05, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(pad, 0, 0, 4000, 0, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 45.0f, 0.9f, 133);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 50.0f, 1.5f, 133);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 45.0f, 0.01f, 133);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 45.0f, 0.005f, 133);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.008, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 520);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 505);

    hist->SetStats(0);
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetTitle("Residue ID");
    hist->GetYaxis()->SetTitle((par_id == 0) ? "Amplitude" : "Width");
    hist->Draw();
}

void DemoPainter::PrintInfoInResidueIDPad(
    TVirtualPad * pad, TPaveText * text, const ModelObject * model_object,
    const std::string & chain_id, int residue_size)
{
    pad->cd();
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, 0.903, 0.003, 0.05, 0.01);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.1);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7, 0.5);
    ROOTHelper::SetTextAttribute(text, 40.0f, 133, 12);
    ROOTHelper::SetLineAttribute(text, 1, 0);
    text->AddText("---------------------");
    text->AddText(("#font[102]{PDB-" + model_object->GetPdbID() +"}").data());
    text->AddText(("#font[102]{"+ model_object->GetEmdID() +"}").data());
    text->AddText(Form("Chain ID = %s", chain_id.data()));
    text->AddText(Form("#Residues = %d", residue_size));
    text->AddText("---------------------");
    pad->Update();
}

#endif