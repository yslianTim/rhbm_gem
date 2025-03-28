#include "ModelPainter.hpp"
#include "ModelObject.hpp"
#include "DataObjectBase.hpp"
#include "PotentialEntryIterator.hpp"
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

ModelPainter::ModelPainter(void) :
    m_folder_path{ "./" }
{

}

ModelPainter::~ModelPainter()
{

}

void ModelPainter::SetFolder(const std::string & folder_path)
{
    m_folder_path = FilePathHelper::EnsureTrailingSlash(folder_path);
}

void ModelPainter::AddDataObject(DataObjectBase * data_object)
{
    m_model_object_list.push_back(dynamic_cast<ModelObject *>(data_object));
}

void ModelPainter::AddReferenceDataObject(DataObjectBase * data_object, const std::string & label)
{
    m_ref_model_object_map[label] = dynamic_cast<ModelObject *>(data_object);
}

void ModelPainter::Painting(void)
{
    std::cout <<"- ModelPainter::Painting"<<std::endl;
    std::cout <<"  Folder path: "<< m_folder_path << std::endl;
    std::cout <<"  Number of model objects to be painted: "<< m_model_object_list.size() << std::endl;
    for (auto model_object : m_model_object_list)
    {
        auto plot_main_chain_name{ "residue_class_group_gaus_main_"+ model_object->GetPdbID() +".pdf" };
        PaintResidueClassGroupGausMainChain(model_object, plot_main_chain_name);
        auto plot_side_chain_name{ "residue_class_group_gaus_side_"+ model_object->GetPdbID() +".pdf" };
        PaintResidueClassGroupGausSideChain(model_object, plot_side_chain_name);
    }
    if (m_model_object_list.size() == 4)
    {
        PaintResidueClassGroupGausMainChain("figure_1_part2.pdf");
    }
}

void ModelPainter::PaintResidueClassGroupGausMainChain(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintResidueClassGroupGausMainChain"<< std::endl;

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
    int marker_element[primary_element_size]  { 25, 24, 26, 32 };
    double marker_size{ 1.5 };

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 3000, 2000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size_x{ 6 };
    const int pad_size_y{ 4 };
    const int pad_extra_size{ 4 };
    std::unique_ptr<TPad> pad[pad_size_x][pad_size_y];
    std::unique_ptr<TPad> pad_extra[pad_extra_size];
    std::unique_ptr<TH2> frame[pad_size_x][pad_size_y];
    std::unique_ptr<TH2> frame_extra[pad_extra_size];
    
    const std::string fsc_list[pad_size_y]{ "1.25", "1.52", "1.71", "1.90" };
    const double x_pos[pad_size_x + 1]{ 0.00, 0.10, 0.40, 0.47, 0.77, 0.84, 1.00 };
    const double y_pos[pad_size_y + 1]{ 0.20, 0.43, 0.62, 0.81, 1.00 };
    for (int i = 0; i < pad_size_x; i++)
    {
        for (int j = 0; j < pad_size_y; j++)
        {
            auto pad_name{ "pad"+ std::to_string(i) + std::to_string(j) };
            pad[i][j] = ROOTHelper::CreatePad(pad_name.data(),"", x_pos[i], y_pos[j], x_pos[i+1], y_pos[j+1]);
        }
    }

    const double x_extra_pos[pad_extra_size + 1]{ 0.00, 0.25, 0.50, 0.75, 1.00 };
    for (int i = 0; i < pad_extra_size; i++)
    {
        auto pad_name{ "pad_extra"+ std::to_string(i) };
        auto frame_name{ "frame_extra"+ std::to_string(i) };
        pad_extra[i] = ROOTHelper::CreatePad(pad_name.data(),"", x_extra_pos[i], 0.00, x_extra_pos[i+1], 0.20);
        frame_extra[i] = ROOTHelper::CreateHist2D(frame_name.data(),"", 100, 0.0, 1.0, 100, 0.0, 1.0);
    }

    std::unique_ptr<TPaveText> info_text[pad_size_y];
    std::unique_ptr<TPaveText> icon_text[pad_size_y];
    std::unique_ptr<TGraphErrors> amplitude_graph[pad_size_y][primary_element_size];
    std::unique_ptr<TGraphErrors> width_graph[pad_size_y][primary_element_size];
    std::unique_ptr<TGraphErrors> correlation_graph[pad_size_y][primary_element_size];
    std::unique_ptr<TH2D> amplitude_hist[pad_size_y][primary_element_size];
    std::unique_ptr<TH2D> width_hist[pad_size_y][primary_element_size];
    std::vector<double> width_array_total;
    for (int j = 0; j < pad_size_y; j++)
    {
        std::vector<double> amplitude_array, width_array;
        info_text[j] = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
        icon_text[j] = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(m_model_object_list.at(j)) };
        for (int k = 0; k < primary_element_size; k++)
        {
            auto element{ element_index[k] };
            auto remoteness{ remoteness_index[k] };
            amplitude_graph[j][k] = ROOTHelper::CreateGraphErrors();
            width_graph[j][k] = ROOTHelper::CreateGraphErrors();
            correlation_graph[j][k] = ROOTHelper::CreateGraphErrors();
            for (auto residue : AtomicInfoHelper::GetStandardResidueList())
            {
                auto group_key{ std::make_tuple(residue, element, remoteness, 0, false) };
                if (entry_iter->IsAvailableGroupKey(group_key) == false) continue;
                auto gaus_estimate{ entry_iter->GetGausEstimatePrior(group_key) };
                auto gaus_variance{ entry_iter->GetGausVariancePrior(group_key) };
                amplitude_array.push_back(std::get<0>(gaus_estimate));
                width_array.push_back(std::get<1>(gaus_estimate));
                width_array_total.push_back(std::get<1>(gaus_estimate));
                amplitude_graph[j][k]->SetPoint(residue, residue, std::get<0>(gaus_estimate));
                amplitude_graph[j][k]->SetPointError(residue, 0.0, std::get<0>(gaus_variance));
                width_graph[j][k]->SetPoint(residue, residue, std::get<1>(gaus_estimate));
                width_graph[j][k]->SetPointError(residue, 0.0, std::get<1>(gaus_variance));
                correlation_graph[j][k]->SetPoint(residue, std::get<1>(gaus_estimate), std::get<0>(gaus_estimate));
                correlation_graph[j][k]->SetPointError(residue, std::get<1>(gaus_variance), std::get<0>(gaus_variance));
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
            ROOTHelper::SetFillAttribute(amplitude_hist[j][k].get(), 1001, color_element[k], 0.3);
            ROOTHelper::SetFillAttribute(width_hist[j][k].get(), 1001, color_element[k], 0.3);

            amplitude_hist[j][k]->SetBarWidth(0.6);
            width_hist[j][k]->SetBarWidth(0.6);
        }
        frame[0][j] = ROOTHelper::CreateHist2D(("frame0"+std::to_string(j)).data(),"", 100, 0.0, 1.0, 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        frame[1][j] = ROOTHelper::CreateHist2D(("frame1"+std::to_string(j)).data(),"", 100, 0.0, 1.0, 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        frame[2][j] = ROOTHelper::CreateHist2D(("frame2"+std::to_string(j)).data(),"", 100, 0.0, 1.0, 100, std::get<0>(width_range), std::get<1>(width_range));
        frame[3][j] = ROOTHelper::CreateHist2D(("frame3"+std::to_string(j)).data(),"", 100, 0.0, 1.0, 100, std::get<0>(width_range), std::get<1>(width_range));
        frame[4][j] = ROOTHelper::CreateHist2D(("frame4"+std::to_string(j)).data(),"", 100, std::get<0>(width_range), std::get<1>(width_range), 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    }

    std::unique_ptr<TPaveText> side_chain_title_text[pad_extra_size];
    std::vector<std::unique_ptr<TGraphErrors>> amplitude_graph_list[2];
    std::vector<std::unique_ptr<TGraphErrors>> width_graph_list[2];
    std::vector<std::string> label_list[2];
    std::vector<double> amplitude_array, width_array;
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(m_model_object_list.at(0)) };
    int residue_list[2]{ 3, 6 };
    for (int i = 0; i < 2; i++)
    {
        auto residue{ residue_list[i] };
        auto count_total{ 0 };
        for (auto element : AtomicInfoHelper::GetStandardElementList())
        {
            auto amplitude_graph{ ROOTHelper::CreateGraphErrors() };
            auto width_graph{ ROOTHelper::CreateGraphErrors() };
            auto count_element{ 0 };
            for (auto remoteness : AtomicInfoHelper::GetStandardRemotenessList())
            {
                for (auto branch : AtomicInfoHelper::GetStandardBranchList())
                {
                    auto group_key{ std::make_tuple(residue, element, remoteness, branch, false) };
                    if (entry_iter->IsAvailableGroupKey(group_key) == false) continue;
                    auto gaus_estimate{ entry_iter->GetGausEstimatePrior(group_key) };
                    auto gaus_variance{ entry_iter->GetGausVariancePrior(group_key) };
                    amplitude_array.push_back(std::get<0>(gaus_estimate));
                    width_array.push_back(std::get<1>(gaus_estimate));
                    amplitude_graph->SetPoint(count_element, count_total, std::get<0>(gaus_estimate));
                    amplitude_graph->SetPointError(count_element, 0.0, std::get<0>(gaus_variance));
                    width_graph->SetPoint(count_element, count_total, std::get<1>(gaus_estimate));
                    width_graph->SetPointError(count_element, 0.0, std::get<1>(gaus_variance));
                    auto element_label{ AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::ElementTag>(element) };
                    auto remoteness_label{ AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::RemotenessTag>(remoteness) };
                    auto branch_label{ AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::BranchTag>(branch) };
                    label_list[i].emplace_back(element_label +"_{"+ remoteness_label + branch_label+"}");
                    count_element++;
                    count_total++;
                }
            }
            auto color{ AtomicInfoHelper::AtomColorMapping<AtomicInfoHelper::ElementTag>(element) };
            ROOTHelper::SetMarkerAttribute(amplitude_graph.get(), kFullCircle, marker_size, color);
            ROOTHelper::SetMarkerAttribute(width_graph.get(), kFullCircle, marker_size, color);
            amplitude_graph_list[i].push_back(std::move(amplitude_graph));
            width_graph_list[i].push_back(std::move(width_graph));
        }
    }
    auto scaling{ 0.3 };
    auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, scaling) };
    auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, scaling) };
    frame_extra[0]->GetYaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    frame_extra[1]->GetYaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    frame_extra[2]->GetYaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));
    frame_extra[3]->GetYaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));

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
    for (int i = 0; i < pad_extra_size; i++)
    {
        side_chain_title_text[i] = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
        ROOTHelper::SetPadDefaultStyle(pad_extra[i].get());
        pad_extra[i]->Draw();
    }

    for (int j = 0; j < pad_size_y; j++)
    {
        auto draw_axis_flag{ (j == 0) ? true : false };

        pad[0][j]->cd();
        PrintInfoMainChainPad(pad[0][j].get(), info_text[j].get(), m_model_object_list.at(j)->GetPdbID(), m_model_object_list.at(j)->GetEmdID(), draw_axis_flag);
        PrintIconMainChainPad(pad[0][j].get(), icon_text[j].get(), fsc_list[j], draw_axis_flag);

        pad[1][j]->cd();
        PrintResultPad(pad[1][j].get(), frame[0][j].get(), draw_axis_flag);
        for (int k = 0; k < primary_element_size; k++) amplitude_graph[j][k]->Draw("PL X0");
        if (j == pad_size_y - 1) PrintTitlePad(pad[1][j].get(), amplitude_title_text.get(), "Amplitude");

        pad[2][j]->cd();
        PrintSummaryPad(pad[2][j].get(), frame[1][j].get(), draw_axis_flag);
        for (int k = 0; k < primary_element_size; k++) amplitude_hist[j][k]->Draw("CANDLE3 SAME");

        pad[3][j]->cd();
        PrintResultPad(pad[3][j].get(), frame[2][j].get(), draw_axis_flag);
        for (int k = 0; k < primary_element_size; k++) width_graph[j][k]->Draw("PL X0");
        if (j == pad_size_y - 1) PrintTitlePad(pad[3][j].get(), width_title_text.get(), "Width");

        pad[4][j]->cd();
        PrintSummaryPad(pad[4][j].get(), frame[3][j].get(), draw_axis_flag);
        for (int k = 0; k < primary_element_size; k++) width_hist[j][k]->Draw("CANDLE3 SAME");

        pad[5][j]->cd();
        PrintCorrelationPad(pad[5][j].get(), frame[4][j].get(), draw_axis_flag);
        for (int k = 0; k < primary_element_size; k++) correlation_graph[j][k]->Draw("P X0");
        if (j == pad_size_y - 1) PrintTitlePad(pad[5][j].get(), correlation_title_text.get(), "Correlation");
    }

    PrintLeftSideChainPad(pad_extra[0].get(), frame_extra[0].get(), residue_list[0],"Amplitude", label_list[0]);
    PrintRightSideChainPad(pad_extra[1].get(), frame_extra[1].get(), residue_list[1], label_list[1]);
    PrintLeftSideChainPad(pad_extra[2].get(), frame_extra[2].get(), residue_list[0],"Width", label_list[0]);
    PrintRightSideChainPad(pad_extra[3].get(), frame_extra[3].get(), residue_list[1], label_list[1]);

    pad_extra[0]->cd();
    for (auto & graph : amplitude_graph_list[0]) graph->Draw("P X0");
    PrintTitleSideChainPad(pad_extra[0].get(), side_chain_title_text[0].get(), AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::ResidueTag>(residue_list[0]));

    pad_extra[1]->cd();
    for (auto & graph : amplitude_graph_list[1]) graph->Draw("P X0");
    PrintTitleSideChainPad(pad_extra[1].get(), side_chain_title_text[1].get(), AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::ResidueTag>(residue_list[1]));

    pad_extra[2]->cd();
    for (auto & graph : width_graph_list[0]) graph->Draw("P X0");
    PrintTitleSideChainPad(pad_extra[2].get(), side_chain_title_text[2].get(), AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::ResidueTag>(residue_list[0]));

    pad_extra[3]->cd();
    for (auto & graph : width_graph_list[1]) graph->Draw("P X0");
    PrintTitleSideChainPad(pad_extra[3].get(), side_chain_title_text[3].get(), AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::ResidueTag>(residue_list[1]));

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ModelPainter::PaintResidueClassGroupGausMainChain(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintResidueClassGroupGausMainChain"<< std::endl;
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    
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
    int marker_element[primary_element_size]  { 54, 53, 55, 59 };

    gStyle->SetLineScalePS(1.0);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1200, 600) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size{ 6 };
    std::unique_ptr<TPad> pad[pad_size];
    std::unique_ptr<TH2> frame[pad_size];
    pad[0] = ROOTHelper::CreatePad("pad0","", 0.00, 0.00, 0.52, 0.55); // The bottom-left pad
    pad[1] = ROOTHelper::CreatePad("pad1","", 0.00, 0.55, 0.52, 1.00); // The top-left pad
    pad[2] = ROOTHelper::CreatePad("pad2","", 0.52, 0.00, 0.65, 0.55); // The bottom-middle pad
    pad[3] = ROOTHelper::CreatePad("pad3","", 0.52, 0.55, 0.65, 1.00); // The top-middle pad
    pad[4] = ROOTHelper::CreatePad("pad4","", 0.65, 0.00, 1.00, 0.80); // The bottom-right pad
    pad[5] = ROOTHelper::CreatePad("pad5","", 0.65, 0.80, 1.00, 1.00); // The top-right pad

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
        for (auto residue : AtomicInfoHelper::GetStandardResidueList())
        {
            auto group_key{ std::make_tuple(residue, element, remoteness, 0, false) };
            if (entry_iter->IsAvailableGroupKey(group_key) == false) continue;
            auto gaus_estimate{ entry_iter->GetGausEstimatePrior(group_key) };
            auto gaus_variance{ entry_iter->GetGausVariancePrior(group_key) };
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

    frame[1] = ROOTHelper::CreateHist2D("hist_1","", 100, 0.0, 1.0, 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    frame[0] = ROOTHelper::CreateHist2D("hist_0","", 100, 0.0, 1.0, 100, std::get<0>(width_range), std::get<1>(width_range));
    frame[3] = ROOTHelper::CreateHist2D("hist_3","", 100, 0.0, 1.0, 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    frame[2] = ROOTHelper::CreateHist2D("hist_2","", 100, 0.0, 1.0, 100, std::get<0>(width_range), std::get<1>(width_range));
    frame[4] = ROOTHelper::CreateHist2D("hist_4","", 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range), 100, std::get<0>(width_range), std::get<1>(width_range));
    auto info_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    auto icon_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };

    PrintAmplitudePad(pad[1].get(), frame[1].get());
    PrintWidthPad(pad[0].get(), frame[0].get());
    PrintAmplitudeSummaryPad(pad[3].get(), frame[3].get());
    PrintWidthSummaryPad(pad[2].get(), frame[2].get());
    PrintGausSummaryPad(pad[4].get(), frame[4].get());
    PrintInfoPad(pad[5].get(), info_text.get(), model_object->GetPdbID() , model_object->GetEmdID());
    PrintIconPad(pad[5].get(), icon_text.get());

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

void ModelPainter::PaintResidueClassGroupGausSideChain(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintResidueClassGroupGausSideChain"<<std::endl;
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    
    #ifdef HAVE_ROOT

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1200, 300) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int pad_size{ 3 };
    std::unique_ptr<TPad> pad[pad_size];
    std::unique_ptr<TH2> frame[pad_size];
    pad[0] = ROOTHelper::CreatePad("pad0","", 0.00, 0.00, 0.50, 1.00); // The left pad
    pad[1] = ROOTHelper::CreatePad("pad1","", 0.50, 0.00, 1.00, 1.00); // The right pad
    pad[2] = ROOTHelper::CreatePad("pad2","", 0.46, 0.80, 0.54, 0.99); // The middle-top pad
    frame[0] = ROOTHelper::CreateHist2D("hist_0","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[1] = ROOTHelper::CreateHist2D("hist_1","", 100, 0.0, 1.0, 100, 0.0, 1.0);

    std::vector<std::unique_ptr<TGraphErrors>> amplitude_graph_list[AtomicInfoHelper::GetStandardResidueCount()];
    std::vector<std::unique_ptr<TGraphErrors>> width_graph_list[AtomicInfoHelper::GetStandardResidueCount()];
    std::unique_ptr<TPaveText> info_text[AtomicInfoHelper::GetStandardResidueCount()];

    for (auto residue : AtomicInfoHelper::GetStandardResidueList())
    {
        info_text[residue] = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
        std::vector<double> amplitude_array, width_array;
        std::vector<std::string> label_list;
        auto count_total{ 0 };
        for (auto element : AtomicInfoHelper::GetStandardElementList())
        {
            auto amplitude_graph{ ROOTHelper::CreateGraphErrors() };
            auto width_graph{ ROOTHelper::CreateGraphErrors() };
            auto count_element{ 0 };
            for (auto remoteness : AtomicInfoHelper::GetStandardRemotenessList())
            {
                for (auto branch : AtomicInfoHelper::GetStandardBranchList())
                {
                    auto group_key{ std::make_tuple(residue, element, remoteness, branch, false) };
                    if (entry_iter->IsAvailableGroupKey(group_key) == false) continue;
                    auto gaus_estimate{ entry_iter->GetGausEstimatePrior(group_key) };
                    auto gaus_variance{ entry_iter->GetGausVariancePrior(group_key) };
                    amplitude_array.push_back(std::get<0>(gaus_estimate));
                    width_array.push_back(std::get<1>(gaus_estimate));
                    amplitude_graph->SetPoint(count_element, count_total, std::get<0>(gaus_estimate));
                    amplitude_graph->SetPointError(count_element, 0.0, std::get<0>(gaus_variance));
                    width_graph->SetPoint(count_element, count_total, std::get<1>(gaus_estimate));
                    width_graph->SetPointError(count_element, 0.0, std::get<1>(gaus_variance));
                    auto element_label{ AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::ElementTag>(element) };
                    auto remoteness_label{ AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::RemotenessTag>(remoteness) };
                    auto branch_label{ AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::BranchTag>(branch) };
                    label_list.emplace_back(element_label +"_{"+ remoteness_label + branch_label+"}");
                    count_element++;
                    count_total++;
                }
            }
            auto color{ AtomicInfoHelper::AtomColorMapping<AtomicInfoHelper::ElementTag>(element) };
            ROOTHelper::SetMarkerAttribute(amplitude_graph.get(), kFullCircle, 1.5, color);
            ROOTHelper::SetMarkerAttribute(width_graph.get(), kFullCircle, 1.5, color);
            amplitude_graph_list[residue].push_back(std::move(amplitude_graph));
            width_graph_list[residue].push_back(std::move(width_graph));
        }

        auto scaling{ 0.3 };
        auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, scaling) };
        auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, scaling) };

        canvas->cd();
        for (int i = 0; i < pad_size; i++)
        {
            ROOTHelper::SetPadDefaultStyle(pad[i].get());
            pad[i]->Draw();
        }
        frame[0]->GetYaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        frame[1]->GetYaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));

        PrintAmplitudeSideChainPad(pad[0].get(), frame[0].get(), residue, label_list);
        PrintWidthSideChainPad(pad[1].get(), frame[1].get(), residue, label_list);
        PrintInfoSideChainPad(pad[2].get(), info_text[residue].get(), AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::ResidueTag>(residue));

        pad[0]->cd();
        for (auto & graph : amplitude_graph_list[residue])
        {
            graph->Draw("P X0");
        }

        pad[1]->cd();
        for (auto & graph : width_graph_list[residue])
        {
            graph->Draw("P X0");
        }

        pad[2]->cd();
        info_text[residue]->Draw();

        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ModelPainter::PaintResidueClassMapValue(ModelObject * model_object, const std::string & name)
{

}

#ifdef HAVE_ROOT

void ModelPainter::PrintIconPad(TPad * pad, TPaveText * text)
{
    pad->cd();
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, 0.26, 0.005, 0.01, 0.02);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.1);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7);
    ROOTHelper::SetLineAttribute(text, 1, 0);
    pad->Update();
}

void ModelPainter::PrintInfoPad(
    TPad * pad, TPaveText * text, const std::string & pdb_id, const std::string & emd_id)
{
    pad->cd();
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, 0.005, 0.1, 0.01, 0.02);
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

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, 21);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 506);
    hist->GetXaxis()->SetLimits(-1.0, 20.0);
    hist->SetStats(0);
    hist->GetYaxis()->SetTitle("Amplitude");
    hist->GetYaxis()->CenterTitle();
    hist->Draw();
}

void ModelPainter::PrintWidthPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.08, 0.005, 0.13, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0, 1.4);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 32.0, 0.11, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0, 0.01);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, 21);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 505);
    hist->GetXaxis()->SetLimits(-1.0, 20.0);
    hist->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    hist->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);
    for (int i = 0; i < AtomicInfoHelper::GetStandardResidueCount(); i++)
    {
        auto label{ AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::ResidueTag>(i) };
        hist->GetXaxis()->ChangeLabel(i+2, 90.0, -1, 12, -1, -1, label.data());
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

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, 5);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 506);
    hist->GetXaxis()->SetLimits(-1.0, 4.0);
    hist->SetStats(0);
    hist->Draw();
}

void ModelPainter::PrintWidthSummaryPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.005, 0.13, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 35.0, 1.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 35.0, 0.005, 103);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 0.0);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, 5);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 505);
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

void ModelPainter::PrintGausSummaryPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.07, 0.13, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 35.0, 1.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0, 1.2);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 30.0);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0, 0.01);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.03, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, 505);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 505);
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("Amplitude");
    hist->GetYaxis()->SetTitle("Width");
    hist->Draw("Y+");
}

void ModelPainter::PrintInfoSideChainPad(
    TPad * pad, TPaveText * text, const std::string & residue_name)
{
    pad->cd();
    pad->Clear();
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, 0.00, 0.00, 0.00, 0.00);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.2);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(text, 45, 103, 22, 0.0, kYellow-10);
    text->AddText(residue_name.data());
    pad->Update();
}

void ModelPainter::PrintAmplitudeSideChainPad(
    TPad * pad, TH2 * hist, int residue, const std::vector<std::string> & label_list)
{
    pad->cd();
    pad->Clear();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.07, 0.005, 0.16, 0.04);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0, 1.2);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 32.0, 0.005, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0, 0.01);
    ModifyAxisLabelSideChain(pad, hist, residue, label_list);
    hist->SetStats(0);
    hist->GetYaxis()->SetTitle("Amplitude");
    hist->GetYaxis()->CenterTitle();
    hist->Draw();
}

void ModelPainter::PrintWidthSideChainPad(
    TPad * pad, TH2 * hist, int residue, const std::vector<std::string> & label_list)
{
    pad->cd();
    pad->Clear();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.07, 0.16, 0.04);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0, 1.2);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 32.0, 0.005, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0, 0.01);
    ModifyAxisLabelSideChain(pad, hist, residue, label_list);
    hist->SetStats(0);
    hist->GetYaxis()->SetTitle("Width");
    hist->GetYaxis()->CenterTitle();
    hist->Draw("Y+");
}

void ModelPainter::ModifyAxisLabelSideChain(
    TPad * pad, TH2 * hist, int residue, const std::vector<std::string> & label_list)
{
    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, label_list.size()+1);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 506);

    hist->GetXaxis()->SetLimits(-1.0, label_list.size());
    hist->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    hist->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);

    for (int i = 0; i < label_list.size(); i++)
    {
        hist->GetXaxis()->ChangeLabel(i+2, 0.0, -1, -1, -1, -1, label_list.at(i).data());
    }
}

void ModelPainter::PrintIconMainChainPad(
    TPad * pad, TPaveText * text, const std::string & fsc, bool is_bottom_pad)
{
    pad->cd();
    auto bottom_margin{ (is_bottom_pad) ? 0.16 : 0.12 };
    auto top_margin{ 0.01 };
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, 0.005, 0.005, bottom_margin, top_margin);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.2);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(text, 80, 133, 22, 0.0, kYellow-10);
    text->AddText((fsc + " #AA").data());
    pad->Update();
}

void ModelPainter::PrintInfoMainChainPad(
    TPad * pad, TPaveText * text,
    const std::string & pdb_id, const std::string & emd_id, bool is_bottom_pad)
{
    pad->cd();
    auto bottom_margin{ (is_bottom_pad) ? 0.045 : 0.005 };
    auto top_margin{ 0.06 };
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

void ModelPainter::PrintTitlePad(TPad * pad, TPaveText * text, const std::string & title)
{
    pad->cd();
    auto left_margin{ (0.15 + pad->GetLeftMargin()) * pad->GetAbsWNDC() };
    auto right_margin{ (0.15 + pad->GetRightMargin()) * pad->GetAbsWNDC() };
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, left_margin, right_margin, 0.165, 0.001);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.2);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(text, 40, 133, 22, 0.0, kYellow-10);
    ROOTHelper::SetLineAttribute(text, 1, 0);
    text->AddText(title.data());
    pad->Update();
}

void ModelPainter::PrintResultPad(TPad * pad, TH2 * hist, bool draw_x_axis)
{
    pad->cd();
    auto bottom_margin{ (draw_x_axis) ? 0.045 : 0.005 };
    auto title_size{ (draw_x_axis) ? 40.0 : 0.0 };
    auto label_size{ (draw_x_axis) ? 40.0 : 0.0 };
    ROOTHelper::SetPadMarginInCanvas(pad, 0.025, 0.005, bottom_margin, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), title_size);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), label_size, 0.1, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 40.0, 0.01);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.008, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, 21);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 505);
    hist->GetXaxis()->SetLimits(-1.0, 20.0);
    hist->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    hist->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);
    for (int i = 0; i < AtomicInfoHelper::GetStandardResidueCount(); i++)
    {
        auto label{ AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::ResidueTag>(i) };
        hist->GetXaxis()->ChangeLabel(i+2, 90.0, -1, 12, -1, -1, label.data());
    }

    hist->SetStats(0);
    hist->Draw();
}

void ModelPainter::PrintSummaryPad(TPad * pad, TH2 * hist, bool draw_x_axis)
{
    pad->cd();
    auto bottom_margin{ (draw_x_axis) ? 0.045 : 0.005 };
    auto title_size{ (draw_x_axis) ? 45.0 : 0.0 };
    auto label_size{ (draw_x_axis) ? 40.0 : 0.0 };
    ROOTHelper::SetPadMarginInCanvas(pad, 0.001, 0.001, bottom_margin, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), title_size, 0.9);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), label_size, 0.005, 103);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 0.0);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.008, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, 5);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 505);
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

void ModelPainter::PrintCorrelationPad(TPad * pad, TH2 * hist, bool draw_x_axis)
{
    pad->cd();
    auto bottom_margin{ (draw_x_axis) ? 0.045 : 0.005 };
    auto title_size{ (draw_x_axis) ? 45.0 : 0.0 };
    auto label_size{ (draw_x_axis) ? 40.0 : 0.0 };
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.04, bottom_margin, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), title_size, 0.9);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 45.0, 1.3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), label_size);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 40.0, 0.02);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.014, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.008, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), x_tick_length, 505);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), y_tick_length, 505);
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("Width");
    hist->GetYaxis()->SetTitle("Amplitude");
    hist->Draw("Y+");
}

void ModelPainter::PrintLeftSideChainPad(
    TPad * pad, TH2 * hist, int residue, const std::string & y_title,
    const std::vector<std::string> & label_list)
{
    pad->cd();
    pad->Clear();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.04, 0.005, 0.03, 0.02);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 45.0, 1.2);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 40.0, 0.005, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 40.0, 0.01);
    ModifyAxisLabelSideChain(pad, hist, residue, label_list);
    hist->SetStats(0);
    hist->GetYaxis()->SetTitle(y_title.data());
    hist->GetYaxis()->CenterTitle();
    hist->Draw();
}

void ModelPainter::PrintRightSideChainPad(
    TPad * pad, TH2 * hist, int residue,
    const std::vector<std::string> & label_list)
{
    pad->cd();
    pad->Clear();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.04, 0.03, 0.02);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 40.0, 0.005, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 0.0);
    ModifyAxisLabelSideChain(pad, hist, residue, label_list);
    hist->SetStats(0);
    hist->Draw();
}

void ModelPainter::PrintTitleSideChainPad(
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