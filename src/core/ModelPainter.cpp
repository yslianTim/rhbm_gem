#include "ModelPainter.hpp"
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

#include <vector>
#include <tuple>

ModelPainter::ModelPainter(void) :
    m_folder_path{ "./" }, m_atom_classifier{ std::make_unique<AtomClassifier>() }
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
    m_ref_model_object_map[label].push_back(dynamic_cast<ModelObject *>(data_object));
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
        model_object->BuildKDTreeRoot();
        PaintResidueClassKNN(model_object, "residue_class_knn_"+ model_object->GetPdbID() +".pdf");
        PaintResidueClassXYPosition(model_object, "residue_class_xy_position_"+ model_object->GetPdbID() +".pdf");
        PaintResidueClassGroupGausScatter(model_object, "amplitude_scatter_"+ model_object->GetPdbID() +".pdf", 0);
        PaintResidueClassGroupGausScatter(model_object, "width_scatter_"+ model_object->GetPdbID() +".pdf", 1);
        PaintAtomGausScatter(model_object, "atom_gaus_scatter_"+ model_object->GetPdbID() +".pdf");
    }

    if (m_ref_model_object_map.find("buried_charge") != m_ref_model_object_map.end())
    {
        auto sim_buried_charge_model_object{ m_ref_model_object_map.at("buried_charge") };
        for (auto model_object : sim_buried_charge_model_object)
        {
            auto plot_main_chain_name{ "residue_class_group_gaus_main_"+ model_object->GetKeyTag() +".pdf" };
            PaintResidueClassGroupGausMainChain(model_object, plot_main_chain_name, true);
            //auto plot_side_chain_name{ "residue_class_group_gaus_side_"+ model_object->GetKeyTag() +".pdf" };
            //PaintResidueClassGroupGausSideChain(model_object, plot_side_chain_name);
        }
    }

    if (m_ref_model_object_map.find("no_charge") != m_ref_model_object_map.end())
    {
        auto sim_no_charge_model_object{ m_ref_model_object_map.at("no_charge") };
        for (auto model_object : sim_no_charge_model_object)
        {
            auto plot_main_chain_name{ "residue_class_group_gaus_main_"+ model_object->GetKeyTag() +".pdf" };
            PaintResidueClassGroupGausMainChain(model_object, plot_main_chain_name, true);
        }
    }

    if (m_model_object_list.size() == 4)
    {
        PaintResidueClassGroupGausMainChain("figure_1_part2.pdf");
        PaintResidueClassMapValue(m_model_object_list.at(3), "figure_1_part3.pdf");
    }

    if (m_model_object_list.size() == 11)
    {
        PaintElementClassGroupGausToFSC("figure6_a.pdf");
    }
}

void ModelPainter::PaintResidueClassGroupGausMainChain(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    auto residue_class{ AtomicInfoHelper::GetResidueClassKey() };
    std::cout <<"- ModelPainter::PaintResidueClassGroupGausMainChain"<< std::endl;

    #ifdef HAVE_ROOT

    const int primary_element_size{ 4 };
    short color_element[primary_element_size] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    short marker_element[primary_element_size]{ 25, 24, 26, 32 };
    float marker_size{ 1.5f };

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

    std::unique_ptr<TPaveText> side_chain_title_text[pad_extra_size];
    std::vector<std::unique_ptr<TGraphErrors>> amplitude_graph_list[2];
    std::vector<std::unique_ptr<TGraphErrors>> width_graph_list[2];
    std::vector<std::string> label_list[2];
    std::vector<double> amplitude_array, width_array;
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(m_model_object_list.at(3)) };
    Residue residue_list[2]{ Residue::ASN, Residue::GLN };
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
        auto pdb_text{ m_model_object_list.at(static_cast<size_t>(j))->GetPdbID() };
        auto emd_text{ m_model_object_list.at(static_cast<size_t>(j))->GetEmdID() };
        auto resolution{ m_model_object_list.at(static_cast<size_t>(j))->GetResolution() };
        PrintInfoMainChainPad(pad[0][j].get(), info_text[j].get(), pdb_text, emd_text, draw_axis_flag);
        PrintIconMainChainPad(pad[0][j].get(), icon_text[j].get(), resolution, draw_axis_flag);

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

void ModelPainter::PaintResidueClassGroupGausMainChain(
    ModelObject * model_object, const std::string & name, bool is_simulation)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintResidueClassGroupGausMainChain"<< std::endl;
    auto residue_class{ AtomicInfoHelper::GetResidueClassKey() };

    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    
    #ifdef HAVE_ROOT

    const int primary_element_size{ 4 };
    short color_element[primary_element_size] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    short marker_element[primary_element_size]{ 54, 53, 55, 59 };

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
    amplitude_array.reserve(static_cast<size_t>(primary_element_size * AtomicInfoHelper::GetStandardResidueCount()));
    width_array.reserve(static_cast<size_t>(primary_element_size * AtomicInfoHelper::GetStandardResidueCount()));
    for (size_t i = 0; i < primary_element_size; i++)
    {
        auto group_key_list{ m_atom_classifier->GetMainChainResidueClassGroupKeyList(i) };
        amplitude_graph[i] = entry_iter->CreateGausEstimateToResidueGraph(group_key_list, residue_class, 0);
        width_graph[i] = entry_iter->CreateGausEstimateToResidueGraph(group_key_list, residue_class, 1);
        correlation_graph[i] = entry_iter->CreateGausEstimateScatterGraph(group_key_list, residue_class);
        for (int p = 0; p < amplitude_graph[i]->GetN(); p++)
        {
            amplitude_array.push_back(amplitude_graph[i]->GetPointY(p));
            width_array.push_back(width_graph[i]->GetPointY(p));
        }
        ROOTHelper::SetMarkerAttribute(amplitude_graph[i].get(), marker_element[i], 1.3f, color_element[i]);
        ROOTHelper::SetMarkerAttribute(width_graph[i].get(), marker_element[i], 1.3f, color_element[i]);
        ROOTHelper::SetMarkerAttribute(correlation_graph[i].get(), marker_element[i], 1.3f, color_element[i]);
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
        ROOTHelper::SetFillAttribute(amplitude_hist[i].get(), 1001, color_element[i], 0.3f);
        ROOTHelper::SetFillAttribute(width_hist[i].get(), 1001, color_element[i], 0.3f);

        amplitude_hist[i]->SetBarWidth(0.6f);
        width_hist[i]->SetBarWidth(0.6f);
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
    if (is_simulation == true)
    {
        PrintInfoPad(pad[5].get(), info_text.get(), model_object->GetPdbID(), "Simulation");
    }
    else
    {
        PrintInfoPad(pad[5].get(), info_text.get(), model_object->GetPdbID(), model_object->GetEmdID());
    }

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
    auto residue_class{ AtomicInfoHelper::GetResidueClassKey() };

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

    size_t residue_size{ static_cast<size_t>(AtomicInfoHelper::GetStandardResidueCount()) };
    std::vector<std::vector<std::unique_ptr<TGraphErrors>>> amplitude_graph_list(residue_size);
    std::vector<std::vector<std::unique_ptr<TGraphErrors>>> width_graph_list(residue_size);
    std::vector<std::unique_ptr<TPaveText>> info_text(residue_size);
    for (auto residue : AtomicInfoHelper::GetStandardResidueList())
    {
        size_t residue_index{ static_cast<size_t>(residue) };
        info_text[residue_index] = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
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
                    auto group_key{ KeyPackerResidueClass::Pack(residue, element, remoteness, branch, false) };
                    if (entry_iter->IsAvailableGroupKey(group_key, residue_class) == false) continue;
                    auto gaus_estimate{ entry_iter->GetGausEstimatePrior(group_key, residue_class) };
                    auto gaus_variance{ entry_iter->GetGausVariancePrior(group_key, residue_class) };
                    amplitude_array.push_back(std::get<0>(gaus_estimate));
                    width_array.push_back(std::get<1>(gaus_estimate));
                    amplitude_graph->SetPoint(count_element, count_total, std::get<0>(gaus_estimate));
                    amplitude_graph->SetPointError(count_element, 0.0, std::get<0>(gaus_variance));
                    width_graph->SetPoint(count_element, count_total, std::get<1>(gaus_estimate));
                    width_graph->SetPointError(count_element, 0.0, std::get<1>(gaus_variance));
                    auto element_label{ AtomicInfoHelper::GetLabel(element) };
                    auto remoteness_label{ AtomicInfoHelper::GetLabel(remoteness) };
                    auto branch_label{ AtomicInfoHelper::GetLabel(branch) };
                    label_list.emplace_back(element_label +"_{"+ remoteness_label + branch_label+"}");
                    count_element++;
                    count_total++;
                }
            }
            auto color{ AtomicInfoHelper::GetDisplayColor(element) };
            ROOTHelper::SetMarkerAttribute(amplitude_graph.get(), kFullCircle, 1.5f, color);
            ROOTHelper::SetMarkerAttribute(width_graph.get(), kFullCircle, 1.5f, color);
            amplitude_graph_list[residue_index].push_back(std::move(amplitude_graph));
            width_graph_list[residue_index].push_back(std::move(width_graph));
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
        PrintInfoSideChainPad(pad[2].get(), info_text[residue_index].get(), AtomicInfoHelper::GetLabel(residue));

        pad[0]->cd();
        for (auto & graph : amplitude_graph_list[residue_index])
        {
            graph->Draw("P X0");
        }

        pad[1]->cd();
        for (auto & graph : width_graph_list[residue_index])
        {
            graph->Draw("P X0");
        }

        pad[2]->cd();
        info_text[residue_index]->Draw();

        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ModelPainter::PaintResidueClassGroupGausScatter(
    ModelObject * model_object, const std::string & name, int par_id)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintResidueClassGroupGausScatter"<< std::endl;
    auto residue_class{ AtomicInfoHelper::GetResidueClassKey() };

    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

    const int column_size{ 3 };
    const int row_size{ 3 };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), column_size, row_size, 0.10f, 0.02f, 0.08f, 0.02f, 0.02f, 0.02f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    
    const size_t col_element_index[column_size]{ 1, 2, 3 };
    const size_t row_element_index[row_size]{ 0, 3, 2 };
    const std::string element_label[4]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };
    const std::string title_label[2]
    {
        "Amplitude",
        "Width"
    };
    short color_element[4] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };

    std::unique_ptr<TH2> frame[column_size][row_size];
    std::vector<std::unique_ptr<TGraphErrors>> graph_list[column_size][row_size];
    std::vector<double> data_array;
    for (size_t i = 0; i < column_size; i++)
    {
        for (size_t j = 0; j < row_size; j++)
        {
            if (col_element_index[i] == row_element_index[j]) continue;
            if (i == 2 && j == 2) continue;
            for (auto residue : AtomicInfoHelper::GetStandardResidueList())
            {
                auto col_group_key{ m_atom_classifier->GetMainChainResidueClassGroupKeyList(col_element_index[i], residue).at(0) };
                if (entry_iter->IsAvailableGroupKey(col_group_key, residue_class) == false) continue;
                auto row_group_key{ m_atom_classifier->GetMainChainResidueClassGroupKeyList(row_element_index[j], residue).at(0) };
                if (entry_iter->IsAvailableGroupKey(row_group_key, residue_class) == false) continue;
                auto graph{ entry_iter->CreateGausEstimateScatterGraph(col_group_key, row_group_key, residue_class, par_id) };
                ROOTHelper::SetMarkerAttribute(graph.get(), AtomicInfoHelper::GetDisplayMarker(residue), 1.3f, AtomicInfoHelper::GetDisplayColor(residue));
                double x_min_tmp, x_max_tmp, y_min_tmp, y_max_tmp;
                graph->ComputeRange(x_min_tmp, y_min_tmp, x_max_tmp, y_max_tmp);
                data_array.emplace_back(x_min_tmp);
                data_array.emplace_back(x_max_tmp);
                data_array.emplace_back(y_min_tmp);
                data_array.emplace_back(y_max_tmp);
                graph_list[i][j].emplace_back(std::move(graph));
            }
        }
    }

    auto scaling{ 0.1 };
    auto data_range{ ArrayStats<double>::ComputeScalingRangeTuple(data_array, scaling) };
    auto min{ std::get<0>(data_range) };
    auto max{ std::get<1>(data_range) };
    auto line_ref{ ROOTHelper::CreateLine(min, min, max, max) };
    ROOTHelper::SetLineAttribute(line_ref.get(), 2, 1, kRed);
    std::unique_ptr<TPaveText> col_title_text[column_size];
    std::unique_ptr<TPaveText> row_title_text[row_size];
    std::unique_ptr<TLegend> legend;
    for (int i = 0; i < column_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, min, max, 500, min, max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 0);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 35, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 0);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->SetStats(0);
            
            bool has_graph{ graph_list[i][j].size() > 0 };
            if (has_graph == true) frame[i][j]->Draw();
            for (auto & graph : graph_list[i][j])
            {
                graph->Draw("P X0");
            }
            if (has_graph == true) line_ref->Draw();

            if (j == 0)
            {
                col_title_text[i] = ROOTHelper::CreatePaveText(0.01, -0.27, 0.99, -0.13, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(col_title_text[i].get());
                ROOTHelper::SetPaveAttribute(col_title_text[i].get(), 0, 0.2);
                ROOTHelper::SetTextAttribute(col_title_text[i].get(), 35, 133, 22);
                ROOTHelper::SetFillAttribute(col_title_text[i].get(), 1001, color_element[col_element_index[i]], 0.5f);
                col_title_text[i]->AddText(element_label[col_element_index[i]].data());
                col_title_text[i]->Draw();
            }
            if (i == 0)
            {
                row_title_text[j] = ROOTHelper::CreatePaveText(-0.35, 0.01, -0.20, 0.99, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(row_title_text[j].get());
                ROOTHelper::SetPaveAttribute(row_title_text[j].get(), 0, 0.1);
                ROOTHelper::SetTextAttribute(row_title_text[j].get(), 35, 133, 22);
                ROOTHelper::SetFillAttribute(row_title_text[j].get(), 1001, color_element[row_element_index[j]], 0.5f);
                row_title_text[j]->AddText(element_label[row_element_index[j]].data());
                row_title_text[j]->Draw();
                auto text{ row_title_text[j]->GetLineWith(element_label[row_element_index[j]].data()) };
                text->SetTextAngle(90.0f);
            }
            if (i == column_size - 1 && j == 1)
            {
                legend = ROOTHelper::CreateLegend(0.00, 0.00, 1.00, 1.00, true);
                ROOTHelper::SetLegendDefaultStyle(legend.get());
                ROOTHelper::SetTextAttribute(legend.get(), 25, 133, 12);
                ROOTHelper::SetFillAttribute(legend.get(), 4000);
                for (auto residue : AtomicInfoHelper::GetStandardResidueList())
                {
                    auto residue_id{ static_cast<size_t>(residue) };
                    auto atom_size{ std::to_string(graph_list[0][0].at(residue_id)->GetN()) };
                    std::string label{ "#font[102]{" + AtomicInfoHelper::GetLabel(residue) + " }("+ atom_size +")" };
                    legend->AddEntry(graph_list[0][0].at(residue_id).get(), label.data(), "p");
                }
                legend->SetNColumns(2);
                legend->Draw();
            }
        }
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ModelPainter::PaintResidueClassMapValue(ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintResidueClassMapValue"<<std::endl;
    auto residue_class{ AtomicInfoHelper::GetResidueClassKey() };

    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    
    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    auto frame{ ROOTHelper::CreateHist2D("hist","", 100, 0.0, 1.0, 100, 0.0, 1.0) };
    for (auto residue : AtomicInfoHelper::GetStandardResidueList())
    {
        for (auto element : AtomicInfoHelper::GetStandardElementList())
        {
            for (auto remoteness : AtomicInfoHelper::GetStandardRemotenessList())
            {
                for (auto branch : AtomicInfoHelper::GetStandardBranchList())
                {
                    auto group_key{ KeyPackerResidueClass::Pack(residue, element, remoteness, branch, false) };
                    if (entry_iter->IsAvailableGroupKey(group_key, residue_class) == false) continue;
                    std::vector<std::unique_ptr<TGraphErrors>> map_value_graph_list;
                    double y_min{ 0.0 };
                    double y_max{ 1.0 };
                    for (auto atom : entry_iter->GetAtomObjectList(group_key, residue_class))
                    {
                        auto atom_iter{ std::make_unique<PotentialEntryIterator>(atom) };
                        auto map_value_graph{ atom_iter->CreateBinnedDistanceToMapValueGraph() };
                        map_value_graph_list.emplace_back(std::move(map_value_graph));

                        auto map_value_range{ atom_iter->GetMapValueRange(0.2) };
                        if (std::get<0>(map_value_range) < y_min) y_min = std::get<0>(map_value_range);
                        if (std::get<1>(map_value_range) > y_max) y_max = std::get<1>(map_value_range);
                    }

                    ROOTHelper::SetPadDefaultStyle(gPad);
                    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
                    ROOTHelper::SetPadMarginInCanvas(gPad, 0.16, 0.02, 0.14, 0.02);
                    ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 60.0f, 1.1f);
                    ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 60.0f, 1.3f);
                    ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), 60.0f, 0.01f);
                    ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 60.0f, 0.02f);

                    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, kWhite, 1, 0);

                    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.04, 0) };
                    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.05, 1) };
                    ROOTHelper::SetAxisTickAttribute(frame->GetXaxis(), static_cast<float>(x_tick_length));
                    ROOTHelper::SetAxisTickAttribute(frame->GetYaxis(), static_cast<float>(y_tick_length));
                    frame->GetXaxis()->SetLimits(0.0, 1.5);
                    frame->GetYaxis()->SetLimits(y_min, y_max);
                    frame->SetStats(0);
                    frame->GetXaxis()->SetTitle("Radial Distance #[]{#AA}");
                    frame->GetYaxis()->SetTitle("Normalized Map Value");
                    frame->Draw();

                    auto ref_line{ ROOTHelper::CreateLine(0.0, 0.0, 1.5, 0.0) };
                    ROOTHelper::SetLineAttribute(ref_line.get(), 2, 1);
                    ref_line->Draw();

                    for (auto & graph : map_value_graph_list)
                    {
                        ROOTHelper::SetLineAttribute(graph.get(), 1, 5, kAzure-7, 0.3f);
                        graph->Draw("LX");
                    }

                    auto gaus_function{ entry_iter->CreateGroupGausFunctionPrior(group_key, residue_class) };
                    ROOTHelper::SetLineAttribute(gaus_function.get(), 9, 5, kRed+1);
                    gaus_function->Draw("SAME");

                    auto legend{ ROOTHelper::CreateLegend(0.0, 0.0, 1.0, 1.0, false) };
                    ROOTHelper::SetLegendMarginInCanvas(gPad, legend.get(), 0.20, 0.02, 0.85, 0.02);
                    ROOTHelper::SetLegendDefaultStyle(legend.get());
                    ROOTHelper::SetFillAttribute(legend.get(), 4000);
                    ROOTHelper::SetTextAttribute(legend.get(), 55, 133, 12, 0.0, kAzure-7);
                    legend->SetMargin(0.15f);
                    legend->AddEntry(gaus_function.get(), "Single Gaus Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})}", "l");
                    legend->AddEntry(map_value_graph_list.at(0).get(), "Atomic Map Value Distribution", "l");
                    legend->Draw();

                    auto result_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC", false) };
                    ROOTHelper::SetPaveTextMarginInCanvas(gPad, result_text.get(), 0.20, 0.02, 0.20, 0.70);
                    ROOTHelper::SetPaveTextDefaultStyle(result_text.get());
                    ROOTHelper::SetPaveAttribute(result_text.get(), 0);
                    ROOTHelper::SetFillAttribute(result_text.get(), 4000);
                    ROOTHelper::SetTextAttribute(result_text.get(), 85, 133, 22, 0.0, kAzure-7);
                    auto gaus_prior{ entry_iter->GetGausEstimatePrior(group_key, residue_class) };
                    auto amplitude{ std::get<0>(gaus_prior) };
                    auto width{ std::get<1>(gaus_prior) };
                    result_text->AddText(Form("#color[633]{#font[1]{A} = %.2f  ;  #font[1]{#tau} = %.2f}", amplitude, width));
/*
                    auto residue_label{ AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::ResidueTag>(residue) };
                    auto element_label{ AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::ElementTag>(element) };
                    auto remoteness_label{ AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::RemotenessTag>(remoteness) };
                    auto branch_label{ AtomicInfoHelper::AtomLabelMapping<AtomicInfoHelper::BranchTag>(branch) };
                    auto info_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
                    ROOTHelper::SetPaveTextMarginInCanvas(gPad, info_text.get(), 0.5, 0.01, 0.85, 0.01);
                    ROOTHelper::SetPaveTextDefaultStyle(info_text.get());
                    ROOTHelper::SetPaveAttribute(info_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(info_text.get(), 1001, kAzure-7);
                    ROOTHelper::SetTextAttribute(info_text.get(), 50, 103, 22, 0.0, kYellow-10);
                    info_text->AddText(("#font[102]{"+ residue_label + " " + element_label + " " + remoteness_label + branch_label + "}").data());
                    //info_text->Draw();
*/

                    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
                }
            }
        }
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ModelPainter::PaintResidueClassKNN(ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintResidueClassKNN"<< std::endl;
    auto residue_class{ AtomicInfoHelper::GetResidueClassKey() };

    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), 4, 2, 0.1f, 0.01f, 0.15f, 0.01f, 0.01f, 0.005f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int column_size{ 4 };
    const int row_size{ 2 };
    const std::string element_label[column_size]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };
    const std::string y_axis_title[row_size]
    {
        "Width #hat{#sigma}",
        "Amplitude #hat{A}"
    };

    short color_element[column_size] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    short marker_element[column_size]{ 54, 53, 55, 59 };
    std::unique_ptr<TH2> frame[column_size][row_size];

    std::vector<std::unique_ptr<TGraphErrors>> graph_list[column_size][row_size];
    double x_min[column_size]{ 0.0 };
    double x_max[column_size]{ 1.0 };
    double y_min[row_size]{ 0.0 };
    double y_max[row_size]{ 1.0 };
    for (size_t i = 0; i < column_size; i++)
    {
        for (auto residue : AtomicInfoHelper::GetStandardResidueList())
        {
            auto group_key{ m_atom_classifier->GetMainChainResidueClassGroupKeyList(i, residue).at(0) };
            if (entry_iter->IsAvailableGroupKey(group_key, residue_class) == false) continue;
            auto amplitude_graph{ entry_iter->CreateInRangeAtomsToGausEstimateGraph(group_key, residue_class, 5.0, 0) };
            auto width_graph{ entry_iter->CreateInRangeAtomsToGausEstimateGraph(group_key, residue_class, 5.0, 1) };
            double knn_min, knn_max, amplitude_min, amplitude_max;
            amplitude_graph->ComputeRange(knn_min, amplitude_min, knn_max, amplitude_max);
            x_min[i] = (knn_min < x_min[i]) ? knn_min : x_min[i];
            x_max[i] = (knn_max > x_max[i]) ? knn_max : x_max[i];
            y_min[1] = (amplitude_min < y_min[1]) ? amplitude_min : y_min[1];
            y_max[1] = (amplitude_max > y_max[1]) ? amplitude_max : y_max[1];
            double width_min, width_max;
            width_graph->ComputeRange(knn_min, width_min, knn_max, width_max);
            x_min[i] = (knn_min < x_min[i]) ? knn_min : x_min[i];
            x_max[i] = (knn_max > x_max[i]) ? knn_max : x_max[i];
            y_min[0] = (width_min < y_min[0]) ? width_min : y_min[0];
            y_max[0] = (width_max > y_max[0]) ? width_max : y_max[0];

            graph_list[i][1].emplace_back(std::move(amplitude_graph));
            graph_list[i][0].emplace_back(std::move(width_graph));
        }
    }

    for (int i = 0; i < column_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min[j], y_max[j]);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 35, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 35, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 35, 1.6f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("#In Range atoms");
            frame[i][j]->GetYaxis()->SetTitle(y_axis_title[j].data());
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            for (auto & graph : graph_list[i][j])
            {
                ROOTHelper::SetMarkerAttribute(graph.get(), marker_element[i], 1.3f, color_element[i]);
                graph->Draw("P X0");
            }
        }
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ModelPainter::PaintResidueClassXYPosition(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintResidueClassXYPosition"<< std::endl;
    auto residue_class{ AtomicInfoHelper::GetResidueClassKey() };

    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    auto x_range_tuple{ model_object->GetModelPositionRange(0, 0.1) };
    auto y_range_tuple{ model_object->GetModelPositionRange(1, 0.1) };
    auto x_min{ std::get<0>(x_range_tuple) };
    auto x_max{ std::get<1>(x_range_tuple) };
    auto y_min{ std::get<0>(y_range_tuple) };
    auto y_max{ std::get<1>(y_range_tuple) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1200, 900) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), 4, 3, 0.1f, 0.05f, 0.15f, 0.01f, 0.01f, 0.005f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int column_size{ 4 };
    const int row_size{ 3 };
    const std::string element_label[column_size]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };
    short color_element[column_size]   { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    short marker_element[column_size]  { 54, 53, 55, 59 };
    std::unique_ptr<TH2> frame[column_size][row_size];

    std::unique_ptr<TGraphErrors> position_graph[column_size][row_size];
    std::unique_ptr<TGraph2DErrors> amplitude_2d_graph[column_size][row_size];
    std::unique_ptr<TGraph2DErrors> width_2d_graph[column_size][row_size];
    std::unique_ptr<TGraphErrors> total_graph[row_size];
    for (size_t i = 0; i < column_size; i++)
    {
        auto group_key_list{ m_atom_classifier->GetMainChainResidueClassGroupKeyList(i) };
        position_graph[i][0] = entry_iter->CreateXYPositionTomographyGraph(group_key_list, residue_class, 0.3, 0.1) ;
        position_graph[i][1] = entry_iter->CreateXYPositionTomographyGraph(group_key_list, residue_class, 0.4, 0.1);
        position_graph[i][2] = entry_iter->CreateXYPositionTomographyGraph(group_key_list, residue_class, 0.5, 0.1);
        amplitude_2d_graph[i][0] = entry_iter->CreateXYPositionTomographyToGausEstimateGraph2D(group_key_list, residue_class, 0.3, 0.1, 0);
        amplitude_2d_graph[i][1] = entry_iter->CreateXYPositionTomographyToGausEstimateGraph2D(group_key_list, residue_class, 0.4, 0.1, 0);
        amplitude_2d_graph[i][2] = entry_iter->CreateXYPositionTomographyToGausEstimateGraph2D(group_key_list, residue_class, 0.5, 0.1, 0);
        width_2d_graph[i][0] = entry_iter->CreateXYPositionTomographyToGausEstimateGraph2D(group_key_list, residue_class, 0.3, 0.1, 1);
        width_2d_graph[i][1] = entry_iter->CreateXYPositionTomographyToGausEstimateGraph2D(group_key_list, residue_class, 0.4, 0.1, 1);
        width_2d_graph[i][2] = entry_iter->CreateXYPositionTomographyToGausEstimateGraph2D(group_key_list, residue_class, 0.5, 0.1, 1);
    }
    total_graph[0] = entry_iter->CreateXYPositionTomographyGraph(0.3, 0.1);
    total_graph[1] = entry_iter->CreateXYPositionTomographyGraph(0.4, 0.1);
    total_graph[2] = entry_iter->CreateXYPositionTomographyGraph(0.5, 0.1);
    for (int j = 0; j < row_size; j++)
    {
        ROOTHelper::SetMarkerAttribute(total_graph[j].get(), 54, 1.3f, kGray, 0.1f);
    }

    for (int i = 0; i < column_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min, x_max, 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 35, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 35, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 35, 1.6f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("Position X");
            frame[i][j]->GetYaxis()->SetTitle("Position Y");
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            total_graph[j]->Draw("P X0");
            ROOTHelper::SetMarkerAttribute(position_graph[i][j].get(), marker_element[i], 1.3f, color_element[i]);
            position_graph[i][j]->Draw("P X0");
        }
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    
    for (int i = 0; i < column_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            gPad->SetTheta(90.0);
            gPad->SetPhi(0.00);
            ROOTHelper::SetMarkerAttribute(amplitude_2d_graph[i][j].get(), 20, 1.3f, color_element[i]);
            amplitude_2d_graph[i][j]->SetTitle("");
            amplitude_2d_graph[i][j]->Draw("PCOLZ");
        }
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    for (int i = 0; i < column_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            gPad->SetTheta(90.0);
            gPad->SetPhi(0.00);
            ROOTHelper::SetMarkerAttribute(width_2d_graph[i][j].get(), 20, 1.3f, color_element[i]);
            width_2d_graph[i][j]->SetTitle("");
            width_2d_graph[i][j]->Draw("PCOLZ");
        }
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ModelPainter::PaintElementClassGroupGausToFSC(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintElementClassGroupGausToFSC"<< std::endl;

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int column_size{ 4 };
    const int row_size{ 1 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 950, 300) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), column_size, row_size, 0.1f, 0.01f, 0.11f, 0.07f, 0.01f, 0.005f);
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
    std::vector<double> x_array, y_array;
    for (size_t i = 0; i < primary_element_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        graph[i] = ROOTHelper::CreateGraphErrors();
        //BuildGausEstimateToFSCGraph(group_key, graph[i].get(), m_model_object_list);

        for (int p = 0; p < graph[i]->GetN(); p++)
        {
            x_array.push_back(graph[i]->GetPointX(p));
            y_array.push_back(graph[i]->GetPointY(p));
        }

        ROOTHelper::SetMarkerAttribute(graph[i].get(), 20, 1.2f, color_element[i]);
        ROOTHelper::SetLineAttribute(graph[i].get(), 1, 2, color_element[i]);
    }

    auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.05) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.3) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };

    std::unique_ptr<TH2> frame[column_size][row_size];
    std::unique_ptr<TPaveText> title_text[column_size];
    for (int i = 0; i < column_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min, x_max, 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 30, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 30, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 35, 1.3f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("FSC");
            frame[i][j]->GetYaxis()->SetTitle("Width");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            graph[i]->Draw("PL X0");
        }
        title_text[i] = ROOTHelper::CreatePaveText(0.01, 1.01, 0.99, 1.16, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
        ROOTHelper::SetPaveAttribute(title_text[i].get(), 0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[i].get(), 25, 133, 22);
        ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, color_element[i], 0.5f);
        title_text[i]->AddText(element_label[i]);
        title_text[i]->Draw();
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ModelPainter::PaintAtomGausScatter(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintAtomGausScatter"<< std::endl;
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    
    const int element_size{ 3 };
    Element element_list[element_size]{ Element::CARBON, Element::NITROGEN, Element::OXYGEN };

    std::unique_ptr<TH2> frame;
    std::vector<std::unique_ptr<TGraphErrors>> graph_list;
    graph_list.reserve(element_size);
    std::vector<double> x_array, y_array;
    x_array.reserve(model_object->GetNumberOfSelectedAtom());
    y_array.reserve(model_object->GetNumberOfSelectedAtom());
    int amplitude_cut_count[element_size]{ 0 };
    int width_cut_count[element_size]{ 0 };
    int atom_count[element_size]{ 0 };
    for (size_t i = 0; i < element_size; i++)
    {
        auto graph{ entry_iter->CreateGausEstimateScatterGraph(element_list[i]) };
        ROOTHelper::SetMarkerAttribute(graph.get(),
            AtomicInfoHelper::GetDisplayMarker(element_list[i]), 1.2f,
            AtomicInfoHelper::GetDisplayColor(element_list[i]));
        for (int p = 0; p < graph->GetN(); p++)
        {
            auto x{ graph->GetPointX(p) };
            auto y{ graph->GetPointY(p) };
            x_array.emplace_back(x);
            y_array.emplace_back(y);
            if (x <= 50.0) amplitude_cut_count[i]++;
            if (y <= 0.6) width_cut_count[i]++;
        }
        atom_count[i] = graph->GetN();
        std::cout << "Ratio of amplitude selected atoms: "<< amplitude_cut_count[i] <<" / "<< atom_count[i] << std::endl;
        std::cout << "Ratio of width selected atoms: "<< width_cut_count[i] <<" / "<< atom_count[i] << std::endl;
        graph_list.emplace_back(std::move(graph));
    }

    auto x_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(x_array, 0.1) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array, 0.1) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };

    frame = ROOTHelper::CreateHist2D("frame","", 500, x_min, x_max, 500, y_min, y_max);
    frame->SetStats(0);
    frame->Draw();
    for (auto & graph : graph_list)
    {
        
        graph->Draw("P X0");
    }

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
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
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, 0.005, 0.07, 0.01, 0.02);
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
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0f, 1.4f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0f, 0.01f);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 21);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 506);
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
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0, 1.4f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 32.0, 0.11f, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0, 0.01f);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
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
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 5);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 506);
    hist->GetXaxis()->SetLimits(-1.0, 4.0);
    hist->SetStats(0);
    hist->Draw();
}

void ModelPainter::PrintWidthSummaryPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.005, 0.13, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 35.0f, 1.0f);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 35.0f, 0.005f, 103);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 0.0f);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
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

void ModelPainter::PrintGausSummaryPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.07, 0.13, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 35.0f, 1.0f);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0f, 1.2f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 30.0f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0f, 0.01f);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.03, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 505);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
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
    TPad * pad, TH2 * hist, Residue residue, const std::vector<std::string> & label_list)
{
    pad->cd();
    pad->Clear();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.07, 0.005, 0.16, 0.04);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0, 1.2f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 32.0, 0.005f, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0, 0.01f);
    ModifyAxisLabelSideChain(pad, hist, residue, label_list);
    hist->SetStats(0);
    hist->GetYaxis()->SetTitle("Amplitude");
    hist->GetYaxis()->CenterTitle();
    hist->Draw();
}

void ModelPainter::PrintWidthSideChainPad(
    TPad * pad, TH2 * hist, Residue residue, const std::vector<std::string> & label_list)
{
    pad->cd();
    pad->Clear();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.07, 0.16, 0.04);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0, 1.2f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 32.0, 0.005f, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 30.0, 0.01f);
    ModifyAxisLabelSideChain(pad, hist, residue, label_list);
    hist->SetStats(0);
    hist->GetYaxis()->SetTitle("Width");
    hist->GetYaxis()->CenterTitle();
    hist->Draw("Y+");
}

void ModelPainter::ModifyAxisLabelSideChain(
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

void ModelPainter::PrintIconMainChainPad(
    TPad * pad, TPaveText * text, double resolution, bool is_bottom_pad)
{
    pad->cd();
    auto bottom_margin{ (is_bottom_pad) ? 0.16 : 0.12 };
    auto top_margin{ 0.01 };
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, 0.005, 0.005, bottom_margin, top_margin);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.2);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(text, 80, 133, 22, 0.0, kYellow-10);
    text->AddText(Form("%.2f #AA", resolution));
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
    auto title_size{ (draw_x_axis) ? 40.0f : 0.0f };
    auto label_size{ (draw_x_axis) ? 40.0f : 0.0f };
    ROOTHelper::SetPadMarginInCanvas(pad, 0.025, 0.005, bottom_margin, 0.01);
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

void ModelPainter::PrintSummaryPad(TPad * pad, TH2 * hist, bool draw_x_axis)
{
    pad->cd();
    auto bottom_margin{ (draw_x_axis) ? 0.045 : 0.005 };
    auto title_size{ (draw_x_axis) ? 45.0f : 0.0f };
    auto label_size{ (draw_x_axis) ? 40.0f : 0.0f };
    ROOTHelper::SetPadMarginInCanvas(pad, 0.001, 0.001, bottom_margin, 0.01);
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

void ModelPainter::PrintCorrelationPad(TPad * pad, TH2 * hist, bool draw_x_axis)
{
    pad->cd();
    auto bottom_margin{ (draw_x_axis) ? 0.045 : 0.005 };
    auto title_size{ (draw_x_axis) ? 45.0f : 0.0f };
    auto label_size{ (draw_x_axis) ? 40.0f : 0.0f };
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.04, bottom_margin, 0.01);
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

void ModelPainter::PrintLeftSideChainPad(
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

void ModelPainter::PrintRightSideChainPad(
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