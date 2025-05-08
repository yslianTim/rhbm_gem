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
#include <TH1.h>
#endif

#include <iostream>
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
    
    if (m_model_object_list.size() > 1) return;

    for (auto model_object : m_model_object_list)
    {
        auto plot_main_chain_name{ "structure_class_group_gaus_main_"+ model_object->GetPdbID() +".pdf" };
        PaintStructureClassGroupGausMainChain(model_object, plot_main_chain_name);
        auto plot_side_chain_name{ "structure_class_group_gaus_side_"+ model_object->GetPdbID() +".pdf" };
        PaintStructureClassGroupGausSideChain(model_object, plot_side_chain_name);
        auto plot_map_value_name{ "residue_class_map_value_"+ model_object->GetPdbID() +".pdf" };
        PaintResidueClassMapValue(model_object, plot_map_value_name);
        model_object->BuildKDTreeRoot();
        PaintResidueClassWidthScatterPlot(model_object, "residue_class_com_"+ model_object->GetPdbID() +".pdf", 0, true);
        //PaintResidueClassWidthScatterPlot(model_object, "residue_class_knn_"+ model_object->GetPdbID() +".pdf", 1, true);
        PaintResidueClassXYPosition(model_object, "residue_class_xy_position_"+ model_object->GetPdbID() +".pdf");
        PaintResidueClassGroupGausScatter(model_object, "amplitude_scatter_"+ model_object->GetPdbID() +".pdf", 0);
        PaintResidueClassGroupGausScatter(model_object, "width_scatter_"+ model_object->GetPdbID() +".pdf", 1);
        PaintAtomXYPosition(model_object, "atom_position_"+ model_object->GetPdbID() +".pdf");
        PaintAtomGausScatter(model_object, "atom_gaus_scatter_"+ model_object->GetPdbID() +".pdf", false);
        PaintAtomGausMainChain(model_object, "atom_gaus_main_chain_"+ model_object->GetPdbID() +".pdf");
    }

    if (m_ref_model_object_map.find("with_charge") != m_ref_model_object_map.end())
    {
        auto sim_model_object_list{ m_ref_model_object_map.at("with_charge") };
        if (sim_model_object_list.size() == 1)
        {
            auto plot_main_chain_name{ "residue_class_group_gaus_main_simulation.pdf" };
            PaintStructureClassGroupGausMainChain(sim_model_object_list.at(0), plot_main_chain_name, true);
        }
    }
}

void ModelPainter::PaintStructureClassGroupGausMainChain(
    ModelObject * model_object, const std::string & name, bool is_simulation)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintStructureClassGroupGausMainChain"<< std::endl;

    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    
    #ifdef HAVE_ROOT

    const int primary_element_size{ 4 };
    const int structure_size{ 3 };
    short color_element[structure_size][primary_element_size]
    {
        { kRed+1, kOrange+1, kGreen+2, kAzure+2 },
        { kGray, kGray, kGray, kGray },
        { kGray+1, kGray+1, kGray+1, kGray+1 }
    };
    short marker_element[structure_size][primary_element_size]
    {
        { 21, 20, 22, 23 },
        { 54, 53, 55, 59 },
        { 54, 53, 55, 59 }
    };
    const std::string class_key_list[structure_size]
    {
        AtomicInfoHelper::GetResidueClassKey(),
        AtomicInfoHelper::GetStructureClassKey(),
        AtomicInfoHelper::GetStructureClassKey()
    };
    const Structure structure_list[structure_size]
    {
        Structure::FREE,
        Structure::FREE,
        Structure::HELX_P
    };

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

    std::vector<uint64_t> group_key_list[primary_element_size][structure_size];
    std::unique_ptr<TGraphErrors> amplitude_graph[primary_element_size][structure_size];
    std::unique_ptr<TGraphErrors> width_graph[primary_element_size][structure_size];
    std::unique_ptr<TGraphErrors> correlation_graph[primary_element_size][structure_size];
    std::unique_ptr<TH2D> amplitude_hist[primary_element_size][structure_size];
    std::unique_ptr<TH2D> width_hist[primary_element_size][structure_size];
    std::unique_ptr<TH1D> count_hist[structure_size];
    std::vector<double> amplitude_array, width_array;
    amplitude_array.reserve(80);
    width_array.reserve(80);
    for (size_t i = 0; i < primary_element_size; i++)
    {
        group_key_list[i][0] = m_atom_classifier->GetMainChainResidueClassGroupKeyList(i);
        group_key_list[i][1] = m_atom_classifier->GetMainChainStructureClassGroupKeyList(i, Structure::FREE);
        group_key_list[i][2] = m_atom_classifier->GetMainChainStructureClassGroupKeyList(i, Structure::HELX_P);
        for (size_t j = 0; j < structure_size; j++)
        {
            amplitude_graph[i][j] = entry_iter->CreateGausEstimateToResidueGraph(group_key_list[i][j], class_key_list[j], 0);
            width_graph[i][j] = entry_iter->CreateGausEstimateToResidueGraph(group_key_list[i][j], class_key_list[j], 1);
            correlation_graph[i][j] = entry_iter->CreateGausEstimateScatterGraph(group_key_list[i][j], class_key_list[j], 0, 1);
            for (int p = 0; p < amplitude_graph[i][j]->GetN(); p++)
            {
                if (j != 0) break;
                amplitude_array.push_back(amplitude_graph[i][j]->GetPointY(p));
                width_array.push_back(width_graph[i][j]->GetPointY(p));
            }
            ROOTHelper::SetMarkerAttribute(amplitude_graph[i][j].get(), marker_element[j][i], 1.3f, color_element[j][i]);
            ROOTHelper::SetMarkerAttribute(width_graph[i][j].get(), marker_element[j][i], 1.3f, color_element[j][i]);
            ROOTHelper::SetMarkerAttribute(correlation_graph[i][j].get(), marker_element[j][i], 1.3f, color_element[j][i]);
            ROOTHelper::SetLineAttribute(amplitude_graph[i][j].get(), 1, 1, color_element[j][i]);
            ROOTHelper::SetLineAttribute(width_graph[i][j].get(), 1, 1, color_element[j][i]);
            ROOTHelper::SetLineAttribute(correlation_graph[i][j].get(), 1, 1, color_element[j][i]);
        }
    }

    auto scaling{ 0.3 };
    auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, scaling) };
    auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, scaling) };
    float width_size_list[structure_size]{ 0.3f, 0.3f, 0.3f };
    float offset_list[structure_size]{ -0.2f, 0.0f, 0.2f };
    for (int i = 0; i < primary_element_size; i++)
    {
        for (size_t j = 0; j < structure_size; j++)
        {
            std::string name_amplitude{ "amplitude_hist_"+ std::to_string(i) + std::to_string(j) };
            std::string name_width{ "width_hist_"+ std::to_string(i) + std::to_string(j) };
            amplitude_hist[i][j] = ROOTHelper::CreateHist2D(name_amplitude.data(),"", 4, -0.5, 3.5, 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
            width_hist[i][j] = ROOTHelper::CreateHist2D(name_width.data(),"", 4, -0.5, 3.5, 100, std::get<0>(width_range), std::get<1>(width_range));
            for (int p = 0; p < amplitude_graph[i][j]->GetN(); p++)
            {
                amplitude_hist[i][j]->Fill(i, amplitude_graph[i][j]->GetPointY(p));
            }
            for (int p = 0; p < width_graph[i][j]->GetN(); p++)
            {
                width_hist[i][j]->Fill(i, width_graph[i][j]->GetPointY(p));
            }
            ROOTHelper::SetLineAttribute(amplitude_hist[i][j].get(), 1, 1, color_element[j][i]);
            ROOTHelper::SetLineAttribute(width_hist[i][j].get(), 1, 1, color_element[j][i]);
            ROOTHelper::SetFillAttribute(amplitude_hist[i][j].get(), 1001, color_element[j][i], 0.3f);
            ROOTHelper::SetFillAttribute(width_hist[i][j].get(), 1001, color_element[j][i], 0.3f);

            amplitude_hist[i][j]->SetBarWidth(width_size_list[j]);
            width_hist[i][j]->SetBarWidth(width_size_list[j]);

            amplitude_hist[i][j]->SetBarOffset(offset_list[j]);
            width_hist[i][j]->SetBarOffset(offset_list[j]);
        }
    }
    auto max_count{ 0 };
    for (size_t j = 0; j < structure_size; j++)
    {
        std::string count_name{ "count_hist_"+ std::to_string(j) };
        count_hist[j] = entry_iter->CreateResidueCountHistogram(class_key_list[j], structure_list[j]);
        count_hist[j]->SetName(count_name.data());
        if (count_hist[j]->GetMaximum() > max_count)
        {
            max_count = static_cast<int>(count_hist[j]->GetMaximum());
        }
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
    frame[5] = ROOTHelper::CreateHist2D("hist_5","", 100, 0.0, 1.0, 500, 0.0, max_count*1.1);
    auto info_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    auto resolution_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };

    PrintAmplitudePad(pad[1].get(), frame[1].get());
    PrintWidthPad(pad[0].get(), frame[0].get());
    PrintAmplitudeSummaryPad(pad[3].get(), frame[3].get());
    PrintWidthSummaryPad(pad[2].get(), frame[2].get());
    PrintGausSummaryPad(pad[4].get(), frame[4].get());
    PrintCountSummaryPad(pad[6].get(), frame[5].get());
    auto emd_id{ (is_simulation == true) ? "Simulation" : model_object->GetEmdID() };
    PrintDataInfoPad(pad[5].get(), info_text.get(), model_object->GetPdbID(), emd_id);
    PrintResolutionInfoPad(pad[5].get(), resolution_text.get(), model_object->GetResolution());

    pad[1]->cd();
    for (int i = 0; i < primary_element_size; i++)
    {
        amplitude_graph[i][0]->Draw("PL X0");
        //amplitude_graph[i][1]->Draw("P X0");
        amplitude_graph[i][2]->Draw("P X0");
    }

    pad[0]->cd();
    for (int i = 0; i < primary_element_size; i++)
    {
        width_graph[i][0]->Draw("PL X0");
        //width_graph[i][1]->Draw("P X0");
        width_graph[i][2]->Draw("P X0");
    }

    pad[3]->cd();
    for (int i = 0; i < primary_element_size; i++)
    {
        amplitude_hist[i][0]->Draw("CANDLE3 SAME");
        //amplitude_hist[i][1]->Draw("CANDLE3 SAME");
        amplitude_hist[i][2]->Draw("CANDLE3 SAME");
    }

    pad[2]->cd();
    for (int i = 0; i < primary_element_size; i++)
    {
        width_hist[i][0]->Draw("CANDLE3 SAME");
        //width_hist[i][1]->Draw("CANDLE3 SAME");
        width_hist[i][2]->Draw("CANDLE3 SAME");
    }

    pad[4]->cd();
    for (int i = 0; i < primary_element_size; i++)
    {
        correlation_graph[i][0]->Draw("P X0");
        //correlation_graph[i][1]->Draw("P X0");
        correlation_graph[i][2]->Draw("P X0");
    }

    pad[6]->cd();
    gStyle->SetTextFont(132);
    ROOTHelper::SetFillAttribute(count_hist[0].get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetLineAttribute(count_hist[0].get(), 1, 1, kAzure-7);
    ROOTHelper::SetMarkerAttribute(count_hist[0].get(), 20, 7.0f, kAzure);
    count_hist[0]->Draw("HIST TEXT0 SAME");

    ROOTHelper::SetFillAttribute(count_hist[2].get(), 1001, kGray);
    ROOTHelper::SetLineAttribute(count_hist[2].get(), 1, 1, kGray);
    count_hist[2]->Draw("HIST SAME");

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ModelPainter::PaintStructureClassGroupGausSideChain(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintStructureClassGroupGausSideChain"<<std::endl;

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
    std::vector<std::vector<std::unique_ptr<TGraphErrors>>> amplitude_mix_graph_list(residue_size);
    std::vector<std::vector<std::unique_ptr<TGraphErrors>>> amplitude_free_graph_list(residue_size);
    std::vector<std::vector<std::unique_ptr<TGraphErrors>>> amplitude_helix_graph_list(residue_size);
    std::vector<std::vector<std::unique_ptr<TGraphErrors>>> width_mix_graph_list(residue_size);
    std::vector<std::vector<std::unique_ptr<TGraphErrors>>> width_free_graph_list(residue_size);
    std::vector<std::vector<std::unique_ptr<TGraphErrors>>> width_helix_graph_list(residue_size);
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
            auto amplitude_mix_graph{ ROOTHelper::CreateGraphErrors() };
            auto amplitude_free_graph{ ROOTHelper::CreateGraphErrors() };
            auto amplitude_helix_graph{ ROOTHelper::CreateGraphErrors() };
            auto width_mix_graph{ ROOTHelper::CreateGraphErrors() };
            auto width_free_graph{ ROOTHelper::CreateGraphErrors() };
            auto width_helix_graph{ ROOTHelper::CreateGraphErrors() };
            auto count_mix{ 0 };
            auto count_free{ 0 };
            auto count_helix{ 0 };
            for (auto remoteness : AtomicInfoHelper::GetStandardRemotenessList())
            {
                for (auto branch : AtomicInfoHelper::GetStandardBranchList())
                {
                    auto mix_group_key{ KeyPackerResidueClass::Pack(residue, element, remoteness, branch, false) };
                    auto free_group_key{ KeyPackerStructureClass::Pack(Structure::FREE, residue, element, remoteness, branch, false) };
                    auto helix_group_key{ KeyPackerStructureClass::Pack(Structure::HELX_P, residue, element, remoteness, branch, false) };
                    bool has_data{ false };
                    if (entry_iter->IsAvailableGroupKey(mix_group_key, AtomicInfoHelper::GetResidueClassKey()) == true)
                    {
                        has_data = true;
                        auto amplitude_value{ entry_iter->GetGausEstimatePrior(mix_group_key, AtomicInfoHelper::GetResidueClassKey(), 0) };
                        auto width_value{ entry_iter->GetGausEstimatePrior(mix_group_key, AtomicInfoHelper::GetResidueClassKey(), 1) };
                        auto amplitude_error{ entry_iter->GetGausVariancePrior(mix_group_key, AtomicInfoHelper::GetResidueClassKey(), 0) };
                        auto width_error{ entry_iter->GetGausVariancePrior(mix_group_key, AtomicInfoHelper::GetResidueClassKey(), 1) };
                        amplitude_array.emplace_back(amplitude_value);
                        width_array.emplace_back(width_value);
                        amplitude_mix_graph->SetPoint(count_mix, count_total, amplitude_value);
                        amplitude_mix_graph->SetPointError(count_mix, 0.0, amplitude_error);
                        width_mix_graph->SetPoint(count_mix, count_total, width_value);
                        width_mix_graph->SetPointError(count_mix, 0.0, width_error);
                        count_mix++;
                    }

                    if (entry_iter->IsAvailableGroupKey(free_group_key, AtomicInfoHelper::GetStructureClassKey()) == true)
                    {
                        auto amplitude_value{ entry_iter->GetGausEstimatePrior(free_group_key, AtomicInfoHelper::GetStructureClassKey(), 0) };
                        auto width_value{ entry_iter->GetGausEstimatePrior(free_group_key, AtomicInfoHelper::GetStructureClassKey(), 1) };
                        auto amplitude_error{ entry_iter->GetGausVariancePrior(free_group_key, AtomicInfoHelper::GetStructureClassKey(), 0) };
                        auto width_error{ entry_iter->GetGausVariancePrior(free_group_key, AtomicInfoHelper::GetStructureClassKey(), 1) };
                        amplitude_array.emplace_back(amplitude_value);
                        width_array.emplace_back(width_value);
                        amplitude_free_graph->SetPoint(count_free, count_total, amplitude_value);
                        amplitude_free_graph->SetPointError(count_free, 0.0, amplitude_error);
                        width_free_graph->SetPoint(count_free, count_total, width_value);
                        width_free_graph->SetPointError(count_free, 0.0, width_error);
                        count_free++;
                    }

                    if (entry_iter->IsAvailableGroupKey(helix_group_key, AtomicInfoHelper::GetStructureClassKey()) == true)
                    {
                        auto amplitude_value{ entry_iter->GetGausEstimatePrior(helix_group_key, AtomicInfoHelper::GetStructureClassKey(), 0) };
                        auto width_value{ entry_iter->GetGausEstimatePrior(helix_group_key, AtomicInfoHelper::GetStructureClassKey(), 1) };
                        auto amplitude_error{ entry_iter->GetGausVariancePrior(helix_group_key, AtomicInfoHelper::GetStructureClassKey(), 0) };
                        auto width_error{ entry_iter->GetGausVariancePrior(helix_group_key, AtomicInfoHelper::GetStructureClassKey(), 1) };
                        amplitude_array.emplace_back(amplitude_value);
                        width_array.emplace_back(width_value);
                        amplitude_helix_graph->SetPoint(count_helix, count_total, amplitude_value);
                        amplitude_helix_graph->SetPointError(count_helix, 0.0, amplitude_error);
                        width_helix_graph->SetPoint(count_helix, count_total, width_value);
                        width_helix_graph->SetPointError(count_helix, 0.0, width_error);
                        count_helix++;
                    }

                    if (has_data)
                    {
                        auto element_label{ AtomicInfoHelper::GetLabel(element) };
                        auto remoteness_label{ AtomicInfoHelper::GetLabel(remoteness) };
                        auto branch_label{ AtomicInfoHelper::GetLabel(branch) };
                        label_list.emplace_back(element_label +"_{"+ remoteness_label + branch_label+"}");
                        count_total++;
                    }
                }
            }
            auto color{ AtomicInfoHelper::GetDisplayColor(element) };
            ROOTHelper::SetMarkerAttribute(amplitude_mix_graph.get(), kFullCircle, 1.5f, color);
            ROOTHelper::SetMarkerAttribute(width_mix_graph.get(), kFullCircle, 1.5f, color);
            ROOTHelper::SetMarkerAttribute(amplitude_free_graph.get(), kOpenTriangleUp, 1.5f, color);
            ROOTHelper::SetMarkerAttribute(width_free_graph.get(), kOpenTriangleUp, 1.5f, color);
            ROOTHelper::SetMarkerAttribute(amplitude_helix_graph.get(), kOpenSquare, 1.5f, color);
            ROOTHelper::SetMarkerAttribute(width_helix_graph.get(), kOpenSquare, 1.5f, color);
            amplitude_mix_graph_list[residue_index].push_back(std::move(amplitude_mix_graph));
            amplitude_free_graph_list[residue_index].push_back(std::move(amplitude_free_graph));
            amplitude_helix_graph_list[residue_index].push_back(std::move(amplitude_helix_graph));
            width_mix_graph_list[residue_index].push_back(std::move(width_mix_graph));
            width_free_graph_list[residue_index].push_back(std::move(width_free_graph));
            width_helix_graph_list[residue_index].push_back(std::move(width_helix_graph));
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
        for (auto & graph : amplitude_mix_graph_list[residue_index]) graph->Draw("P X0");
        for (auto & graph : amplitude_free_graph_list[residue_index]) graph->Draw("P X0");
        for (auto & graph : amplitude_helix_graph_list[residue_index]) graph->Draw("P X0");

        pad[1]->cd();
        for (auto & graph : width_mix_graph_list[residue_index]) graph->Draw("P X0");
        for (auto & graph : width_free_graph_list[residue_index]) graph->Draw("P X0");
        for (auto & graph : width_helix_graph_list[residue_index]) graph->Draw("P X0");

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

    const int col_size{ 3 };
    const int row_size{ 3 };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.10f, 0.02f, 0.08f, 0.02f, 0.02f, 0.02f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    
    const size_t col_element_index[col_size]{ 2, 1, 0 };
    const size_t row_element_index[row_size]{ 3, 0, 1 };
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

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::vector<std::unique_ptr<TGraphErrors>> graph_list[col_size][row_size];
    std::vector<double> data_array;
    for (size_t i = 0; i < col_size; i++)
    {
        for (size_t j = 0; j < row_size; j++)
        {
            if (col_element_index[i] == row_element_index[j]) continue;
            if (i == 2 && j == 2) continue;
            for (auto residue : AtomicInfoHelper::GetStandardResidueList())
            {
                auto col_group_key{ m_atom_classifier->GetMainChainResidueClassGroupKey(col_element_index[i], residue) };
                if (entry_iter->IsAvailableGroupKey(col_group_key, residue_class) == false) continue;
                auto row_group_key{ m_atom_classifier->GetMainChainResidueClassGroupKey(row_element_index[j], residue) };
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

    auto data_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(data_array, 0.1, 0.01, 0.99) };
    auto min{ std::get<0>(data_range) };
    auto max{ std::get<1>(data_range) };
    auto line_ref{ ROOTHelper::CreateLine(min, min, max, max) };
    ROOTHelper::SetLineAttribute(line_ref.get(), 2, 1, kRed);
    std::unique_ptr<TPaveText> col_title_text[col_size];
    std::unique_ptr<TPaveText> row_title_text[row_size];
    std::unique_ptr<TLegend> legend;
    for (int i = 0; i < col_size; i++)
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
            if (i == col_size - 1 && j == 1)
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
    canvas->cd();

    auto pad_title{ ROOTHelper::CreatePad("pad_title","", 0.40, 0.69, 1.00, 1.00) };
    ROOTHelper::SetPadDefaultStyle(pad_title.get());
    ROOTHelper::SetFillAttribute(pad_title.get(), 4000);
    pad_title->Draw();
    pad_title->cd();
    auto title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextMarginInCanvas(pad_title.get(), title_text.get(), 0.01, 0.01, 0.17, 0.02);
    ROOTHelper::SetPaveTextDefaultStyle(title_text.get());
    ROOTHelper::SetPaveAttribute(title_text.get(), 0, 0.2);
    ROOTHelper::SetLineAttribute(title_text.get(), 1, 0);
    ROOTHelper::SetFillAttribute(title_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(title_text.get(), 38.0f, 23, 22, 0.0, kYellow-10);
    title_text->AddText("Residue-Based On-Site Main-Chain");
    title_text->AddText(Form("%s Estimate Scatter Plot", title_label[par_id].data()));
    title_text->Draw();
    
    auto subtitle1_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextMarginInCanvas(pad_title.get(), subtitle1_text.get(), 0.01, 0.43, 0.10, 0.15);
    ROOTHelper::SetPaveTextDefaultStyle(subtitle1_text.get());
    ROOTHelper::SetPaveAttribute(subtitle1_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle1_text.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetTextAttribute(subtitle1_text.get(), 40, 133, 22, 0.0, kYellow-10);
    subtitle1_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
    subtitle1_text->Draw();

    auto subtitle2_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextMarginInCanvas(pad_title.get(), subtitle2_text.get(), 0.18, 0.23, 0.10, 0.15);
    ROOTHelper::SetPaveTextDefaultStyle(subtitle2_text.get());
    ROOTHelper::SetPaveAttribute(subtitle2_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle2_text.get(), 1001, kAzure-7, 0.3f);
    ROOTHelper::SetTextAttribute(subtitle2_text.get(), 40, 103, 22);
    subtitle2_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
    subtitle2_text->Draw();

    auto subtitle3_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextMarginInCanvas(pad_title.get(), subtitle3_text.get(), 0.38, 0.01, 0.10, 0.15);
    ROOTHelper::SetPaveTextDefaultStyle(subtitle3_text.get());
    ROOTHelper::SetPaveAttribute(subtitle3_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle3_text.get(), 1001, kAzure-7, 0.3f);
    ROOTHelper::SetTextAttribute(subtitle3_text.get(), 40, 103, 22);
    subtitle3_text->AddText(model_object->GetEmdID().data());
    subtitle3_text->Draw();

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
                    auto amplitude{ entry_iter->GetGausEstimatePrior(group_key, residue_class, 0) };
                    auto width{ entry_iter->GetGausEstimatePrior(group_key, residue_class, 1) };
                    result_text->AddText(Form("#color[633]{#font[1]{A} = %.2f  ;  #font[1]{#tau} = %.2f}", amplitude, width));

                    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
                }
            }
        }
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ModelPainter::PaintResidueClassWidthScatterPlot(
    ModelObject * model_object, const std::string & name, int par_id, bool draw_box_plot)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintResidueClassWidthScatterPlot"<< std::endl;
    auto residue_class{ AtomicInfoHelper::GetResidueClassKey() };

    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    const int col_size{ 2 };
    const int row_size{ 2 };

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.15f, 0.02f, 0.12f, 0.02f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int main_chain_element_size{ 4 };
    const std::string element_label[main_chain_element_size]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };

    short color_element[main_chain_element_size] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    
    std::vector<std::unique_ptr<TGraphErrors>> graph_list[col_size][row_size];
    std::vector<double> x_array[col_size][row_size];
    std::vector<double> y_array[col_size][row_size];
    std::vector<double> global_x_array[col_size];
    std::vector<double> global_y_array[row_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto element_id{ static_cast<size_t>(j * col_size + i) };
            for (auto residue : AtomicInfoHelper::GetStandardResidueList())
            {
                auto group_key{ m_atom_classifier->GetMainChainResidueClassGroupKey(element_id, residue) };
                if (entry_iter->IsAvailableGroupKey(group_key, residue_class) == false) continue;
                auto gaus_graph
                {
                    (par_id == 0) ?
                    entry_iter->CreateCOMDistanceToGausEstimateGraph(group_key, residue_class, 1) :
                    entry_iter->CreateInRangeAtomsToGausEstimateGraph(group_key, residue_class, 5.0, 1)
                };
                for (int p = 0; p < gaus_graph->GetN(); p++)
                {
                    x_array[i][j].emplace_back(gaus_graph->GetPointX(p));
                    y_array[i][j].emplace_back(gaus_graph->GetPointY(p));
                    global_x_array[i].emplace_back(gaus_graph->GetPointX(p));
                    global_y_array[j].emplace_back(gaus_graph->GetPointY(p));
                }
                graph_list[i][j].emplace_back(std::move(gaus_graph));
            }
        }
    }

    double x_min[col_size];
    double x_max[col_size];
    for (int i = 0; i < col_size; i++)
    {
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(global_x_array[i], 0.05) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);
    }
    double y_min[row_size];
    double y_max[row_size];
    for (int j = 0; j < row_size; j++)
    {
        auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(global_y_array[j], 0.1) };
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }

    std::unique_ptr<TH2D> summary_hist[col_size][row_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array[i][j], 0.0) };
            auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array[i][j], 0.0) };
            summary_hist[i][j] = ROOTHelper::CreateHist2D(
                Form("hist_%d_%d", i, j), "",
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
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto element_id{ j * col_size + i };
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min[j], y_max[j]);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 45.0f, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 45.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 55.0f, 1.4f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 45.0f, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            auto title{ (par_id == 0)? "Distance to C.M. #[]{#AA}" : "In Range Atoms" };
            frame[i][j]->GetXaxis()->SetTitle(title);
            frame[i][j]->GetYaxis()->SetTitle("Width");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            if (draw_box_plot == true)
            {
                ROOTHelper::SetLineAttribute(summary_hist[i][j].get(), 1, 1, color_element[element_id]);
                ROOTHelper::SetFillAttribute(summary_hist[i][j].get(), 1001, color_element[element_id], 0.3f);
                summary_hist[i][j]->Draw("CANDLE3 SAME");
            }
            else
            {
                for (auto & graph : graph_list[i][j])
                {
                    ROOTHelper::SetMarkerAttribute(graph.get(), 24, 1.0f, color_element[element_id]);
                    graph->Draw("P");
                }
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
    auto x_range_tuple{ model_object->GetModelPositionRange(0) };
    auto y_range_tuple{ model_object->GetModelPositionRange(1) };
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

    const int col_size{ 4 };
    const int row_size{ 3 };
    const std::string element_label[col_size]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };
    short color_element[col_size]   { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    short marker_element[col_size]  { 54, 53, 55, 59 };
    std::unique_ptr<TH2> frame[col_size][row_size];

    std::unique_ptr<TGraphErrors> position_graph[col_size][row_size];
    std::unique_ptr<TGraph2DErrors> amplitude_2d_graph[col_size][row_size];
    std::unique_ptr<TGraph2DErrors> width_2d_graph[col_size][row_size];
    std::unique_ptr<TGraphErrors> total_graph[row_size];
    for (size_t i = 0; i < col_size; i++)
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

    for (int i = 0; i < col_size; i++)
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
    
    for (int i = 0; i < col_size; i++)
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
    for (int i = 0; i < col_size; i++)
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

void ModelPainter::PaintAtomXYPosition(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintAtomXYPosition"<< std::endl;
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    auto normalized_z_position{ 0.5 };
    auto z_ratio_window{ 0.1 };
    auto graph{ entry_iter->CreateXYPositionTomographyGraph(normalized_z_position, z_ratio_window, true) };
    auto z_position{ model_object->GetModelPosition(2, normalized_z_position) };
    auto z_window{ model_object->GetModelLength(2) * z_ratio_window };
    auto com_pos{ model_object->GetCenterOfMassPosition() };
    std::vector<double> x_array, y_array;
    x_array.reserve(static_cast<size_t>(graph->GetN()));
    y_array.reserve(static_cast<size_t>(graph->GetN()));
    for (int p = 0; p < graph->GetN(); p++)
    {
        x_array.emplace_back(graph->GetPointX(p));
        y_array.emplace_back(graph->GetPointY(p));
    }

    auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.1) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.1) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };

    canvas->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.15, 0.05, 0.12, 0.12);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);

    auto frame{ ROOTHelper::CreateHist2D("frame","", 500, x_min, x_max, 500, y_min, y_max) };
    ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 50.0f, 1.1f);
    ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 50.0f, 1.4f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), 45.0f, 0.01f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 45.0f, 0.01f);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("Atom Position #font[2]{x} #[]{#AA}");
    frame->GetYaxis()->SetTitle("Atom Position #font[2]{y} #[]{#AA}");
    frame->Draw();
    ROOTHelper::SetMarkerAttribute(graph.get(), 6, 1.0f, kAzure-7);
    graph->Draw("P X0");

    auto pos_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextMarginInCanvas(gPad, pos_text.get(), 0.45, 0.06, 0.83, 0.12);
    ROOTHelper::SetPaveTextDefaultStyle(pos_text.get());
    ROOTHelper::SetFillAttribute(pos_text.get(), 4000);
    ROOTHelper::SetTextAttribute(pos_text.get(), 35, 133, 32);
    pos_text->AddText(
        Form("Atom Position #font[2]{z} = %.1f #pm %.1f #[]{#AA}",
        z_position - com_pos.at(2), z_window*0.5));
    pos_text->Draw();

    auto title_text{ ROOTHelper::CreatePaveText(0.02, 0.89, 0.67, 0.98, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(title_text.get());
    ROOTHelper::SetPaveAttribute(title_text.get(), 0, 0.2);
    ROOTHelper::SetLineAttribute(title_text.get(), 1, 0);
    ROOTHelper::SetFillAttribute(title_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(title_text.get(), 47, 23, 22, 0.0, kYellow-10);
    title_text->AddText("Tomography of Atomic Model");
    title_text->Draw();

    auto model_text{ ROOTHelper::CreatePaveText(0.69, 0.89, 0.95, 0.98, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(model_text.get());
    ROOTHelper::SetPaveAttribute(model_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(model_text.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetTextAttribute(model_text.get(), 50, 103, 22);
    model_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
    model_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ModelPainter::PaintAtomGausScatter(
    ModelObject * model_object, const std::string & name, bool do_normalize)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintAtomGausScatter"<< std::endl;
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    auto amplitude_min{ entry_iter->GetGausEstimateMinimum(0, Element::OXYGEN) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TH2> frame;
    std::vector<std::unique_ptr<TGraphErrors>> graph_list;
    std::vector<int> atom_count_list;
    std::vector<std::string> element_name_list;
    graph_list.reserve(AtomicInfoHelper::GetElementCount());
    atom_count_list.reserve(AtomicInfoHelper::GetElementCount());
    element_name_list.reserve(AtomicInfoHelper::GetElementCount());
    std::vector<double> x_array, y_array;
    x_array.reserve(model_object->GetNumberOfSelectedAtom());
    y_array.reserve(model_object->GetNumberOfSelectedAtom());
    for (auto & [element_type, element_name] : AtomicInfoHelper::GetElementLabelMap())
    {
        auto graph
        {
            (do_normalize == true) ?
            entry_iter->CreateNormalizedGausEstimateScatterGraph(element_type, amplitude_min) :
            entry_iter->CreateGausEstimateScatterGraph(element_type)
        };
        auto marker_size{ (AtomicInfoHelper::IsStandardElement(element_type)) ? 1.5f : 2.0f };
        ROOTHelper::SetMarkerAttribute(graph.get(),
            AtomicInfoHelper::GetDisplayMarker(element_type), marker_size,
            AtomicInfoHelper::GetDisplayColor(element_type));
        for (int p = 0; p < graph->GetN(); p++)
        {
            auto x{ graph->GetPointX(p) };
            auto y{ graph->GetPointY(p) };
            x_array.emplace_back(x);
            y_array.emplace_back(y);
            //if (x <= 50.0) amplitude_cut_count[i]++;
            //if (y <= 0.6) width_cut_count[i]++;
        }
        //std::cout << "Ratio of amplitude selected atoms: "<< amplitude_cut_count[i] <<" / "<< atom_count[i] << std::endl;
        //std::cout << "Ratio of width selected atoms: "<< width_cut_count[i] <<" / "<< atom_count[i] << std::endl;
        if (graph->GetN() == 0) continue;
        atom_count_list.emplace_back(graph->GetN());
        element_name_list.emplace_back(element_name);
        graph_list.emplace_back(std::move(graph));
    }

    auto x_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(x_array, 0.2) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array, 0.2) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };

    canvas->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.15, 0.05, 0.12, 0.20);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);

    frame = ROOTHelper::CreateHist2D("frame","", 500, x_min, x_max, 500, y_min, y_max);
    ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 50.0f, 1.1f);
    ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 50.0f, 1.4f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), 45.0f, 0.01f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 45.0f, 0.01f);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("Amplitude");
    frame->GetYaxis()->SetTitle("Width");
    frame->Draw();
    for (auto & graph : graph_list)
    {
        graph->Draw("P X0");
    }

    auto legend{ ROOTHelper::CreateLegend(0.15, 0.8, 0.95, 0.9, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 40.0f, 133, 12);
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    legend->SetNColumns(3);
    for (size_t i = 0; i < graph_list.size(); i++)
    {
        legend->AddEntry(graph_list.at(i).get(),
                         Form("#font[102]{%s} (%d)", element_name_list.at(i).data(), atom_count_list.at(i)), "p");
    }
    legend->Draw();

    auto title_text{ ROOTHelper::CreatePaveText(0.02, 0.91, 0.98, 0.98, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(title_text.get());
    ROOTHelper::SetPaveAttribute(title_text.get(), 0, 0.2);
    ROOTHelper::SetLineAttribute(title_text.get(), 1, 0);
    ROOTHelper::SetFillAttribute(title_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(title_text.get(), 38, 23, 22, 0.0, kYellow-10);
    title_text->AddText("Atom-Based Main/Side-Chain Gaussian Estimats Scatter Plot");
    title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ModelPainter::PaintAtomGausMainChain(ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ModelPainter::PaintAtomGausMainChain"<< std::endl;
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 1 };
    const int row_size{ 3 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 900) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.08f, 0.02f, 0.10f, 0.10f, 0.02f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int main_chain_element_count{ 4 };
    const std::string element_label[main_chain_element_count]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };
    const std::string y_title[row_size]
    {
        "Intensity",
        "Width",
        "Amplitude"
    };
    short color_element[main_chain_element_count] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TGraphErrors> gaus_graph[row_size][main_chain_element_count];
    std::string chain_id;
    for (size_t k = 0; k < main_chain_element_count; k++)
    {
        gaus_graph[2][k] = entry_iter->CreateGausEstimateToResidueIDGraph(k, chain_id, 0);
        gaus_graph[1][k] = entry_iter->CreateGausEstimateToResidueIDGraph(k, chain_id, 1);
        gaus_graph[0][k] = entry_iter->CreateGausEstimateToResidueIDGraph(k, chain_id, 2);
    }

    std::vector<double> y_array[row_size];
    double y_min[row_size], y_max[row_size];
    for (int j = 0; j < row_size; j++)
    {
        y_array[j].reserve(static_cast<size_t>(gaus_graph[j][0]->GetN() * 4));
        for (size_t k = 0; k < main_chain_element_count; k++)
        {
            for (int p = 0; p < gaus_graph[j][k]->GetN(); p++)
            {
                y_array[j].emplace_back(gaus_graph[j][k]->GetPointY(p));
            }
        }
        auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array[j], 0.2) };
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }
    
    double x_min, x_max, y_min_tmp, y_max_tmp;
    gaus_graph[0][0]->ComputeRange(x_min, y_min_tmp, x_max, y_max_tmp);

    std::unique_ptr<TPaveText> subtitle1_text;
    std::unique_ptr<TPaveText> subtitle2_text;
    std::unique_ptr<TPaveText> subtitle3_text;
    std::unique_ptr<TLegend> legend;
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 500, 0, x_max+5, 500, y_min[j], y_max[j]);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 45.0f, 0.9f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 510);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 45.0f, 1.3f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 40.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.02/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("Residue ID");
            frame[i][j]->GetYaxis()->SetTitle(y_title[j].data());
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            for (size_t k = 0; k < main_chain_element_count; k++)
            {
                ROOTHelper::SetMarkerAttribute(gaus_graph[j][k].get(), 20, 0.7f, color_element[k]);
                ROOTHelper::SetLineAttribute(gaus_graph[j][k].get(), 1, 1, color_element[k]);
                gaus_graph[j][k]->Draw("PL X0");
            }

            if (i == 0 && j == row_size - 1)
            {
                subtitle1_text = ROOTHelper::CreatePaveText(0.00, 1.01, 0.15, 1.32, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(subtitle1_text.get());
                ROOTHelper::SetPaveAttribute(subtitle1_text.get(), 0, 0.2);
                ROOTHelper::SetFillAttribute(subtitle1_text.get(), 1001, kAzure-7);
                ROOTHelper::SetTextAttribute(subtitle1_text.get(), 60.0f, 133, 22, 0.0, kYellow-10);
                subtitle1_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
                subtitle1_text->Draw();

                subtitle2_text = ROOTHelper::CreatePaveText(0.16, 1.01, 0.35, 1.32, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(subtitle2_text.get());
                ROOTHelper::SetPaveAttribute(subtitle2_text.get(), 0, 0.2);
                ROOTHelper::SetFillAttribute(subtitle2_text.get(), 1001, kAzure-7, 0.5f);
                ROOTHelper::SetTextAttribute(subtitle2_text.get(), 50.0f, 103, 22);
                subtitle2_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
                subtitle2_text->Draw();

                subtitle3_text = ROOTHelper::CreatePaveText(0.36, 1.01, 0.57, 1.32, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(subtitle3_text.get());
                ROOTHelper::SetPaveAttribute(subtitle3_text.get(), 0, 0.2);
                ROOTHelper::SetFillAttribute(subtitle3_text.get(), 1001, kAzure-7, 0.5f);
                ROOTHelper::SetTextAttribute(subtitle3_text.get(), 50.0f, 103, 22);
                subtitle3_text->AddText(model_object->GetEmdID().data());
                subtitle3_text->Draw();

                legend = ROOTHelper::CreateLegend(0.58, 1.01, 0.98, 1.32, true);
                ROOTHelper::SetLegendDefaultStyle(legend.get());
                ROOTHelper::SetTextAttribute(legend.get(), 30.0f, 133, 12);
                ROOTHelper::SetFillAttribute(legend.get(), 4000);
                for (size_t k = 0; k < main_chain_element_count; k++)
                {
                    legend->AddEntry(gaus_graph[j][k].get(), element_label[k].data(), "pl");
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

#ifdef HAVE_ROOT

void ModelPainter::PrintResolutionInfoPad(TPad * pad, TPaveText * text, double resolution)
{
    pad->cd();
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, 0.25, 0.07, 0.01, 0.02);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.1);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7);
    ROOTHelper::SetLineAttribute(text, 1, 0);
    ROOTHelper::SetTextAttribute(text, 80, 133, 22, 0.0, kYellow-10);
    text->AddText(Form("%.2f #AA", resolution));
    pad->Update();
}

void ModelPainter::PrintDataInfoPad(
    TPad * pad, TPaveText * text, const std::string & pdb_id, const std::string & emd_id)
{
    pad->cd();
    ROOTHelper::SetPaveTextMarginInCanvas(pad, text, 0.005, 0.24, 0.01, 0.02);
    ROOTHelper::SetPaveTextDefaultStyle(text);
    ROOTHelper::SetPaveAttribute(text, 0, 0.1);
    ROOTHelper::SetFillAttribute(text, 1001, kAzure-7, 0.5);
    ROOTHelper::SetTextAttribute(text, 55, 133, 12);
    text->AddText(("#font[102]{PDB-" + pdb_id +"}").data());
    text->AddText(("#font[102]{"+ emd_id +"}").data());
    pad->Update();
}

void ModelPainter::PrintAmplitudePad(TPad * pad, TH2 * hist)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.07, 0.005, 0.01, 0.01);
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
    ROOTHelper::SetPadMarginInCanvas(pad, 0.07, 0.005, 0.11, 0.01);
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
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.005, 0.01, 0.01);
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
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.005, 0.11, 0.01);
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
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.07, 0.11, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 35.0f, 1.0f);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 35.0f, 1.5f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 35.0f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 35.0f, 0.01f);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.03, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 505);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("Amplitude Estimate");
    hist->GetYaxis()->SetTitle("Width Estimate");
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->Draw("Y+");
}

void ModelPainter::PrintCountSummaryPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.07, 0.005, 0.01, 0.02);
    ROOTHelper::SetPadLayout(pad, 1, 0, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(pad, 0, 0, 4000, 0, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 32.0f, 1.2f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), 0.0f, 21);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), 0.0f);
    hist->GetXaxis()->SetLimits(-1.0, 20.0);
    hist->SetMinimum(0.5);
    hist->SetStats(0);
    hist->GetYaxis()->SetTitle("#splitline{Residue}{Counts}");
    hist->GetYaxis()->CenterTitle();
    hist->Draw();
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

#endif