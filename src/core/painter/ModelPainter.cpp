#include <rhbm_gem/core/painter/ModelPainter.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <detail/LocalPotentialAccess.hpp>
#include <detail/ModelPotentialView.hpp>
#include "PotentialPlotBuilder.hpp"
#include <detail/PotentialSeriesOps.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include "data/object/AtomStyleCatalog.hpp"
#include "data/object/BondStyleCatalog.hpp"
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/BondClassifier.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include "detail/PainterSupport.hpp"

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
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

namespace rhbm_gem {

ModelPainter::ModelPainter() = default;

ModelPainter::~ModelPainter()
{

}

void ModelPainter::AddModel(ModelObject & data_object)
{
    m_model_object_list.push_back(&data_object);
}

void ModelPainter::Painting()
{
    Logger::Log(LogLevel::Info, "ModelPainter::Painting() called.");
    Logger::Log(LogLevel::Info, "Folder path: " + m_folder_path);
    Logger::Log(LogLevel::Info, "Number of atom objects to be painted: "
                + std::to_string(m_model_object_list.size()));

    for (auto model_object : m_model_object_list)
    {
        auto label{ painter_internal::BuildPainterOutputLabel(*model_object) };
        
        PaintAtomGroupGausMainChain(model_object, "atom_group_gaus_main_chain_"+ label);
        //PaintBondGroupGausMainChain(model_object, "bond_group_gaus_main_chain_"+ label);

        PaintAtomGroupGausNucleotideMainChain(model_object, "atom_group_gaus_nucleotide_main_chain_"+ label);
//        model_object->BuildKDTreeRoot();
        //PaintGroupWidthScatterPlot(model_object, "group_gaus_com_"+ label, 0, true);
        //PaintGroupWidthScatterPlot(model_object, "group_gaus_knn_"+ label, 1, true);
        //PaintAtomXYPosition(model_object, "atom_position_"+ label);
        PaintAtomGausScatterPlot(model_object, "atom_gaus_scatter_"+ label, false);
        //PaintBondGausScatterPlot(model_object, "bond_gaus_scatter_"+ label, false);
        PaintAtomGausMainChain(model_object, "atom_gaus_main_chain_"+ label);
        //PaintBondGausMainChain(model_object, "bond_gaus_main_chain_"+ label);
        PaintAtomMapValueMainChain(model_object, "atom_map_value_main_chain_"+ label);
        //PaintBondMapValueMainChain(model_object, "bond_map_value_main_chain_"+ label);
        PaintAtomRankMainChain(model_object, "atom_rank_main_chain_"+ label);
    }
}

void ModelPainter::PaintAtomGroupGausMainChain(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomGroupGausMainChain");

    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    
    #ifdef HAVE_ROOT

    const int primary_element_size{ 4 };
    const int structure_size{ 3 };

    const Structure structure_list[structure_size+1]
    {
        Structure::FREE,
        Structure::FREE,
        Structure::HELX_P,
        Structure::SHEET
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
    frame[1] = ROOTHelper::CreateHist2D("hist_1","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[0] = ROOTHelper::CreateHist2D("hist_0","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[3] = ROOTHelper::CreateHist2D("hist_3","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[2] = ROOTHelper::CreateHist2D("hist_2","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[4] = ROOTHelper::CreateHist2D("hist_4","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[5] = ROOTHelper::CreateHist2D("hist_5","", 100, 0.0, 1.0, 100, 0.0, 1.0);

    std::vector<GroupKey> group_key_list[primary_element_size][structure_size+1];
    for (size_t i = 0; i < primary_element_size; i++)
    {
        group_key_list[i][0] = AtomClassifier::GetMainChainComponentAtomClassGroupKeyList(i);
        group_key_list[i][1] = AtomClassifier::GetMainChainStructureAtomClassGroupKeyList(i, structure_list[1]);
        group_key_list[i][2] = AtomClassifier::GetMainChainStructureAtomClassGroupKeyList(i, structure_list[2]);
        group_key_list[i][3] = AtomClassifier::GetMainChainStructureAtomClassGroupKeyList(i, structure_list[3]);
    }

    for (size_t j = 0; j < structure_size + 1; j++)
    {
        std::unique_ptr<TGraphErrors> amplitude_graph[primary_element_size];
        std::unique_ptr<TGraphErrors> width_graph[primary_element_size];
        std::unique_ptr<TGraphErrors> correlation_graph[primary_element_size];
        std::unique_ptr<TH2D> amplitude_hist[primary_element_size];
        std::unique_ptr<TH2D> width_hist[primary_element_size];
        std::vector<double> amplitude_array, width_array;
        amplitude_array.reserve(80);
        width_array.reserve(80);
        auto class_key{ (j == 0) ? ChemicalDataHelper::GetComponentAtomClassKey() : ChemicalDataHelper::GetStructureAtomClassKey() };
        for (size_t i = 0; i < primary_element_size; i++)
        {
            amplitude_graph[i] = plot_builder->CreateAtomGausEstimateToResidueGraph(group_key_list[i][j], class_key, 0);
            width_graph[i] = plot_builder->CreateAtomGausEstimateToResidueGraph(group_key_list[i][j], class_key, 1);
            correlation_graph[i] = plot_builder->CreateAtomGausEstimateScatterGraph(group_key_list[i][j], class_key, 0, 1);
            for (int p = 0; p < amplitude_graph[i]->GetN(); p++)
            {
                amplitude_array.push_back(amplitude_graph[i]->GetPointY(p));
                width_array.push_back(width_graph[i]->GetPointY(p));
            }
            auto element_color{ AtomStyleCatalog::GetMainChainElementColor(i) };
            auto element_marker{ AtomStyleCatalog::GetMainChainElementSolidMarker(i) };
            ROOTHelper::SetMarkerAttribute(amplitude_graph[i].get(), element_marker, 1.3f, element_color);
            ROOTHelper::SetMarkerAttribute(width_graph[i].get(), element_marker, 1.3f, element_color);
            ROOTHelper::SetMarkerAttribute(correlation_graph[i].get(), element_marker, 1.3f, element_color);
            ROOTHelper::SetLineAttribute(amplitude_graph[i].get(), 1, 1, element_color);
            ROOTHelper::SetLineAttribute(width_graph[i].get(), 1, 1, element_color);
            ROOTHelper::SetLineAttribute(correlation_graph[i].get(), 1, 1, element_color);
        }
        if (amplitude_graph[0]->GetN() == 0) continue;

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
            auto element_color{ AtomStyleCatalog::GetMainChainElementColor(static_cast<size_t>(i)) };
            ROOTHelper::SetLineAttribute(amplitude_hist[i].get(), 1, 1, element_color);
            ROOTHelper::SetLineAttribute(width_hist[i].get(), 1, 1, element_color);
            ROOTHelper::SetFillAttribute(amplitude_hist[i].get(), 1001, element_color, 0.3f);
            ROOTHelper::SetFillAttribute(width_hist[i].get(), 1001, element_color, 0.3f);
            amplitude_hist[i]->SetBarWidth(0.5f);
            width_hist[i]->SetBarWidth(0.5f);
        }

        auto count_hist{ plot_builder->CreateAtomResidueCountHistogram(class_key, structure_list[j]) };
        auto max_count{ static_cast<int>(count_hist->GetMaximum()) };

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
        frame[5]->GetYaxis()->SetLimits(0.0, max_count*1.1);

        pad[1]->cd();
        PrintAmplitudePad(pad[1].get(), frame[1].get());
        for (int i = 0; i < primary_element_size; i++) amplitude_graph[i]->Draw("PL X0");

        pad[0]->cd();
        PrintWidthPad(pad[0].get(), frame[0].get());
        for (int i = 0; i < primary_element_size; i++) width_graph[i]->Draw("PL X0");

        pad[3]->cd();
        PrintAmplitudeSummaryPad(pad[3].get(), frame[3].get());
        for (int i = 0; i < primary_element_size; i++) amplitude_hist[i]->Draw("CANDLE3 SAME");

        pad[2]->cd();
        PrintAtomWidthSummaryPad(pad[2].get(), frame[2].get());
        for (int i = 0; i < primary_element_size; i++) width_hist[i]->Draw("CANDLE3 SAME");

        pad[4]->cd();
        PrintGausSummaryPad(pad[4].get(), frame[4].get());
        for (int i = 0; i < primary_element_size; i++) correlation_graph[i]->Draw("P X0");

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
        gStyle->SetTextFont(132);
        ROOTHelper::SetFillAttribute(count_hist.get(), 1001, kAzure-7, 0.5f);
        ROOTHelper::SetLineAttribute(count_hist.get(), 1, 1, kAzure-7);
        ROOTHelper::SetMarkerAttribute(count_hist.get(), 20, 7.0f, kAzure);
        count_hist->Draw("HIST TEXT0 SAME");

        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintBondGroupGausMainChain(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintBondGroupGausMainChain");

    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    
    #ifdef HAVE_ROOT

    const int primary_element_size{ 4 };

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

    std::vector<GroupKey> group_key_list[primary_element_size];
    for (size_t i = 0; i < primary_element_size; i++)
    {
        group_key_list[i] = BondClassifier::GetMainChainComponentBondClassGroupKeyList(i);
    }

    std::unique_ptr<TGraphErrors> amplitude_graph[primary_element_size];
    std::unique_ptr<TGraphErrors> width_graph[primary_element_size];
    std::unique_ptr<TGraphErrors> correlation_graph[primary_element_size];
    std::unique_ptr<TH2D> amplitude_hist[primary_element_size];
    std::unique_ptr<TH2D> width_hist[primary_element_size];
    std::vector<double> amplitude_array, width_array;
    amplitude_array.reserve(80);
    width_array.reserve(80);
    auto class_key{ ChemicalDataHelper::GetComponentBondClassKey() };
    for (size_t i = 0; i < primary_element_size; i++)
    {
        amplitude_graph[i] = plot_builder->CreateBondGausEstimateToResidueGraph(group_key_list[i], class_key, 0);
        width_graph[i] = plot_builder->CreateBondGausEstimateToResidueGraph(group_key_list[i], class_key, 1);
        correlation_graph[i] = plot_builder->CreateBondGausEstimateScatterGraph(group_key_list[i], class_key, 0, 1);
        for (int p = 0; p < amplitude_graph[i]->GetN(); p++)
        {
            amplitude_array.push_back(amplitude_graph[i]->GetPointY(p));
            width_array.push_back(width_graph[i]->GetPointY(p));
        }
        auto element_color{ BondStyleCatalog::GetMainChainMemberColor(i) };
        auto element_marker{ BondStyleCatalog::GetMainChainMemberSolidMarker(i) };
        ROOTHelper::SetMarkerAttribute(amplitude_graph[i].get(), element_marker, 1.3f, element_color);
        ROOTHelper::SetMarkerAttribute(width_graph[i].get(), element_marker, 1.3f, element_color);
        ROOTHelper::SetMarkerAttribute(correlation_graph[i].get(), element_marker, 1.3f, element_color);
        ROOTHelper::SetLineAttribute(amplitude_graph[i].get(), 1, 1, element_color);
        ROOTHelper::SetLineAttribute(width_graph[i].get(), 1, 1, element_color);
        ROOTHelper::SetLineAttribute(correlation_graph[i].get(), 1, 1, element_color);
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
        auto element_color{ BondStyleCatalog::GetMainChainMemberColor(static_cast<size_t>(i)) };
        ROOTHelper::SetLineAttribute(amplitude_hist[i].get(), 1, 1, element_color);
        ROOTHelper::SetLineAttribute(width_hist[i].get(), 1, 1, element_color);
        ROOTHelper::SetFillAttribute(amplitude_hist[i].get(), 1001, element_color, 0.3f);
        ROOTHelper::SetFillAttribute(width_hist[i].get(), 1001, element_color, 0.3f);
        amplitude_hist[i]->SetBarWidth(0.5f);
        width_hist[i]->SetBarWidth(0.5f);
    }

    auto count_hist{ plot_builder->CreateBondResidueCountHistogram(class_key) };
    auto max_count{ static_cast<int>(count_hist->GetMaximum()) };

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
    frame[5]->GetYaxis()->SetLimits(0.0, max_count*1.1);

    pad[1]->cd();
    PrintAmplitudePad(pad[1].get(), frame[1].get());
    for (int i = 0; i < primary_element_size; i++) amplitude_graph[i]->Draw("PL X0");

    pad[0]->cd();
    PrintWidthPad(pad[0].get(), frame[0].get());
    for (int i = 0; i < primary_element_size; i++) width_graph[i]->Draw("PL X0");

    pad[3]->cd();
    PrintAmplitudeSummaryPad(pad[3].get(), frame[3].get());
    for (int i = 0; i < primary_element_size; i++) amplitude_hist[i]->Draw("CANDLE3 SAME");

    pad[2]->cd();
    PrintBondWidthSummaryPad(pad[2].get(), frame[2].get());
    for (int i = 0; i < primary_element_size; i++) width_hist[i]->Draw("CANDLE3 SAME");

    pad[4]->cd();
    PrintGausSummaryPad(pad[4].get(), frame[4].get());
    for (int i = 0; i < primary_element_size; i++) correlation_graph[i]->Draw("P X0");

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
    gStyle->SetTextFont(132);
    ROOTHelper::SetFillAttribute(count_hist.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetLineAttribute(count_hist.get(), 1, 1, kAzure-7);
    ROOTHelper::SetMarkerAttribute(count_hist.get(), 20, 7.0f, kAzure);
    count_hist->Draw("HIST TEXT0 SAME");

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintAtomGroupGausNucleotideMainChain(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomGroupGausNucleotideMainChain");

    if (model_object->HasStandardRNAComponent() == false)
    {
        Logger::Log(LogLevel::Warning,
            " The model does not have standard RNA components. Skipped.");
        return;
    }

    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    
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
            AtomClassifier::GetNucleotideMainChainComponentAtomClassGroupKeyList(component_id[i]);
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
        amplitude_graph[i] = plot_builder->CreateAtomGausEstimateToSpotGraph(group_key_list[i], class_key, 0);
        width_graph[i] = plot_builder->CreateAtomGausEstimateToSpotGraph(group_key_list[i], class_key, 1);
        correlation_graph[i] = plot_builder->CreateAtomGausEstimateScatterGraph(group_key_list[i], class_key, 0, 1);
        for (int p = 0; p < amplitude_graph[i]->GetN(); p++)
        {
            amplitude_array.push_back(amplitude_graph[i]->GetPointY(p));
            width_array.push_back(width_graph[i]->GetPointY(p));
        }
        auto component_color{ AtomStyleCatalog::GetMainChainElementColor(i) }; // TODO: nucleotide color
        auto component_marker{ AtomStyleCatalog::GetMainChainElementSolidMarker(i) }; // TODO: nucleotide marker
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
        auto element_color{ AtomStyleCatalog::GetMainChainElementColor(static_cast<size_t>(i)) };
        ROOTHelper::SetLineAttribute(amplitude_hist[i].get(), 1, 1, element_color);
        ROOTHelper::SetLineAttribute(width_hist[i].get(), 1, 1, element_color);
        ROOTHelper::SetFillAttribute(amplitude_hist[i].get(), 1001, element_color, 0.3f);
        ROOTHelper::SetFillAttribute(width_hist[i].get(), 1001, element_color, 0.3f);
        amplitude_hist[i]->SetBarWidth(0.5f);
        width_hist[i]->SetBarWidth(0.5f);
    }

    //auto count_hist{ plot_builder->CreateAtomResidueCountHistogram(class_key, Structure::FREE) };
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
    PrintAmplitudeSummaryPad(pad[3].get(), frame[3].get());
    for (int i = 0; i < component_size; i++) amplitude_hist[i]->Draw("CANDLE3 SAME");

    pad[2]->cd();
    PrintAtomWidthSummaryPad(pad[2].get(), frame[2].get());
    for (int i = 0; i < component_size; i++) width_hist[i]->Draw("CANDLE3 SAME");

    pad[4]->cd();
    PrintGausSummaryPad(pad[4].get(), frame[4].get());
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

} // namespace rhbm_gem
#include <rhbm_gem/core/painter/ModelPainter.hpp>
#include "PotentialPlotBuilder.hpp"
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/BondClassifier.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <detail/LocalPotentialAccess.hpp>
#include <detail/ModelPotentialView.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#endif

namespace rhbm_gem {

void ModelPainter::PaintAtomMapValueMainChain(ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomMapValueMainChain");
    (void)model_object;

    #ifdef HAVE_ROOT
    const auto & class_key{ ChemicalDataHelper::GetSimpleAtomClassKey() };
    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    const int col_size{ 4 };
    const int row_size{ 1 };
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 350) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.10f, 0.02f, 0.20f, 0.20f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const size_t main_chain_element_size{ 4 };

    std::vector<std::unique_ptr<TGraphErrors>> map_value_graph_list[main_chain_element_size];
    std::unique_ptr<TF1> gaus_function[main_chain_element_size];
    double amplitude_prior[main_chain_element_size];
    double width_prior[main_chain_element_size];
    std::vector<double> y_array;
    y_array.reserve(model_object->GetNumberOfSelectedAtom());
        for (size_t k = 0; k < main_chain_element_size; k++)
    {
        auto group_key{ AtomClassifier::GetMainChainSimpleAtomClassGroupKey(k) };
        if (!entry_iter->HasAtomGroup(group_key, class_key)) continue;
        for (auto atom : entry_iter->GetAtomObjectList(group_key, class_key))
        {
            auto atom_plot_builder{ std::make_unique<PotentialPlotBuilder>(atom) };
            auto graph{ atom_plot_builder->CreateBinnedDistanceToMapValueGraph() };
            ROOTHelper::SetLineAttribute(graph.get(), 1, 2, static_cast<short>(kAzure-7), 0.3f);
            map_value_graph_list[k].emplace_back(std::move(graph));
            auto map_value_range{ series_ops::ComputeMapValueRange(
                RequireLocalPotentialEntry(*atom), 0.0) };
            y_array.emplace_back(std::get<0>(map_value_range));
            y_array.emplace_back(std::get<1>(map_value_range));
        }
        gaus_function[k] = plot_builder->CreateAtomGroupGausFunctionPrior(group_key, class_key);
        amplitude_prior[k] = entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 0);
        width_prior[k] = entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 1);
    }

    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.15) };
    auto x_min{ 0.01 };
    auto x_max{ 1.49 };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };
    auto line_ref{ ROOTHelper::CreateLine(x_min, 0.0, x_max, 0.0) };

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> element_text[col_size];
    std::unique_ptr<TPaveText> result_text[col_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 100, x_min, x_max, 100, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 35.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 505);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 40.0f, 1.2f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35.0f, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->SetStats(0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("Map Value");
            frame[i][j]->Draw();

            for (auto & graph : map_value_graph_list[i])
            {
                graph->Draw("L X0");
            }

            ROOTHelper::SetLineAttribute(line_ref.get(), 2, 1, kBlack);
            line_ref->Draw();

            ROOTHelper::SetLineAttribute(gaus_function[i].get(), 2, 2, kRed);
            gaus_function[i]->Draw("SAME");

            element_text[i] = ROOTHelper::CreatePaveText(0.70, 0.75, 0.99, 0.99, "nbNDC ARC", true);
            ROOTHelper::SetPaveTextDefaultStyle(element_text[i].get());
            ROOTHelper::SetPaveAttribute(element_text[i].get(), 0, 0.2);
            ROOTHelper::SetTextAttribute(element_text[i].get(), 40.0f, 103, 22);
            ROOTHelper::SetFillAttribute(element_text[i].get(), 1001, AtomStyleCatalog::GetMainChainElementColor(static_cast<size_t>(i)), 0.5f);
            element_text[i]->AddText(AtomStyleCatalog::GetMainChainElementLabel(static_cast<size_t>(i)).data());
            element_text[i]->Draw();

            result_text[i] = ROOTHelper::CreatePaveText(0.05, 0.08, 0.95, 0.25, "nbNDC", true);
            ROOTHelper::SetPaveTextDefaultStyle(result_text[i].get());
            ROOTHelper::SetTextAttribute(result_text[i].get(), 20.0f, 133, 22, 0.0, kRed);
            ROOTHelper::SetFillAttribute(result_text[i].get(), 4000);
            result_text[i]->AddText(Form("#font[2]{A} = %.2f  ;  #tau = %.2f", amplitude_prior[i], width_prior[i]));
            result_text[i]->Draw();
        }
    }

    canvas->cd();
    auto pad_title{ ROOTHelper::CreatePad("pad_title","", 0.00, 0.80, 1.00, 1.00) };
    ROOTHelper::SetPadDefaultStyle(pad_title.get());
    pad_title->Draw();
    pad_title->cd();

    auto subtitle1_text{ ROOTHelper::CreatePaveText(0.10, 0.10, 0.25, 0.90, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(subtitle1_text.get());
    ROOTHelper::SetPaveAttribute(subtitle1_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle1_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(subtitle1_text.get(), 40.0f, 133, 22, 0.0, kYellow-10);
    subtitle1_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
    subtitle1_text->Draw();

    auto subtitle2_text{ ROOTHelper::CreatePaveText(0.26, 0.10, 0.45, 0.90, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(subtitle2_text.get());
    ROOTHelper::SetPaveAttribute(subtitle2_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle2_text.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetTextAttribute(subtitle2_text.get(), 35.0f, 103, 22);
    subtitle2_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
    subtitle2_text->Draw();

    auto subtitle3_text{ ROOTHelper::CreatePaveText(0.46, 0.10, 0.68, 0.90, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(subtitle3_text.get());
    ROOTHelper::SetPaveAttribute(subtitle3_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle3_text.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetTextAttribute(subtitle3_text.get(), 35.0f, 103, 22);
    subtitle3_text->AddText(model_object->GetEmdID().data());
    subtitle3_text->Draw();

    auto legend{ ROOTHelper::CreateLegend(0.69, 0.10, 0.99, 0.90, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    ROOTHelper::SetTextAttribute(legend.get(), 25.0f, 133, 12, 0.0);
    legend->SetMargin(0.15f);
    legend->AddEntry(gaus_function[0].get(),
        "Gaussian Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})}", "l");
    legend->AddEntry(map_value_graph_list[0].at(0).get(),
        "Map Value", "l");
    legend->Draw();

    canvas->cd();
    auto pad_bottom{ ROOTHelper::CreatePad("pad_bottom","", 0.10, 0.00, 0.98, 0.11) };
    ROOTHelper::SetPadDefaultStyle(pad_bottom.get());
    pad_bottom->Draw();
    pad_bottom->cd();
    auto x_title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(x_title_text.get());
    ROOTHelper::SetFillAttribute(x_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(x_title_text.get(), 35.0f, 133, 22);
    x_title_text->AddText("Radial Distance #[]{#AA}");
    x_title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintBondMapValueMainChain(ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintBondMapValueMainChain");
    (void)model_object;

    #ifdef HAVE_ROOT
    const auto & class_key{ ChemicalDataHelper::GetSimpleBondClassKey() };
    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    const int col_size{ 4 };
    const int row_size{ 1 };
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 350) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.10f, 0.02f, 0.20f, 0.20f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const size_t main_chain_element_size{ 4 };

    std::vector<std::unique_ptr<TGraphErrors>> map_value_graph_list[main_chain_element_size];
    std::unique_ptr<TF1> gaus_function[main_chain_element_size];
    double amplitude_prior[main_chain_element_size];
    double width_prior[main_chain_element_size];
    int member_entries[main_chain_element_size];
    std::vector<double> y_array;
    y_array.reserve(model_object->GetNumberOfSelectedBond());
        for (size_t k = 0; k < main_chain_element_size; k++)
    {
        auto group_key{ BondClassifier::GetMainChainSimpleBondClassGroupKey(k) };
        if (!entry_iter->HasBondGroup(group_key, class_key, true)) continue;

        for (auto bond : entry_iter->GetBondObjectList(group_key, class_key))
        {
            auto bond_plot_builder{ std::make_unique<PotentialPlotBuilder>(bond) };
            auto graph{ bond_plot_builder->CreateBinnedDistanceToMapValueGraph() };
            ROOTHelper::SetLineAttribute(graph.get(), 1, 2, static_cast<short>(kAzure-7), 0.3f);
            map_value_graph_list[k].emplace_back(std::move(graph));
            auto map_value_range{ series_ops::ComputeMapValueRange(
                RequireLocalPotentialEntry(*bond), 0.0) };
            y_array.emplace_back(std::get<0>(map_value_range));
            y_array.emplace_back(std::get<1>(map_value_range));
        }
        gaus_function[k] = plot_builder->CreateBondGroupGausFunctionPrior(group_key, class_key);
        amplitude_prior[k] = entry_iter->GetBondGausEstimatePrior(group_key, class_key, 0);
        width_prior[k] = entry_iter->GetBondGausEstimatePrior(group_key, class_key, 1);
        member_entries[k] = static_cast<int>(entry_iter->GetBondObjectList(group_key, class_key).size());
    }

    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.15) };
    auto x_min{ 0.01 };
    auto x_max{ 1.49 };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };
    auto line_ref{ ROOTHelper::CreateLine(x_min, 0.0, x_max, 0.0) };

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> element_text[col_size];
    std::unique_ptr<TPaveText> result_text[col_size];
    std::unique_ptr<TPaveText> count_text[col_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 100, x_min, x_max, 100, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 35.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 505);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 40.0f, 1.2f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35.0f, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->SetStats(0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("Map Value");
            frame[i][j]->Draw();

            for (auto & graph : map_value_graph_list[i])
            {
                graph->Draw("L X0");
            }

            ROOTHelper::SetLineAttribute(line_ref.get(), 2, 1, kBlack);
            line_ref->Draw();

            ROOTHelper::SetLineAttribute(gaus_function[i].get(), 2, 2, kRed);
            gaus_function[i]->Draw("SAME");

            element_text[i] = ROOTHelper::CreatePaveText(0.55, 0.75, 0.99, 0.99, "nbNDC ARC", true);
            ROOTHelper::SetPaveTextDefaultStyle(element_text[i].get());
            ROOTHelper::SetPaveAttribute(element_text[i].get(), 0, 0.2);
            ROOTHelper::SetTextAttribute(element_text[i].get(), 40.0f, 103, 22);
            ROOTHelper::SetFillAttribute(element_text[i].get(), 1001, BondStyleCatalog::GetMainChainMemberColor(static_cast<size_t>(i)), 0.5f);
            element_text[i]->AddText(BondStyleCatalog::GetMainChainMemberLabel(static_cast<size_t>(i)).data());
            element_text[i]->Draw();

            count_text[i] = ROOTHelper::CreatePaveText(0.03, 0.90, 0.54, 0.99, "nbNDC", true);
            ROOTHelper::SetPaveTextDefaultStyle(count_text[i].get());
            ROOTHelper::SetTextAttribute(count_text[i].get(), 15.0f, 133, 12, 0.0, kGray+2);
            ROOTHelper::SetFillAttribute(count_text[i].get(), 4000);
            count_text[i]->AddText(Form("#Entries = %d", member_entries[i]));
            count_text[i]->Draw();

            result_text[i] = ROOTHelper::CreatePaveText(0.05, 0.08, 0.95, 0.25, "nbNDC", true);
            ROOTHelper::SetPaveTextDefaultStyle(result_text[i].get());
            ROOTHelper::SetTextAttribute(result_text[i].get(), 20.0f, 133, 22, 0.0, kRed);
            ROOTHelper::SetFillAttribute(result_text[i].get(), 4000);
            result_text[i]->AddText(Form("#font[2]{A} = %.2f  ;  #tau = %.2f", amplitude_prior[i], width_prior[i]));
            result_text[i]->Draw();
        }
    }

    canvas->cd();
    auto pad_title{ ROOTHelper::CreatePad("pad_title","", 0.00, 0.80, 1.00, 1.00) };
    ROOTHelper::SetPadDefaultStyle(pad_title.get());
    pad_title->Draw();
    pad_title->cd();

    auto subtitle1_text{ ROOTHelper::CreatePaveText(0.10, 0.10, 0.25, 0.90, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(subtitle1_text.get());
    ROOTHelper::SetPaveAttribute(subtitle1_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle1_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(subtitle1_text.get(), 40.0f, 133, 22, 0.0, kYellow-10);
    subtitle1_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
    subtitle1_text->Draw();

    auto subtitle2_text{ ROOTHelper::CreatePaveText(0.26, 0.10, 0.45, 0.90, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(subtitle2_text.get());
    ROOTHelper::SetPaveAttribute(subtitle2_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle2_text.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetTextAttribute(subtitle2_text.get(), 35.0f, 103, 22);
    subtitle2_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
    subtitle2_text->Draw();

    auto subtitle3_text{ ROOTHelper::CreatePaveText(0.46, 0.10, 0.68, 0.90, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(subtitle3_text.get());
    ROOTHelper::SetPaveAttribute(subtitle3_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle3_text.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetTextAttribute(subtitle3_text.get(), 35.0f, 103, 22);
    subtitle3_text->AddText(model_object->GetEmdID().data());
    subtitle3_text->Draw();

    auto legend{ ROOTHelper::CreateLegend(0.69, 0.10, 0.99, 0.90, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    ROOTHelper::SetTextAttribute(legend.get(), 25.0f, 133, 12, 0.0);
    legend->SetMargin(0.15f);
    legend->AddEntry(gaus_function[0].get(),
        "Gaussian Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})}", "l");
    legend->AddEntry(map_value_graph_list[0].at(0).get(),
        "Map Value", "l");
    legend->Draw();

    canvas->cd();
    auto pad_bottom{ ROOTHelper::CreatePad("pad_bottom","", 0.10, 0.00, 0.98, 0.11) };
    ROOTHelper::SetPadDefaultStyle(pad_bottom.get());
    pad_bottom->Draw();
    pad_bottom->cd();
    auto x_title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(x_title_text.get());
    ROOTHelper::SetFillAttribute(x_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(x_title_text.get(), 35.0f, 133, 22);
    x_title_text->AddText("Radial Distance #[]{#AA}");
    x_title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

} // namespace rhbm_gem
#include <rhbm_gem/core/painter/ModelPainter.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <detail/LocalPotentialAccess.hpp>
#include <detail/ModelPotentialView.hpp>
#include "PotentialPlotBuilder.hpp"
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/BondClassifier.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
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

namespace rhbm_gem {
void ModelPainter::PaintGroupWidthScatterPlot(
    ModelObject * model_object, const std::string & name, int par_id, bool draw_box_plot)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintGroupWidthScatterPlot");
    auto residue_class{ ChemicalDataHelper::GetComponentAtomClassKey() };
    (void)par_id;
    (void)draw_box_plot;

    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    const int col_size{ 4 };
    const int row_size{ 1 };

    auto canvas{ ROOTHelper::CreateCanvas("test","", 2000, 500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.15f, 0.02f, 0.12f, 0.02f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    
    std::vector<std::unique_ptr<TGraphErrors>> graph_list[col_size][row_size];
    std::vector<double> x_array[col_size][row_size];
    std::vector<double> y_array[col_size][row_size];
    std::vector<double> global_x_array[col_size];
    std::vector<double> global_y_array;
    for (size_t i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            for (auto residue : ChemicalDataHelper::GetStandardAminoAcidList())
            {
                auto group_key{ AtomClassifier::GetMainChainComponentAtomClassGroupKey(i, residue) };
                if (!entry_iter->HasAtomGroup(group_key, residue_class)) continue;
                auto gaus_graph
                {
                    (par_id == 0) ?
                    plot_builder->CreateCOMDistanceToGausEstimateGraph(group_key, residue_class, 1) :
                    plot_builder->CreateInRangeAtomsToGausEstimateGraph(group_key, residue_class, 5.0, 1)
                };
                for (int p = 0; p < gaus_graph->GetN(); p++)
                {
                    x_array[i][j].emplace_back(gaus_graph->GetPointX(p));
                    y_array[i][j].emplace_back(gaus_graph->GetPointY(p));
                    global_x_array[i].emplace_back(gaus_graph->GetPointX(p));
                    global_y_array.emplace_back(gaus_graph->GetPointY(p));
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

    auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(global_y_array, 0.1) };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };

    std::unique_ptr<TH2D> summary_hist[col_size][row_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto x_range_tmp{ ArrayStats<double>::ComputeScalingRangeTuple(x_array[i][j], 0.0) };
            auto y_range_tmp{ ArrayStats<double>::ComputeScalingRangeTuple(y_array[i][j], 0.0) };
            summary_hist[i][j] = ROOTHelper::CreateHist2D(
                Form("hist_%d_%d", i, j), "",
                5, std::get<0>(x_range_tmp), std::get<1>(x_range_tmp),
                100, std::get<0>(y_range_tmp), std::get<1>(y_range_tmp));
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
    for (size_t i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto element_color{ AtomStyleCatalog::GetMainChainElementColor(i) };
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), static_cast<int>(i), j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", static_cast<int>(i), j),"", 500, x_min[i], x_max[i], 500, y_min, y_max);
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
                ROOTHelper::SetLineAttribute(summary_hist[i][j].get(), 1, 1, element_color);
                ROOTHelper::SetFillAttribute(summary_hist[i][j].get(), 1001, element_color, 0.3f);
                summary_hist[i][j]->Draw("CANDLE3 SAME");
            }
            else
            {
                for (auto & graph : graph_list[i][j])
                {
                    ROOTHelper::SetMarkerAttribute(graph.get(), 24, 1.0f, element_color);
                    graph->Draw("P");
                }
            }
        }
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintAtomXYPosition(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomXYPosition");
    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    auto normalized_z_position{ 0.5 };
    auto z_ratio_window{ 0.1 };
    auto graph{ plot_builder->CreateAtomXYPositionTomographyGraph(
        normalized_z_position, z_ratio_window, true)
    };
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
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintAtomGausScatterPlot(
    ModelObject * model_object, const std::string & name, bool do_normalize)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomGausScatterPlot");
    (void)model_object;
    (void)do_normalize;

    #ifdef HAVE_ROOT
    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    auto amplitude_min{ entry_iter->GetAtomGausEstimateMinimum(0, Element::OXYGEN) };

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1300, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TH2> frame;
    std::unordered_map<Element, std::unique_ptr<TGraphErrors>> graph_map;
    std::vector<double> x_array, y_array;
    x_array.reserve(model_object->GetNumberOfSelectedAtom());
    y_array.reserve(model_object->GetNumberOfSelectedAtom());
    for (auto & [element_type, element_name] : ChemicalDataHelper::GetElementLabelMap())
    {
        auto graph
        {
            (do_normalize == true) ?
            plot_builder->CreateNormalizedAtomGausEstimateScatterGraph(element_type, amplitude_min) :
            plot_builder->CreateAtomGausEstimateScatterGraph(element_type)
        };
        auto atomic_number{ ChemicalDataHelper::GetAtomicNumber(element_type) };
        auto marker_size{ (atomic_number <= 8) ? 1.2f : 2.0f };
        short marker_color{ (ChemicalDataHelper::IsStandardElement(element_type)) ?
            ChemicalDataHelper::GetDisplayColor(element_type) : static_cast<short>(kRed)
        };
        auto transparency{ (atomic_number <= 8) ? 0.2f : 1.0f };
        ROOTHelper::SetMarkerAttribute(graph.get(),
            ChemicalDataHelper::GetDisplayMarker(element_type), marker_size,
            marker_color, transparency);
        for (int p = 0; p < graph->GetN(); p++)
        {
            auto x{ graph->GetPointX(p) };
            auto y{ graph->GetPointY(p) };
            x_array.emplace_back(x);
            y_array.emplace_back(y);
        }
        if (graph->GetN() == 0) continue;
        graph_map.emplace(element_type, std::move(graph));
    }

    auto x_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(x_array, 0.2, 0.005, 0.995) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array, 0.2, 0.005, 0.995) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };

    canvas->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.12, 0.20, 0.12, 0.10);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);

    frame = ROOTHelper::CreateHist2D("frame","", 500, x_min, x_max, 500, y_min, y_max);
    ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 50.0f, 1.1f);
    ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 50.0f, 1.4f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), 45.0f, 0.01f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 45.0f, 0.01f);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("Amplitude #font[2]{A}");
    frame->GetYaxis()->SetTitle("Width #font[2]{#tau}");
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->Draw();
    auto legend{ ROOTHelper::CreateLegend(0.81, 0.12, 0.98, 0.90, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 40.0f, 133, 12);
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    for (int i = 1; i <= 90; i++)
    {
        auto element{ ChemicalDataHelper::GetElementFromAtomicNumber(i) };
        if (graph_map.find(element) == graph_map.end()) continue;
        auto & graph{ graph_map.at(element) };
        graph->Draw("P X0");
        legend->AddEntry(graph.get(),
                         Form("#font[102]{%s} (%d)",
                            ChemicalDataHelper::GetElementLabelMap().at(element).data(),
                            graph->GetN()), "p");
    }
    legend->Draw();

    auto title_text{ ROOTHelper::CreatePaveText(0.12, 0.91, 0.98, 0.99, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(title_text.get());
    ROOTHelper::SetPaveAttribute(title_text.get(), 0, 0.2);
    ROOTHelper::SetLineAttribute(title_text.get(), 1, 0);
    ROOTHelper::SetFillAttribute(title_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(title_text.get(), 50.0f, 23, 22, 0.0, kYellow-10);
    title_text->AddText("Atom-Based Main/Side-Chain #font[2]{A}#minus#font[2]{#tau} Scatter Plot");
    title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintBondGausScatterPlot(
    ModelObject * model_object, const std::string & name, bool do_normalize)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintBondGausScatterPlot");
    (void)model_object;
    (void)do_normalize;

    #ifdef HAVE_ROOT
    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    auto amplitude_min{ entry_iter->GetBondGausEstimateMinimum(0) };

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1300, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TH2> frame;
    std::unordered_map<Element, std::unique_ptr<TGraphErrors>> graph_map;
    std::vector<double> x_array, y_array;
    x_array.reserve(model_object->GetNumberOfSelectedBond());
    y_array.reserve(model_object->GetNumberOfSelectedBond());
    for (auto & [element_type, element_name] : ChemicalDataHelper::GetElementLabelMap())
    {
        auto graph
        {
            (do_normalize == true) ?
            plot_builder->CreateNormalizedBondGausEstimateScatterGraph(element_type, amplitude_min) :
            plot_builder->CreateBondGausEstimateScatterGraph(element_type)
        };
        auto marker_size{ 2.0f };
        auto marker_color{ (ChemicalDataHelper::IsStandardElement(element_type)) ?
            ChemicalDataHelper::GetDisplayColor(element_type) : static_cast<short>(kRed)
        };
        ROOTHelper::SetMarkerAttribute(graph.get(),
            ChemicalDataHelper::GetDisplayMarker(element_type), marker_size, marker_color);
        for (int p = 0; p < graph->GetN(); p++)
        {
            auto x{ graph->GetPointX(p) };
            auto y{ graph->GetPointY(p) };
            x_array.emplace_back(x);
            y_array.emplace_back(y);
        }
        if (graph->GetN() == 0) continue;
        graph_map.emplace(element_type, std::move(graph));
    }

    auto x_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(x_array, 0.2, 0.005, 0.995) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array, 0.2, 0.005, 0.995) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };

    canvas->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.12, 0.20, 0.12, 0.10);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);

    frame = ROOTHelper::CreateHist2D("frame","", 500, x_min, x_max, 500, y_min, y_max);
    ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 50.0f, 1.1f);
    ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 50.0f, 1.4f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), 45.0f, 0.01f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 45.0f, 0.01f);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("Amplitude #font[2]{A}");
    frame->GetYaxis()->SetTitle("Width #font[2]{#tau}");
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->Draw();
    auto legend{ ROOTHelper::CreateLegend(0.81, 0.12, 0.98, 0.90, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 40.0f, 133, 12);
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    for (int i = 1; i <= 90; i++)
    {
        auto element{ ChemicalDataHelper::GetElementFromAtomicNumber(i) };
        if (graph_map.find(element) == graph_map.end()) continue;
        auto & graph{ graph_map.at(element) };
        graph->Draw("P X0");
        legend->AddEntry(graph.get(),
                         Form("#font[102]{%s} (%d)",
                            ChemicalDataHelper::GetElementLabelMap().at(element).data(),
                            graph->GetN()), "p");
    }
    legend->Draw();

    auto title_text{ ROOTHelper::CreatePaveText(0.12, 0.91, 0.98, 0.99, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(title_text.get());
    ROOTHelper::SetPaveAttribute(title_text.get(), 0, 0.2);
    ROOTHelper::SetLineAttribute(title_text.get(), 1, 0);
    ROOTHelper::SetFillAttribute(title_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(title_text.get(), 50.0f, 23, 22, 0.0, kYellow-10);
    title_text->AddText("Bond-Based Main/Side-Chain #font[2]{A}#minus#font[2]{#tau} Scatter Plot");
    title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintAtomGausMainChain(ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomGausMainChain");
    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };

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
    const std::string y_title[row_size]
    {
        "Intensity",
        "Width",
        "Amplitude"
    };
    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> gaus_graph_map[row_size][main_chain_element_count];
    for (size_t k = 0; k < main_chain_element_count; k++)
    {
        gaus_graph_map[2][k] = plot_builder->CreateAtomGausEstimateToSequenceIDGraphMap(k, 0);
        gaus_graph_map[1][k] = plot_builder->CreateAtomGausEstimateToSequenceIDGraphMap(k, 1);
        gaus_graph_map[0][k] = plot_builder->CreateAtomGausEstimateToSequenceIDGraphMap(k, 2);
    }

    for (auto & [chain_id, gaus_graph] : gaus_graph_map[0][0])
    {
        std::vector<double> x_array;
        std::vector<double> y_array[row_size];
        double y_min[row_size], y_max[row_size];
        for (int j = 0; j < row_size; j++)
        {
            y_array[j].reserve(static_cast<size_t>(gaus_graph->GetN() * 4));
            for (size_t k = 0; k < main_chain_element_count; k++)
            {
                if (gaus_graph_map[j][k].find(chain_id) == gaus_graph_map[j][k].end()) continue;
                for (int p = 0; p < gaus_graph_map[j][k].at(chain_id)->GetN(); p++)
                {
                    if (j == 0) x_array.emplace_back(gaus_graph_map[j][k].at(chain_id)->GetPointX(p));
                    y_array[j].emplace_back(gaus_graph_map[j][k].at(chain_id)->GetPointY(p));
                }
            }
            auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array[j], 0.2) };
            y_min[j] = std::get<0>(y_range);
            y_max[j] = std::get<1>(y_range);
        }
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.05) };
        auto x_min{ std::get<0>(x_range) };
        auto x_max{ std::get<1>(x_range) };

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
                if (frame[i][j] == nullptr)
                {
                    frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 500, 0.0, 1.0, 500, 0.0, 1.0);
                    ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 45.0f, 0.9f, 133);
                    ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.01f, 133);
                    ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 510);
                    ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 45.0f, 1.3f, 133);
                    ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 40.0f, 0.005f, 133);
                    ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.02/y_factor), 506);
                    ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
                    frame[i][j]->GetXaxis()->CenterTitle();
                    frame[i][j]->GetYaxis()->CenterTitle();
                    frame[i][j]->SetStats(0);
                }
                frame[i][j]->GetXaxis()->SetTitle(Form("Residue ID #[]{Chain %s}", chain_id.data()));
                frame[i][j]->GetYaxis()->SetTitle(y_title[j].data());
                frame[i][j]->GetXaxis()->SetLimits(x_min, x_max);
                frame[i][j]->GetYaxis()->SetLimits(y_min[j], y_max[j]);
                frame[i][j]->Draw();
                for (size_t k = 0; k < main_chain_element_count; k++)
                {
                    if (gaus_graph_map[j][k].find(chain_id) == gaus_graph_map[j][k].end()) continue;
                    auto element_color{ AtomStyleCatalog::GetMainChainElementColor(static_cast<size_t>(k)) };
                    ROOTHelper::SetMarkerAttribute(gaus_graph_map[j][k].at(chain_id).get(), 20, 0.7f, element_color);
                    ROOTHelper::SetLineAttribute(gaus_graph_map[j][k].at(chain_id).get(), 1, 1, element_color);
                    gaus_graph_map[j][k].at(chain_id)->Draw("PL X0");
                }

                if (i == 0 && j == row_size - 1)
                {
                    subtitle1_text = ROOTHelper::CreatePaveText(0.00, 1.02, 0.15, 1.37, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle1_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle1_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle1_text.get(), 1001, kAzure-7);
                    ROOTHelper::SetTextAttribute(subtitle1_text.get(), 60.0f, 133, 22, 0.0, kYellow-10);
                    subtitle1_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
                    subtitle1_text->Draw();

                    subtitle2_text = ROOTHelper::CreatePaveText(0.16, 1.02, 0.35, 1.37, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle2_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle2_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle2_text.get(), 1001, kAzure-7, 0.5f);
                    ROOTHelper::SetTextAttribute(subtitle2_text.get(), 50.0f, 103, 22);
                    subtitle2_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
                    subtitle2_text->Draw();

                    subtitle3_text = ROOTHelper::CreatePaveText(0.36, 1.02, 0.57, 1.37, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle3_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle3_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle3_text.get(), 1001, kAzure-7, 0.5f);
                    ROOTHelper::SetTextAttribute(subtitle3_text.get(), 50.0f, 103, 22);
                    subtitle3_text->AddText(model_object->GetEmdID().data());
                    subtitle3_text->Draw();

                    legend = ROOTHelper::CreateLegend(0.58, 1.02, 0.98, 1.37, true);
                    ROOTHelper::SetLegendDefaultStyle(legend.get());
                    ROOTHelper::SetTextAttribute(legend.get(), 70.0f, 133, 12);
                    ROOTHelper::SetFillAttribute(legend.get(), 4000);
                    for (size_t k = 0; k < main_chain_element_count; k++)
                    {
                        if (gaus_graph_map[j][k].find(chain_id) == gaus_graph_map[j][k].end()) continue;
                        auto label{ "#font[102]{" + AtomStyleCatalog::GetMainChainElementLabel(k) +"}" };
                        auto color{ AtomStyleCatalog::GetMainChainElementColor(k) };
                        auto result{ Form("#color[%d]{%s}", color, label.data()) };
                        legend->AddEntry(gaus_graph_map[j][k].at(chain_id).get(), result, "pl");
                    }
                    legend->SetNColumns(4);
                    legend->Draw();
                }
            }
        }
        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintBondGausMainChain(ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintBondGausMainChain");
    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };

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
    const std::string y_title[row_size]
    {
        "Intensity",
        "Width",
        "Amplitude"
    };
    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> gaus_graph_map[row_size][main_chain_element_count];
    for (size_t k = 0; k < main_chain_element_count; k++)
    {
        gaus_graph_map[2][k] = plot_builder->CreateBondGausEstimateToSequenceIDGraphMap(k, 0);
        gaus_graph_map[1][k] = plot_builder->CreateBondGausEstimateToSequenceIDGraphMap(k, 1);
        gaus_graph_map[0][k] = plot_builder->CreateBondGausEstimateToSequenceIDGraphMap(k, 2);
    }

    for (auto & [chain_id, gaus_graph] : gaus_graph_map[0][0])
    {
        std::vector<double> x_array;
        std::vector<double> y_array[row_size];
        double y_min[row_size], y_max[row_size];
        for (int j = 0; j < row_size; j++)
        {
            y_array[j].reserve(static_cast<size_t>(gaus_graph->GetN() * 4));
            for (size_t k = 0; k < main_chain_element_count; k++)
            {
                if (gaus_graph_map[j][k].find(chain_id) == gaus_graph_map[j][k].end()) continue;
                for (int p = 0; p < gaus_graph_map[j][k].at(chain_id)->GetN(); p++)
                {
                    if (j == 0) x_array.emplace_back(gaus_graph_map[j][k].at(chain_id)->GetPointX(p));
                    y_array[j].emplace_back(gaus_graph_map[j][k].at(chain_id)->GetPointY(p));
                }
            }
            auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array[j], 0.2) };
            y_min[j] = std::get<0>(y_range);
            y_max[j] = std::get<1>(y_range);
        }
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.05) };
        auto x_min{ std::get<0>(x_range) };
        auto x_max{ std::get<1>(x_range) };

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
                if (frame[i][j] == nullptr)
                {
                    frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 500, 0.0, 1.0, 500, 0.0, 1.0);
                    ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 45.0f, 0.9f, 133);
                    ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.01f, 133);
                    ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 510);
                    ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 45.0f, 1.3f, 133);
                    ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 40.0f, 0.005f, 133);
                    ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.02/y_factor), 506);
                    ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
                    frame[i][j]->GetXaxis()->CenterTitle();
                    frame[i][j]->GetYaxis()->CenterTitle();
                    frame[i][j]->SetStats(0);
                }
                frame[i][j]->GetXaxis()->SetTitle(Form("Residue ID #[]{Chain %s}", chain_id.data()));
                frame[i][j]->GetYaxis()->SetTitle(y_title[j].data());
                frame[i][j]->GetXaxis()->SetLimits(x_min, x_max);
                frame[i][j]->GetYaxis()->SetLimits(y_min[j], y_max[j]);
                frame[i][j]->Draw();
                for (size_t k = 0; k < main_chain_element_count; k++)
                {
                    if (gaus_graph_map[j][k].find(chain_id) == gaus_graph_map[j][k].end()) continue;
                    auto element_color{ BondStyleCatalog::GetMainChainMemberColor(static_cast<size_t>(k)) };
                    ROOTHelper::SetMarkerAttribute(gaus_graph_map[j][k].at(chain_id).get(), 20, 0.7f, element_color);
                    ROOTHelper::SetLineAttribute(gaus_graph_map[j][k].at(chain_id).get(), 1, 1, element_color);
                    gaus_graph_map[j][k].at(chain_id)->Draw("PL X0");
                }

                if (i == 0 && j == row_size - 1)
                {
                    subtitle1_text = ROOTHelper::CreatePaveText(0.00, 1.02, 0.15, 1.37, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle1_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle1_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle1_text.get(), 1001, kAzure-7);
                    ROOTHelper::SetTextAttribute(subtitle1_text.get(), 60.0f, 133, 22, 0.0, kYellow-10);
                    subtitle1_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
                    subtitle1_text->Draw();

                    subtitle2_text = ROOTHelper::CreatePaveText(0.16, 1.02, 0.35, 1.37, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle2_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle2_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle2_text.get(), 1001, kAzure-7, 0.5f);
                    ROOTHelper::SetTextAttribute(subtitle2_text.get(), 50.0f, 103, 22);
                    subtitle2_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
                    subtitle2_text->Draw();

                    subtitle3_text = ROOTHelper::CreatePaveText(0.36, 1.02, 0.57, 1.37, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle3_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle3_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle3_text.get(), 1001, kAzure-7, 0.5f);
                    ROOTHelper::SetTextAttribute(subtitle3_text.get(), 50.0f, 103, 22);
                    subtitle3_text->AddText(model_object->GetEmdID().data());
                    subtitle3_text->Draw();

                    legend = ROOTHelper::CreateLegend(0.58, 1.02, 0.98, 1.37, true);
                    ROOTHelper::SetLegendDefaultStyle(legend.get());
                    ROOTHelper::SetTextAttribute(legend.get(), 70.0f, 133, 12);
                    ROOTHelper::SetFillAttribute(legend.get(), 4000);
                    for (size_t k = 0; k < main_chain_element_count; k++)
                    {
                        if (gaus_graph_map[j][k].find(chain_id) == gaus_graph_map[j][k].end()) continue;
                        auto label{ "#font[102]{" + BondStyleCatalog::GetMainChainMemberLabel(k) +"}" };
                        auto color{ BondStyleCatalog::GetMainChainMemberColor(k) };
                        auto result{ Form("#color[%d]{%s}", color, label.data()) };
                        legend->AddEntry(gaus_graph_map[j][k].at(chain_id).get(), result, "pl");
                    }
                    legend->SetNColumns(4);
                    legend->Draw();
                }
            }
        }
        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintAtomRankMainChain(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomRankMainChain");
    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 4 };
    const int row_size{ 3 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 800) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.08f, 0.05f, 0.10f, 0.10f, 0.02f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const std::string y_title[row_size]
    {
        "Intensity",
        "Width",
        "Amplitude"
    };

    int chain_size;
    std::unique_ptr<TH2> frame[col_size][row_size];
    std::vector<std::unique_ptr<TH1D>> rank_hist_list[row_size][3];
    std::unique_ptr<TH1D> rank_reference[row_size];
    std::vector<Residue> veto_residues_list{  };
    rank_hist_list[2][0] = plot_builder->CreateMainChainAtomGausRankHistogram(0, chain_size, Residue::UNK, 2, veto_residues_list);
    rank_hist_list[1][0] = plot_builder->CreateMainChainAtomGausRankHistogram(1, chain_size, Residue::UNK, 1, veto_residues_list);
    rank_hist_list[0][0] = plot_builder->CreateMainChainAtomGausRankHistogram(2, chain_size, Residue::UNK, 0, veto_residues_list);
    rank_hist_list[2][1] = plot_builder->CreateMainChainAtomGausRankHistogram(0, chain_size, Residue::GLY);
    rank_hist_list[1][1] = plot_builder->CreateMainChainAtomGausRankHistogram(1, chain_size, Residue::GLY);
    rank_hist_list[0][1] = plot_builder->CreateMainChainAtomGausRankHistogram(2, chain_size, Residue::GLY);
    rank_hist_list[2][2] = plot_builder->CreateMainChainAtomGausRankHistogram(0, chain_size, Residue::PRO);
    rank_hist_list[1][2] = plot_builder->CreateMainChainAtomGausRankHistogram(1, chain_size, Residue::PRO);
    rank_hist_list[0][2] = plot_builder->CreateMainChainAtomGausRankHistogram(2, chain_size, Residue::PRO);
    
    for (int j = 0; j < row_size; j++)
    {
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
    std::unique_ptr<TPaveText> element_text[col_size];
    std::unique_ptr<TPaveText> par_text[row_size];
    std::unique_ptr<TLegend> legend;
    for (int i = 0; i < col_size; i++)
    {
        auto element_color{ AtomStyleCatalog::GetMainChainElementColor(static_cast<size_t>(i)) };
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 0, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 100, 0.1, 4.9, 100, 0.01, 110.0);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 45.0f, 0.8f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.01f, 103, kCyan+3);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 5);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 28.0f, 1.3f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 25.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("Rank");
            frame[i][j]->GetYaxis()->SetTitle("Count Ratio (%)");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();

            ROOTHelper::SetFillAttribute(rank_reference[j].get(), 1001, kGray, 0.2f);
            rank_reference[j]->SetBarWidth(0.9f);
            rank_reference[j]->SetBarOffset(0.05f);
            rank_reference[j]->Draw("BAR SAME");

            auto text_size{ 0.2f * 4.0f/static_cast<float>(gPad->GetAbsHNDC()) };
            ROOTHelper::SetFillAttribute(rank_hist_list[j][0].at(static_cast<size_t>(i)).get(), 1001, kAzure-7, 0.5f);
            ROOTHelper::SetLineAttribute(rank_hist_list[j][0].at(static_cast<size_t>(i)).get(), 1, 0);
            ROOTHelper::SetMarkerAttribute(rank_hist_list[j][0].at(static_cast<size_t>(i)).get(), 1, text_size, kAzure-7);
            rank_hist_list[j][0].at(static_cast<size_t>(i))->SetBarWidth(0.25f);
            rank_hist_list[j][0].at(static_cast<size_t>(i))->SetBarOffset(0.125f);
            rank_hist_list[j][0].at(static_cast<size_t>(i))->Draw("BAR TEXT0 SAME");

            ROOTHelper::SetFillAttribute(rank_hist_list[j][1].at(static_cast<size_t>(i)).get(), 1001, kOrange-6, 0.5f);
            ROOTHelper::SetLineAttribute(rank_hist_list[j][1].at(static_cast<size_t>(i)).get(), 1, 0);
            ROOTHelper::SetMarkerAttribute(rank_hist_list[j][1].at(static_cast<size_t>(i)).get(), 1, text_size, kOrange-6);
            rank_hist_list[j][1].at(static_cast<size_t>(i))->SetBarWidth(0.25f);
            rank_hist_list[j][1].at(static_cast<size_t>(i))->SetBarOffset(0.375f);
            rank_hist_list[j][1].at(static_cast<size_t>(i))->Draw("BAR TEXT0 SAME");

            ROOTHelper::SetFillAttribute(rank_hist_list[j][2].at(static_cast<size_t>(i)).get(), 1001, kMagenta-2, 0.5f);
            ROOTHelper::SetLineAttribute(rank_hist_list[j][2].at(static_cast<size_t>(i)).get(), 1, 0);
            ROOTHelper::SetMarkerAttribute(rank_hist_list[j][2].at(static_cast<size_t>(i)).get(), 1, text_size, kMagenta-2);
            rank_hist_list[j][2].at(static_cast<size_t>(i))->SetBarWidth(0.25f);
            rank_hist_list[j][2].at(static_cast<size_t>(i))->SetBarOffset(0.625f);
            rank_hist_list[j][2].at(static_cast<size_t>(i))->Draw("BAR TEXT0 SAME");

            if (j == row_size - 1)
            {
                element_text[i] = ROOTHelper::CreatePaveText(0.01, 1.01, 0.99, 1.20, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(element_text[i].get());
                ROOTHelper::SetPaveAttribute(element_text[i].get(), 0, 0.2);
                ROOTHelper::SetTextAttribute(element_text[i].get(), 30.0f, 133, 22);
                ROOTHelper::SetFillAttribute(element_text[i].get(), 1001, element_color, 0.5f);
                element_text[i]->AddText(AtomStyleCatalog::GetMainChainElementLabel(static_cast<size_t>(i)).data());
                element_text[i]->Draw();
            }
            if (i == col_size - 1)
            {
                par_text[j] = ROOTHelper::CreatePaveText(1.01, 0.01, 1.22, 0.99, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(par_text[j].get());
                ROOTHelper::SetPaveAttribute(par_text[j].get(), 0, 0.2);
                ROOTHelper::SetTextAttribute(par_text[j].get(), 30.0f, 133, 22, 0.0f, kYellow-10);
                ROOTHelper::SetLineAttribute(par_text[j].get(), 1, 0);
                ROOTHelper::SetFillAttribute(par_text[j].get(), 1001, kAzure-7);
                par_text[j]->AddText(y_title[j].data());
                auto text{ par_text[j]->GetLineWith(y_title[j].data()) };
                text->SetTextAngle(90.0f);
                par_text[j]->Draw();
            }
        }
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}


} // namespace rhbm_gem
#include <rhbm_gem/core/painter/ModelPainter.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>

#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/BondClassifier.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
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

namespace rhbm_gem {

#ifdef HAVE_ROOT

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
    for (size_t i = 0; i < ChemicalDataHelper::GetStandardAminoAcidCount(); i++)
    {
        auto residue{ ChemicalDataHelper::GetStandardAminoAcidList().at(i) };
        auto label{ ChemicalDataHelper::GetLabel(residue) };
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

void ModelPainter::PrintAtomWidthSummaryPad(TPad * pad, TH2 * hist)
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

    for (size_t i = 0; i < 4; i++)
    {
        auto index{ static_cast<int>(i) + 2 };
        auto color{ AtomStyleCatalog::GetMainChainElementColor(i) };
        auto label{ AtomStyleCatalog::GetMainChainElementLabel(i) };
        hist->GetXaxis()->ChangeLabel(index, 0.0, -1, -1, color, -1, label.data());
    }
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("Element");
    hist->GetXaxis()->CenterTitle();
    hist->Draw();
}

void ModelPainter::PrintBondWidthSummaryPad(TPad * pad, TH2 * hist)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.005, 0.11, 0.01);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 35.0f, 1.0f);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 35.0f, 0.055f, 103);
    ROOTHelper::SetAxisLabelAttribute(hist->GetYaxis(), 0.0f);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.0, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad, 0.015, 1) };
    ROOTHelper::SetAxisTickAttribute(hist->GetXaxis(), static_cast<float>(x_tick_length), 5);
    ROOTHelper::SetAxisTickAttribute(hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
    hist->GetXaxis()->SetLimits(-1.0, 4.0);
    hist->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    hist->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);

    for (size_t i = 0; i < 4; i++)
    {
        auto index{ static_cast<int>(i) + 2 };
        auto color{ BondStyleCatalog::GetMainChainMemberColor(i) };
        auto label{ BondStyleCatalog::GetMainChainMemberLabel(i) };
        hist->GetXaxis()->ChangeLabel(index, 45.0, -1, 21, color, -1, label.data());
    }
    hist->SetStats(0);
    hist->GetXaxis()->SetTitle("");
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

void ModelPainter::ModifyAxisLabelSideChain(
    TPad * pad, TH2 * hist, Residue residue, const std::vector<std::string> & label_list)
{
    if (ChemicalDataHelper::IsStandardResidue(residue) == false) return;

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


} // namespace rhbm_gem
