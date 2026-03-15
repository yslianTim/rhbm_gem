#include <rhbm_gem/core/painter/DemoPainter.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/PotentialEntryQuery.hpp>
#include <rhbm_gem/core/painter/PotentialPlotBuilder.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/core/painter/GausPainter.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
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
#include <TF1.h>
#include <TLine.h>
#include <TLatex.h>
#include <TEllipse.h>
#endif

#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <sstream>

#include <boost/math/distributions/students_t.hpp>

namespace rhbm_gem {
void DemoPainter::PaintGroupWidthScatterPlot(
    const std::vector<ModelObject *> & model_list, const std::string & name,
    int par_id, bool draw_box_plot)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " DemoPainter::PaintGroupWidthScatterPlot");
    auto residue_class{ ChemicalDataHelper::GetComponentAtomClassKey() };
    (void)par_id;
    (void)draw_box_plot;

    for (auto & model : model_list) if (model == nullptr) return;

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 4 };
    const int row_size{ 4 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 800) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.17f, 0.08f, 0.10f, 0.05f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::vector<std::unique_ptr<TGraphErrors>> graph_list[col_size][row_size];
    std::vector<double> x_array[col_size][row_size];
    std::vector<double> y_array[col_size][row_size];
    std::vector<double> global_x_array[col_size];
    std::vector<double> global_y_array[row_size];
    for (size_t i = 0; i < col_size; i++)
    {
        auto element_id{ i };
        for (size_t j = 0; j < row_size; j++)
        {
            auto model_id{ j };
            auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_list.at(model_id)) };
            auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_list.at(model_id)) };
            for (auto residue : ChemicalDataHelper::GetStandardAminoAcidList())
            {
                auto group_key{ m_atom_classifier->GetMainChainComponentAtomClassGroupKey(element_id, residue) };
                if (entry_iter->IsAvailableAtomGroupKey(group_key, residue_class) == false) continue;
                auto graph{ (par_id == 0) ?
                    plot_builder->CreateCOMDistanceToGausEstimateGraph(group_key, residue_class, 1) :
                    plot_builder->CreateInRangeAtomsToGausEstimateGraph(group_key, residue_class, 5.0, 1) };
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
        auto x_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(global_x_array[i], 0.12, 0.005, 0.995) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);
    }
    for (size_t j = 0; j < row_size; j++)
    {
        auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(global_y_array[j], 0.38, 0.005, 0.995) };
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }

    std::unique_ptr<TH2D> summary_hist[col_size][row_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto x_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(x_array[i][j], 0.0, 0.005, 0.995) };
            //auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array[i][j], 0.15) };
            auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array[i][j], 0.5, 0.005, 0.995) };
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
    std::unique_ptr<TPaveText> resolution_text[col_size];
    std::unique_ptr<TPaveText> title_x_text[col_size];
    std::unique_ptr<TPaveText> title_y_text[row_size];
    for (int i = 0; i < col_size; i++)
    {
        auto element_id{ i };
        auto element_color{ AtomClassifier::GetMainChainElementColor(static_cast<size_t>(element_id)) };
        auto element_label{ AtomClassifier::GetMainChainElementLabel(static_cast<size_t>(element_id)) };
        for (int j = 0; j < row_size; j++)
        {
            auto model_id{ j };
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min[j], y_max[j]);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 50.0f, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.08/x_factor), 505);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 50.0f, 1.2f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 45.0f, 0.02f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw("Y+");

            
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
                    ROOTHelper::SetLineAttribute(graph.get(), 1, 2, element_color);
                    graph->Draw("P");
                }
            }

            if (i == 0)
            {
                title_y_text[j] = ROOTHelper::CreatePaveText(-0.90, 0.20, 0.02, 0.80, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(title_y_text[j].get());
                ROOTHelper::SetPaveAttribute(title_y_text[j].get(), 0, 0.1);
                ROOTHelper::SetLineAttribute(title_y_text[j].get(), 1, 0);
                ROOTHelper::SetTextAttribute(title_y_text[j].get(), 40.0f, 103, 32);
                ROOTHelper::SetFillAttribute(title_y_text[j].get(), 1001, kAzure-7, 0.5f);
                title_y_text[j]->AddText(model_list.at(static_cast<size_t>(model_id))->GetPdbID().data());
                title_y_text[j]->AddText(model_list.at(static_cast<size_t>(model_id))->GetEmdID().data());
                title_y_text[j]->Draw();

                resolution_text[j] = ROOTHelper::CreatePaveText(-0.92, 0.50,-0.40, 0.90, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(resolution_text[j].get());
                ROOTHelper::SetPaveAttribute(resolution_text[j].get(), 0, 0.2);
                ROOTHelper::SetFillAttribute(resolution_text[j].get(), 1001, kAzure-7);
                ROOTHelper::SetTextAttribute(resolution_text[j].get(), 45.0f, 133, 22, 0.0, kYellow-10);
                resolution_text[j]->AddText(Form("%.2f #AA", model_list.at(static_cast<size_t>(model_id))->GetResolution()));
                resolution_text[j]->Draw();
            }

            if (j == row_size - 1)
            {
                title_x_text[i] = ROOTHelper::CreatePaveText(0.02, 0.95, 0.35, 1.23, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(title_x_text[i].get());
                ROOTHelper::SetPaveAttribute(title_x_text[i].get(), 0, 0.2);
                ROOTHelper::SetTextAttribute(title_x_text[i].get(), 50.0f, 103, 22);
                ROOTHelper::SetFillAttribute(title_x_text[i].get(), 1001, element_color, 0.5f);
                title_x_text[i]->AddText(element_label.data());
                //title_x_text[i]->Draw();
            }
        }
    }

    canvas->cd();
    auto pad_extra1{ ROOTHelper::CreatePad("pad_extra1","", 0.07, 0.00, 0.92, 0.06) };
    pad_extra1->Draw();
    pad_extra1->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra1.get());
    ROOTHelper::SetFillAttribute(pad_extra1.get(), 4000);
    auto bottom_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(bottom_title_text.get());
    ROOTHelper::SetFillAttribute(bottom_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(bottom_title_text.get(), 40.0f, 133, 22);
    bottom_title_text->AddText("Number of In-Range Atoms");
    bottom_title_text->Draw();
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);

    canvas->cd();
    auto pad_extra2{ ROOTHelper::CreatePad("pad_extra2","", 0.96, 0.10, 1.00, 0.90) };
    pad_extra2->Draw();
    pad_extra2->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra2.get());
    ROOTHelper::SetFillAttribute(pad_extra2.get(), 4000);
    auto right_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(right_title_text.get());
    ROOTHelper::SetFillAttribute(right_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(right_title_text.get(), 50.0f, 133, 22);
    right_title_text->AddText("Width #tau");
    auto text{ right_title_text->GetLineWith("Width") };
    text->SetTextAngle(90.0f);
    right_title_text->Draw();
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void DemoPainter::PaintAtomGausMainChainDemo(
    ModelObject * model_object1, ModelObject * model_object2, const std::string & name, int par_id)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " DemoPainter::PaintAtomGausMainChainDemo");
    (void)par_id;
    
    if (model_object1 == nullptr || model_object2 == nullptr) return;
    auto entry_iter1{ std::make_unique<PotentialEntryQuery>(model_object1) };
    auto plot_builder1{ std::make_unique<PotentialPlotBuilder>(model_object1) };
    auto entry_iter2{ std::make_unique<PotentialEntryQuery>(model_object2) };
    auto plot_builder2{ std::make_unique<PotentialPlotBuilder>(model_object2) };
    
    #ifdef HAVE_ROOT

    const int main_chain_element_count{ 4 };

    gStyle->SetLineScalePS(1.0);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 600) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TGraphErrors> gaus_graph[main_chain_element_count];
    std::unique_ptr<TGraphErrors> gaus_gly_graph;
    std::unique_ptr<TGraphErrors> gaus_pro_graph;
    std::vector<double> x_array, y_array;
    for (size_t i = 0; i < main_chain_element_count; i++)
    {
        gaus_graph[i] = std::move(plot_builder1->CreateAtomGausEstimateToSequenceIDGraphMap(i, par_id).at("A"));
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

    gaus_gly_graph = std::move(plot_builder1->CreateAtomGausEstimateToSequenceIDGraphMap(0, par_id, Residue::GLY).at("A"));
    gaus_pro_graph = std::move(plot_builder1->CreateAtomGausEstimateToSequenceIDGraphMap(0, par_id, Residue::PRO).at("A"));
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
    frame->GetXaxis()->SetTitle(Form("Residue ID (Chain #font[102]{A} in #font[102]{%s})", model_object1->GetPdbID().data()));
    frame->GetYaxis()->SetTitle((par_id == 0) ? "Amplitude #font[2]{A}" : "Width #font[2]{#tau}");
    frame->Draw();
    for (size_t i = 0; i < main_chain_element_count; i++)
    {
        auto element_color{ AtomClassifier::GetMainChainElementColor(i) };
        ROOTHelper::SetMarkerAttribute(gaus_graph[i].get(), 20, 1.0f, element_color);
        ROOTHelper::SetLineAttribute(gaus_graph[i].get(), 1, 2, element_color);
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
    legend->AddEntry(gaus_gly_graph.get(), "GLY", "p");
    legend->AddEntry(gaus_pro_graph.get(), "PRO", "p");
    for (size_t i = 0; i < main_chain_element_count; i++)
    {
        auto element_label{ AtomClassifier::GetMainChainElementLabel(i) };
        legend->AddEntry(gaus_graph[i].get(), element_label.data(), "pl");
    }
    legend->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void DemoPainter::PaintAtomGausMainChainDemoSingle(
    ModelObject * model_object, const std::string & name, int par_id)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " DemoPainter::PaintAtomGausMainChainDemoSingle");
    (void)par_id;
    
    if (model_object == nullptr) return;
    auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    
    #ifdef HAVE_ROOT

    const int main_chain_element_count{ 4 };

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
        gaus_graph[i] = std::move(plot_builder->CreateAtomGausEstimateToSequenceIDGraphMap(i, par_id).at("A"));
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

    gaus_gly_graph = std::move(plot_builder->CreateAtomGausEstimateToSequenceIDGraphMap(0, par_id, Residue::GLY).at("A"));
    gaus_pro_graph = std::move(plot_builder->CreateAtomGausEstimateToSequenceIDGraphMap(0, par_id, Residue::PRO).at("A"));
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
        auto element_color{ AtomClassifier::GetMainChainElementColor(i) };
        ROOTHelper::SetMarkerAttribute(gaus_graph[i].get(), 20, 1.0f, element_color);
        ROOTHelper::SetLineAttribute(gaus_graph[i].get(), 1, 2, element_color);
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
    legend->AddEntry(gaus_gly_graph.get(), "GLY", "p");
    legend->AddEntry(gaus_pro_graph.get(), "PRO", "p");
    for (size_t i = 0; i < main_chain_element_count; i++)
    {
        auto element_label{ AtomClassifier::GetMainChainElementLabel(i) };
        legend->AddEntry(gaus_graph[i].get(), element_label.data(), "pl");
    }
    legend->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void DemoPainter::PaintGroupWidthAlphaCarbonDemo(
    const std::vector<ModelObject *> & model_list, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    auto residue_class{ ChemicalDataHelper::GetComponentAtomClassKey() };
    Logger::Log(LogLevel::Info, " DemoPainter::PaintGroupWidthAlphaCarbonDemo");

    for (auto & model : model_list) if (model == nullptr) return;

    #ifdef HAVE_ROOT
    float marker_size{ 1.5f };

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    auto pad{ ROOTHelper::CreatePad("pad_0","", 0.0, 0.0, 1.0, 1.0) };
    std::unique_ptr<TH2> frame;

    const size_t model_size{ 6 };
    const short color_list[model_size] = { kRed-4, kBlue-4, kGreen-3, kMagenta-4, kOrange-3, kCyan+1 };
    const short marker_list[model_size] = { 20, 21, 22, 23, 33, 34 };
    const short line_style_list[model_size] = { 1, 2, 3, 4, 5, 6 };
    const int member_id{ 0 }; // Alpha Carbon
    std::unique_ptr<TGraphErrors> width_graph[model_size];
    std::vector<double> width_array;
    for (size_t j = 0; j < model_size; j++)
    {
        auto model_object{ model_list.at(j) };
        auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
        auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };

        auto group_key_list{ m_atom_classifier->GetMainChainComponentAtomClassGroupKeyList(member_id) };
        width_graph[j] = plot_builder->CreateAtomGausEstimateToResidueGraph(group_key_list, residue_class, 1);
        for (int p = 0; p < width_graph[j]->GetN(); p++)
        {
            width_array.push_back(width_graph[j]->GetPointY(p));
        }
        ROOTHelper::SetMarkerAttribute(width_graph[j].get(), marker_list[j], marker_size, color_list[j]);
        ROOTHelper::SetLineAttribute(width_graph[j].get(), line_style_list[j], 1, color_list[j]);
    }

    auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, 0.2) };
    frame = ROOTHelper::CreateHist2D("frame_total","", 100, 0.0, 1.0, 100, std::get<0>(width_range), std::get<1>(width_range));
    frame->GetYaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));
    frame->GetYaxis()->SetTitle("Width #font[2]{#tau}");
    frame->GetYaxis()->CenterTitle();

    auto legend{ ROOTHelper::CreateLegend(0.05, 0.85, 0.95, 1.00, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 40.0f, 133, 12);
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    for (size_t j = 0; j < model_size; j++)
    {
        auto model_object{ model_list.at(j) };
        auto text{ Form("#font[102]{%s} (FSC = %.2f #AA)", model_object->GetPdbID().data(), model_object->GetResolution()) };
        legend->AddEntry(width_graph[j].get(), text, "pl");
    }
    legend->SetNColumns(2);

    canvas->cd();
    ROOTHelper::SetPadDefaultStyle(pad.get());
    pad->Draw();
    pad->cd();
    auto left_margin{ 0.15 };
    auto right_margin{ 0.05 };
    auto bottom_margin{ 0.12 };
    auto top_margin{ 0.15 };
    auto label_size{ 55.0f };
    ROOTHelper::SetPadMarginInCanvas(gPad, left_margin, right_margin, bottom_margin, top_margin);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 60.0f, 1.2f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), label_size, 0.06f, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 50.0f, 0.01f);

    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.00, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.02, 1) };
    ROOTHelper::SetAxisTickAttribute(frame->GetXaxis(), static_cast<float>(x_tick_length), 21);
    ROOTHelper::SetAxisTickAttribute(frame->GetYaxis(), static_cast<float>(y_tick_length), 505);
    frame->GetXaxis()->SetLimits(-1.0, 20.0);
    frame->GetXaxis()->ChangeLabel(1, -1.0, 0.0);
    frame->GetXaxis()->ChangeLabel(-1, -1.0, 0.0);
    for (size_t i = 0; i < ChemicalDataHelper::GetStandardAminoAcidCount(); i++)
    {
        auto residue{ ChemicalDataHelper::GetStandardAminoAcidList().at(i) };
        auto label{ ChemicalDataHelper::GetLabel(residue) };
        auto label_index{ static_cast<int>(i) + 2 };
        frame->GetXaxis()->ChangeLabel(label_index, 90.0, -1, 12, -1, -1, label.data());
    }
    frame->SetStats(0);
    frame->Draw();
    legend->Draw();

    for (size_t j = 0; j < model_size; j++) width_graph[j]->Draw("PL X0");

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void DemoPainter::PaintGroupGausMergeResidueDemo(
    const std::vector<ModelObject *> & model_list, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    auto residue_class{ ChemicalDataHelper::GetComponentAtomClassKey() };
    Logger::Log(LogLevel::Info, " DemoPainter::PaintGroupGausMergeResidueDemo");

    for (auto & model : model_list) if (model == nullptr) return;

    const int col_size{ 3 };
    const int row_size{ 2 };
    auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    (void)row_size;
    const std::vector<Spot> spot_list{ Spot::CA, Spot::C, Spot::N };
    std::map<Spot, std::vector<GroupKey>> group_key_list_map[col_size];
    for (size_t i = 0; i < model_list.size(); i++)
    {
        auto model_object{ model_list.at(i) };
        auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
        auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
        for (size_t p = 0; p < spot_list.size(); p++)
        {
            auto spot{ spot_list.at(p) };
            group_key_list_map[i].emplace(
                spot, m_atom_classifier->GetMainChainComponentAtomClassGroupKeyList(p));
            auto & group_key_list{ group_key_list_map[i].at(spot) };
            for (auto it = group_key_list.begin(); it != group_key_list.end(); )
            {
                if (entry_iter->IsAvailableAtomGroupKey(*it, class_key) == false)
                {
                    it = group_key_list.erase(it);
                }
                else
                {
                    ++it;
                }
            }
        }
    }

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1200, 500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.09f, 0.07f, 0.09f, 0.12f, 0.02f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::map<Spot, std::unique_ptr<TGraphErrors>> graph_map[col_size];
    std::map<Spot, std::unique_ptr<TH2D>> hist_map[col_size][row_size];
    std::vector<double> y_array[row_size];
    std::vector<double> y_array_last[row_size];
    std::map<Spot, std::vector<double>> width_array_map[col_size];
    
    for (size_t i = 0; i < model_list.size(); i++)
    {
        auto model_object{ model_list.at(i) };
        auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
        auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
        for (auto & [spot, group_key_list] : group_key_list_map[i])
        {
            graph_map[i][spot] = plot_builder->CreateAtomGausEstimateScatterGraph(group_key_list, class_key, 0, 1);
            for (int p = 0; p < graph_map[i][spot]->GetN(); p++)
            {
                if (i == model_list.size() - 1)
                {
                    y_array_last[1].push_back(graph_map[i][spot]->GetPointX(p));
                    y_array_last[0].push_back(graph_map[i][spot]->GetPointY(p));
                }
                else
                {
                    y_array[1].push_back(graph_map[i][spot]->GetPointX(p));
                    y_array[0].push_back(graph_map[i][spot]->GetPointY(p));
                }
                width_array_map[i][spot].push_back(graph_map[i][spot]->GetPointY(p));
            }
        }
    }

    auto scaling{ 0.3 };
    std::tuple<double, double> range[row_size];
    std::tuple<double, double> range_last[row_size];
    range[0] = ArrayStats<double>::ComputeScalingRangeTuple(y_array[0], scaling);
    range[1] = ArrayStats<double>::ComputeScalingRangeTuple(y_array[1], scaling);
    range_last[0] = ArrayStats<double>::ComputeScalingRangeTuple(y_array_last[0], scaling);
    range_last[1] = ArrayStats<double>::ComputeScalingRangeTuple(y_array_last[1], scaling);
    auto spot_count{ spot_list.size() };
    std::vector<std::string> spot_label_list;
    spot_label_list.reserve(spot_count);
    for (size_t i = 0; i < model_list.size(); i++)
    {
        for (size_t p = 0; p < spot_count; p++)
        {
            auto spot{ spot_list.at(p) };
            auto spot_label{ AtomClassifier::GetMainChainSpotLabel(spot) };
            auto x_value{ static_cast<double>(p) };
            if (i == 0) spot_label_list.emplace_back(spot_label);
            std::string name_amplitude{ "amplitude_hist_"+ spot_label + std::to_string(i) };
            std::string name_width{ "width_hist_"+ spot_label + std::to_string(i) };
            auto amplitude_range_tmp{ (i == model_list.size() - 1) ? range_last[1] : range[1] };
            auto width_range_tmp{ (i == model_list.size() - 1) ? range_last[0] : range[0] };
            auto amplitude_hist{
                ROOTHelper::CreateHist2D(
                    name_amplitude.data(),"",
                    static_cast<int>(spot_count), -0.5, static_cast<int>(spot_count)-0.5,
                    100, std::get<0>(amplitude_range_tmp), std::get<1>(amplitude_range_tmp))
            };
            auto width_hist{
                ROOTHelper::CreateHist2D(
                    name_width.data(),"",
                    static_cast<int>(spot_count), -0.5, static_cast<int>(spot_count)-0.5,
                    100, std::get<0>(width_range_tmp), std::get<1>(width_range_tmp))
            };
            for (int k = 0; k < graph_map[i].at(spot)->GetN(); k++)
            {
                amplitude_hist->Fill(x_value, graph_map[i].at(spot)->GetPointX(k));
                width_hist->Fill(x_value, graph_map[i].at(spot)->GetPointY(k));
            }
            hist_map[i][1].emplace(spot, std::move(amplitude_hist));
            hist_map[i][0].emplace(spot, std::move(width_hist));
        }
    }


    for (size_t i = 0; i < model_list.size(); i++)
    {
        auto model_object{ model_list.at(i) };
        auto data_carbon{ width_array_map[i].at(Spot::CA) };
        data_carbon.insert(
            data_carbon.end(),
            width_array_map[i].at(Spot::C).begin(),
            width_array_map[i].at(Spot::C).end()
        );
        auto data_nitrogen{ width_array_map[i].at(Spot::N) };


        auto width_carbon_mean{
            ArrayStats<double>::ComputeMean(data_carbon.data(), data_carbon.size())
        };
        auto width_carbon_std{
            ArrayStats<double>::ComputeStandardDeviation(
                data_carbon.data(), data_carbon.size(), width_carbon_mean)
        };
        auto width_nitrogen_mean{
            ArrayStats<double>::ComputeMean(data_nitrogen.data(), data_nitrogen.size())
        };
        auto width_nitrogen_std{
            ArrayStats<double>::ComputeStandardDeviation(
                data_nitrogen.data(), data_nitrogen.size(), width_nitrogen_mean)
        };

        auto n1{ static_cast<double>(data_carbon.size()) };
        auto n2{ static_cast<double>(data_nitrogen.size()) };
        //auto s1_square{ width_carbon_std * width_carbon_std };
        //auto s2_square{ width_nitrogen_std * width_nitrogen_std };
        
        auto std_square_total{ 0.0 };
        for (auto & value : data_carbon)
        {
            std_square_total += (value - width_carbon_mean) * (value - width_carbon_mean);
        }
        for (auto & value : data_nitrogen)
        {
            std_square_total += (value - width_nitrogen_mean) * (value - width_nitrogen_mean);
        }
        std_square_total /= (n1 + n2 - 2);

        auto t_value{ (width_carbon_mean - width_nitrogen_mean) /
            std::sqrt(std_square_total / n1 + std_square_total / n2)
        };

        auto dof{ n1 + n2 - 2 };
        //auto dof{ std::pow(s1_square / n1 + s2_square / n2, 2) /
        //    ((s1_square * s1_square) / (n1 * n1 * (n1 - 1)) + (s2_square * s2_square) / (n2 * n2 * (n2 - 1)))
        //};

        auto p_value{ std::numeric_limits<double>::quiet_NaN() };
        if (std::isfinite(t_value) && dof > 0.0)
        {
            boost::math::students_t_distribution<double> dist(dof);
            p_value = 2.0 * boost::math::cdf(
                boost::math::complement(dist, std::fabs(t_value))
            );
        }

        std::ostringstream width_summary;
        width_summary << "Width estimate summary for model " << model_object->GetPdbID()
                      << ": t-value=" << t_value
                      << ", dof=" << dof
                      << ", p-value=" << p_value
                      << ", carbon mean=" << width_carbon_mean
                      << ", carbon sd=" << width_carbon_std
                      << ", carbon n=" << n1
                      << ", nitrogen mean=" << width_nitrogen_mean
                      << ", nitrogen sd=" << width_nitrogen_std
                      << ", nitrogen n=" << n2
                      << ", total sd=" << std::sqrt(std_square_total);
        Logger::Log(LogLevel::Debug, width_summary.str());
        
    }

    
    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> title_text[col_size];
    for (size_t i = 0; i < col_size; i++)
    {
        auto model_object{ model_list.at(i) };
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), static_cast<int>(i), j);
            ROOTHelper::SetPadLayout(gPad, 0, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto range_tmp{ (i == model_list.size() - 1) ? range_last[j] : range[j] };
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", static_cast<int>(i), j),"", 100, 0.0, 3.0, 500, std::get<0>(range_tmp), std::get<1>(range_tmp));
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 35.0f, 0.005f, 103, kCyan+3);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            if (i == col_size - 1) ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 0.0);
            else ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 35.0f, 1.5f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35.0f, 0.02f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle((j == 0) ? "Width" : "Amplitude");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->GetXaxis()->SetLimits(0.0, 3.0);
            frame[i][j]->GetYaxis()->SetLimits(std::get<0>(range_tmp), std::get<1>(range_tmp));
            GausPainter::RemodelAxisLabels(frame[i][j]->GetXaxis(), spot_label_list, 0.0, -1);
            frame[i][j]->Draw((i == col_size - 1) ? "Y+" : "");
            for (auto & [spot, hist] : hist_map[i][j])
            {
                auto spot_color{ AtomClassifier::GetMainChainSpotColor(spot) };
                auto spot_marker{ AtomClassifier::GetMainChainSpotOpenMarker(spot) };
                ROOTHelper::SetLineAttribute(hist.get(), 1, 1, spot_color);
                ROOTHelper::SetFillAttribute(hist.get(), 1001, spot_color, 0.3f);
                ROOTHelper::SetMarkerAttribute(hist.get(), spot_marker, 1.0f, spot_color);
                hist->SetBarWidth(0.6f);
                hist->Draw("VIOLIN SAME");
            }

            if (j == row_size - 1)
            {
                title_text[i] = ROOTHelper::CreatePaveText(0.02, 1.02, 0.98, 1.30, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
                ROOTHelper::SetPaveAttribute(title_text[i].get(), 0, 0.2);
                ROOTHelper::SetTextAttribute(title_text[i].get(), 35.0f, 133, 22);
                ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, kAzure-7, 0.5f);
                ROOTHelper::SetLineAttribute(title_text[i].get(), 1, 0);
                title_text[i]->AddText(Form("#font[102]{%s} (FSC = %.2f #AA)",
                    model_object->GetPdbID().data(), model_object->GetResolution()));
                title_text[i]->Draw();
            }
        }
    }
/*
    canvas->cd();
    auto pad_extra1{ ROOTHelper::CreatePad("pad_extra1","", 0.00, 0.00, 0.50, 0.35) };
    pad_extra1->Draw();
    pad_extra1->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra1.get());
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.09, 0.06, 0.11, 0.03);
    ROOTHelper::SetFillAttribute(pad_extra1.get(), 4000);
    auto frame_extra1{ ROOTHelper::CreateHist2D("frame_extra1","", 100, 0.0, 1.0, 100, 0.0, 1.0) };
    ROOTHelper::SetAxisTitleAttribute(frame_extra1->GetXaxis(), 35.0f, 0.9f, 133);
    ROOTHelper::SetAxisLabelAttribute(frame_extra1->GetXaxis(), 35.0f, 0.005f, 133);
    ROOTHelper::SetAxisTickAttribute(frame_extra1->GetXaxis(), 0.06f, 506);
    ROOTHelper::SetAxisTitleAttribute(frame_extra1->GetYaxis(), 35.0f, 1.5f, 133);
    ROOTHelper::SetAxisLabelAttribute(frame_extra1->GetYaxis(), 35.0f, 0.02f, 133);
    ROOTHelper::SetAxisTickAttribute(frame_extra1->GetYaxis(), 0.03f, 506);
    ROOTHelper::SetLineAttribute(frame_extra1.get(), 1, 0);
    frame_extra1->GetXaxis()->SetTitle("FSC Resolution (#AA)");
    frame_extra1->GetYaxis()->SetTitle("Width Mean");
    frame_extra1->GetXaxis()->CenterTitle();
    frame_extra1->GetYaxis()->CenterTitle();
    frame_extra1->SetStats(0);
    frame_extra1->GetXaxis()->SetLimits(1.0, 3.5);
    frame_extra1->GetYaxis()->SetLimits(0.5, 2.5);
    frame_extra1->Draw();
    for (auto & [spot, graph] : summary_graph1_map)
    {
        auto spot_color{ AtomClassifier::GetMainChainSpotColor(spot) };
        auto spot_marker{ AtomClassifier::GetMainChainSpotOpenMarker(spot) };
        auto spot_line_style{ AtomClassifier::GetMainChainSpotLineStyle(spot) };
        ROOTHelper::SetLineAttribute(graph.get(), spot_line_style, 1, spot_color);
        ROOTHelper::SetMarkerAttribute(graph.get(), spot_marker, 1.0f, spot_color);
        graph->Draw("PL");
    }

    canvas->cd();
    auto pad_extra2{ ROOTHelper::CreatePad("pad_extra2","", 0.50, 0.00, 1.00, 0.35) };
    pad_extra2->Draw();
    pad_extra2->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra2.get());
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.09, 0.06, 0.11, 0.03);
    ROOTHelper::SetFillAttribute(pad_extra2.get(), 4000);
    auto frame_extra2{ ROOTHelper::CreateHist2D("frame_extra2","", 100, 0.0, 1.0, 100, 0.0, 1.0) };
    ROOTHelper::SetAxisTitleAttribute(frame_extra2->GetXaxis(), 35.0f, 0.9f, 133);
    ROOTHelper::SetAxisLabelAttribute(frame_extra2->GetXaxis(), 35.0f, 0.005f, 133);
    ROOTHelper::SetAxisTickAttribute(frame_extra2->GetXaxis(), 0.06f, 506);
    ROOTHelper::SetAxisTitleAttribute(frame_extra2->GetYaxis(), 35.0f, 1.7f, 133);
    ROOTHelper::SetAxisLabelAttribute(frame_extra2->GetYaxis(), 35.0f, 0.02f, 133);
    ROOTHelper::SetAxisTickAttribute(frame_extra2->GetYaxis(), 0.03f, 506);
    ROOTHelper::SetLineAttribute(frame_extra2.get(), 1, 0);
    frame_extra2->GetXaxis()->SetTitle("FSC Resolution (#AA)");
    frame_extra2->GetYaxis()->SetTitle("Width S.D.");
    frame_extra2->GetXaxis()->CenterTitle();
    frame_extra2->GetYaxis()->CenterTitle();
    frame_extra2->SetStats(0);
    frame_extra2->GetXaxis()->SetLimits(1.0, 3.5);
    frame_extra2->GetYaxis()->SetLimits(0.0, 0.15);
    frame_extra2->Draw();
    for (auto & [spot, graph] : summary_graph2_map)
    {
        auto spot_color{ AtomClassifier::GetMainChainSpotColor(spot) };
        auto spot_marker{ AtomClassifier::GetMainChainSpotOpenMarker(spot) };
        auto spot_line_style{ AtomClassifier::GetMainChainSpotLineStyle(spot) };
        ROOTHelper::SetLineAttribute(graph.get(), spot_line_style, 1, spot_color);
        ROOTHelper::SetMarkerAttribute(graph.get(), spot_marker, 1.0f, spot_color);
        graph->Draw("PL");
    }

    auto legend{ ROOTHelper::CreateLegend(0.88, 0.35, 1.00, 0.93, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 25.0f, 133, 12);
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    legend->AddEntry(summary_graph1_map.at(Spot::CA).get(), "C_{#alpha}", "lp");
    legend->AddEntry(summary_graph1_map.at(Spot::C).get(), "C", "lp");
    legend->AddEntry(summary_graph1_map.at(Spot::N).get(), "N", "lp");
    legend->SetMargin(0.55f);
    legend->Draw();
*/

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

#ifdef HAVE_ROOT
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

void DemoPainter::PrintGausResultGlobalPad(
    TPad * pad, TH2 * hist,
    double left_margin, double right_margin, double bottom_margin, double top_margin,
    bool is_right_side_pad)
{
    pad->cd();
    ROOTHelper::SetPadMarginInCanvas(pad, left_margin, right_margin, bottom_margin, top_margin);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 55.0f, 0.11f, 103, kCyan+3);
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
    ROOTHelper::SetPadFrameAttribute(pad, 0, 0, 4000, 0, 0, 0, 0);
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

void DemoPainter::PrintGausCorrelationPad(TPad * pad, TH2 * hist, bool draw_x_axis, bool draw_title_label)
{
    pad->cd();
    auto bottom_margin{ (draw_x_axis) ? 0.105 : 0.010 };
    auto top_margin{ (draw_title_label) ? 0.105 : 0.010 };
    auto title_size{ (draw_x_axis) ? 80.0f : 0.0f };
    //auto label_size{ (draw_x_axis) ? 50.0f : 12.0f };
    ROOTHelper::SetPadMarginInCanvas(pad, 0.005, 0.005, bottom_margin, top_margin);
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(pad, 0, 0, 4000, 0, 0, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist->GetXaxis(), title_size, 0.5f);
    ROOTHelper::SetAxisTitleAttribute(hist->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(hist->GetXaxis(), 0.0f);
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

void DemoPainter::BuildMapValueScatterGraph(
    GroupKey group_key, TGraphErrors * graph, ModelObject * model1, ModelObject * model2,
    int bin_size, double x_min, double x_max)
{
    auto entry1_iter{ std::make_unique<PotentialEntryQuery>(model1) };
    auto plot_builder1{ std::make_unique<PotentialPlotBuilder>(model1) };
    auto entry2_iter{ std::make_unique<PotentialEntryQuery>(model2) };
    auto plot_builder2{ std::make_unique<PotentialPlotBuilder>(model2) };
    if (entry1_iter->IsAvailableAtomGroupKey(group_key, ChemicalDataHelper::GetSimpleAtomClassKey()) == false) return;
    if (entry2_iter->IsAvailableAtomGroupKey(group_key, ChemicalDataHelper::GetSimpleAtomClassKey()) == false) return;
    auto model1_atom_map{ entry1_iter->GetAtomObjectMap(group_key, ChemicalDataHelper::GetSimpleAtomClassKey()) };
    auto model2_atom_map{ entry2_iter->GetAtomObjectMap(group_key, ChemicalDataHelper::GetSimpleAtomClassKey()) };
    auto count{ 0 };
    for (auto & [atom_id, atom_object1] : model1_atom_map)
    {
        if (model2_atom_map.find(atom_id) == model2_atom_map.end()) continue;
        auto atom_object2{ model2_atom_map.at(atom_id) };
        auto atom1_iter{ std::make_unique<PotentialEntryQuery>(atom_object1) };
        auto atom_plot_builder1{ std::make_unique<PotentialPlotBuilder>(atom_object1) };
        auto atom2_iter{ std::make_unique<PotentialEntryQuery>(atom_object2) };
        auto atom_plot_builder2{ std::make_unique<PotentialPlotBuilder>(atom_object2) };
        auto data1_array{ atom1_iter->GetBinnedDistanceAndMapValueList(bin_size, x_min, x_max) };
        auto data2_array{ atom2_iter->GetBinnedDistanceAndMapValueList(bin_size, x_min, x_max) };
        for (size_t i = 0; i < static_cast<size_t>(bin_size); i++)
        {
            auto x_value{ std::get<1>(data1_array.at(i)) };
            auto y_value{ std::get<1>(data2_array.at(i)) };
            graph->SetPoint(count, x_value, y_value);
            count++;
        }
    }
}

#endif


} // namespace rhbm_gem
