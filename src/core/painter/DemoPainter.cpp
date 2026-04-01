#include <rhbm_gem/core/painter/DemoPainter.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/LocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelPotentialView.hpp>
#include "PotentialPlotBuilder.hpp"
#include <detail/PotentialSeriesOps.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/core/painter/GausPainter.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include "detail/PainterTypeCheck.hpp"
#include "detail/PainterSupport.hpp"

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

DemoPainter::DemoPainter() :
    m_atom_classifier{ std::make_unique<AtomClassifier>() }
{

}

DemoPainter::~DemoPainter()
{

}

void DemoPainter::AddDataObject(DataObjectBase * data_object)
{
    painter_internal::AppendPainterObject<ModelObject>(
        data_object, "DemoPainter", "AddDataObject", m_model_object_list);
}

void DemoPainter::AddReferenceDataObject(DataObjectBase * data_object, const std::string & label)
{
    painter_internal::AppendPainterReferenceObject<ModelObject>(
        data_object,
        label,
        "DemoPainter",
        "AddReferenceDataObject",
        m_ref_model_object_list_map);
}

void DemoPainter::Painting()
{
    Logger::Log(LogLevel::Info, "DemoPainter::Painting() called.");
    Logger::Log(LogLevel::Info, "Folder path: " + m_folder_path);

    for (auto & model : m_model_object_list)
    {
        model->BuildKDTreeRoot();
    }

    ModelObject * demo_model_object{ nullptr };
    std::vector<ModelObject *> demo_model_list{ nullptr, nullptr, nullptr, nullptr };
    std::vector<ModelObject *> demo_fsc_model_list
    {
        nullptr, nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr, nullptr,
        nullptr
    };
    std::vector<ModelObject *> demo_alpha_carbon_list
    {
        nullptr, nullptr, nullptr
    };
    for (auto model : m_model_object_list)
    {
        if (model->GetPdbID() == "6Z6U" && model->GetEmdID() == "EMD-11103")
        {
            demo_model_object = model;
            demo_model_list[3] = model;
            demo_fsc_model_list[0] = model;
            demo_alpha_carbon_list[0] = model;
        }
        else if (model->GetPdbID() == "8DQV" && model->GetEmdID() == "EMD-27661")
        {
            demo_model_list[2] = model;
            demo_fsc_model_list[1] = model;
        }
        else if (model->GetPdbID() == "9EVX" && model->GetEmdID() == "EMD-50019")
        {
            demo_model_list[1] = model;
            demo_fsc_model_list[2] = model;
        }
        else if (model->GetPdbID() == "6CVM" && model->GetEmdID() == "EMD-7770" )
        {
            demo_model_list[0] = model;
            demo_fsc_model_list[3] = model;
        }
        else if (model->GetPdbID() == "8RQB" && model->GetEmdID() == "EMD-19436")
        {
            demo_fsc_model_list[4] = model;
        }
        else if (model->GetPdbID() == "7A6A" && model->GetEmdID() == "EMD-11668")
        {
            demo_fsc_model_list[5] = model;
        }
        else if (model->GetPdbID() == "7A6B" && model->GetEmdID() == "EMD-11669")
        {
            demo_fsc_model_list[6] = model;
        }
        else if (model->GetPdbID() == "6Z9E" && model->GetEmdID() == "EMD-11121")
        {
            demo_fsc_model_list[7] = model;
        }
        else if (model->GetPdbID() == "6Z9F" && model->GetEmdID() == "EMD-11122")
        {
            demo_fsc_model_list[8] = model;
        }
        else if (model->GetPdbID() == "7KOD" && model->GetEmdID() == "EMD-22972")
        {
            demo_fsc_model_list[9] = model;
        }
        else if (model->GetPdbID() == "6DRV" && model->GetEmdID() == "EMD-8908" )
        {
            demo_fsc_model_list[10] = model;
            demo_alpha_carbon_list[1] = model; 
        }
        else if (model->GetPdbID() == "9GXM" && model->GetEmdID() == "EMD-51667" )
        {
            demo_alpha_carbon_list[2] = model;
        }
    };

    //PaintGroupGausMergeResidueDemo(demo_alpha_carbon_list, "figure_gaus_backbone_boxplot.pdf");

    PaintAtomMapValueExample(demo_model_object, "figure_1_a.pdf");
    PaintGroupGausMainChainSummary(demo_model_list, "figure_1_b.pdf");
    PaintAtomGausMainChainDemoSingle(demo_model_object, "figure_2_c1.pdf", 0);

    PaintGroupGausToFSC(demo_fsc_model_list, "figure_3_a.pdf");
    PaintAtomWidthScatterPlotSingle(demo_model_object, "figure_3_b.pdf", true);
    PaintGroupWidthScatterPlot(demo_model_list, "figure_3_c.pdf", 1, true);
    
    for (auto & [key, model_list] : m_ref_model_object_list_map)
    {
        if (key != "with_charge") continue;
        for (auto & model : model_list)
        {
            PaintAtomGausMainChainDemoSingle(model, "figure_2_c2.pdf", 0);
            PaintGroupGausMainChainSingle(model, "figure_2_d.pdf");
            PainMapValueComparisonSingle("figure_5.pdf", demo_model_object, model);
        }
    }
}

void DemoPainter::PainMapValueComparisonSingle(
    const std::string & name,
    ModelObject * model_object,
    ModelObject * ref_model_object)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " DemoPainter::PainMapValueComparisonSingle");
    (void)model_object;
    (void)ref_model_object;

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 4 };
    const int row_size{ 1 };

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 600) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.08f, 0.02f, 0.20f, 0.12f, 0.01f, 0.005f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::vector<std::unique_ptr<TGraphErrors>> scatter_graph_list;
    std::vector<std::unique_ptr<TF1>> fit_function_list;

    TGraphErrors * scatter_graph[col_size];
    TF1 * fit_function[col_size];
    double r_square[col_size]{ 0.0 };
    double slope[col_size]{ 0.0 };
    double intercept[col_size]{ 0.0 };
    std::vector<double> x_array, y_array;
    for (size_t i = 0; i < col_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainSimpleAtomClassGroupKey(i) };
        auto graph{ ROOTHelper::CreateGraphErrors() };
        painter_internal::BuildMapValueScatterGraph(
            group_key, graph.get(), ref_model_object, model_object, 15, 0.0, 1.5);
        r_square[i] = ROOTHelper::PerformLinearRegression(graph.get(), slope[i], intercept[i]);
        auto function{ ROOTHelper::CreateFunction1D(Form("fit_%d", static_cast<int>(i)), "x*[1]+[0]") };
        function->SetParameters(intercept[i], slope[i]);
        scatter_graph[i] = graph.get();
        fit_function[i] = function.get();
        scatter_graph_list.emplace_back(std::move(graph));
        fit_function_list.emplace_back(std::move(function));

        for (int p = 0; p < scatter_graph[i]->GetN(); p++)
        {
            x_array.emplace_back(scatter_graph[i]->GetPointX(p));
            y_array.emplace_back(scatter_graph[i]->GetPointY(p));
        }
    }

    auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.05) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.2) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> title_text[col_size];
    std::unique_ptr<TPaveText> r_square_text[col_size];
    std::unique_ptr<TPaveText> fit_info_text[col_size];
    for (size_t i = 0; i < col_size; i++)
    {
        auto element_color{ AtomClassifier::GetMainChainElementColor(i) };
        auto element_marker{ AtomClassifier::GetMainChainElementOpenMarker(i) };
        auto element_label{ AtomClassifier::GetMainChainElementLabel(i) };
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), static_cast<int>(i), j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", static_cast<int>(i), j),"", 500, x_min, x_max, 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 50.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 60.0f, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 50.0f, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.08/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("Real Map Value");
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            ROOTHelper::SetMarkerAttribute(scatter_graph[i], element_marker, 1.0f, element_color);
            ROOTHelper::SetLineAttribute(scatter_graph[i], 1, 2, element_color);
            scatter_graph[i]->Draw("P");
            ROOTHelper::SetLineAttribute(fit_function[i], 2, 2, kRed);
            fit_function[i]->SetRange(x_min, x_max);
            fit_function[i]->Draw("SAME");
        }
        title_text[i] = ROOTHelper::CreatePaveText(0.30, 1.01, 0.70, 1.16, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
        ROOTHelper::SetPaveAttribute(title_text[i].get(), 0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[i].get(), 60.0f, 103, 22);
        ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, element_color, 0.5f);
        title_text[i]->AddText(element_label.data());
        title_text[i]->Draw();

        r_square_text[i] = ROOTHelper::CreatePaveText(0.30, 0.05, 0.95, 0.18, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(r_square_text[i].get());
        ROOTHelper::SetPaveAttribute(r_square_text[i].get(), 0, 0.5);
        ROOTHelper::SetLineAttribute(r_square_text[i].get(), 1, 0);
        ROOTHelper::SetTextAttribute(r_square_text[i].get(), 45.0f, 133, 22);
        ROOTHelper::SetFillAttribute(r_square_text[i].get(), 1001, kAzure-7, 0.20f);
        r_square_text[i]->AddText(Form("R^{2} = %.2f", r_square[i]));
        r_square_text[i]->Draw();

        fit_info_text[i] = ROOTHelper::CreatePaveText(0.05, 0.88, 0.95, 0.99, "nbNDC", true);
        ROOTHelper::SetPaveTextDefaultStyle(fit_info_text[i].get());
        ROOTHelper::SetTextAttribute(fit_info_text[i].get(), 45.0f, 133, 22, 0.0, kGray+2);
        if (intercept[i] > 0.0)
        {
            fit_info_text[i]->AddText(Form("#font[1]{y} = %.2f#font[1]{x}+ %.2f", slope[i], intercept[i]));
        }
        else
        {
            fit_info_text[i]->AddText(Form("#font[1]{y} = %.2f#font[1]{x}#minus %.2f", slope[i], std::fabs(intercept[i])));
        }
        fit_info_text[i]->Draw();
    }
    canvas->cd();
    auto pad_extra{ ROOTHelper::CreatePad("pad_extra","", 0.08, 0.00, 0.98, 0.12) };
    pad_extra->Draw();
    pad_extra->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra.get());
    ROOTHelper::SetFillAttribute(pad_extra.get(), 4000);
    auto bottom_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(bottom_title_text.get());
    ROOTHelper::SetFillAttribute(bottom_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(bottom_title_text.get(), 50.0f, 133, 22);
    bottom_title_text->AddText("Simulated Map Value");
    bottom_title_text->Draw();
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);

    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void DemoPainter::PaintAtomMapValueExample(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    Logger::Log(LogLevel::Info, " DemoPainter::PaintAtomMapValueExample");

    if (model_object == nullptr) return;
    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };

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
    AtomClassifier classifier;
    auto group_key{ classifier.GetMainChainComponentAtomClassGroupKey(0, Residue::ALA) };
    if (!entry_iter->HasAtomGroup(group_key, class_key)) return;
    for (auto atom : entry_iter->GetAtomObjectList(group_key, class_key))
    {
        auto atom_iter{ std::make_unique<LocalPotentialView>(*atom) };
        auto atom_plot_builder{ std::make_unique<PotentialPlotBuilder>(atom) };
        auto graph{ atom_plot_builder->CreateBinnedDistanceToMapValueGraph() };
        ROOTHelper::SetLineAttribute(graph.get(), 1, 5, static_cast<short>(kAzure-7), 0.3f);
        map_value_graph_list.emplace_back(std::move(graph));
        auto map_value_range{ series_ops::ComputeMapValueRange(*atom_iter, 0.0) };
        y_array.emplace_back(std::get<0>(map_value_range));
        y_array.emplace_back(std::get<1>(map_value_range));
    }
    gaus_function = plot_builder->CreateAtomGroupGausFunctionPrior(group_key, class_key);
    amplitude_prior = entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 0);
    width_prior = entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 1);


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
    frame->GetYaxis()->SetTitle("Map Value");
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
    ROOTHelper::SetTextAttribute(legend.get(), 65.0f, 133, 12, 0.0);
    legend->SetMargin(0.15f);
    legend->AddEntry(gaus_function.get(),
        "Gaussian Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})}", "l");
    legend->AddEntry(map_value_graph_list.at(0).get(),
        "Map Value", "l");
    legend->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void DemoPainter::PaintGroupGausMainChainSummary(
    const std::vector<ModelObject *> & model_list, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    auto residue_class{ ChemicalDataHelper::GetComponentAtomClassKey() };
    Logger::Log(LogLevel::Info, " DemoPainter::PaintGroupGausMainChainSummary");

    for (auto & model : model_list) if (model == nullptr) return;

    #ifdef HAVE_ROOT
    const int main_chain_element_count{ 4 };
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
        auto model_object{ model_list.at(j) };
        std::vector<double> amplitude_array, width_array;
        info_text[j] = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
        resolution_text[j] = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
        auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
        auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
        for (size_t k = 0; k < main_chain_element_count; k++)
        {
            auto group_key_list{ m_atom_classifier->GetMainChainComponentAtomClassGroupKeyList(k) };
            amplitude_graph[j][k] = plot_builder->CreateAtomGausEstimateToResidueGraph(group_key_list, residue_class, 0);
            width_graph[j][k] = plot_builder->CreateAtomGausEstimateToResidueGraph(group_key_list, residue_class, 1);
            correlation_graph[j][k] = plot_builder->CreateAtomGausEstimateScatterGraph(group_key_list, residue_class, 1, 0);
            for (int p = 0; p < amplitude_graph[j][k]->GetN(); p++)
            {
                amplitude_array.push_back(amplitude_graph[j][k]->GetPointY(p));
                width_array.push_back(width_graph[j][k]->GetPointY(p));
                width_array_total.push_back(width_graph[j][k]->GetPointY(p));
            }
            auto element_color{ AtomClassifier::GetMainChainElementColor(k) };
            auto element_marker{ AtomClassifier::GetMainChainElementOpenMarker(k) };
            ROOTHelper::SetMarkerAttribute(amplitude_graph[j][k].get(), element_marker, marker_size, element_color);
            ROOTHelper::SetMarkerAttribute(width_graph[j][k].get(), element_marker, marker_size, element_color);
            ROOTHelper::SetMarkerAttribute(correlation_graph[j][k].get(), element_marker, marker_size, element_color);
            ROOTHelper::SetLineAttribute(amplitude_graph[j][k].get(), 1, 1, element_color);
            ROOTHelper::SetLineAttribute(width_graph[j][k].get(), 1, 1, element_color);
            ROOTHelper::SetLineAttribute(correlation_graph[j][k].get(), 1, 1, element_color);
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

    for (int j = 0; j < pad_size_y; j++)
    {
        auto model_object{ model_list.at(static_cast<size_t>(j)) };
        auto draw_axis_flag{ (j == 0) ? true : false };
        auto draw_title_flag{ (j == pad_size_y - 1) ? true : false };

        pad[0][j]->cd();
        auto pdb_id{ model_object->GetPdbID() };
        auto emd_id{ model_object->GetEmdID() };
        auto resolution{ model_object->GetResolution() };
    
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
        if (j == pad_size_y - 1) painter_internal::PrintGausTitlePad(
            pad[1][j].get(), amplitude_title_text.get(), "Amplitude #font[2]{A}", 80.0f);

        pad[2][j]->cd();
        PrintGausCorrelationPad(pad[2][j].get(), frame[1][j].get(), draw_axis_flag, draw_title_flag);
        for (int k = 0; k < main_chain_element_count; k++) correlation_graph[j][k]->Draw("P X0");
        if (j == pad_size_y - 1) painter_internal::PrintGausTitlePad(
            pad[2][j].get(), correlation_title_text.get(), "#tau#minus#font[2]{A} Plot", 80.0f);

        pad[3][j]->cd();
        PrintGausResultPad(pad[3][j].get(), frame[2][j].get(), draw_axis_flag, draw_title_flag, true);
        for (int k = 0; k < main_chain_element_count; k++) width_graph[j][k]->Draw("PL X0");
        if (j == pad_size_y - 1) painter_internal::PrintGausTitlePad(
            pad[3][j].get(), width_title_text.get(), "Width #font[2]{#tau}", 80.0f);
    }

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void DemoPainter::PaintGroupGausMainChainSingle(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    auto residue_class{ ChemicalDataHelper::GetComponentAtomClassKey() };
    Logger::Log(LogLevel::Info, " DemoPainter::PaintGroupGausMainChainSingle");

    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };

    #ifdef HAVE_ROOT

    const int main_chain_element_count{ 4 };
    float marker_size{ 1.5f };

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 3000, 500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size{ 4 };
    std::unique_ptr<TPad> pad[pad_size];
    std::unique_ptr<TH2> frame[pad_size];
    
    const double x_pos[pad_size + 1]{ 0.00, 0.16, 0.52, 0.64, 1.00 };
    for (int i = 0; i < pad_size; i++)
    {
        auto pad_name{ "pad"+ std::to_string(i) };
        pad[i] = ROOTHelper::CreatePad(pad_name.data(),"", x_pos[i], 0.0, x_pos[i+1], 1.0);
    }

    std::unique_ptr<TPaveText> info_text;
    std::unique_ptr<TPaveText> resolution_text;
    std::unique_ptr<TGraphErrors> amplitude_graph[main_chain_element_count];
    std::unique_ptr<TGraphErrors> width_graph[main_chain_element_count];
    std::unique_ptr<TGraphErrors> correlation_graph[main_chain_element_count];
    std::vector<double> amplitude_array, width_array;
    info_text = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
    resolution_text = ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false);
    for (size_t k = 0; k < main_chain_element_count; k++)
    {
        auto group_key_list{ m_atom_classifier->GetMainChainComponentAtomClassGroupKeyList(k) };
        amplitude_graph[k] = plot_builder->CreateAtomGausEstimateToResidueGraph(group_key_list, residue_class, 0);
        width_graph[k] = plot_builder->CreateAtomGausEstimateToResidueGraph(group_key_list, residue_class, 1);
        correlation_graph[k] = plot_builder->CreateAtomGausEstimateScatterGraph(group_key_list, residue_class, 1, 0);
        for (int p = 0; p < amplitude_graph[k]->GetN(); p++)
        {
            amplitude_array.push_back(amplitude_graph[k]->GetPointY(p));
            width_array.push_back(width_graph[k]->GetPointY(p));
        }
        auto element_color{ AtomClassifier::GetMainChainElementColor(k) };
        auto element_marker{ AtomClassifier::GetMainChainElementOpenMarker(k) };
        ROOTHelper::SetMarkerAttribute(amplitude_graph[k].get(), element_marker, marker_size, element_color);
        ROOTHelper::SetMarkerAttribute(width_graph[k].get(), element_marker, marker_size, element_color);
        ROOTHelper::SetMarkerAttribute(correlation_graph[k].get(), element_marker, marker_size, element_color);
        ROOTHelper::SetLineAttribute(amplitude_graph[k].get(), 1, 1, element_color);
        ROOTHelper::SetLineAttribute(width_graph[k].get(), 1, 1, element_color);
        ROOTHelper::SetLineAttribute(correlation_graph[k].get(), 1, 1, element_color);
    }
    auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, 0.2) };
    auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, 0.1) };

    frame[0] = ROOTHelper::CreateHist2D("frame0","", 100, 0.0, 1.0, 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    frame[1] = ROOTHelper::CreateHist2D("frame1","", 100, std::get<0>(width_range), std::get<1>(width_range), 100, std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    frame[2] = ROOTHelper::CreateHist2D("frame2","", 100, 0.0, 1.0, 100, std::get<0>(width_range), std::get<1>(width_range));

    auto amplitude_title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    auto width_title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    auto correlation_title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };

    canvas->cd();
    for (int i = 0; i < pad_size; i++)
    {
        ROOTHelper::SetPadDefaultStyle(pad[i].get());
        pad[i]->Draw();
    }

    pad[0]->cd();
    auto pdb_id{ model_object->GetPdbID() };
    auto emd_id{ model_object->GetEmdID() };

    // PDB & EMD Text
    ROOTHelper::SetPaveTextMarginInCanvas(pad[0].get(), info_text.get(), 0.003, 0.00, 0.25, 0.25);
    ROOTHelper::SetPaveTextDefaultStyle(info_text.get());
    ROOTHelper::SetPaveAttribute(info_text.get(), 0, 0.1);
    ROOTHelper::SetFillAttribute(info_text.get(), 1001, kAzure-7, 0.3f);
    ROOTHelper::SetTextAttribute(info_text.get(), 75.0f, 103, 32);
    ROOTHelper::SetLineAttribute(info_text.get(), 1, 0);
    info_text->AddText(pdb_id.data());
    info_text->AddText(emd_id.data());
    info_text->Draw();
    
    // Resolution Text
    ROOTHelper::SetPaveTextMarginInCanvas(pad[0].get(), resolution_text.get(), 0.001, 0.07, 0.50, 0.20);
    ROOTHelper::SetPaveTextDefaultStyle(resolution_text.get());
    ROOTHelper::SetPaveAttribute(resolution_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(resolution_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(resolution_text.get(), 80.0f, 133, 22, 0.0, kYellow-10);
    auto resolution{ model_object->GetResolution() };
    if (emd_id == "Simulation")
    {
        resolution_text->AddText(Form("#sigma_{B} = %.1f", resolution));
    }
    else
    {
        resolution_text->AddText(Form("%.2f #AA", resolution));
    }
    resolution_text->Draw();

    auto legend{ ROOTHelper::CreateLegend(0.15, 0.00, 1.00, 0.20, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 70.0f, 103, 12);
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    for (size_t k = 0; k < main_chain_element_count; k++)
    {
        auto element_color{ AtomClassifier::GetMainChainElementColor(k) };
        auto element_label{ AtomClassifier::GetMainChainElementLabel(k) };
        auto label{ Form("#color[%d]{%s}", element_color, element_label.data()) };
        legend->AddEntry(correlation_graph[k].get(), label, "p");
    }
    legend->SetNColumns(4);
    //legend->Draw();

    pad[1]->cd();
    painter_internal::PrintGausResultGlobalPad(
        pad[1].get(), frame[0].get(), 0.030, 0.005, 0.205, 0.205, false);
    for (int k = 0; k < main_chain_element_count; k++) amplitude_graph[k]->Draw("PL X0");
    painter_internal::PrintGausTitlePad(
        pad[1].get(), amplitude_title_text.get(), "Amplitude #font[2]{A}", 80.0f);

    pad[2]->cd();
    ROOTHelper::SetPadMarginInCanvas(pad[2].get(), 0.005, 0.005, 0.205, 0.205);
    ROOTHelper::SetPadLayout(pad[2].get(), 1, 1, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(frame[1]->GetXaxis(), 60.0f, 0.75f);
    ROOTHelper::SetAxisTitleAttribute(frame[1]->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(frame[1]->GetXaxis(), 47.0f);
    ROOTHelper::SetAxisLabelAttribute(frame[1]->GetYaxis(), 0.0f);
    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad[2].get(), 0.040, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(pad[2].get(), 0.008, 1) };
    ROOTHelper::SetAxisTickAttribute(frame[1]->GetXaxis(), static_cast<float>(x_tick_length), 505);
    ROOTHelper::SetAxisTickAttribute(frame[1]->GetYaxis(), static_cast<float>(y_tick_length), 505);
    frame[1]->SetStats(0);
    frame[1]->GetXaxis()->SetTitle("Width #font[2]{#tau}");
    frame[1]->GetXaxis()->CenterTitle();
    frame[1]->Draw("");
    for (int k = 0; k < main_chain_element_count; k++) correlation_graph[k]->Draw("P X0");
    painter_internal::PrintGausTitlePad(
        pad[2].get(), correlation_title_text.get(), "#tau#minus#font[2]{A} Plot", 80.0f);

    pad[3]->cd();
    painter_internal::PrintGausResultGlobalPad(
        pad[3].get(), frame[2].get(), 0.005, 0.030, 0.205, 0.205, true);
    for (int k = 0; k < main_chain_element_count; k++) width_graph[k]->Draw("PL X0");
    painter_internal::PrintGausTitlePad(
        pad[3].get(), width_title_text.get(), "Width #font[2]{#tau}", 80.0f);

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

} // namespace rhbm_gem
#include <rhbm_gem/core/painter/DemoPainter.hpp>
#include "PotentialPlotBuilder.hpp"
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/LocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelPotentialView.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#endif

#include <vector>

namespace rhbm_gem {

void DemoPainter::PaintGroupGausToFSC(
    const std::vector<ModelObject *> & model_list, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " DemoPainter::PaintGroupGausToFSC");

    for (auto & model : model_list)
    {
        if (model == nullptr) return;
        Logger::Log(LogLevel::Debug, "Find Model : "+ model->GetPdbID());
        
    }

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 4 };
    const int row_size{ 1 };

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 300) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.08f, 0.01f, 0.30f, 0.02f, 0.01f, 0.005f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TGraphErrors> graph[col_size];
    std::unique_ptr<TF1> fit_function[col_size];
    std::vector<double> x_array, y_array;
    double r_square[col_size]{ 0.0 };
    double slope[col_size]{ 0.0 };
    double intercept[col_size]{ 0.0 };
    for (size_t i = 0; i < AtomClassifier::GetMainChainMemberCount(); i++)
    {
        if (i >= col_size) continue;
        auto group_key{ m_atom_classifier->GetMainChainSimpleAtomClassGroupKey(i) };
        graph[i] = ROOTHelper::CreateGraphErrors();
        auto count{ 0 };
        for (auto model : model_list)
        {
            auto entry_iter{ std::make_unique<ModelPotentialView>(*model) };
            auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model) };
            auto width_value{ entry_iter->GetAtomGausEstimatePrior(group_key, ChemicalDataHelper::GetSimpleAtomClassKey(), 1) };
            auto width_error{ entry_iter->GetAtomGausVariancePrior(group_key, ChemicalDataHelper::GetSimpleAtomClassKey(), 1) };
            graph[i]->SetPoint(count, model->GetResolution(), width_value);
            graph[i]->SetPointError(count, 0.0, width_error);
            count++;
        }
        for (int p = 0; p < graph[i]->GetN(); p++)
        {
            x_array.push_back(graph[i]->GetPointX(p));
            y_array.push_back(graph[i]->GetPointY(p));
        }
        auto element_color{ AtomClassifier::GetMainChainElementColor(i) };
        ROOTHelper::SetMarkerAttribute(graph[i].get(), 20, 1.2f, element_color);
        ROOTHelper::SetLineAttribute(graph[i].get(), 1, 2, element_color);

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
    for (size_t i = 0; i < col_size; i++)
    {
        auto element_color{ AtomClassifier::GetMainChainElementColor(i) };
        auto element_label{ AtomClassifier::GetMainChainElementLabel(i) };
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), static_cast<int>(i), j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", static_cast<int>(i), j),"", 500, x_min, x_max, 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 45.0f, 1.3f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 45.0f, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.08/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("Width #tau");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            graph[i]->Draw("P X0");
            fit_function[i]->SetRange(x_min, x_max);
            fit_function[i]->Draw("SAME");
        }
        title_text[i] = ROOTHelper::CreatePaveText(0.10, 0.75, 0.40, 0.99, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
        ROOTHelper::SetPaveAttribute(title_text[i].get(), 0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[i].get(), 60.0f, 103, 22);
        ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, element_color, 0.5f);
        title_text[i]->AddText(element_label.data());
        title_text[i]->Draw();

        r_square_text[i] = ROOTHelper::CreatePaveText(0.48, 0.05, 0.99, 0.25, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(r_square_text[i].get());
        ROOTHelper::SetPaveAttribute(r_square_text[i].get(), 0, 0.5);
        ROOTHelper::SetLineAttribute(r_square_text[i].get(), 1, 0);
        ROOTHelper::SetTextAttribute(r_square_text[i].get(), 40.0f, 133, 22);
        ROOTHelper::SetFillAttribute(r_square_text[i].get(), 1001, kAzure-7, 0.20f);
        r_square_text[i]->AddText(Form("R^{2} = %.2f", r_square[i]));
        r_square_text[i]->Draw();
    }
    canvas->cd();
    auto pad_extra{ ROOTHelper::CreatePad("pad_extra","", 0.08, 0.00, 0.98, 0.19) };
    pad_extra->Draw();
    pad_extra->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra.get());
    ROOTHelper::SetFillAttribute(pad_extra.get(), 4000);
    auto bottom_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(bottom_title_text.get());
    ROOTHelper::SetFillAttribute(bottom_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(bottom_title_text.get(), 40.0f, 133, 22);
    bottom_title_text->AddText("FSC Resolution #[]{#AA}");
    bottom_title_text->Draw();
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void DemoPainter::PaintAtomWidthScatterPlotSingle(
    ModelObject * model_object, const std::string & name, bool draw_box_plot)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " DemoPainter::PaintAtomWidthScatterPlotSingle");
    auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    (void)draw_box_plot;

    if (model_object == nullptr) return;
    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    const int pad_size{ 5 };

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 400) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), pad_size, 1, 0.08f, 0.08f, 0.25f, 0.13f, 0.011f, 0.011f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    auto normalized_z_position{ 0.5 };
    auto z_ratio_window{ 0.1 };
    auto graph_position{
        plot_builder->CreateAtomXYPositionTomographyGraph(normalized_z_position, z_ratio_window, true)
    };
    auto z_position{ model_object->GetModelPosition(2, normalized_z_position) };
    auto z_window{ model_object->GetModelLength(2) * z_ratio_window };
    auto com_pos{ model_object->GetCenterOfMassPosition() };
    std::vector<double> x_pos_array, y_pos_array;
    x_pos_array.reserve(static_cast<size_t>(graph_position->GetN()));
    y_pos_array.reserve(static_cast<size_t>(graph_position->GetN()));
    for (int p = 0; p < graph_position->GetN(); p++)
    {
        x_pos_array.emplace_back(graph_position->GetPointX(p));
        y_pos_array.emplace_back(graph_position->GetPointY(p));
    }

    auto x_pos_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_pos_array, 0.1) };
    double x_pos_min{ std::get<0>(x_pos_range) };
    double x_pos_max{ std::get<1>(x_pos_range) };

    auto y_pos_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_pos_array, 0.1) };
    double y_pos_min{ std::get<0>(y_pos_range) };
    double y_pos_max{ std::get<1>(y_pos_range) };

    std::vector<std::unique_ptr<TGraphErrors>> graph_list[pad_size - 1];
    std::vector<double> x_array[pad_size - 1];
    std::vector<double> y_array;
    for (size_t i = 0; i < pad_size - 1; i++)
    {
        for (auto residue : ChemicalDataHelper::GetStandardAminoAcidList())
        {
            auto group_key{ m_atom_classifier->GetMainChainComponentAtomClassGroupKey(i, residue) };
            if (!entry_iter->HasAtomGroup(group_key, class_key)) continue;
            auto gaus_graph{ plot_builder->CreateCOMDistanceToGausEstimateGraph(group_key, class_key, 1) };
            for (int p = 0; p < gaus_graph->GetN(); p++)
            {
                x_array[i].emplace_back(gaus_graph->GetPointX(p));
                y_array.emplace_back(gaus_graph->GetPointY(p));
            }
            graph_list[i].emplace_back(std::move(gaus_graph));
        }
    }

    double x_min[pad_size - 1];
    double x_max[pad_size - 1];
    for (int i = 0; i < pad_size - 1; i++)
    {
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array[i], 0.1) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);
    }

    auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array, 0.1) };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };

    std::unique_ptr<TH2D> summary_hist[pad_size - 1];
    for (int i = 0; i < pad_size - 1; i++)
    {
        summary_hist[i] = ROOTHelper::CreateHist2D(
            Form("hist_%d", i), "", 5, 35.0, 65.0, 100, y_min, y_max);
        for (auto & graph : graph_list[i])
        {
            for (int p = 0; p < graph->GetN(); p++)
            {
                summary_hist[i]->Fill(graph->GetPointX(p), graph->GetPointY(p));
            }
        }
    }

    std::unique_ptr<TH2> frame[pad_size];
    ROOTHelper::FindPadInCanvasPartition(canvas.get(), 0, 0);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    auto x_factor_tmp{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
    auto y_factor_tmp{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
    frame[0] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", 0, 0),"", 500, x_pos_min, x_pos_max, 500, y_pos_min, y_pos_max);
    ROOTHelper::SetAxisTitleAttribute(frame[0]->GetXaxis(), 50.0f, 0.8f, 133);
    ROOTHelper::SetAxisLabelAttribute(frame[0]->GetXaxis(), 45.0f, 0.005f, 133);
    ROOTHelper::SetAxisTickAttribute(frame[0]->GetXaxis(), static_cast<float>(y_factor_tmp*0.05/x_factor_tmp), 506);
    ROOTHelper::SetAxisTitleAttribute(frame[0]->GetYaxis(), 50.0f, 1.1f, 133);
    ROOTHelper::SetAxisLabelAttribute(frame[0]->GetYaxis(), 45.0f, 0.01f, 133);
    ROOTHelper::SetAxisTickAttribute(frame[0]->GetYaxis(), static_cast<float>(x_factor_tmp*0.05/y_factor_tmp), 506);
    ROOTHelper::SetLineAttribute(frame[0].get(), 1, 0);
    frame[0]->GetXaxis()->SetTitle("#font[2]{x} #[]{#AA}");
    frame[0]->GetYaxis()->SetTitle("#font[2]{y} #[]{#AA}");
    frame[0]->GetXaxis()->CenterTitle();
    frame[0]->GetYaxis()->CenterTitle();
    frame[0]->SetStats(0);
    frame[0]->Draw();
    ROOTHelper::SetMarkerAttribute(graph_position.get(), 20, 0.2f, kAzure-7);
    graph_position->Draw("P X0");

    auto circle_inner{ new TEllipse(0.0, 0.0, 40.0) };
    ROOTHelper::SetFillAttribute(circle_inner, 0);
    ROOTHelper::SetLineAttribute(circle_inner, 1, 1, kRed);
    circle_inner->Draw();
    auto circle_outer{ new TEllipse(0.0, 0.0, 60.0) };
    ROOTHelper::SetFillAttribute(circle_outer, 0);
    ROOTHelper::SetLineAttribute(circle_outer, 1, 1, kRed);
    circle_outer->Draw();

    auto pos_text{ ROOTHelper::CreatePaveText(0.08, 1.0, 1.0, 1.2, "nbNDC", true) };
    ROOTHelper::SetPaveTextDefaultStyle(pos_text.get());
    ROOTHelper::SetFillAttribute(pos_text.get(), 4000);
    ROOTHelper::SetTextAttribute(pos_text.get(), 35, 133, 32, 0.0, kAzure-7);
    auto slicing_z_center{ z_position - com_pos.at(2) };
    pos_text->AddText(
        Form("#font[2]{z} #in [%.1f, %.1f] #AA",
        slicing_z_center - z_window*0.5, slicing_z_center + z_window*0.5));
    pos_text->Draw();

    std::unique_ptr<TPaveText> title_text[pad_size - 1];
    for (int i = 1; i < pad_size; i++)
    {
        auto element_id{ i - 1 };
        ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, 0);
        ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
        auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
        frame[i] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, 0),"", 500, x_min[element_id], x_max[element_id], 500, y_min, y_max);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetXaxis(), 45.0f, 0.005f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 505);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetYaxis(), 55.0f, 1.0f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetYaxis(), 45.0f, 0.01f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
        ROOTHelper::SetLineAttribute(frame[i].get(), 1, 0);
        frame[i]->GetXaxis()->SetTitle("");
        frame[i]->GetYaxis()->SetTitle("Width #tau");
        frame[i]->GetXaxis()->CenterTitle();
        frame[i]->GetYaxis()->CenterTitle();
        frame[i]->SetStats(0);
        frame[i]->GetXaxis()->SetLimits(35.0, 65.0);
        frame[i]->Draw("Y+");

        auto element_color{ AtomClassifier::GetMainChainElementColor(static_cast<size_t>(element_id)) };
        auto element_label{ AtomClassifier::GetMainChainElementLabel(static_cast<size_t>(element_id)) };
        if (draw_box_plot == true)
        {
            ROOTHelper::SetLineAttribute(summary_hist[element_id].get(), 1, 1, element_color);
            ROOTHelper::SetFillAttribute(summary_hist[element_id].get(), 1001, element_color, 0.3f);
            summary_hist[element_id]->Draw("CANDLE3 SAME");
        }
        else
        {
            for (auto & graph : graph_list[element_id])
            {
                ROOTHelper::SetMarkerAttribute(graph.get(), 24, 0.1f, element_color);
                graph->Draw("P");
            }
        }
        title_text[element_id] = ROOTHelper::CreatePaveText(0.01, 0.95, 0.39, 1.20, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[element_id].get());
        ROOTHelper::SetPaveAttribute(title_text[element_id].get(), 0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[element_id].get(), 60.0f, 103, 22);
        ROOTHelper::SetFillAttribute(title_text[element_id].get(), 1001, element_color, 0.5f);
        title_text[element_id]->AddText(element_label.data());
        title_text[element_id]->Draw();
    }

    canvas->cd();
    auto pad_extra{ ROOTHelper::CreatePad("pad_extra","", 0.25, 0.00, 0.92, 0.15) };
    pad_extra->Draw();
    pad_extra->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra.get());
    ROOTHelper::SetFillAttribute(pad_extra.get(), 4000);
    auto bottom_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(bottom_title_text.get());
    ROOTHelper::SetFillAttribute(bottom_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(bottom_title_text.get(), 50.0f, 133, 22);
    bottom_title_text->AddText("Distance to Center of Mass #font[2]{d} #[]{#AA}");
    bottom_title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

} // namespace rhbm_gem
#include <rhbm_gem/core/painter/DemoPainter.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/LocalPotentialView.hpp>
#include <rhbm_gem/data/object/ModelPotentialView.hpp>
#include "PotentialPlotBuilder.hpp"
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
            auto entry_iter{ std::make_unique<ModelPotentialView>(*model_list.at(model_id)) };
            auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_list.at(model_id)) };
            for (auto residue : ChemicalDataHelper::GetStandardAminoAcidList())
            {
                auto group_key{ m_atom_classifier->GetMainChainComponentAtomClassGroupKey(element_id, residue) };
                if (!entry_iter->HasAtomGroup(group_key, residue_class)) continue;
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
    auto entry_iter1{ std::make_unique<ModelPotentialView>(*model_object1) };
    auto plot_builder1{ std::make_unique<PotentialPlotBuilder>(model_object1) };
    auto entry_iter2{ std::make_unique<ModelPotentialView>(*model_object2) };
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
    auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
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
        auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
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
        auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
        auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
        for (size_t p = 0; p < spot_list.size(); p++)
        {
            auto spot{ spot_list.at(p) };
            group_key_list_map[i].emplace(
                spot, m_atom_classifier->GetMainChainComponentAtomClassGroupKeyList(p));
            auto & group_key_list{ group_key_list_map[i].at(spot) };
            for (auto it = group_key_list.begin(); it != group_key_list.end(); )
            {
                if (!entry_iter->HasAtomGroup(*it, class_key))
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
        auto entry_iter{ std::make_unique<ModelPotentialView>(*model_object) };
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

#endif


} // namespace rhbm_gem
