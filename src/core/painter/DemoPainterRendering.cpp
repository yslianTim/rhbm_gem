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
#include "internal/PainterTypeCheck.hpp"

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
    auto & typed_data_object{
        painter_internal::RequirePainterObject<ModelObject>(
            data_object, "DemoPainter", "AddDataObject") };
    AppendModelObject(typed_data_object);
}

void DemoPainter::AddReferenceDataObject(DataObjectBase * data_object, const std::string & label)
{
    auto & typed_data_object{
        painter_internal::RequirePainterObject<ModelObject>(
            data_object, "DemoPainter", "AddReferenceDataObject") };
    AppendReferenceModelObject(typed_data_object, label);
}

void DemoPainter::AppendModelObject(ModelObject & data_object)
{
    m_model_object_list.push_back(&data_object);
}

void DemoPainter::AppendReferenceModelObject(
    ModelObject & data_object,
    const std::string & label)
{
    m_ref_model_object_list_map[label].push_back(&data_object);
}

void DemoPainter::Painting()
{
    Logger::Log(LogLevel::Info, "DemoPainter::Painting() called.");
    Logger::Log(LogLevel::Info, "Folder path: " + m_folder_path);

    for (auto & model : m_model_object_list)
    {
        model->BuildKDTreeRoot();
    }

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

    PaintGroupGausMergeResidueDemo(demo_alpha_carbon_list, "figure_gaus_backbone_boxplot.pdf");
/*
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
    }*/
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
        BuildMapValueScatterGraph(group_key, graph.get(), ref_model_object, model_object, 15, 0.0, 1.5);
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
    auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
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
    if (entry_iter->IsAvailableAtomGroupKey(group_key, class_key) == false) return;
    for (auto atom : entry_iter->GetAtomObjectList(group_key, class_key))
    {
        auto atom_iter{ std::make_unique<PotentialEntryQuery>(atom) };
        auto atom_plot_builder{ std::make_unique<PotentialPlotBuilder>(atom) };
        auto graph{ atom_plot_builder->CreateBinnedDistanceToMapValueGraph() };
        ROOTHelper::SetLineAttribute(graph.get(), 1, 5, static_cast<short>(kAzure-7), 0.3f);
        map_value_graph_list.emplace_back(std::move(graph));
        auto map_value_range{ atom_iter->GetMapValueRange(0.0) };
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
        auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
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
        if (j == pad_size_y - 1) PrintGausTitlePad(pad[1][j].get(), amplitude_title_text.get(), "Amplitude #font[2]{A}", 80.0f);

        pad[2][j]->cd();
        PrintGausCorrelationPad(pad[2][j].get(), frame[1][j].get(), draw_axis_flag, draw_title_flag);
        for (int k = 0; k < main_chain_element_count; k++) correlation_graph[j][k]->Draw("P X0");
        if (j == pad_size_y - 1) PrintGausTitlePad(pad[2][j].get(), correlation_title_text.get(), "#tau#minus#font[2]{A} Plot", 80.0f);

        pad[3][j]->cd();
        PrintGausResultPad(pad[3][j].get(), frame[2][j].get(), draw_axis_flag, draw_title_flag, true);
        for (int k = 0; k < main_chain_element_count; k++) width_graph[j][k]->Draw("PL X0");
        if (j == pad_size_y - 1) PrintGausTitlePad(pad[3][j].get(), width_title_text.get(), "Width #font[2]{#tau}", 80.0f);
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

    auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
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
    PrintGausResultGlobalPad(pad[1].get(), frame[0].get(), 0.030, 0.005, 0.205, 0.205, false);
    for (int k = 0; k < main_chain_element_count; k++) amplitude_graph[k]->Draw("PL X0");
    PrintGausTitlePad(pad[1].get(), amplitude_title_text.get(), "Amplitude #font[2]{A}", 80.0f);

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
    PrintGausTitlePad(pad[2].get(), correlation_title_text.get(), "#tau#minus#font[2]{A} Plot", 80.0f);

    pad[3]->cd();
    PrintGausResultGlobalPad(pad[3].get(), frame[2].get(), 0.005, 0.030, 0.205, 0.205, true);
    for (int k = 0; k < main_chain_element_count; k++) width_graph[k]->Draw("PL X0");
    PrintGausTitlePad(pad[3].get(), width_title_text.get(), "Width #font[2]{#tau}", 80.0f);

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

} // namespace rhbm_gem
