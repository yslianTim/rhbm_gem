#include <rhbm_gem/core/painter/ComparisonPainter.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/LocalPotentialEntry.hpp>
#include <rhbm_gem/data/object/ModelPotentialView.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include "detail/PainterSupport.hpp"

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TColor.h>
#include <TMarker.h>
#include <TAxis.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLine.h>
#endif

#include <tuple>

namespace rhbm_gem {

ComparisonPainter::ComparisonPainter() :
    m_atom_classifier{ std::make_unique<AtomClassifier>() }
{

}

ComparisonPainter::~ComparisonPainter()
{

}

void ComparisonPainter::AddModel(ModelObject & data_object)
{
    m_model_object_list.emplace_back(&data_object);
    m_resolution_list.emplace_back(data_object.GetResolution());
}

void ComparisonPainter::AddReferenceModel(ModelObject & data_object, std::string_view label)
{
    m_ref_model_object_list_map[std::string(label)].push_back(&data_object);
}

void ComparisonPainter::Painting()
{
    Logger::Log(LogLevel::Info, "ComparisonPainter::Painting() called.");
    Logger::Log(LogLevel::Info, "Folder path: " + m_folder_path);
    Logger::Log(LogLevel::Info, "Number of atom objects to be painted: "
                + std::to_string(m_model_object_list.size()));

    if (m_ref_model_object_list_map.find("with_charge") != m_ref_model_object_list_map.end() &&
        m_ref_model_object_list_map.find("no_charge") != m_ref_model_object_list_map.end())
    {
        PaintGroupGausEstimateComparison("figure_4_ab.pdf");
        PaintGausEstimateResidueClassDenseComparison("figure_4_sup.pdf");
    }

    for (auto & model_object : m_model_object_list)
    {
        auto plot_name{ "figure_method_charge_"+ model_object->GetPdbID() +".pdf" };
        auto & ref_model_object_list{ m_ref_model_object_list_map.at("with_charge") };
        if (model_object->GetPdbID() == ref_model_object_list.at(0)->GetPdbID())
        {
            PainMapValueComparison(plot_name, model_object, ref_model_object_list);
        }
    }

}

void ComparisonPainter::PaintGroupGausEstimateComparison(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, "- ComparisonPainter::PaintGroupGausEstimateComparison");

    auto sim_with_charge_model_object_list{ m_ref_model_object_list_map.at("with_charge")};
    auto sim_no_charge_model_object_list{ m_ref_model_object_list_map.at("no_charge")};
    auto sim_amber95_model_object_list{ m_ref_model_object_list_map.at("amber95")};

    #ifdef HAVE_ROOT
    const char * data_index[10]{"A","B","C","D","E","F","G","H","I","J"};
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int extra_pad_size{ 5 };
    const int col_size{ 3 };
    const int row_size{ 4 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 2000, 1600) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.20f, 0.30f, 0.10f, 0.06f, 0.02f, 0.02f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TPad> pad_extra[extra_pad_size];
    std::unique_ptr<TH2> frame_extra[extra_pad_size];
    pad_extra[0] = ROOTHelper::CreatePad("pad_0","", 0.70, 0.00, 1.00, 0.35); // Right bottom
    pad_extra[1] = ROOTHelper::CreatePad("pad_1","", 0.70, 0.35, 1.00, 0.70); // Right middle
    pad_extra[2] = ROOTHelper::CreatePad("pad_2","", 0.70, 0.70, 1.00, 1.00); // Right top
    pad_extra[3] = ROOTHelper::CreatePad("pad_3","", 0.20, 0.00, 0.70, 0.06); // Left bottom
    pad_extra[4] = ROOTHelper::CreatePad("pad_4","", 0.13, 0.12, 0.17, 0.94); // Left middle

    size_t ref_id[col_size]{ 3, 3, 3 };
    Residue residue_id[row_size]{ Residue::TRP, Residue::ILE, Residue::GLU, Residue::UNK };

    std::unique_ptr<TGraphErrors> data_graph[col_size][row_size];
    std::unique_ptr<TGraphErrors> sim_with_charge_graph[col_size][row_size];
    std::unique_ptr<TGraphErrors> sim_no_charge_graph[col_size][row_size];
    std::unique_ptr<TGraphErrors> sim_additional_graph[col_size][row_size];
    double x_min[col_size]{ 0.0 };
    double x_max[col_size]{ 1.0 };
    double y_min[row_size]{ 0.0 };
    double y_max[row_size]{ 1.0 };
    std::vector<double> x_array[col_size], y_array[row_size];
    for (size_t i = 0; i < col_size; i++)
    {
        for (size_t j = 0; j < row_size; j++)
        {
            data_graph[i][j] = ROOTHelper::CreateGraphErrors();
            sim_with_charge_graph[i][j] = ROOTHelper::CreateGraphErrors();
            sim_no_charge_graph[i][j] = ROOTHelper::CreateGraphErrors();
            sim_additional_graph[i][j] = ROOTHelper::CreateGraphErrors();

            auto class_key{ (j == row_size - 1) ? ChemicalDataHelper::GetSimpleAtomClassKey() : ChemicalDataHelper::GetComponentAtomClassKey() };
            BuildAmplitudeRatioToWidthGraph(i, ref_id[i], data_graph[i][j].get(), m_model_object_list, class_key, true, residue_id[j]);
            BuildAmplitudeRatioToWidthGraph(i, ref_id[i], sim_with_charge_graph[i][j].get(), sim_with_charge_model_object_list, class_key, false, residue_id[j]);
            BuildAmplitudeRatioToWidthGraph(i, ref_id[i], sim_no_charge_graph[i][j].get(), sim_no_charge_model_object_list, class_key, false, residue_id[j]);
            BuildAmplitudeRatioToWidthGraph(i, ref_id[i], sim_additional_graph[i][j].get(), sim_amber95_model_object_list, class_key, false, residue_id[j]);

            for (int p = 0; p < data_graph[i][j]->GetN(); p++)
            {
                x_array[i].emplace_back(data_graph[i][j]->GetPointX(p));
                y_array[j].emplace_back(data_graph[i][j]->GetPointY(p));
            }
        }
    }

    for (size_t i = 0; i < col_size; i++)
    {
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array[i], 0.16) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);
    }
    for (size_t j = 0; j < row_size; j++)
    {
        auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array[j], 0.24) };
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }

    // Simulation plots
    size_t target_id{ 0 };
    size_t reference_id{ 1 };
    const int par_id[2]{ 1, 0 };
    std::unique_ptr<TGraphErrors> simulation1_graph[2];
    std::unique_ptr<TGraphErrors> simulation2_graph[2];
    std::unique_ptr<TGraphErrors> simulation3_graph[2];
    std::vector<double> x_array_sim, y_array_sim[2];
    for (int i = 0; i < 2; i++)
    {
        simulation1_graph[i] = ROOTHelper::CreateGraphErrors();
        simulation2_graph[i] = ROOTHelper::CreateGraphErrors();
        simulation3_graph[i] = ROOTHelper::CreateGraphErrors();
        BuildGausRatioToResolutionGraph(par_id[i], target_id, reference_id, simulation1_graph[i].get(), sim_with_charge_model_object_list, ChemicalDataHelper::GetSimpleAtomClassKey());
        BuildGausRatioToResolutionGraph(par_id[i], target_id, reference_id, simulation2_graph[i].get(), sim_amber95_model_object_list, ChemicalDataHelper::GetSimpleAtomClassKey());
        BuildGausRatioToResolutionGraph(par_id[i], target_id, reference_id, simulation3_graph[i].get(), sim_no_charge_model_object_list, ChemicalDataHelper::GetSimpleAtomClassKey());
        for (int p = 0; p < simulation1_graph[i]->GetN(); p++)
        {
            x_array_sim.push_back(simulation1_graph[i]->GetPointX(p));
            y_array_sim[i].push_back(simulation1_graph[i]->GetPointY(p));
        }
        for (int p = 0; p < simulation2_graph[i]->GetN(); p++)
        {
            x_array_sim.push_back(simulation2_graph[i]->GetPointX(p));
            y_array_sim[i].push_back(simulation2_graph[i]->GetPointY(p));
        }
        for (int p = 0; p < simulation3_graph[i]->GetN(); p++)
        {
            x_array_sim.push_back(simulation3_graph[i]->GetPointX(p));
            y_array_sim[i].push_back(simulation3_graph[i]->GetPointY(p));
        }
    }

    auto x_range_sim{ ArrayStats<double>::ComputeScalingRangeTuple(x_array_sim, 0.1) };
    auto x_min_sim{ std::get<0>(x_range_sim) };
    auto x_max_sim{ std::get<1>(x_range_sim) };
    double y_min_sim[2]{ 0.0 };
    double y_max_sim[2]{ 0.0 };
    
    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> info_text[row_size];
    std::unique_ptr<TPaveText> info_x_text[col_size];
    for (size_t i = 0; i < col_size; i++)
    {
        for (size_t j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), static_cast<int>(i), static_cast<int>(j));
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 0, 0, 4000, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", static_cast<int>(i), static_cast<int>(j)),"", 500, x_min[i], x_max[i], 500, y_min[j], y_max[j]);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 50.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.08/x_factor), 504);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 50.0f, 0.02f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            auto element_label_x{ AtomClassifier::GetMainChainElementLabel(i) };
            auto ref_element_label_x{ AtomClassifier::GetMainChainElementLabel(ref_id[i]) };
            auto element_color_x{  AtomClassifier::GetMainChainElementColor(i)};
            auto element_label_y{ AtomClassifier::GetMainChainElementLabel(j) };
            auto ref_element_label_y{ AtomClassifier::GetMainChainElementLabel(3) };
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("");
            frame[i][j]->Draw();

            auto element_marker{ AtomClassifier::GetMainChainElementSolidMarker(i) };
            ROOTHelper::SetMarkerAttribute(data_graph[i][j].get(), element_marker, 2.0f, element_color_x);
            ROOTHelper::SetLineAttribute(sim_with_charge_graph[i][j].get(), 1, 3, kRed);
            ROOTHelper::SetLineAttribute(sim_no_charge_graph[i][j].get(),   2, 3, kGray+2);
            ROOTHelper::SetLineAttribute(sim_additional_graph[i][j].get(),  3, 3, kBlue);
            sim_with_charge_graph[i][j]->Draw("C X0");
            sim_no_charge_graph[i][j]->Draw("C X0");
            sim_additional_graph[i][j]->Draw("C X0");
            data_graph[i][j]->Draw("P X0");

            if (i == 0)
            {
                info_text[j] = ROOTHelper::CreatePaveText(-1.28, 0.25, -0.50, 0.75, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(info_text[j].get());
                ROOTHelper::SetPaveAttribute(info_text[j].get(), 0, 0.1);
                ROOTHelper::SetFillAttribute(info_text[j].get(), 1001, kAzure-7);
                ROOTHelper::SetTextAttribute(info_text[j].get(), 70.0f, 133, 22, 0.0, kYellow-10);
                if (j == row_size - 1) info_text[j]->AddText("#splitline{Whole}{Protein}");
                else info_text[j]->AddText(Form("#font[102]{%s}", ChemicalDataHelper::GetLabel(residue_id[j]).data()));
                info_text[j]->Draw();
            }

            if (j == row_size - 1)
            {
                info_x_text[i] = ROOTHelper::CreatePaveText( 0.10, 1.01, 0.90, 1.27, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(info_x_text[i].get());
                ROOTHelper::SetPaveAttribute(info_x_text[i].get(), 0, 0.2);
                ROOTHelper::SetFillAttribute(info_x_text[i].get(), 1001, element_color_x, 0.50f);
                ROOTHelper::SetTextAttribute(info_x_text[i].get(), 70.0f, 133, 22);
                ROOTHelper::SetLineAttribute(info_x_text[i].get(), 1, 0);
                info_x_text[i]->AddText(Form("#font[2]{t_{i}} = #font[102]{%s}", element_label_x.data()));
                info_x_text[i]->Draw();
            }
        }
    }

   
    canvas->cd();
    for (int i = 0; i < extra_pad_size; i++)
    {
        ROOTHelper::SetPadDefaultStyle(pad_extra[i].get());
        pad_extra[i]->Draw();
    }

    std::unique_ptr<TH2> frame_sim[2];
    for (int i = 0; i < 2; i++)
    {
        auto y_range_sim{ ArrayStats<double>::ComputeScalingRangeTuple(y_array_sim[i], 0.2) };
        y_min_sim[i] = std::get<0>(y_range_sim);
        y_max_sim[i] = std::get<1>(y_range_sim);
    }

    const std::string y_axis_title[2]
    {
        Form("#tau_{#color[%d]{#font[102]{%s}}} / #tau_{#color[%d]{#font[102]{%s}}}",
            AtomClassifier::GetMainChainElementColor(0),
            AtomClassifier::GetMainChainElementLabel(0).data(),
            AtomClassifier::GetMainChainElementColor(1),
            AtomClassifier::GetMainChainElementLabel(1).data()),
        Form("#font[2]{A}_{#color[%d]{#font[102]{%s}}} / #font[2]{A}_{#color[%d]{#font[102]{%s}}}",
            AtomClassifier::GetMainChainElementColor(0),
            AtomClassifier::GetMainChainElementLabel(0).data(),
            AtomClassifier::GetMainChainElementColor(1),
            AtomClassifier::GetMainChainElementLabel(1).data())
    };

    pad_extra[0]->cd(); // right bottom
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.08, 0.01, 0.10, 0.01);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 0, 0, 4000, 0);
    frame_sim[0] = ROOTHelper::CreateHist2D("hist_sim_0","", 500, x_min_sim, x_max_sim, 500, y_min_sim[0], y_max_sim[0]);
    ROOTHelper::SetAxisTitleAttribute(frame_sim[0]->GetXaxis(), 70.0f, 1.0f, 133);
    ROOTHelper::SetAxisLabelAttribute(frame_sim[0]->GetXaxis(), 50.0f, 0.005f, 133);
    ROOTHelper::SetAxisTickAttribute(frame_sim[0]->GetXaxis(), 0.05f, 506);
    ROOTHelper::SetAxisTitleAttribute(frame_sim[0]->GetYaxis(), 75.0f, 0.9f, 133);
    ROOTHelper::SetAxisLabelAttribute(frame_sim[0]->GetYaxis(), 50.0f, 0.02f, 133);
    ROOTHelper::SetAxisTickAttribute(frame_sim[0]->GetYaxis(), 0.05f, 506);
    ROOTHelper::SetLineAttribute(frame_sim[0].get(), 1, 0);
    frame_sim[0]->GetXaxis()->SetTitle("Blurring Width #sigma_{B}");
    frame_sim[0]->GetYaxis()->SetTitle(y_axis_title[0].data());
    frame_sim[0]->GetYaxis()->CenterTitle();
    frame_sim[0]->SetStats(0);
    frame_sim[0]->Draw();

    pad_extra[1]->cd(); // right middle
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.08, 0.01, 0.01, 0.10);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 0, 0, 4000, 0);
    frame_sim[1] = ROOTHelper::CreateHist2D("hist_sim_1","", 500, x_min_sim, x_max_sim, 500, y_min_sim[1], y_max_sim[1]);
    ROOTHelper::SetAxisTitleAttribute(frame_sim[1]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(frame_sim[1]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTickAttribute(frame_sim[1]->GetXaxis(), 0.05f, 506);
    ROOTHelper::SetAxisTitleAttribute(frame_sim[1]->GetYaxis(), 75.0f, 0.9f, 133);
    ROOTHelper::SetAxisLabelAttribute(frame_sim[1]->GetYaxis(), 50.0f, 0.02f, 133);
    ROOTHelper::SetAxisTickAttribute(frame_sim[1]->GetYaxis(), 0.05f, 506);
    ROOTHelper::SetLineAttribute(frame_sim[1].get(), 1, 0);
    frame_sim[1]->GetXaxis()->SetTitle("");
    frame_sim[1]->GetYaxis()->SetTitle(y_axis_title[1].data());
    frame_sim[1]->GetYaxis()->CenterTitle();
    frame_sim[1]->SetStats(0);
    frame_sim[1]->Draw();
    auto title_sim_text{ ROOTHelper::CreatePaveText( 0.27, 0.73, 0.95, 0.90, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(title_sim_text.get());
    ROOTHelper::SetPaveAttribute(title_sim_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(title_sim_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(title_sim_text.get(), 70.0f, 133, 22, 0.0, kYellow-10);
    ROOTHelper::SetLineAttribute(title_sim_text.get(), 1, 0);
    title_sim_text->AddText("Simulation");
    title_sim_text->Draw();
    
    for (int i = 0; i < 2; i++)
    {
        pad_extra[i]->cd();
        ROOTHelper::SetMarkerAttribute(simulation1_graph[i].get(), 20, 1.5f, kRed);
        ROOTHelper::SetLineAttribute(simulation1_graph[i].get(), 1, 3, kRed);
        ROOTHelper::SetMarkerAttribute(simulation2_graph[i].get(), 20, 1.5f, kBlue);
        ROOTHelper::SetLineAttribute(simulation2_graph[i].get(), 2, 3, kBlue);
        ROOTHelper::SetMarkerAttribute(simulation3_graph[i].get(), 20, 1.5f, kGray+2);
        ROOTHelper::SetLineAttribute(simulation3_graph[i].get(), 3, 3, kGray+2);
        simulation3_graph[i]->Draw("PL X0");
        simulation2_graph[i]->Draw("PL X0");
        simulation1_graph[i]->Draw("PL X0");
    }

    pad_extra[2]->cd(); // right top legend
    std::unique_ptr<TPaveText> data_text[2];
    double x_pos[4]{ 0.02, 0.50, 1.00 };
    for (int i = 0; i < 2; i++)
    {
        data_text[i] = ROOTHelper::CreatePaveText(x_pos[i], 0.50, x_pos[i+1], 0.97, "nbNDC", false);
        ROOTHelper::SetPaveTextDefaultStyle(data_text[i].get());
        ROOTHelper::SetFillAttribute(data_text[i].get(), 4000);
        ROOTHelper::SetTextAttribute(data_text[i].get(), 60.0f, 133, 12);
    }
    auto count{ 0 };
    for (auto & resolution : m_resolution_list)
    {
        auto id{ (count <= 3) ? 0 : 1 };
        data_text[id]->AddText(Form("#font[102]{[#color[853]{%s}]} %.2f #AA", data_index[count], resolution));
        count++;
        if (count >= 10) break;
    }
    data_text[1]->AddText("");
    data_text[0]->Draw();
    data_text[1]->Draw();

    auto legend{ ROOTHelper::CreateLegend(0.02, 0.00, 1.00, 0.45, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 55.0f, 133, 12);
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    legend->AddEntry(sim_with_charge_graph[0][0].get(), "Thomas #font[2]{et al.}", "l");
    legend->AddEntry(sim_additional_graph[0][0].get(), "Cornell #font[2]{et al.}", "l");
    legend->AddEntry(sim_no_charge_graph[0][0].get(), "All Neutral #alpha = 0", "l");
    legend->Draw();

    pad_extra[3]->cd(); // bottom X title
    auto bottom_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(bottom_title_text.get());
    ROOTHelper::SetFillAttribute(bottom_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(bottom_title_text.get(), 80.0f, 133, 22);
    bottom_title_text->AddText("Width #tau_{#font[2]{t_{i}}}");
    bottom_title_text->Draw();

    pad_extra[4]->cd(); // left Y title
    auto left_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(left_title_text.get());
    ROOTHelper::SetFillAttribute(left_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(left_title_text.get(), 80.0f, 133, 22);
    left_title_text->AddText(Form("Amplitude Ratio #font[1]{A}_{#font[2]{t_{i}}} / #font[2]{A}_{#color[%d]{#font[102]{%s}}}",
                    AtomClassifier::GetMainChainElementColor(3), AtomClassifier::GetMainChainElementLabel(3).data()));
    left_title_text->Draw();
    auto text{ left_title_text->GetLineWith("Amplitude") };
    text->SetTextAngle(90.0f);
    left_title_text->Draw();
    

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ComparisonPainter::PaintGausEstimateResidueClassDenseComparison(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ComparisonPainter::PaintGausEstimateResidueClassDenseComparison");

    auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    
    auto sim_with_charge_model_object_list{ m_ref_model_object_list_map.at("with_charge")};
    auto sim_no_charge_model_object_list{ m_ref_model_object_list_map.at("no_charge")};
    auto sim_amber95_model_object_list{ m_ref_model_object_list_map.at("amber95")};

    #ifdef HAVE_ROOT
    const char * data_index[10]{"A","B","C","D","E","F","G","H","I","J"};
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 5 };
    const int row_size{ 4 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 2000, 1800) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.09f, 0.00f, 0.09f, 0.15f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int primary_element_size{ 3 };
    const int standard_residue_size{ 20 };
    size_t ref_id[primary_element_size]{ 3, 3, 3 };

    std::unordered_map<Residue, std::unique_ptr<TGraphErrors>> data_graph_map[primary_element_size];
    std::unordered_map<Residue, std::unique_ptr<TGraphErrors>> sim_with_charge_graph_map[primary_element_size];
    std::unordered_map<Residue, std::unique_ptr<TGraphErrors>> sim_no_charge_graph_map[primary_element_size];
    std::unordered_map<Residue, std::unique_ptr<TGraphErrors>> sim_additional_graph_map[primary_element_size];
    for (size_t k = 0; k < primary_element_size; k++)
    {
        auto element_marker{ AtomClassifier::GetMainChainElementSolidMarker(k) };
        auto element_color{ AtomClassifier::GetMainChainElementColor(k) };
        auto element_label{ AtomClassifier::GetMainChainElementLabel(k) };
        auto ref_element_color{ AtomClassifier::GetMainChainElementColor(ref_id[k]) };
        auto ref_element_label{ AtomClassifier::GetMainChainElementLabel(ref_id[k]) };

        std::vector<double> x_array, y_array;
        for (int i = 0; i < standard_residue_size; i++)
        {
            auto residue{ static_cast<Residue>(i+1) };
            auto data_graph{ ROOTHelper::CreateGraphErrors() };
            auto sim_with_charge_graph{ ROOTHelper::CreateGraphErrors() };
            auto sim_no_charge_graph{ ROOTHelper::CreateGraphErrors() };
            auto sim_additional_graph{ ROOTHelper::CreateGraphErrors() };
            BuildAmplitudeRatioToWidthGraph(k, ref_id[k], data_graph.get(), m_model_object_list, class_key, true, residue);
            BuildAmplitudeRatioToWidthGraph(k, ref_id[k], sim_with_charge_graph.get(), sim_with_charge_model_object_list, class_key, false, residue);
            BuildAmplitudeRatioToWidthGraph(k, ref_id[k], sim_no_charge_graph.get(), sim_no_charge_model_object_list, class_key, false, residue);
            BuildAmplitudeRatioToWidthGraph(k, ref_id[k], sim_additional_graph.get(), sim_amber95_model_object_list, class_key, false, residue);

            for (int p = 0; p < data_graph->GetN(); p++)
            {
                x_array.emplace_back(data_graph->GetPointX(p));
                y_array.emplace_back(data_graph->GetPointY(p));
            }
            data_graph_map[k][residue] = std::move(data_graph);
            sim_with_charge_graph_map[k][residue] = std::move(sim_with_charge_graph);
            sim_no_charge_graph_map[k][residue] = std::move(sim_no_charge_graph);
            sim_additional_graph_map[k][residue] = std::move(sim_additional_graph);
        }

        auto x_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(x_array, 0.30) };
        auto x_min{ std::get<0>(x_range) };
        auto x_max{ std::get<1>(x_range) };
        x_min = (x_min <= 0.4) ? 0.45 : x_min;

        auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.45) };
        auto y_min{ std::get<0>(y_range) };
        auto y_max{ std::get<1>(y_range) };
        y_min = (y_min <= 0.0) ? 0.1 : y_min;

        std::unique_ptr<TH2> frame[standard_residue_size];
        std::unique_ptr<TPaveText> title_text[standard_residue_size];
        for (int i = 0; i < standard_residue_size; i++)
        {
            auto residue{ static_cast<Residue>(i+1) };
            auto col_id{ i % col_size };
            auto row_id{ row_size - (i / col_size) - 1 };
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), col_id, row_id);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 0, 0, 4000, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i] = ROOTHelper::CreateHist2D(Form("hist_%d", i),"", 500, x_min, x_max, 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i]->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(frame[i]->GetXaxis(), 50.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 504);
            ROOTHelper::SetAxisTitleAttribute(frame[i]->GetYaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(frame[i]->GetYaxis(), 50.0f, 0.02f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
            ROOTHelper::SetLineAttribute(frame[i].get(), 1, 0);
            frame[i]->GetXaxis()->SetTitle("");
            frame[i]->GetYaxis()->SetTitle("");
            frame[i]->SetStats(0);
            frame[i]->Draw();

            ROOTHelper::SetMarkerAttribute(data_graph_map[k].at(residue).get(), element_marker, 2.0f, element_color);
            ROOTHelper::SetLineAttribute(sim_with_charge_graph_map[k].at(residue).get(), 1, 2, kRed);
            ROOTHelper::SetLineAttribute(sim_no_charge_graph_map[k].at(residue).get(), 2, 2, kGray+2);
            ROOTHelper::SetLineAttribute(sim_additional_graph_map[k].at(residue).get(), 3, 2, kBlue);
            sim_with_charge_graph_map[k].at(residue)->Draw("C X0");
            sim_no_charge_graph_map[k].at(residue)->Draw("C X0");
            sim_additional_graph_map[k].at(residue)->Draw("C X0");
            data_graph_map[k].at(residue)->Draw("P X0");

            auto x1_pos{ (k == 3) ? 0.60 : 0.10 };
            auto x2_pos{ (k == 3) ? 0.95 : 0.45 };
            title_text[i] = ROOTHelper::CreatePaveText(x1_pos, 0.80, x2_pos, 1.00, "nbNDC ARC", true);
            ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
            ROOTHelper::SetPaveAttribute(title_text[i].get(), 0, 0.2);
            ROOTHelper::SetTextAttribute(title_text[i].get(), 60.0f, 103, 22, 0.0, kYellow-10);
            ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, kAzure-7);
            title_text[i]->AddText(ChemicalDataHelper::GetLabel(residue).data());
            title_text[i]->Draw();
        }
        canvas->cd();
        auto pad_extra_t{ ROOTHelper::CreatePad("pad_extra_t","", 0.09, 0.86, 1.00, 1.00) };
        pad_extra_t->Draw();
        pad_extra_t->cd();
        ROOTHelper::SetPadDefaultStyle(pad_extra_t.get());
        ROOTHelper::SetFillAttribute(pad_extra_t.get(), 4000);

        auto legend{ ROOTHelper::CreateLegend(0.00, 0.00, 0.40, 1.00, false) };
        ROOTHelper::SetLegendDefaultStyle(legend.get());
        ROOTHelper::SetTextAttribute(legend.get(), 55.0f, 133, 12);
        ROOTHelper::SetFillAttribute(legend.get(), 4000);
        legend->AddEntry(sim_with_charge_graph_map[k].at(Residue::ALA).get(), "Thomas #font[2]{et al.}", "l");
        legend->AddEntry(sim_additional_graph_map[k].at(Residue::ALA).get(), "Cornell #font[2]{et al.}", "l");
        legend->AddEntry(sim_no_charge_graph_map[k].at(Residue::ALA).get(), "All Neutral #alpha = 0", "l");
        legend->Draw();

        std::unique_ptr<TPaveText> info_text[3];
        double x_pos[4]{ 0.40, 0.60, 0.80, 1.00 };
        for (int i = 0; i < 3; i++)
        {
            info_text[i] = ROOTHelper::CreatePaveText(x_pos[i], 0.00, x_pos[i+1], 1.00, "nbNDC", false);
            ROOTHelper::SetPaveTextDefaultStyle(info_text[i].get());
            ROOTHelper::SetFillAttribute(info_text[i].get(), 4000);
            ROOTHelper::SetTextAttribute(info_text[i].get(), 55.0f, 133, 12);
        }
        auto count{ 0 };
        for (auto & resolution : m_resolution_list)
        {
            auto id{ count / 3 };
            info_text[id]->AddText(Form("#font[102]{[#color[853]{%s}]} %.2f #AA", data_index[count], resolution));
            count++;
            if (count >= 10) break;
        }
        info_text[2]->AddText("");
        info_text[2]->AddText("");
        info_text[0]->Draw();
        info_text[1]->Draw();
        info_text[2]->Draw();

        canvas->cd();
        auto pad_extra_x{ ROOTHelper::CreatePad("pad_extra_x","", 0.10, 0.00, 1.00, 0.05) };
        pad_extra_x->Draw();
        pad_extra_x->cd();
        ROOTHelper::SetPadDefaultStyle(pad_extra_x.get());
        ROOTHelper::SetFillAttribute(pad_extra_x.get(), 4000);
        auto x_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
        ROOTHelper::SetPaveTextDefaultStyle(x_title_text.get());
        ROOTHelper::SetFillAttribute(x_title_text.get(), 4000);
        ROOTHelper::SetTextAttribute(x_title_text.get(), 70.0f, 133, 22);
        x_title_text->AddText(Form("Width #tau_{#font[102]{#color[%d]{%s}}}", element_color, element_label.data()));
        x_title_text->Draw();

        canvas->cd();
        auto pad_extra_y{ ROOTHelper::CreatePad("pad_extra_y","", 0.00, 0.10, 0.05, 0.85) };
        pad_extra_y->Draw();
        pad_extra_y->cd();
        ROOTHelper::SetPadDefaultStyle(pad_extra_y.get());
        ROOTHelper::SetFillAttribute(pad_extra_y.get(), 4000);
        auto y_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
        ROOTHelper::SetPaveTextDefaultStyle(y_title_text.get());
        ROOTHelper::SetFillAttribute(y_title_text.get(), 4000);
        ROOTHelper::SetTextAttribute(y_title_text.get(), 80.0f, 133, 22);
        y_title_text->AddText(
            Form("Amplitude Ratio #font[1]{A}_{#font[102]{#color[%d]{%s}}} / #font[1]{A}_{#font[102]{#color[%d]{%s}}}",
                element_color, element_label.data(), ref_element_color, ref_element_label.data()));
        auto text{ y_title_text->GetLineWith("Amplitude") };
        text->SetTextAngle(90.0f);
        y_title_text->Draw();
        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ComparisonPainter::PainMapValueComparison(
    const std::string & name,
    ModelObject * model_object,
    const std::vector<ModelObject *> & ref_model_object_list)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ComparisonPainter::PainMapValueComparison");
    (void)model_object;
    (void)ref_model_object_list;

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
    for (auto & ref_model_object : ref_model_object_list)
    {
        if (model_object->GetPdbID() != ref_model_object->GetPdbID()) continue;
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
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

#ifdef HAVE_ROOT
void ComparisonPainter::BuildGausRatioToResolutionGraph(
    int par_id, size_t target_id, size_t reference_id, TGraphErrors * graph,
    const std::vector<ModelObject *> & model_list, const std::string & class_key, Residue residue)
{
    auto group_key
    {
        (class_key == ChemicalDataHelper::GetSimpleAtomClassKey()) ?
        m_atom_classifier->GetMainChainSimpleAtomClassGroupKey(target_id) :
        m_atom_classifier->GetMainChainComponentAtomClassGroupKey(target_id, residue)
    };
    auto ref_group_key
    {
        (class_key == ChemicalDataHelper::GetSimpleAtomClassKey()) ?
        m_atom_classifier->GetMainChainSimpleAtomClassGroupKey(reference_id) :
        m_atom_classifier->GetMainChainComponentAtomClassGroupKey(reference_id, residue)
    };
    auto count{ 0 };
    for (auto model_object : model_list)
    {
        const ModelPotentialView entry_view{ *model_object };
        if (!entry_view.HasAtomGroup(group_key, class_key)) continue;
        if (!entry_view.HasAtomGroup(ref_group_key, class_key)) continue;
        auto x_value{ model_object->GetResolution() };
        auto y_value{ entry_view.GetAtomGausEstimatePrior(group_key, class_key, par_id) };
        auto y_error{ entry_view.GetAtomGausVariancePrior(group_key, class_key, par_id) };
        auto ref_y_value{ entry_view.GetAtomGausEstimatePrior(ref_group_key, class_key, par_id) };
        auto ref_y_error{ entry_view.GetAtomGausVariancePrior(ref_group_key, class_key, par_id) };
        if (x_value == 0.0 || ref_y_value == 0.0) continue;
        auto ratio{ y_value/ref_y_value };
        auto error{ ratio * std::sqrt(std::pow(y_error/y_value, 2) + std::pow(ref_y_error/ref_y_value, 2)) };
        graph->SetPoint(count, x_value, ratio);
        graph->SetPointError(count, 0.0, error);
        count++;
    }
}

void ComparisonPainter::BuildAmplitudeRatioToWidthGraph(
    size_t target_id, size_t reference_id, TGraphErrors * graph,
    const std::vector<ModelObject *> & model_list, const std::string & class_key,
    bool draw_index, Residue residue)
{
    const char * data_index[11]{"A","B","C","D","E","F","G","H","I","J","K"};
    auto group_key
    {
        (class_key == ChemicalDataHelper::GetSimpleAtomClassKey()) ?
        m_atom_classifier->GetMainChainSimpleAtomClassGroupKey(target_id) :
        m_atom_classifier->GetMainChainComponentAtomClassGroupKey(target_id, residue)
    };
    auto ref_group_key
    {
        (class_key == ChemicalDataHelper::GetSimpleAtomClassKey()) ?
        m_atom_classifier->GetMainChainSimpleAtomClassGroupKey(reference_id) :
        m_atom_classifier->GetMainChainComponentAtomClassGroupKey(reference_id, residue)
    };
    auto count{ 0 };
    auto model_count{ 0 };
    for (auto model_object : model_list)
    {
        model_count++;
        const ModelPotentialView entry_view{ *model_object };
        if (!entry_view.HasAtomGroup(group_key, class_key)) continue;
        if (!entry_view.HasAtomGroup(ref_group_key, class_key)) continue;
        auto x_value{ entry_view.GetAtomGausEstimatePrior(group_key, class_key, 1) };
        auto y_value{ entry_view.GetAtomGausEstimatePrior(group_key, class_key, 0) };
        auto x_error{ std::sqrt(entry_view.GetAtomGausVariancePrior(group_key, class_key, 1)) };
        auto y_error{ std::sqrt(entry_view.GetAtomGausVariancePrior(group_key, class_key, 0)) };
        auto ref_y_value{ entry_view.GetAtomGausEstimatePrior(ref_group_key, class_key, 0) };
        auto ref_y_error{ std::sqrt(entry_view.GetAtomGausVariancePrior(ref_group_key, class_key, 0)) };
        if (x_value == 0.0 || ref_y_value == 0.0) continue;
        auto ratio{ y_value/ref_y_value };
        auto error{ ratio * std::sqrt(std::pow(y_error/y_value, 2) + std::pow(ref_y_error/ref_y_value, 2)) };
        graph->SetPoint(count, x_value, ratio);
        graph->SetPointError(count, x_error, error);
        count++;
        if (draw_index == true)
        {
            std::string label;
            if (model_count > 11)
            {
                Logger::Log(LogLevel::Warning, "Data label size exceeds 12, label switch to numbers.");
                label = std::to_string(model_count);
            }
            else
            {
                label = "#splitline{" + std::string(data_index[model_count-1]) + "}{ }";
            }
            TLatex * latex = new TLatex(x_value, ratio, label.data());
            ROOTHelper::SetTextAttribute(latex, 35.0f, 103, 21, 0.0, kAzure-7);
            graph->GetListOfFunctions()->Add(latex);
        }
    }
}

#endif

} // namespace rhbm_gem
