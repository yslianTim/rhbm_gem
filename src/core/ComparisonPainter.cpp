#include "ComparisonPainter.hpp"
#include "ModelObject.hpp"
#include "AtomObject.hpp"
#include "DataObjectBase.hpp"
#include "PotentialEntryIterator.hpp"
#include "AtomicPotentialEntry.hpp"
#include "FilePathHelper.hpp"
#include "ArrayStats.hpp"
#include "AtomClassifier.hpp"
#include "AtomicInfoHelper.hpp"
#include "GlobalEnumClass.hpp"
#include "DataObjectManager.hpp"

#ifdef HAVE_ROOT
#include "ROOTHelper.hpp"
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

ComparisonPainter::ComparisonPainter(void) :
    m_folder_path{ "./" }, m_atom_classifier{ std::make_unique<AtomClassifier>() }
{

}

ComparisonPainter::~ComparisonPainter()
{

}

void ComparisonPainter::SetFolder(const std::string & folder_path)
{
    m_folder_path = FilePathHelper::EnsureTrailingSlash(folder_path);
}

void ComparisonPainter::AddDataObject(DataObjectBase * data_object)
{
    auto model_object{ dynamic_cast<ModelObject *>(data_object) };
    m_model_object_list.emplace_back(model_object);
    m_resolution_list.emplace_back(model_object->GetResolution());
}

void ComparisonPainter::AddReferenceDataObject(DataObjectBase * data_object, const std::string & label)
{
    m_ref_model_object_list_map[label].emplace_back(dynamic_cast<ModelObject *>(data_object));
}

void ComparisonPainter::Painting(void)
{
    std::cout <<"- ComparisonPainter::Painting"<<std::endl;
    std::cout <<"  Folder path: "<< m_folder_path << std::endl;
    std::cout <<"  Number of model objects to be painted: "<< m_model_object_list.size() << std::endl;

    if (m_ref_model_object_list_map.find("with_charge") == m_ref_model_object_list_map.end() ||
        m_ref_model_object_list_map.find("no_charge") == m_ref_model_object_list_map.end())
    {
        std::cout <<"[ComparisonPainter::Painting] Missing simulation data [with_charge] and [no_charge] class. Skip..."<< std::endl;
        return;
    }

    PaintGausEstimateElementClassComparison("figure_3_a.pdf");
    PaintGausEstimateResidueClassComparison("figure_3_b.pdf");
    PaintGausEstimateResidueClassDenseComparison("figure_3_sup.pdf");

    for (auto & model_object : m_model_object_list)
    {
        auto plot_name{ "figure_2_charge_"+ model_object->GetPdbID() +".pdf" };
        auto & ref_model_object_list{ m_ref_model_object_list_map.at("with_charge") };
        if (model_object->GetPdbID() == ref_model_object_list.at(0)->GetPdbID())
        {
            PainMapValueComparison(plot_name, model_object, ref_model_object_list);
        }
    }

}

void ComparisonPainter::PaintGausEstimateElementClassComparison(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PaintGausEstimateElementClassComparison"<< std::endl;

    auto class_key{ AtomicInfoHelper::GetElementClassKey() };

    auto sim_with_charge_model_object_list{ m_ref_model_object_list_map.at("with_charge")};
    auto sim_no_charge_model_object_list{ m_ref_model_object_list_map.at("no_charge")};
    auto sim_amber95_model_object_list{ m_ref_model_object_list_map.at("amber95")};
    auto sim_test_model_object_list{ m_ref_model_object_list_map.at("sim_test")};

    const char * data_index[10]{"A","B","C","D","E","F","G","H","I","J"};

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int pad_size{ 5 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 2000, 800) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TPad> pad[pad_size];
    std::unique_ptr<TH2> frame[pad_size];
    pad[0] = ROOTHelper::CreatePad("pad_0","", 0.00, 0.70, 1.00, 1.00);
    pad[1] = ROOTHelper::CreatePad("pad_1","", 0.00, 0.00, 0.25, 0.70);
    pad[2] = ROOTHelper::CreatePad("pad_2","", 0.25, 0.00, 0.50, 0.70);
    pad[3] = ROOTHelper::CreatePad("pad_3","", 0.50, 0.00, 0.75, 0.70);
    pad[4] = ROOTHelper::CreatePad("pad_4","", 0.75, 0.00, 1.00, 0.70);

    const auto primary_element_size{ 4 };
    size_t ref_id[primary_element_size]{ 3, 3, 3, 1 };

    std::unique_ptr<TGraphErrors> data_graph[primary_element_size];
    std::unique_ptr<TGraphErrors> sim_with_charge_graph[primary_element_size];
    std::unique_ptr<TGraphErrors> sim_no_charge_graph[primary_element_size];
    std::unique_ptr<TGraphErrors> sim_additional_graph[primary_element_size];
    std::unique_ptr<TGraphErrors> sim_test_graph[primary_element_size];
    double x_min[primary_element_size]{ 0.0 };
    double x_max[primary_element_size]{ 1.0 };
    double y_min[primary_element_size]{ 0.0 };
    double y_max[primary_element_size]{ 1.0 };
    for (size_t i = 0; i < primary_element_size; i++)
    {
        data_graph[i] = ROOTHelper::CreateGraphErrors();
        sim_with_charge_graph[i] = ROOTHelper::CreateGraphErrors();
        sim_no_charge_graph[i] = ROOTHelper::CreateGraphErrors();
        sim_additional_graph[i] = ROOTHelper::CreateGraphErrors();
        sim_test_graph[i] = ROOTHelper::CreateGraphErrors();

        BuildAmplitudeRatioToWidthGraph(i, ref_id[i], data_graph[i].get(), m_model_object_list, class_key, true);
        BuildAmplitudeRatioToWidthGraph(i, ref_id[i], sim_with_charge_graph[i].get(), sim_with_charge_model_object_list, class_key);
        BuildAmplitudeRatioToWidthGraph(i, ref_id[i], sim_no_charge_graph[i].get(), sim_no_charge_model_object_list, class_key);
        BuildAmplitudeRatioToWidthGraph(i, ref_id[i], sim_additional_graph[i].get(), sim_amber95_model_object_list, class_key);
        BuildAmplitudeRatioToWidthGraph(i, ref_id[i], sim_test_graph[i].get(), sim_test_model_object_list, class_key);

        std::vector<double> x_array, y_array;
        for (int p = 0; p < data_graph[i]->GetN(); p++)
        {
            x_array.emplace_back(data_graph[i]->GetPointX(p));
            y_array.emplace_back(data_graph[i]->GetPointY(p));
        }
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.10) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);

        auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.30) };
        y_min[i] = std::get<0>(y_range);
        y_max[i] = std::get<1>(y_range);
    }

    canvas->cd();
    for (int i = 0; i < pad_size; i++)
    {
        ROOTHelper::SetPadDefaultStyle(pad[i].get());
        pad[i]->Draw();
    }

    pad[0]->cd();
    
    auto legend{ ROOTHelper::CreateLegend(0.06, 0.43, 0.60, 1.00, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 45.0f, 133, 12);
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    legend->SetNColumns(2);
    legend->AddEntry(sim_with_charge_graph[0].get(), "A. Thomas #font[2]{et al.}", "l");
    legend->AddEntry(sim_additional_graph[0].get(), "W. D. Cornell #font[2]{et al.}", "l");
    legend->AddEntry(sim_no_charge_graph[0].get(), "All Neutral #alpha = 0", "l");
    legend->AddEntry(sim_test_graph[0].get(), "Partial Charge #alpha_{#color[418]{#font[102]{N}}} = #minus0.1", "l");
    legend->Draw();

    std::unique_ptr<TPaveText> info_text[3];
    double x_pos[4]{ 0.60, 0.73, 0.86, 0.99 };
    for (int i = 0; i < 3; i++)
    {
        info_text[i] = ROOTHelper::CreatePaveText(x_pos[i], 0.43, x_pos[i+1], 1.00, "nbNDC", false);
        ROOTHelper::SetPaveTextDefaultStyle(info_text[i].get());
        ROOTHelper::SetFillAttribute(info_text[i].get(), 4000);
        ROOTHelper::SetTextAttribute(info_text[i].get(), 45.0f, 133, 12);
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

    auto top_title_text{ ROOTHelper::CreatePaveText(0.06, 0.01, 0.99, 0.40, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(top_title_text.get());
    ROOTHelper::SetPaveAttribute(top_title_text.get(), 0, 0.1);
    ROOTHelper::SetFillAttribute(top_title_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(top_title_text.get(), 75.0f, 133, 22, 0.0, kYellow-10);
    ROOTHelper::SetLineAttribute(top_title_text.get(), 1, 0);
    top_title_text->AddText("Integrated with All Residues");
    top_title_text->Draw();

    for (size_t i = 0; i < primary_element_size; i++)
    {
        pad[i+1]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.06, 0.01, 0.14, 0.05);
        ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 0, 0, 4000, 0);
        y_min[i] = (y_min[i] < 0.0) ? 0.0 : y_min[i];

        frame[i] = ROOTHelper::CreateHist2D(Form("hist_%d", static_cast<int>(i)),"", 500, x_min[i], x_max[i], 500, y_min[i], y_max[i]);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetXaxis(), 50.0f, 0.95f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetXaxis(), 45.0f, 0.005f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetXaxis(), 0.05f, 505);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetYaxis(), 45.0f, 1.2f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetYaxis(), 45.0f, 0.02f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetYaxis(), 0.05f, 505);
        ROOTHelper::SetLineAttribute(frame[i].get(), 1, 0);
        frame[i]->GetXaxis()->SetTitle("");
        frame[i]->GetXaxis()->CenterTitle();
        frame[i]->GetYaxis()->CenterTitle();
        frame[i]->SetStats(0);
        auto element_label{ AtomClassifier::GetMainChainElementLabel(i) };
        auto ref_element_label{ AtomClassifier::GetMainChainElementLabel(ref_id[i]) };
        auto element_color{  AtomClassifier::GetMainChainElementColor(i)};
        auto ref_element_color{  AtomClassifier::GetMainChainElementColor(ref_id[i]) };
        frame[i]->GetXaxis()->SetTitle(
            Form("Width #tau_{#font[102]{#color[%d]{%s}}}",
                element_color, element_label.data()));
        frame[i]->GetYaxis()->SetTitle(
            Form("Amplitude Ratio #font[1]{A}_{#font[102]{#color[%d]{%s}}} / #font[1]{A}_{#font[102]{#color[%d]{%s}}}",
                element_color, element_label.data(), ref_element_color, ref_element_label.data()));
        frame[i]->Draw();

        auto element_marker{ AtomClassifier::GetMainChainElementSolidMarker(i) };
        ROOTHelper::SetMarkerAttribute(data_graph[i].get(), element_marker, 2.0f, element_color);
        ROOTHelper::SetLineAttribute(sim_with_charge_graph[i].get(), 1, 2, kRed);
        ROOTHelper::SetLineAttribute(sim_no_charge_graph[i].get(),   2, 2, kGray+2);
        ROOTHelper::SetLineAttribute(sim_additional_graph[i].get(),  3, 2, kBlue);
        ROOTHelper::SetLineAttribute(sim_test_graph[i].get(),  4, 2, kViolet);
        sim_with_charge_graph[i]->Draw("C X0");
        sim_no_charge_graph[i]->Draw("C X0");
        sim_additional_graph[i]->Draw("C X0");
        sim_test_graph[i]->Draw("L X0");
        data_graph[i]->Draw("P X0");
    }

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ComparisonPainter::PaintGausEstimateResidueClassComparison(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PaintGausEstimateResidueClassComparison"<< std::endl;

    auto class_key{ AtomicInfoHelper::GetResidueClassKey() };
    
    auto sim_with_charge_model_object_list{ m_ref_model_object_list_map.at("with_charge")};
    auto sim_no_charge_model_object_list{ m_ref_model_object_list_map.at("no_charge")};
    auto sim_amber95_model_object_list{ m_ref_model_object_list_map.at("amber95")};

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int pad_size{ 5 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 2000, 800) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TPad> pad[pad_size];
    std::unique_ptr<TH2> frame[pad_size];
    pad[0] = ROOTHelper::CreatePad("pad_0","", 0.00, 0.70, 1.00, 1.00);
    pad[1] = ROOTHelper::CreatePad("pad_1","", 0.00, 0.00, 0.25, 0.70);
    pad[2] = ROOTHelper::CreatePad("pad_2","", 0.25, 0.00, 0.50, 0.70);
    pad[3] = ROOTHelper::CreatePad("pad_3","", 0.50, 0.00, 0.75, 0.70);
    pad[4] = ROOTHelper::CreatePad("pad_4","", 0.75, 0.00, 1.00, 0.70);

    const int primary_element_size{ 4 };
    size_t ref_id[primary_element_size]{ 3, 3, 3, 1 };

    std::vector<std::unique_ptr<TGraphErrors>> data_graph_list[primary_element_size][20];
    std::unique_ptr<TGraphErrors> data_graph[primary_element_size][20];
    std::unique_ptr<TGraphErrors> sim_with_charge_graph[primary_element_size][20];
    std::unique_ptr<TGraphErrors> sim_no_charge_graph[primary_element_size][20];
    std::unique_ptr<TGraphErrors> sim_additional_graph[primary_element_size][20];
    for (int k = 0; k < 20; k++)
    {
        auto residue{ static_cast<Residue>(k) };
        double x_min[primary_element_size]{ 0.0 };
        double x_max[primary_element_size]{ 1.0 };
        double y_min[primary_element_size]{ 0.0 };
        double y_max[primary_element_size]{ 1.0 };
        for (size_t i = 0; i < primary_element_size; i++)
        {
            data_graph[i][k] = ROOTHelper::CreateGraphErrors();
            sim_with_charge_graph[i][k] = ROOTHelper::CreateGraphErrors();
            sim_no_charge_graph[i][k] = ROOTHelper::CreateGraphErrors();
            sim_additional_graph[i][k] = ROOTHelper::CreateGraphErrors();
            BuildAmplitudeRatioToWidthGraph(i, ref_id[i], data_graph[i][k].get(), m_model_object_list, class_key, true, residue);
            BuildAmplitudeRatioToWidthGraph(i, ref_id[i], sim_with_charge_graph[i][k].get(), sim_with_charge_model_object_list, class_key, false, residue);
            BuildAmplitudeRatioToWidthGraph(i, ref_id[i], sim_no_charge_graph[i][k].get(), sim_no_charge_model_object_list, class_key, false, residue);
            BuildAmplitudeRatioToWidthGraph(i, ref_id[i], sim_additional_graph[i][k].get(), sim_amber95_model_object_list, class_key, false, residue);
            for (auto & model : m_model_object_list)
            {
                auto entry_iter{ std::make_unique<PotentialEntryIterator>(model) };
                data_graph_list[i][k].emplace_back(
                    entry_iter->CreateAmplitudeRatioToWidthScatterGraph(i, ref_id[i], residue));
            }

            std::vector<double> x_array, y_array;
            for (int p = 0; p < data_graph[i][k]->GetN(); p++)
            {
                x_array.emplace_back(data_graph[i][k]->GetPointX(p));
                y_array.emplace_back(data_graph[i][k]->GetPointY(p));
            }
            auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.05) };
            x_min[i] = std::get<0>(x_range);
            x_max[i] = std::get<1>(x_range);

            auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.30) };
            y_min[i] = std::get<0>(y_range);
            y_max[i] = std::get<1>(y_range);
        }

        canvas->cd();
        for (int i = 0; i < pad_size; i++)
        {
            ROOTHelper::SetPadDefaultStyle(pad[i].get());
            pad[i]->Draw();
        }

        pad[0]->cd();
        auto info_text{ ROOTHelper::CreatePaveText(0.06, 0.01, 0.99, 0.40, "nbNDC ARC", false) };
        ROOTHelper::SetPaveTextDefaultStyle(info_text.get());
        ROOTHelper::SetPaveAttribute(info_text.get(), 0, 0.1);
        ROOTHelper::SetFillAttribute(info_text.get(), 1001, kAzure-7);
        ROOTHelper::SetTextAttribute(info_text.get(), 100.0f, 103, 22, 0.0, kYellow-10);
        ROOTHelper::SetLineAttribute(info_text.get(), 1, 0);
        info_text->AddText(AtomicInfoHelper::GetLabel(residue).data());
        info_text->Draw();

        for (size_t i = 0; i < primary_element_size; i++)
        {
            pad[i+1]->cd();
            ROOTHelper::SetPadMarginInCanvas(gPad, 0.06, 0.01, 0.14, 0.05);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 0, 0, 4000, 0);
            y_min[i] = (y_min[i] < 0.0) ? 0.0 : y_min[i];

            frame[i] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", static_cast<int>(i), k),"", 500, x_min[i], x_max[i], 500, y_min[i], y_max[i]);
            ROOTHelper::SetAxisTitleAttribute(frame[i]->GetXaxis(), 50.0f, 0.95f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i]->GetXaxis(), 45.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i]->GetXaxis(), 0.05f, 505);
            ROOTHelper::SetAxisTitleAttribute(frame[i]->GetYaxis(), 45.0f, 1.3f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i]->GetYaxis(), 45.0f, 0.02f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i]->GetYaxis(), 0.05f, 505);
            ROOTHelper::SetLineAttribute(frame[i].get(), 1, 0);
            frame[i]->GetXaxis()->SetTitle("");
            frame[i]->GetXaxis()->CenterTitle();
            frame[i]->GetYaxis()->CenterTitle();
            frame[i]->SetStats(0);
            auto element_label{ AtomClassifier::GetMainChainElementLabel(i) };
            auto ref_element_label{ AtomClassifier::GetMainChainElementLabel(ref_id[i]) };
            auto element_color{  AtomClassifier::GetMainChainElementColor(i)};
            auto ref_element_color{  AtomClassifier::GetMainChainElementColor(ref_id[i]) }; 
            frame[i]->GetXaxis()->SetTitle(
                Form("Width #tau_{#font[102]{#color[%d]{%s}}}",
                    element_color, element_label.data()));
            frame[i]->GetYaxis()->SetTitle(
                Form("Amplitude Ratio #font[1]{A}_{#font[102]{#color[%d]{%s}}} / #font[1]{A}_{#font[102]{#color[%d]{%s}}}",
                    element_color, element_label.data(), ref_element_color, ref_element_label.data()));
            frame[i]->Draw();

            auto element_marker{ AtomClassifier::GetMainChainElementSolidMarker(i)};
            ROOTHelper::SetMarkerAttribute(data_graph[i][k].get(), element_marker, 2.0f, element_color);
            ROOTHelper::SetLineAttribute(sim_with_charge_graph[i][k].get(), 1, 2, kRed);
            ROOTHelper::SetLineAttribute(sim_no_charge_graph[i][k].get(),   2, 2, kGray+2);
            ROOTHelper::SetLineAttribute(sim_additional_graph[i][k].get(),  3, 2, kBlue);
            sim_with_charge_graph[i][k]->Draw("C X0");
            sim_no_charge_graph[i][k]->Draw("C X0");
            sim_additional_graph[i][k]->Draw("C X0");
            data_graph[i][k]->Draw("P X0");
            for (auto & graph : data_graph_list[i][k])
            {
                ROOTHelper::SetMarkerAttribute(graph.get(), element_marker, 1.0f, element_color);
                //graph->Draw("P");
            }
        }

        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ComparisonPainter::PaintGausEstimateResidueClassDenseComparison(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PaintGausEstimateResidueClassDenseComparison"<< std::endl;

    auto class_key{ AtomicInfoHelper::GetResidueClassKey() };
    
    auto sim_with_charge_model_object_list{ m_ref_model_object_list_map.at("with_charge")};
    auto sim_no_charge_model_object_list{ m_ref_model_object_list_map.at("no_charge")};
    auto sim_amber95_model_object_list{ m_ref_model_object_list_map.at("amber95")};

    const char * data_index[10]{"A","B","C","D","E","F","G","H","I","J"};

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 5 };
    const int row_size{ 4 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 2000, 1800) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.09f, 0.00f, 0.09f, 0.15f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int primary_element_size{ 4 };
    const int standard_residue_size{ 20 };
    size_t ref_id[primary_element_size]{ 3, 3, 3, 1 };

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
            auto residue{ static_cast<Residue>(i) };
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
            auto residue{ static_cast<Residue>(i) };
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
            title_text[i]->AddText(AtomicInfoHelper::GetLabel(residue).data());
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
        legend->AddEntry(sim_with_charge_graph_map[k].at(Residue::ALA).get(), "A. Thomas #font[2]{et al.}", "l");
        legend->AddEntry(sim_additional_graph_map[k].at(Residue::ALA).get(), "W. D. Cornell #font[2]{et al.}", "l");
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
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ComparisonPainter::PainMapValueComparison(
    const std::string & name,
    ModelObject * model_object,
    const std::vector<ModelObject *> & ref_model_object_list)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PainMapValueComparison"<< std::endl;

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 4 };
    const int row_size{ 1 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1400, 500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.07f, 0.01f, 0.15f, 0.11f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    for (auto & ref_model_object : ref_model_object_list)
    {
        if (model_object->GetPdbID() != ref_model_object->GetPdbID()) continue;
        std::unique_ptr<TGraphErrors> scatter_graph[col_size];
        std::unique_ptr<TF1> fit_function[col_size];
        double r_square[col_size]{ 0.0 };
        double slope[col_size]{ 0.0 };
        double intercept[col_size]{ 0.0 };
        for (size_t i = 0; i < col_size; i++)
        {
            auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
            scatter_graph[i] = ROOTHelper::CreateGraphErrors();
            BuildMapValueScatterGraph(group_key, scatter_graph[i].get(), ref_model_object, model_object, 15, 0.0, 1.5);
            r_square[i] = ROOTHelper::PerformLinearRegression(scatter_graph[i].get(), slope[i], intercept[i]);
            fit_function[i] = ROOTHelper::CreateFunction1D(Form("fit_%d", static_cast<int>(i)), "x*[1]+[0]");
            fit_function[i]->SetParameters(intercept[i], slope[i]);
        }

        double x_min[col_size]{ 0.0 };
        double x_max[col_size]{ 1.0 };
        std::vector<double> y_array;
        for (int i = 0; i < col_size; i++)
        {
            std::vector<double> x_array;
            for (int p = 0; p < scatter_graph[i]->GetN(); p++)
            {
                x_array.emplace_back(scatter_graph[i]->GetPointX(p));
                y_array.emplace_back(scatter_graph[i]->GetPointY(p));
            }
            auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.24) };
            x_min[i] = std::get<0>(x_range);
            x_max[i] = std::get<1>(x_range);
        }

        auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.15) };
        auto y_min{ std::get<0>(y_range) };
        auto y_max{ std::get<1>(y_range) };

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
                ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
                auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
                auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
                frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", static_cast<int>(i), j),"", 500, x_min[i], x_max[i], 500, y_min, y_max);
                ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 35, 0.9f, 133);
                ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 35, 0.005f, 133);
                ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
                ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 40, 1.1f, 133);
                ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 40, 0.02f, 133);
                ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
                ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
                frame[i][j]->GetXaxis()->SetTitle("Simulated Map Value");
                frame[i][j]->GetYaxis()->SetTitle("Real Map Value");
                frame[i][j]->GetXaxis()->CenterTitle();
                frame[i][j]->GetYaxis()->CenterTitle();
                frame[i][j]->SetStats(0);
                frame[i][j]->Draw();
                ROOTHelper::SetMarkerAttribute(scatter_graph[i].get(), element_marker, 1.0f, element_color);
                ROOTHelper::SetLineAttribute(scatter_graph[i].get(), 1, 2, element_color);
                scatter_graph[i]->Draw("P");
                ROOTHelper::SetLineAttribute(fit_function[i].get(), 1, 2, kRed);
                fit_function[i]->SetRange(x_min[i], x_max[i]);
                fit_function[i]->Draw("SAME");
            }
            title_text[i] = ROOTHelper::CreatePaveText(0.01, 1.01, 0.99, 1.14, "nbNDC ARC", true);
            ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
            ROOTHelper::SetPaveAttribute(title_text[i].get(), 0, 0.2);
            ROOTHelper::SetTextAttribute(title_text[i].get(), 40, 133, 22);
            ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, element_color, 0.5f);
            title_text[i]->AddText(element_label.data());
            title_text[i]->Draw();

            r_square_text[i] = ROOTHelper::CreatePaveText(0.50, 0.05, 0.95, 0.18, "nbNDC ARC", true);
            ROOTHelper::SetPaveTextDefaultStyle(r_square_text[i].get());
            ROOTHelper::SetPaveAttribute(r_square_text[i].get(), 0, 0.5);
            ROOTHelper::SetLineAttribute(r_square_text[i].get(), 1, 0);
            ROOTHelper::SetTextAttribute(r_square_text[i].get(), 35, 133, 22);
            ROOTHelper::SetFillAttribute(r_square_text[i].get(), 1001, kAzure-7, 0.20f);
            r_square_text[i]->AddText(Form("R^{2} = %.2f", r_square[i]));
            r_square_text[i]->Draw();

            fit_info_text[i] = ROOTHelper::CreatePaveText(0.12, 0.88, 0.70, 0.99, "nbNDC", true);
            ROOTHelper::SetPaveTextDefaultStyle(fit_info_text[i].get());
            ROOTHelper::SetTextAttribute(fit_info_text[i].get(), 35, 133, 22, 0.0, kGray);
            fit_info_text[i]->AddText(Form("#font[1]{y} = %.2f#font[1]{x}%+.2f", slope[i], intercept[i]));
            fit_info_text[i]->Draw();
        }
        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

#ifdef HAVE_ROOT
void ComparisonPainter::BuildAmplitudeRatioToWidthGraph(
    size_t target_id, size_t reference_id, TGraphErrors * graph,
    const std::vector<ModelObject *> & model_list, const std::string & class_key,
    bool draw_index, Residue residue)
{
    const char * data_index[11]{"A","B","C","D","E","F","G","H","I","J","K"};
    auto group_key
    {
        (class_key == AtomicInfoHelper::GetElementClassKey()) ?
        m_atom_classifier->GetMainChainElementClassGroupKey(target_id) :
        m_atom_classifier->GetMainChainResidueClassGroupKey(target_id, residue)
    };
    auto ref_group_key
    {
        (class_key == AtomicInfoHelper::GetElementClassKey()) ?
        m_atom_classifier->GetMainChainElementClassGroupKey(reference_id) :
        m_atom_classifier->GetMainChainResidueClassGroupKey(reference_id, residue)
    };
    auto count{ 0 };
    auto model_count{ 0 };
    for (auto model_object : model_list)
    {
        model_count++;
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
        if (entry_iter->IsAvailableGroupKey(group_key, class_key) == false) continue;
        if (entry_iter->IsAvailableGroupKey(ref_group_key, class_key) == false) continue;
        auto x_value{ entry_iter->GetGausEstimatePrior(group_key, class_key, 1) };
        auto y_value{ entry_iter->GetGausEstimatePrior(group_key, class_key, 0) };
        auto x_error{ entry_iter->GetGausVariancePrior(group_key, class_key, 1) };
        auto y_error{ entry_iter->GetGausVariancePrior(group_key, class_key, 0) };
        auto ref_y_value{ entry_iter->GetGausEstimatePrior(ref_group_key, class_key, 0) };
        auto ref_y_error{ entry_iter->GetGausVariancePrior(ref_group_key, class_key, 0) };
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
                std::cout <<"Warning: Data label size exceeds 12, label switch to numbers."<< std::endl;
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

void ComparisonPainter::BuildMapValueScatterGraph(
    uint64_t group_key, TGraphErrors * graph, ModelObject * model1, ModelObject * model2,
    int bin_size, double x_min, double x_max)
{
    auto entry1_iter{ std::make_unique<PotentialEntryIterator>(model1) };
    auto entry2_iter{ std::make_unique<PotentialEntryIterator>(model2) };
    if (entry1_iter->IsAvailableGroupKey(group_key, AtomicInfoHelper::GetElementClassKey()) == false) return;
    if (entry2_iter->IsAvailableGroupKey(group_key, AtomicInfoHelper::GetElementClassKey()) == false) return;
    auto model1_atom_map{ entry1_iter->GetAtomObjectMap(group_key, AtomicInfoHelper::GetElementClassKey()) };
    auto model2_atom_map{ entry2_iter->GetAtomObjectMap(group_key, AtomicInfoHelper::GetElementClassKey()) };
    auto count{ 0 };
    for (auto & [atom_id, atom_object1] : model1_atom_map)
    {
        if (model2_atom_map.find(atom_id) == model2_atom_map.end()) continue;
        auto atom_object2{ model2_atom_map.at(atom_id) };
        auto atom1_iter{ std::make_unique<PotentialEntryIterator>(atom_object1) };
        auto atom2_iter{ std::make_unique<PotentialEntryIterator>(atom_object2) };
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