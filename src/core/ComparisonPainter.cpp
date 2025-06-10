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
    m_model_object_list.push_back(dynamic_cast<ModelObject *>(data_object));
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
    std::cout <<"  Number of reference types to be painted: "<< m_ref_model_object_list_map.size() << std::endl;

    if (m_ref_model_object_list_map.find("with_charge") != m_ref_model_object_list_map.end())
    {
        auto sim_with_charge_model_object_list{ m_ref_model_object_list_map.at("with_charge") };

        if (m_model_object_list.size() != 0 && m_ref_model_object_list_map.size() >= 2 && sim_with_charge_model_object_list.size() > 1)
        {
            //PaintSimulationGaus("figure_3_a.pdf");
            PaintSimulationGausEstimate("figure_3_a.pdf");
            PaintGausEstimateComparison("figure_3_b.pdf");
            PaintGausEstimateResidueComparison("figure_3_d.pdf");
            PaintGausEstimateScatterComparison("figure_3_c.pdf");
        }
        
        if (sim_with_charge_model_object_list.size() > 1)
        {
            //PaintSimulationGausRatio("figure_2_c.pdf", sim_with_charge_model_object_list);
        }

        if (m_model_object_list.size() == 1 && sim_with_charge_model_object_list.size() == 1)
        {
            PainMapValueComparison("figure_2_a.pdf", m_model_object_list.at(0), sim_with_charge_model_object_list.at(0));
            //PainResidueClassGausComparison("figure_width_comparison.pdf", m_model_object_list.at(0), sim_with_charge_model_object_list.at(0), 1);
        }
    }
}

void ComparisonPainter::PaintSimulationGaus(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PaintSimulationGaus"<< std::endl;

    auto sim_with_charge_model_object_list{ m_ref_model_object_list_map.at("with_charge")};
    auto sim_no_charge_model_object_list{ m_ref_model_object_list_map.at("no_charge")};

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 4 };
    const int row_size{ 2 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 950, 600) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.1f, 0.01f, 0.11f, 0.07f, 0.01f, 0.005f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    
    const int primary_element_size{ 4 };
    const char * element_label[primary_element_size]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };
    const std::string y_axis_title[row_size]
    {
        "Width",
        "Amplitude"
    };
    short color_element[primary_element_size] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };

    std::unique_ptr<TGraphErrors> sim_with_charge_graph[primary_element_size][row_size];
    std::unique_ptr<TGraphErrors> sim_no_charge_graph[primary_element_size][row_size];
    std::vector<double> x_array, y_array[row_size];
    for (size_t i = 0; i < primary_element_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        for (size_t j = 0; j < row_size; j++)
        {
            sim_with_charge_graph[i][j] = ROOTHelper::CreateGraphErrors();
            sim_no_charge_graph[i][j] = ROOTHelper::CreateGraphErrors();

            int reverse_index{ row_size - static_cast<int>(j) - 1 };
            BuildGausEstimateToBlurringWidthGraph(group_key, sim_with_charge_graph[i][j].get(), sim_with_charge_model_object_list, reverse_index);
            BuildGausEstimateToBlurringWidthGraph(group_key, sim_no_charge_graph[i][j].get(), sim_no_charge_model_object_list, reverse_index);

            for (int p = 0; p < sim_with_charge_graph[i][j]->GetN(); p++)
            {
                x_array.push_back(sim_with_charge_graph[i][j]->GetPointX(p));
                y_array[j].push_back(sim_with_charge_graph[i][j]->GetPointY(p));
            }
            for (int p = 0; p < sim_no_charge_graph[i][j]->GetN(); p++)
            {
                x_array.push_back(sim_no_charge_graph[i][j]->GetPointX(p));
                y_array[j].push_back(sim_no_charge_graph[i][j]->GetPointY(p));
            }

            ROOTHelper::SetMarkerAttribute(sim_with_charge_graph[i][j].get(), 20, 1.2f, color_element[i]);
            ROOTHelper::SetMarkerAttribute(sim_no_charge_graph[i][j].get(), 53, 1.2f, color_element[i]);
            ROOTHelper::SetLineAttribute(sim_with_charge_graph[i][j].get(), 1, 2, color_element[i]);
            ROOTHelper::SetLineAttribute(sim_no_charge_graph[i][j].get(), 2, 2, color_element[i]);
        }
    }

    auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.05) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    double y_min[row_size], y_max[row_size];
    for (int j = 0; j < row_size; j++)
    {
        auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array[j], 0.3) };
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> title_text[col_size];
    std::unique_ptr<TLegend> legend;
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min, x_max, 500, y_min[j], y_max[j]);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 30, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 30, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 35, 1.3f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("Blurring Width");
            frame[i][j]->GetYaxis()->SetTitle(y_axis_title[j].data());
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            sim_with_charge_graph[i][j]->Draw("PL X0");
            sim_no_charge_graph[i][j]->Draw("PL X0");
            if (i == col_size - 1 && j == row_size - 1)
            {
                legend = ROOTHelper::CreateLegend(0.1, 0.7, 0.9, 0.99, true);
                ROOTHelper::SetLegendDefaultStyle(legend.get());
                ROOTHelper::SetTextAttribute(legend.get(), 20, 133, 12);
                ROOTHelper::SetFillAttribute(legend.get(), 1001, kWhite, 0.5);
                legend->AddEntry(sim_with_charge_graph[i][1].get(), "Partial Charge", "pl");
                legend->AddEntry(sim_no_charge_graph[i][1].get(), "Neutral", "pl");
                legend->Draw();
            }
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

void ComparisonPainter::PaintSimulationGausRatio(
    const std::string & name, const std::vector<ModelObject *> & model_list)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PaintSimulationGausRatio"<< std::endl;

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 2 };
    const int row_size{ 1 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 600) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.11f, 0.11f, 0.15f, 0.05f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    
    const int primary_element_size{ 4 };
    const std::string y_axis_title[col_size]
    {
        "Amplitude Ratio A_{C#alpha} / A_{C}",
        "Width Ratio #tau_{C#alpha} / #tau_{C}"
    };

    std::unique_ptr<TGraphErrors> simulation_graph[primary_element_size][2];
    for (size_t i = 0; i < primary_element_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        for (int j = 0; j < 2; j++)
        {
            simulation_graph[i][j] = ROOTHelper::CreateGraphErrors();
            BuildGausEstimateToBlurringWidthGraph(group_key, simulation_graph[i][j].get(), model_list, j);
        }
    }

    std::unique_ptr<TGraphErrors> ratio_graph[col_size];
    std::vector<double> x_array, y_array;
    for (int i = 0; i < col_size; i++)
    {
        ratio_graph[i] = ROOTHelper::CreateGraphErrors();
        BuildRatioGraph(ratio_graph[i].get(), simulation_graph[0][i].get(), simulation_graph[1][i].get());

        for (int p = 0; p < ratio_graph[i]->GetN(); p++)
        {
            x_array.push_back(ratio_graph[i]->GetPointX(p));
            y_array.push_back(ratio_graph[i]->GetPointY(p));
        }
        ROOTHelper::SetMarkerAttribute(ratio_graph[i].get(), 20, 1.2f, kRed+1);
        ROOTHelper::SetLineAttribute(ratio_graph[i].get(), 1, 2, kRed+1);
    }

    auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.01) };
    auto x_min{ std::get<0>(x_range) };
    auto x_max{ std::get<1>(x_range) };
    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.3) };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };

    std::unique_ptr<TH2> frame[col_size][row_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min, x_max, 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 40, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 40, 1.3f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 40, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("Blurring Width");
            frame[i][j]->GetYaxis()->SetTitle(y_axis_title[i].data());
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            
            if (i == 0)
            {
                frame[i][j]->Draw();
                ratio_graph[i]->Draw("PL X0");
            }
            else
            {
                frame[i][j]->Draw("Y+");
                ratio_graph[i]->Draw("PL X0");
            }
        }
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ComparisonPainter::PaintSimulationGausEstimate(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PaintSimulationGausEstimate"<< std::endl;

    auto sim_with_charge_model_object_list{ m_ref_model_object_list_map.at("with_charge")};
    auto sim_no_charge_model_object_list{ m_ref_model_object_list_map.at("no_charge")};

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int pad_size{ 4 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 2000, 500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TPad> pad[pad_size];
    std::unique_ptr<TH2> frame[pad_size];
    pad[0] = ROOTHelper::CreatePad("pad_0","", 0.00, 0.00, 0.25, 1.00);
    pad[1] = ROOTHelper::CreatePad("pad_1","", 0.25, 0.00, 0.50, 1.00);
    pad[2] = ROOTHelper::CreatePad("pad_2","", 0.50, 0.00, 0.75, 1.00);
    pad[3] = ROOTHelper::CreatePad("pad_3","", 0.75, 0.00, 1.00, 1.00);


    const int primary_element_size{ 4 };
    const char * element_label[primary_element_size]{ "C_{#alpha}", "C", "N", "O" };

    short color_element[primary_element_size] { kRed+1, kViolet+1, kGreen+2, kAzure+2 };
    std::vector<std::string> element_list{"Ap", "Cp", "Nn", "On"};
    std::vector<std::string> charge_list{"1", "3", "5", "7"};

    std::unique_ptr<TGraphErrors> sim_graph[primary_element_size][4+2];
    for (size_t i = 0; i < primary_element_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        sim_graph[i][4] = ROOTHelper::CreateGraphErrors();
        sim_graph[i][5] = ROOTHelper::CreateGraphErrors();

        BuildAmplitudeToWidthGraph(group_key, sim_graph[i][4].get(), sim_with_charge_model_object_list);
        BuildAmplitudeToWidthGraph(group_key, sim_graph[i][5].get(), sim_no_charge_model_object_list);
        
        for (size_t j = 0; j < 4; j++)
        {
            sim_graph[i][j] = ROOTHelper::CreateGraphErrors();
            auto tmp_name{ element_list.at(i) + charge_list.at(j) };
            if (m_ref_model_object_list_map.find(tmp_name) == m_ref_model_object_list_map.end())
            {
                std::cout << tmp_name << std::endl;
                continue;
            }
            BuildAmplitudeToWidthGraph(group_key, sim_graph[i][j].get(), m_ref_model_object_list_map.at(tmp_name));
        }
    }

    double x_min[primary_element_size]{ 0.0 };
    double x_max[primary_element_size]{ 1.0 };
    double y_min[primary_element_size]{ 0.0 };
    double y_max[primary_element_size]{ 1.0 };
    for (int i = 0; i < primary_element_size; i++)
    {
        std::vector<double> x_array, y_array;
        for (size_t j = 0; j < 6; j++)
        {
            for (int p = 0; p < sim_graph[i][j]->GetN(); p++)
            {
                x_array.emplace_back(sim_graph[i][j]->GetPointX(p));
                y_array.emplace_back(sim_graph[i][j]->GetPointY(p));
            }
        }
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.05) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);

        auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.05) };
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
    auto line_ref1{ ROOTHelper::CreateLine(0.0, 0.0, 1.0, 1.0) };
    auto line_ref2{ ROOTHelper::CreateLine(0.0, 0.0, 1.0, 1.0) };
    ROOTHelper::SetLineAttribute(line_ref1.get(), 1, 5, kGray+2);
    ROOTHelper::SetLineAttribute(line_ref2.get(), 2, 5, kGray+2);

    std::unique_ptr<TPaveText> title_text[primary_element_size];
    short line_style[6]{ 1, 1, 1, 1, 1, 2 };
    short line_color[6]{ kGray+2, kGray+2, kGray+2, kGray+2, kRed, kRed };
    float transparent[6]{ 0.2f, 0.4f, 0.6f, 0.8f, 1.0f, 1.0f };
    for (int i = 0; i < primary_element_size; i++)
    {
        pad[i]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.07, 0.02, 0.22, 0.05);
        ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 0, 0, 4000, 0);

        x_min[i] = (x_min[i] < 0.0) ? 0.0 : x_min[i];
        x_max[i] = (i == 3) ? 0.7 : 1.0;
        y_min[i] = 10.0;
        y_max[i] = (i == 3) ? 60.0 : 150.0;

        frame[i] = ROOTHelper::CreateHist2D(Form("hist_%d", i),"", 500, x_min[i], x_max[i], 500, y_min[i], y_max[i]);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetXaxis(), 50.0f, 0.95f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetXaxis(), 45.0f, 0.005f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetXaxis(), 0.05f, 506);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetYaxis(), 50.0f, 1.3f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetYaxis(), 45.0f, 0.02f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetYaxis(), 0.05f, 506);
        ROOTHelper::SetLineAttribute(frame[i].get(), 1, 0);
        frame[i]->GetXaxis()->SetTitle("");
        frame[i]->GetXaxis()->CenterTitle();
        frame[i]->GetYaxis()->CenterTitle();
        frame[i]->SetStats(0);
        frame[i]->GetXaxis()->SetTitle(Form("Width #tau_{#font[102]{#color[%d]{%s}}}", color_element[i], element_label[i]));
        frame[i]->GetYaxis()->SetTitle(Form("Amplitude #font[1]{A}_{#font[102]{#color[%d]{%s}}}", color_element[i], element_label[i]));
        frame[i]->Draw();

        for (size_t j = 0; j < 6; j++)
        {
            ROOTHelper::SetLineAttribute(sim_graph[i][j].get(), line_style[j], 2, line_color[j], transparent[j]);
            sim_graph[i][j]->Draw("L X0");
        }
    }

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ComparisonPainter::PaintGausEstimateComparison(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PaintGausEstimateComparison"<< std::endl;

    auto sim_with_charge_model_object_list{ m_ref_model_object_list_map.at("with_charge")};
    auto sim_no_charge_model_object_list{ m_ref_model_object_list_map.at("no_charge")};
    auto sim_amber95_model_object_list{ m_ref_model_object_list_map.at("amber95")};
    auto sim_test_model_object_list{ m_ref_model_object_list_map.at("sim_test")};

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int pad_size{ 5 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 2000, 2000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TPad> pad[pad_size];
    std::unique_ptr<TH2> frame[pad_size];
    pad[0] = ROOTHelper::CreatePad("pad_0","", 0.00, 0.80, 1.00, 1.00);
    pad[1] = ROOTHelper::CreatePad("pad_1","", 0.00, 0.40, 0.50, 0.80);
    pad[2] = ROOTHelper::CreatePad("pad_2","", 0.50, 0.40, 1.00, 0.80);
    pad[3] = ROOTHelper::CreatePad("pad_3","", 0.00, 0.00, 0.50, 0.40);
    pad[4] = ROOTHelper::CreatePad("pad_4","", 0.50, 0.00, 1.00, 0.40);


    const int primary_element_size{ 4 };
    const char * element_label[primary_element_size]{ "C_{#alpha}", "C", "N", "O" };

    short color_element[primary_element_size] { kRed+1, kViolet+1, kGreen+2, kAzure+2 };
    short marker_element[primary_element_size]{ 21, 20, 22, 23 };

    std::unique_ptr<TGraphErrors> data_graph[primary_element_size][2];
    std::unique_ptr<TGraphErrors> sim_with_charge_graph[primary_element_size][2];
    std::unique_ptr<TGraphErrors> sim_no_charge_graph[primary_element_size][2];
    std::unique_ptr<TGraphErrors> sim_additional_graph[primary_element_size][2];
    std::unique_ptr<TGraphErrors> sim_test_graph[primary_element_size][2];
    for (size_t i = 0; i < primary_element_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        data_graph[i][0] = ROOTHelper::CreateGraphErrors();
        sim_with_charge_graph[i][0] = ROOTHelper::CreateGraphErrors();
        sim_no_charge_graph[i][0] = ROOTHelper::CreateGraphErrors();
        sim_additional_graph[i][0] = ROOTHelper::CreateGraphErrors();
        sim_test_graph[i][0] = ROOTHelper::CreateGraphErrors();
        BuildAmplitudeToWidthGraph(group_key, data_graph[i][0].get(), m_model_object_list);
        BuildAmplitudeToWidthGraph(group_key, sim_with_charge_graph[i][0].get(), sim_with_charge_model_object_list);
        BuildAmplitudeToWidthGraph(group_key, sim_no_charge_graph[i][0].get(), sim_no_charge_model_object_list);
        BuildAmplitudeToWidthGraph(group_key, sim_additional_graph[i][0].get(), sim_amber95_model_object_list);
        BuildAmplitudeToWidthGraph(group_key, sim_test_graph[i][0].get(), sim_test_model_object_list);
    }

    const char * data_index[10]{"A","B","C","D","E","F","G","H","I","J"};
    std::vector<double> resolution_list;
    resolution_list.reserve(m_model_object_list.size());
    for (auto & model : m_model_object_list)
    {
        resolution_list.emplace_back(model->GetResolution());
    }

    int ref_id[primary_element_size]{ 3, 3, 3, 1 };
    double x_min[primary_element_size]{ 0.0 };
    double x_max[primary_element_size]{ 1.0 };
    double y_min[primary_element_size]{ 0.0 };
    double y_max[primary_element_size]{ 1.0 };
    for (int i = 0; i < primary_element_size; i++)
    {
        data_graph[i][1] = ROOTHelper::CreateGraphErrors();
        sim_with_charge_graph[i][1] = ROOTHelper::CreateGraphErrors();
        sim_no_charge_graph[i][1] = ROOTHelper::CreateGraphErrors();
        sim_additional_graph[i][1] = ROOTHelper::CreateGraphErrors();
        sim_test_graph[i][1] = ROOTHelper::CreateGraphErrors();
        BuildRatioGraph(data_graph[i][1].get(), data_graph[i][0].get(), data_graph[ref_id[i]][0].get(), true);
        BuildRatioGraph(sim_with_charge_graph[i][1].get(), sim_with_charge_graph[i][0].get(), sim_with_charge_graph[ref_id[i]][0].get());
        BuildRatioGraph(sim_no_charge_graph[i][1].get(), sim_no_charge_graph[i][0].get(), sim_no_charge_graph[ref_id[i]][0].get());
        BuildRatioGraph(sim_additional_graph[i][1].get(), sim_additional_graph[i][0].get(), sim_additional_graph[ref_id[i]][0].get());
        BuildRatioGraph(sim_test_graph[i][1].get(), sim_test_graph[i][0].get(), sim_test_graph[ref_id[i]][0].get());

        std::vector<double> x_array, y_array;
        for (int p = 0; p < data_graph[i][1]->GetN(); p++)
        {
            x_array.emplace_back(data_graph[i][1]->GetPointX(p));
            y_array.emplace_back(data_graph[i][1]->GetPointY(p));
        }
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.20) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);

        auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.70) };
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
    
    auto legend{ ROOTHelper::CreateLegend(0.05, 0.00, 0.60, 1.00, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 70.0f, 133, 12);
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    legend->AddEntry(sim_with_charge_graph[0][1].get(), "A. Thomas #font[2]{et al.}", "l");
    legend->AddEntry(sim_additional_graph[0][1].get(), "Amber (W. D. Cornell #font[2]{et al.})", "l");
    legend->AddEntry(sim_no_charge_graph[0][1].get(), "Neutral", "l");
    legend->AddEntry(sim_test_graph[0][1].get(), "#alpha_{#font[102]{N}} = #minus0.1", "l");
    legend->Draw();

    auto info_text1{ ROOTHelper::CreatePaveText(0.60, 0.00, 0.80, 1.00, "nbNDC", false) };
    auto info_text2{ ROOTHelper::CreatePaveText(0.80, 0.00, 1.00, 1.00, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(info_text1.get());
    ROOTHelper::SetPaveTextDefaultStyle(info_text2.get());
    ROOTHelper::SetFillAttribute(info_text1.get(), 4000);
    ROOTHelper::SetFillAttribute(info_text2.get(), 4000);
    ROOTHelper::SetTextAttribute(info_text1.get(), 70.0f, 133, 12);
    ROOTHelper::SetTextAttribute(info_text2.get(), 70.0f, 133, 12);
    auto count{ 0 };
    for (auto & line : resolution_list)
    {
        if (count < 4) info_text1->AddText(Form("#font[102]{[#color[853]{%s}]} %.2f #AA", data_index[count], line));
        else info_text2->AddText(Form("#font[102]{[#color[853]{%s}]} %.2f #AA", data_index[count], line));
        count++;
        if (count >= 10) break;
    }
    info_text2->AddText("");
    info_text1->Draw();
    info_text2->Draw();

    std::unique_ptr<TPaveText> title_text[primary_element_size];
    for (int i = 0; i < primary_element_size; i++)
    {
        pad[i+1]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.10, 0.03, 0.08, 0.02);
        ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 0, 0, 4000, 0);
        y_min[i] = (y_min[i] < 0.0) ? 0.0 : y_min[i];
        y_max[i] = (y_max[i] < 1.0) ? 1.0 : y_max[i];

        frame[i] = ROOTHelper::CreateHist2D(Form("hist_%d", i),"", 500, x_min[i], x_max[i], 500, y_min[i], y_max[i]);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetXaxis(), 70.0f, 1.0f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetXaxis(), 65.0f, 0.005f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetXaxis(), 0.05f, 506);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetYaxis(), 70.0f, 1.3f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetYaxis(), 65.0f, 0.02f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetYaxis(), 0.05f, 506);
        ROOTHelper::SetLineAttribute(frame[i].get(), 1, 0);
        frame[i]->GetXaxis()->SetTitle("");
        frame[i]->GetXaxis()->CenterTitle();
        frame[i]->GetYaxis()->CenterTitle();
        frame[i]->SetStats(0);
        frame[i]->GetXaxis()->SetTitle(Form("Width #tau_{#font[102]{#color[%d]{%s}}}", color_element[i], element_label[i]));
        frame[i]->GetYaxis()->SetTitle(Form("Amplitude Ratio #font[1]{A}_{#font[102]{#color[%d]{%s}}} / #font[1]{A}_{#font[102]{#color[%d]{%s}}}", color_element[i], element_label[i], color_element[ref_id[i]], element_label[ref_id[i]]));
        frame[i]->Draw();

        ROOTHelper::SetMarkerAttribute(data_graph[i][1].get(), marker_element[i], 3.0f, color_element[i]);
        ROOTHelper::SetLineAttribute(sim_with_charge_graph[i][1].get(), 1, 4, kRed);
        ROOTHelper::SetLineAttribute(sim_no_charge_graph[i][1].get(),   2, 4, kGray+2);
        ROOTHelper::SetLineAttribute(sim_additional_graph[i][1].get(),  3, 4, kBlue);
        ROOTHelper::SetLineAttribute(sim_test_graph[i][1].get(),  4, 4, kViolet);
        sim_with_charge_graph[i][1]->Draw("L X0");
        sim_no_charge_graph[i][1]->Draw("L X0");
        sim_additional_graph[i][1]->Draw("L X0");
        sim_test_graph[i][1]->Draw("L X0");
        data_graph[i][1]->Draw("P X0");
    }

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ComparisonPainter::PaintGausEstimateResidueComparison(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PaintGausEstimateResidueComparison"<< std::endl;

    auto sim_with_charge_model_object_list{ m_ref_model_object_list_map.at("with_charge")};
    auto sim_no_charge_model_object_list{ m_ref_model_object_list_map.at("no_charge")};
    auto sim_amber95_model_object_list{ m_ref_model_object_list_map.at("amber95")};

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int pad_size{ 5 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 2000, 2000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TPad> pad[pad_size];
    std::unique_ptr<TH2> frame[pad_size];
    pad[0] = ROOTHelper::CreatePad("pad_0","", 0.00, 0.80, 1.00, 1.00);
    pad[1] = ROOTHelper::CreatePad("pad_1","", 0.00, 0.40, 0.50, 0.80);
    pad[2] = ROOTHelper::CreatePad("pad_2","", 0.50, 0.40, 1.00, 0.80);
    pad[3] = ROOTHelper::CreatePad("pad_3","", 0.00, 0.00, 0.50, 0.40);
    pad[4] = ROOTHelper::CreatePad("pad_4","", 0.50, 0.00, 1.00, 0.40);


    const int primary_element_size{ 4 };
    const char * element_label[primary_element_size]{ "C_{#alpha}", "C", "N", "O" };
    const std::string & class_key{ AtomicInfoHelper::GetResidueClassKey() };

    short color_element[primary_element_size] { kRed+1, kViolet+1, kGreen+2, kAzure+2 };
    short marker_element[primary_element_size]{ 21, 20, 22, 23 };

    for (int k = 0; k < 20; k++)
    {
        auto residue{ static_cast<Residue>(k) };
        std::unique_ptr<TGraphErrors> data_graph[primary_element_size][2];
        std::unique_ptr<TGraphErrors> sim_with_charge_graph[primary_element_size][2];
        std::unique_ptr<TGraphErrors> sim_no_charge_graph[primary_element_size][2];
        std::unique_ptr<TGraphErrors> sim_additional_graph[primary_element_size][2];
        for (size_t i = 0; i < primary_element_size; i++)
        {
            auto group_key{ m_atom_classifier->GetMainChainResidueClassGroupKey(i, residue) };
            data_graph[i][0] = ROOTHelper::CreateGraphErrors();
            sim_with_charge_graph[i][0] = ROOTHelper::CreateGraphErrors();
            sim_no_charge_graph[i][0] = ROOTHelper::CreateGraphErrors();
            sim_additional_graph[i][0] = ROOTHelper::CreateGraphErrors();
            BuildAmplitudeToWidthGraph(group_key, data_graph[i][0].get(), m_model_object_list, class_key);
            BuildAmplitudeToWidthGraph(group_key, sim_with_charge_graph[i][0].get(), sim_with_charge_model_object_list, class_key);
            BuildAmplitudeToWidthGraph(group_key, sim_no_charge_graph[i][0].get(), sim_no_charge_model_object_list, class_key);
            BuildAmplitudeToWidthGraph(group_key, sim_additional_graph[i][0].get(), sim_amber95_model_object_list, class_key);
        }

        const char * data_index[10]{"A","B","C","D","E","F","G","H","I","J"};
        std::vector<double> resolution_list;
        resolution_list.reserve(m_model_object_list.size());
        for (auto & model : m_model_object_list)
        {
            resolution_list.emplace_back(model->GetResolution());
        }

        int ref_id[primary_element_size]{ 3, 3, 3, 1 };
        double x_min[primary_element_size]{ 0.0 };
        double x_max[primary_element_size]{ 1.0 };
        double y_min[primary_element_size]{ 0.0 };
        double y_max[primary_element_size]{ 1.0 };
        for (int i = 0; i < primary_element_size; i++)
        {
            data_graph[i][1] = ROOTHelper::CreateGraphErrors();
            sim_with_charge_graph[i][1] = ROOTHelper::CreateGraphErrors();
            sim_no_charge_graph[i][1] = ROOTHelper::CreateGraphErrors();
            sim_additional_graph[i][1] = ROOTHelper::CreateGraphErrors();
            BuildRatioGraph(data_graph[i][1].get(), data_graph[i][0].get(), data_graph[ref_id[i]][0].get(), true);
            BuildRatioGraph(sim_with_charge_graph[i][1].get(), sim_with_charge_graph[i][0].get(), sim_with_charge_graph[ref_id[i]][0].get());
            BuildRatioGraph(sim_no_charge_graph[i][1].get(), sim_no_charge_graph[i][0].get(), sim_no_charge_graph[ref_id[i]][0].get());
            BuildRatioGraph(sim_additional_graph[i][1].get(), sim_additional_graph[i][0].get(), sim_additional_graph[ref_id[i]][0].get());

            std::vector<double> x_array, y_array;
            for (int p = 0; p < data_graph[i][1]->GetN(); p++)
            {
                x_array.emplace_back(data_graph[i][1]->GetPointX(p));
                y_array.emplace_back(data_graph[i][1]->GetPointY(p));
            }
            auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.20) };
            x_min[i] = std::get<0>(x_range);
            x_max[i] = std::get<1>(x_range);

            auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.70) };
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
        
        auto legend{ ROOTHelper::CreateLegend(0.05, 0.00, 0.60, 1.00, false) };
        ROOTHelper::SetLegendDefaultStyle(legend.get());
        ROOTHelper::SetTextAttribute(legend.get(), 70.0f, 133, 12);
        ROOTHelper::SetFillAttribute(legend.get(), 4000);
        legend->AddEntry(sim_with_charge_graph[0][1].get(), "A. Thomas #font[2]{et al.}", "l");
        legend->AddEntry(sim_additional_graph[0][1].get(), "Amber (W. D. Cornell #font[2]{et al.})", "l");
        legend->AddEntry(sim_no_charge_graph[0][1].get(), "Neutral", "l");
        legend->Draw();

        auto info_text1{ ROOTHelper::CreatePaveText(0.60, 0.00, 0.80, 1.00, "nbNDC", false) };
        auto info_text2{ ROOTHelper::CreatePaveText(0.80, 0.00, 1.00, 1.00, "nbNDC", false) };
        ROOTHelper::SetPaveTextDefaultStyle(info_text1.get());
        ROOTHelper::SetPaveTextDefaultStyle(info_text2.get());
        ROOTHelper::SetFillAttribute(info_text1.get(), 4000);
        ROOTHelper::SetFillAttribute(info_text2.get(), 4000);
        ROOTHelper::SetTextAttribute(info_text1.get(), 70.0f, 133, 12);
        ROOTHelper::SetTextAttribute(info_text2.get(), 70.0f, 133, 12);
        auto count{ 0 };
        for (auto & line : resolution_list)
        {
            if (count < 4) info_text1->AddText(Form("#font[102]{[#color[853]{%s}]} %.2f #AA", data_index[count], line));
            else info_text2->AddText(Form("#font[102]{[#color[853]{%s}]} %.2f #AA", data_index[count], line));
            count++;
            if (count >= 10) break;
        }
        info_text2->AddText("");
        info_text1->Draw();
        info_text2->Draw();

        std::unique_ptr<TPaveText> title_text[primary_element_size];
        for (int i = 0; i < primary_element_size; i++)
        {
            pad[i+1]->cd();
            ROOTHelper::SetPadMarginInCanvas(gPad, 0.10, 0.03, 0.08, 0.02);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 0, 0, 4000, 0);
            y_min[i] = (y_min[i] < 0.0) ? 0.0 : y_min[i];
            y_max[i] = (y_max[i] < 1.0) ? 1.0 : y_max[i];

            frame[i] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, k),"", 500, x_min[i], x_max[i], 500, y_min[i], y_max[i]);
            ROOTHelper::SetAxisTitleAttribute(frame[i]->GetXaxis(), 70.0f, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i]->GetXaxis(), 65.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i]->GetXaxis(), 0.05f, 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i]->GetYaxis(), 70.0f, 1.3f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i]->GetYaxis(), 65.0f, 0.02f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i]->GetYaxis(), 0.05f, 506);
            ROOTHelper::SetLineAttribute(frame[i].get(), 1, 0);
            frame[i]->GetXaxis()->SetTitle("");
            frame[i]->GetXaxis()->CenterTitle();
            frame[i]->GetYaxis()->CenterTitle();
            frame[i]->SetStats(0);
            frame[i]->GetXaxis()->SetTitle(Form("Width #tau_{#font[102]{#color[%d]{%s}}}", color_element[i], element_label[i]));
            frame[i]->GetYaxis()->SetTitle(Form("Amplitude Ratio #font[1]{A}_{#font[102]{#color[%d]{%s}}} / #font[1]{A}_{#font[102]{#color[%d]{%s}}}", color_element[i], element_label[i], color_element[ref_id[i]], element_label[ref_id[i]]));
            frame[i]->Draw();

            ROOTHelper::SetMarkerAttribute(data_graph[i][1].get(), marker_element[i], 3.0f, color_element[i]);
            ROOTHelper::SetLineAttribute(sim_with_charge_graph[i][1].get(), 1, 4, kRed);
            ROOTHelper::SetLineAttribute(sim_no_charge_graph[i][1].get(),   2, 4, kGray+2);
            ROOTHelper::SetLineAttribute(sim_additional_graph[i][1].get(),  3, 4, kBlue);
            sim_with_charge_graph[i][1]->Draw("L X0");
            sim_no_charge_graph[i][1]->Draw("L X0");
            sim_additional_graph[i][1]->Draw("L X0");
            data_graph[i][1]->Draw("P X0");
        }

        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ComparisonPainter::PaintGausEstimateComparisonOld(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PaintGausEstimateComparisonOld"<< std::endl;

    auto sim_with_charge_model_object_list{ m_ref_model_object_list_map.at("with_charge")};
    auto sim_no_charge_model_object_list{ m_ref_model_object_list_map.at("no_charge")};

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 4 };
    const int row_size{ 1 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 2000, 500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.06f, 0.20f, 0.20f, 0.12f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int primary_element_size{ 4 };
    const char * element_label[col_size]{ "C_{#alpha}", "C", "N", "O" };

    short color_element[primary_element_size] { kRed+1, kViolet+1, kGreen+2, kAzure+2 };
    short marker_element[primary_element_size]{ 54, 53, 55, 59 };

    std::unique_ptr<TGraphErrors> data_graph[primary_element_size][2];
    std::unique_ptr<TGraphErrors> sim_with_charge_graph[primary_element_size][2];
    std::unique_ptr<TGraphErrors> sim_no_charge_graph[primary_element_size][2];
    for (size_t i = 0; i < primary_element_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        data_graph[i][0] = ROOTHelper::CreateGraphErrors();
        sim_with_charge_graph[i][0] = ROOTHelper::CreateGraphErrors();
        sim_no_charge_graph[i][0] = ROOTHelper::CreateGraphErrors();

        BuildAmplitudeToWidthGraph(group_key, data_graph[i][0].get(), m_model_object_list);
        BuildAmplitudeToWidthGraph(group_key, sim_with_charge_graph[i][0].get(), sim_with_charge_model_object_list);
        BuildAmplitudeToWidthGraph(group_key, sim_no_charge_graph[i][0].get(), sim_no_charge_model_object_list);
    }

    const char * data_index[10]{"A","B","C","D","E","F","G","H","I","J"};
    std::vector<double> resolution_list;
    resolution_list.reserve(m_model_object_list.size());
    for (auto & model : m_model_object_list)
    {
        resolution_list.emplace_back(model->GetResolution());
    }

    int ref_id[col_size]{ 3, 3, 3, 1 };
    double x_min[col_size]{ 0.0 };
    double x_max[col_size]{ 1.0 };
    std::vector<double> y_array;
    std::vector<double> y_array_oxygen;
    for (int i = 0; i < col_size; i++)
    {
        data_graph[i][1] = ROOTHelper::CreateGraphErrors();
        sim_with_charge_graph[i][1] = ROOTHelper::CreateGraphErrors();
        sim_no_charge_graph[i][1] = ROOTHelper::CreateGraphErrors();

        BuildRatioGraph(data_graph[i][1].get(), data_graph[i][0].get(), data_graph[ref_id[i]][0].get(), true);
        BuildRatioGraph(sim_with_charge_graph[i][1].get(), sim_with_charge_graph[i][0].get(), sim_with_charge_graph[ref_id[i]][0].get());
        BuildRatioGraph(sim_no_charge_graph[i][1].get(), sim_no_charge_graph[i][0].get(), sim_no_charge_graph[ref_id[i]][0].get());

        std::vector<double> x_array;
        for (int p = 0; p < data_graph[i][1]->GetN(); p++)
        {
            x_array.emplace_back(data_graph[i][1]->GetPointX(p));
            if (i == 3) y_array_oxygen.emplace_back(data_graph[i][1]->GetPointY(p));
            else y_array.emplace_back(data_graph[i][1]->GetPointY(p));
        }
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.15) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);

        ROOTHelper::SetMarkerAttribute(data_graph[i][1].get(), marker_element[i], 1.3f, color_element[i]);
        ROOTHelper::SetMarkerAttribute(sim_with_charge_graph[i][1].get(), marker_element[i], 1.3f, color_element[i]);
        ROOTHelper::SetMarkerAttribute(sim_no_charge_graph[i][1].get(), marker_element[i], 1.3f, color_element[i]);
        ROOTHelper::SetLineAttribute(data_graph[i][1].get(), 1, 2, color_element[i]);
        ROOTHelper::SetLineAttribute(sim_with_charge_graph[i][1].get(), 1, 2, color_element[i]);
        ROOTHelper::SetLineAttribute(sim_no_charge_graph[i][1].get(), 2, 2, color_element[i]);
    }

    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.35) };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };

    auto y_range_oxygen{ ArrayStats<double>::ComputeScalingRangeTuple(y_array_oxygen, 0.35) };
    auto y_min_oxygen{ std::get<0>(y_range_oxygen) };
    auto y_max_oxygen{ std::get<1>(y_range_oxygen) };

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> title_text[col_size];
    std::unique_ptr<TPaveText> info_text;
    std::unique_ptr<TLegend> legend;
    std::unique_ptr<TLine> line_ref1;
    std::unique_ptr<TLine> line_ref2;
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 0, 0, 4000, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            if (i == 3) frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min_oxygen, y_max_oxygen);
            else frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 50.0f, 0.9f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 45.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 45.0f, 0.9f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 45.0f, 0.02f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            if (i == 3)
            {
                ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 45.0f, 1.5f, 133);
                frame[i][j]->GetYaxis()->SetTitle("Amplitude Ratio #font[1]{A}/#font[1]{A}_{ #font[102]{C}}");
                frame[i][j]->Draw("Y+");
            }
            else
            {
                frame[i][j]->GetYaxis()->SetTitle("Amplitude Ratio #font[1]{A}/#font[1]{A}_{ #font[102]{O}}");
                frame[i][j]->Draw();
            }
            sim_with_charge_graph[i][1]->Draw("L X0");
            sim_no_charge_graph[i][1]->Draw("L X0");
            data_graph[i][1]->Draw("P");
        }

        auto shift_oxygen{ (i == 3) ? 0.5 : 0.0 };
        title_text[i] = ROOTHelper::CreatePaveText(0.10+shift_oxygen, 0.96, 0.40+shift_oxygen, 1.17, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
        ROOTHelper::SetPaveAttribute(title_text[i].get(), 0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[i].get(), 60.0f, 103, 22);
        ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, color_element[i], 0.5f);
        title_text[i]->AddText(element_label[i]);
        title_text[i]->Draw();

        if (i == col_size - 1)
        {
            line_ref1 = ROOTHelper::CreateLine(0.0, 0.0, 1.0, 1.0);
            line_ref2 = ROOTHelper::CreateLine(0.0, 0.0, 1.0, 1.0);
            ROOTHelper::SetLineAttribute(line_ref1.get(), 1, 5, kGray+2);
            ROOTHelper::SetLineAttribute(line_ref2.get(), 2, 5, kGray+2);
            
            legend = ROOTHelper::CreateLegend(1.40, 0.75, 2.10, 1.17, true);
            ROOTHelper::SetLegendDefaultStyle(legend.get());
            ROOTHelper::SetTextAttribute(legend.get(), 35.0f, 133, 12);
            ROOTHelper::SetFillAttribute(legend.get(), 4000);
            legend->AddEntry(line_ref1.get(), "Partial Charge", "l");
            legend->AddEntry(line_ref2.get(), "Neutral", "l");
            legend->Draw();

            info_text = ROOTHelper::CreatePaveText(1.41, -0.28, 2.10, 0.75, "nbNDC", true);
            ROOTHelper::SetPaveTextDefaultStyle(info_text.get());
            ROOTHelper::SetFillAttribute(info_text.get(), 4000);
            ROOTHelper::SetTextAttribute(info_text.get(), 50.0f, 133, 12);
            auto count{ 0 };
            for (auto & line : resolution_list)
            {
                info_text->AddText(Form("#font[102]{[#color[922]{%s}]} %.2f #AA", data_index[count], line));
                count++;
                if (count >= 10) break;
            }
            info_text->Draw();
        }
    }
    canvas->cd();
    auto pad_extra{ ROOTHelper::CreatePad("pad_extra","", 0.10, 0.00, 0.80, 0.12) };
    pad_extra->Draw();
    pad_extra->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra.get());
    ROOTHelper::SetFillAttribute(pad_extra.get(), 4000);
    auto bottom_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(bottom_title_text.get());
    ROOTHelper::SetFillAttribute(bottom_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(bottom_title_text.get(), 50.0f, 133, 22);
    bottom_title_text->AddText("Cubic Width #tau^{3}");
    bottom_title_text->Draw();
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ComparisonPainter::PaintGausEstimateScatterComparison(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PaintGausEstimateScatterComparison"<< std::endl;

    auto sim_with_charge_model_object_list{ m_ref_model_object_list_map.at("with_charge")};
    auto sim_no_charge_model_object_list{ m_ref_model_object_list_map.at("no_charge")};

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 3 };
    const int row_size{ 1 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.06f, 0.20f, 0.20f, 0.12f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int primary_element_size{ 4 };
    const char * element_label[col_size]{ "C_{#alpha}", "C", "N" };

    short color_element[primary_element_size] { kRed+1, kViolet+1, kGreen+2, kAzure+2 };
    short marker_element[primary_element_size]{ 54, 53, 55, 59 };

    std::vector<std::unique_ptr<TGraphErrors>> data_graph_list[primary_element_size];
    std::unique_ptr<TGraphErrors> sim_with_charge_graph[primary_element_size][2];
    std::unique_ptr<TGraphErrors> sim_no_charge_graph[primary_element_size][2];
    size_t reference_id{ 3 };
    for (size_t i = 0; i < primary_element_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        sim_with_charge_graph[i][0] = ROOTHelper::CreateGraphErrors();
        sim_no_charge_graph[i][0] = ROOTHelper::CreateGraphErrors();
        BuildAmplitudeToWidthGraph(group_key, sim_with_charge_graph[i][0].get(), sim_with_charge_model_object_list);
        BuildAmplitudeToWidthGraph(group_key, sim_no_charge_graph[i][0].get(), sim_no_charge_model_object_list);

        for (auto & model : m_model_object_list)
        {
            auto graph{ ROOTHelper::CreateGraphErrors() };
            BuildAmplitudeRatioToWidthScatterGraph(i, reference_id, graph.get(), model);
            data_graph_list[i].emplace_back(std::move(graph));
        }
    }

    const char * data_index[10]{"A","B","C","D","E","F","G","H","I","J"};
    std::vector<double> resolution_list;
    resolution_list.reserve(m_model_object_list.size());
    for (auto & model : m_model_object_list)
    {
        resolution_list.emplace_back(model->GetResolution());
    }

    double x_max[col_size]{ 1.0 };
    std::vector<double> y_array;
    for (int i = 0; i < col_size; i++)
    {
        sim_with_charge_graph[i][1] = ROOTHelper::CreateGraphErrors();
        sim_no_charge_graph[i][1] = ROOTHelper::CreateGraphErrors();

        BuildRatioGraph(sim_with_charge_graph[i][1].get(), sim_with_charge_graph[i][0].get(), sim_with_charge_graph[reference_id][0].get());
        BuildRatioGraph(sim_no_charge_graph[i][1].get(), sim_no_charge_graph[i][0].get(), sim_no_charge_graph[reference_id][0].get());

        std::vector<double> x_array;
        for (auto & graph : data_graph_list[i])
        {
            for (int p = 0; p < graph->GetN(); p++)
            {
                x_array.emplace_back(graph->GetPointX(p));
                y_array.emplace_back(graph->GetPointY(p));
            }
            ROOTHelper::SetMarkerAttribute(graph.get(), marker_element[i], 1.0f, color_element[i]);
        }
        auto x_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(x_array, 0.15) };
        x_max[i] = std::get<1>(x_range);

        ROOTHelper::SetMarkerAttribute(sim_with_charge_graph[i][1].get(), marker_element[i], 1.3f, color_element[i]);
        ROOTHelper::SetMarkerAttribute(sim_no_charge_graph[i][1].get(), marker_element[i], 1.3f, color_element[i]);
        ROOTHelper::SetLineAttribute(sim_with_charge_graph[i][1].get(), 1, 2, color_element[i]);
        ROOTHelper::SetLineAttribute(sim_no_charge_graph[i][1].get(), 2, 2, color_element[i]);
    }

    auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array, 0.15) };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> title_text[col_size];
    std::unique_ptr<TPaveText> info_text;
    std::unique_ptr<TLegend> legend;
    std::unique_ptr<TLine> line_ref1;
    std::unique_ptr<TLine> line_ref2;
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 0, 0, 4000, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, 0.1, x_max[i], 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 50.0f, 0.9f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 45.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 45.0f, 0.9f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 45.0f, 0.02f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("Amplitude Ratio #font[1]{A}/#font[1]{A}_{ #font[102]{O}}");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            sim_with_charge_graph[i][1]->Draw("L X0");
            sim_no_charge_graph[i][1]->Draw("L X0");
            for (auto & graph : data_graph_list[i]) graph->Draw("P");
        }
        title_text[i] = ROOTHelper::CreatePaveText(0.10, 0.96, 0.40, 1.17, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
        ROOTHelper::SetPaveAttribute(title_text[i].get(), 0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[i].get(), 60.0f, 103, 22);
        ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, color_element[i], 0.5f);
        title_text[i]->AddText(element_label[i]);
        title_text[i]->Draw();

        if (i == col_size - 1)
        {
            line_ref1 = ROOTHelper::CreateLine(0.0, 0.0, 1.0, 1.0);
            line_ref2 = ROOTHelper::CreateLine(0.0, 0.0, 1.0, 1.0);
            ROOTHelper::SetLineAttribute(line_ref1.get(), 1, 5, kGray+2);
            ROOTHelper::SetLineAttribute(line_ref2.get(), 2, 5, kGray+2);
            
            legend = ROOTHelper::CreateLegend(0.60, 0.75, 1.80, 1.17, true);
            ROOTHelper::SetLegendDefaultStyle(legend.get());
            ROOTHelper::SetTextAttribute(legend.get(), 50.0f, 133, 12);
            ROOTHelper::SetFillAttribute(legend.get(), 4000);
            legend->AddEntry(line_ref1.get(), "Partial Charge", "l");
            legend->AddEntry(line_ref2.get(), "Neutral", "l");
            legend->Draw();

            info_text = ROOTHelper::CreatePaveText(1.01, -0.28, 1.80, 0.75, "nbNDC", true);
            ROOTHelper::SetPaveTextDefaultStyle(info_text.get());
            ROOTHelper::SetFillAttribute(info_text.get(), 4000);
            ROOTHelper::SetTextAttribute(info_text.get(), 50.0f, 133, 12);
            auto count{ 0 };
            for (auto & line : resolution_list)
            {
                info_text->AddText(Form("#font[102]{[#color[922]{%s}]} %.2f #AA", data_index[count], line));
                count++;
                if (count >= 10) break;
            }
            info_text->Draw();
        }
    }
    canvas->cd();
    auto pad_extra{ ROOTHelper::CreatePad("pad_extra","", 0.10, 0.00, 0.80, 0.12) };
    pad_extra->Draw();
    pad_extra->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra.get());
    ROOTHelper::SetFillAttribute(pad_extra.get(), 4000);
    auto bottom_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(bottom_title_text.get());
    ROOTHelper::SetFillAttribute(bottom_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(bottom_title_text.get(), 50.0f, 133, 22);
    bottom_title_text->AddText("Cubic Width #tau^{3}");
    bottom_title_text->Draw();
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ComparisonPainter::PainMapValueComparison(
    const std::string & name, ModelObject * model_data, ModelObject * model_sim)
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

    const int primary_element_size{ 4 };
    const char * element_label[col_size]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };

    short color_element[primary_element_size] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    short marker_element[primary_element_size]{ 54, 53, 55, 59 };

    std::unique_ptr<TGraphErrors> scatter_graph[col_size];
    std::unique_ptr<TF1> fit_function[col_size];
    double r_square[col_size]{ 0.0 };
    double slope[col_size]{ 0.0 };
    double intercept[col_size]{ 0.0 };
    for (size_t i = 0; i < col_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        scatter_graph[i] = ROOTHelper::CreateGraphErrors();
        BuildMapValueScatterGraph(group_key, scatter_graph[i].get(), model_sim, model_data, 15, 0.0, 1.5);
        r_square[i] = ROOTHelper::PerformLinearRegression(scatter_graph[i].get(), slope[i], intercept[i]);
        ROOTHelper::SetMarkerAttribute(scatter_graph[i].get(), marker_element[i], 1.0f, color_element[i]);
        ROOTHelper::SetLineAttribute(scatter_graph[i].get(), 1, 2, color_element[i]);

        fit_function[i] = ROOTHelper::CreateFunction1D(Form("fit_%d", static_cast<int>(i)), "x*[1]+[0]");
        fit_function[i]->SetParameters(intercept[i], slope[i]);
        ROOTHelper::SetLineAttribute(fit_function[i].get(), 1, 2, kRed);
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
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min, y_max);
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
            scatter_graph[i]->Draw("P");
            fit_function[i]->SetRange(x_min[i], x_max[i]);
            fit_function[i]->Draw("SAME");
        }
        title_text[i] = ROOTHelper::CreatePaveText(0.01, 1.01, 0.99, 1.14, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
        ROOTHelper::SetPaveAttribute(title_text[i].get(), 0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[i].get(), 40, 133, 22);
        ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, color_element[i], 0.5f);
        title_text[i]->AddText(element_label[i]);
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
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

void ComparisonPainter::PainResidueClassGausComparison(
    const std::string & name, ModelObject * model_data, ModelObject * model_sim, int par_id)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PainResidueClassGausComparison"<< std::endl;

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 4 };
    const int row_size{ 1 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1400, 500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.07f, 0.01f, 0.15f, 0.11f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int primary_element_size{ 4 };
    const char * element_label[col_size]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };

    short color_element[primary_element_size] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    short marker_element[primary_element_size]{ 54, 53, 55, 59 };

    std::unique_ptr<TGraphErrors> scatter_graph[col_size][row_size];
    std::unique_ptr<TF1> fit_function[col_size][row_size];
    std::vector<double> x_array[col_size];
    std::vector<double> y_array[row_size];
    double r_square[col_size][row_size]{ {0.0} };
    double slope[col_size][row_size]{ {0.0} };
    double intercept[col_size][row_size]{ {0.0} };
    for (size_t i = 0; i < col_size; i++)
    {
        auto group_key_list{ m_atom_classifier->GetMainChainResidueClassGroupKeyList(i) };
        for (size_t j = 0; j < row_size; j++)
        {
            scatter_graph[i][j] = ROOTHelper::CreateGraphErrors();
            BuildGausScatterGraph(group_key_list, scatter_graph[i][j].get(), model_sim, model_data, AtomicInfoHelper::GetResidueClassKey(), par_id);
            for (int p = 0; p < scatter_graph[i][j]->GetN(); p++)
            {
                x_array[i].emplace_back(scatter_graph[i][j]->GetPointX(p));
                y_array[j].emplace_back(scatter_graph[i][j]->GetPointY(p));
            }
            r_square[i][j] = ROOTHelper::PerformLinearRegression(scatter_graph[i][j].get(), slope[i][j], intercept[i][j]);
            fit_function[i][j] = ROOTHelper::CreateFunction1D(Form("fit_%d_%d", static_cast<int>(i), static_cast<int>(j)), "x*[1]+[0]");
            fit_function[i][j]->SetParameters(intercept[i][j], slope[i][j]);
        }
    }

    double x_min[col_size]{ 0.0 };
    double x_max[col_size]{ 1.0 };
    double y_min[row_size]{ 0.0 };
    double y_max[row_size]{ 1.0 };
    for (size_t i = 0; i < col_size; i++)
    {        
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array[i], 0.1) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);
    }

    for (size_t j = 0; j < row_size; j++)
    {
        auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array[j], 0.1) };
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
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 35, 0.9f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 35, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 40, 1.1f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 40, 0.02f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("Simulated Map");
            frame[i][j]->GetYaxis()->SetTitle("Real Map");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            ROOTHelper::SetMarkerAttribute(scatter_graph[i][j].get(), marker_element[i], 1.0f, color_element[i]);
            ROOTHelper::SetLineAttribute(scatter_graph[i][j].get(), 1, 2, color_element[i]);
            scatter_graph[i][j]->Draw("P X0");

            ROOTHelper::SetLineAttribute(fit_function[i][j].get(), 2, 2, kRed);
            fit_function[i][j]->SetRange(x_min[i], x_max[i]);
            fit_function[i][j]->Draw("SAME");

            r_square_text[i][j] = ROOTHelper::CreatePaveText(0.50, 0.05, 0.95, 0.18, "nbNDC ARC", true);
            ROOTHelper::SetPaveTextDefaultStyle(r_square_text[i][j].get());
            ROOTHelper::SetPaveAttribute(r_square_text[i][j].get(), 0, 0.5);
            ROOTHelper::SetLineAttribute(r_square_text[i][j].get(), 1, 0);
            ROOTHelper::SetTextAttribute(r_square_text[i][j].get(), 35, 133, 22);
            ROOTHelper::SetFillAttribute(r_square_text[i][j].get(), 1001, kAzure-7, 0.20f);
            r_square_text[i][j]->AddText(Form("R^{2} = %.2f", r_square[i][j]));
            r_square_text[i][j]->Draw();

            fit_info_text[i][j] = ROOTHelper::CreatePaveText(0.12, 0.88, 0.70, 0.99, "nbNDC", true);
            ROOTHelper::SetPaveTextDefaultStyle(fit_info_text[i][j].get());
            ROOTHelper::SetTextAttribute(fit_info_text[i][j].get(), 35, 133, 22, 0.0, kGray);
            fit_info_text[i][j]->AddText(Form("#font[1]{y} = %.2f#font[1]{x}%+.2f", slope[i][j], intercept[i][j]));
            fit_info_text[i][j]->Draw();
        }
        title_text[i] = ROOTHelper::CreatePaveText(0.01, 1.01, 0.99, 1.14, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
        ROOTHelper::SetPaveAttribute(title_text[i].get(), 0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[i].get(), 40, 133, 22);
        ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, color_element[i], 0.5f);
        title_text[i]->AddText(element_label[i]);
        title_text[i]->Draw();
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    std::cout <<"  Output file: "<< file_path << std::endl;
    #endif
}

double ComparisonPainter::CalculateErrorPropagation(
    double target_value, double reference_value, double target_error, double reference_error)
{
    double ratio{ target_value / reference_value };
    double error{ ratio * std::sqrt(std::pow(target_error/target_value, 2) + std::pow(reference_error/reference_value, 2)) };
    return error;
}

#ifdef HAVE_ROOT
void ComparisonPainter::BuildRatioGraph(
    TGraphErrors * ratio_graph,
    const TGraphErrors * target_graph, std::vector<TGraphErrors *> reference_graph_list, bool draw_index)
{
    const char * data_index[11]{"A","B","C","D","E","F","G","H","I","J","K"};
    for (int i = 0; i < target_graph->GetN(); i++)
    {
        double x_value_target, y_value_target, x_value_reference, y_value_reference;
        target_graph->GetPoint(i, x_value_target, y_value_target);
        double y_value_tmp;
        for (auto reference_graph : reference_graph_list)
        {
            reference_graph->GetPoint(i, x_value_reference, y_value_tmp);
            y_value_reference += y_value_tmp * y_value_tmp;
        }
        y_value_reference = std::sqrt(y_value_reference);
        auto ratio{ y_value_target / y_value_reference };
        ratio_graph->SetPoint(i, x_value_target, ratio);
        if (draw_index == true)
        {
            std::string label;
            if (i > 11)
            {
                std::cout <<"Warning: Data label size exceeds 12, label switch to numbers."<< std::endl;
                label = std::to_string(i);
            }
            else
            {
                label = "#splitline{" + std::string(data_index[i]) + "}{ }";
            }
            TLatex * latex = new TLatex(x_value_target, ratio, label.data());
            ROOTHelper::SetTextAttribute(latex, 60.0f, 103, 21, 0.0, kGray+2);
            ratio_graph->GetListOfFunctions()->Add(latex);
        }
    }
}

void ComparisonPainter::BuildRatioGraph(
    TGraphErrors * ratio_graph,
    const TGraphErrors * target_graph, const TGraphErrors * reference_graph, bool draw_index)
{
    if (target_graph->GetN() != reference_graph->GetN())
    {
        std::cerr << "Error: target and reference graph have different number of points." << std::endl;
        return;
    }
    const char * data_index[11]{"A","B","C","D","E","F","G","H","I","J","K"};
    for (int i = 0; i < target_graph->GetN(); i++)
    {
        double x_value_target, y_value_target, x_value_reference, y_value_reference;
        double x_error_target, y_error_target, y_error_reference;
        target_graph->GetPoint(i, x_value_target, y_value_target);
        x_error_target = target_graph->GetErrorX(i);
        y_error_target = target_graph->GetErrorY(i);
        reference_graph->GetPoint(i, x_value_reference, y_value_reference);
        y_error_reference = reference_graph->GetErrorY(i);
        auto ratio{ y_value_target / y_value_reference };
        auto error{ CalculateErrorPropagation(y_value_target, y_value_reference, y_error_target, y_error_reference) };
        ratio_graph->SetPoint(i, x_value_target, ratio);
        ratio_graph->SetPointError(i, x_error_target, error);
        if (draw_index == true)
        {
            std::string label;
            if (i > 11)
            {
                std::cout <<"Warning: Data label size exceeds 12, label switch to numbers."<< std::endl;
                label = std::to_string(i);
            }
            else
            {
                label = "#splitline{" + std::string(data_index[i]) + "}{ }";
            }
            TLatex * latex = new TLatex(x_value_target, ratio, label.data());
            ROOTHelper::SetTextAttribute(latex, 60.0f, 103, 21, 0.0, kAzure-7);
            ratio_graph->GetListOfFunctions()->Add(latex);
        }
    }
}

void ComparisonPainter::BuildGausEstimateToBlurringWidthGraph(
    uint64_t group_key, TGraphErrors * graph,
    const std::vector<ModelObject *> & model_list, int par_id)
{
    auto count{ 0 };
    for (auto model_object : model_list)
    {
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
        auto gaus_estimate{ entry_iter->GetGausEstimatePrior(group_key, AtomicInfoHelper::GetElementClassKey(), par_id) };
        auto gaus_variance{ entry_iter->GetGausVariancePrior(group_key, AtomicInfoHelper::GetElementClassKey(), par_id) };
        auto x_value{ count * 0.1 + 0.05 };
        graph->SetPoint(count, x_value, gaus_estimate);
        graph->SetPointError(count, 0.0, gaus_variance);
        count++;
    }
}

void ComparisonPainter::BuildAmplitudeToWidthGraph(
    uint64_t group_key, TGraphErrors * graph,
    const std::vector<ModelObject *> & model_list, std::string class_key)
{
    auto count{ 0 };
    for (auto model_object : model_list)
    {
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
        if (entry_iter->IsAvailableGroupKey(group_key, class_key) == false) continue;
        auto x_value{ entry_iter->GetGausEstimatePrior(group_key, class_key, 1) };
        auto y_value{ entry_iter->GetGausEstimatePrior(group_key, class_key, 0) };
        auto x_error{ entry_iter->GetGausVariancePrior(group_key, class_key, 1) };
        auto y_error{ entry_iter->GetGausVariancePrior(group_key, class_key, 0) };
        graph->SetPoint(count, x_value, y_value);
        //auto x_error_new{ std::pow(x_value, 3) * 3.0 * x_error / x_value };
        graph->SetPointError(count, x_error, y_error);
        count++;
    }
}

void ComparisonPainter::BuildAmplitudeRatioToWidthScatterGraph(
    size_t target_id, size_t reference_id, TGraphErrors * graph, ModelObject * model)
{
    auto element_size{ AtomClassifier::GetMainChainMemberCount() };
    if (target_id >= element_size || reference_id >= element_size)
    {
        std::cerr << "Error: target or reference ID exceeds the number of main chain elements." << std::endl;
        return;
    }
    std::unordered_map<int, std::tuple<double, double>> gaus_estimate_map[4];
    size_t current_id{ 0 };
    for (auto & atom : model->GetSelectedAtomList())
    {
        if (atom->GetSpecialAtomFlag() == true) continue;
        if (AtomClassifier::IsMainChainMember(atom->GetElement(), atom->GetRemoteness(), current_id) == false) continue;
        auto entry{ atom->GetAtomicPotentialEntry() };
        auto residue_id{ atom->GetResidueID() };
        auto amplitude_estimate{ entry->GetAmplitudeEstimateMDPDE() };
        auto width_estimate{ entry->GetWidthEstimateMDPDE() };
        (gaus_estimate_map[current_id])[residue_id] = std::make_tuple(amplitude_estimate, width_estimate);
    }

    auto count{ 0 };
    for (auto & [residue_id, gaus_estimate] : gaus_estimate_map[target_id])
    {
        if (gaus_estimate_map[reference_id].find(residue_id) == gaus_estimate_map[reference_id].end()) continue;
        auto target_amplitude{ std::get<0>(gaus_estimate) };
        auto target_width{ std::get<1>(gaus_estimate) };
        auto reference_amplitude{ std::get<0>(gaus_estimate_map[reference_id].at(residue_id)) };
        if (reference_amplitude <= 0.0) continue;
        auto ratio{ target_amplitude / reference_amplitude };
        graph->SetPoint(count, std::pow(target_width, 3), ratio);
        graph->SetPointError(count, 0.0, 0.0);
        count++;
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

void ComparisonPainter::BuildGausScatterGraph(
    const std::vector<uint64_t> & group_key_list, TGraphErrors * graph,
    ModelObject * model1, ModelObject * model2, const std::string & class_key, int par_id)
{
    auto entry1_iter{ std::make_unique<PotentialEntryIterator>(model1) };
    auto entry2_iter{ std::make_unique<PotentialEntryIterator>(model2) };
    auto count{ 0 };
    for (auto & group_key : group_key_list)
    {
        if (entry1_iter->IsAvailableGroupKey(group_key, class_key) == false) continue;
        if (entry2_iter->IsAvailableGroupKey(group_key, class_key) == false) continue;
        auto x_value{ entry1_iter->GetGausEstimatePrior(group_key, class_key, par_id) };
        auto y_value{ entry2_iter->GetGausEstimatePrior(group_key, class_key, par_id) };
        auto x_error{ entry1_iter->GetGausVariancePrior(group_key, class_key, par_id) };
        auto y_error{ entry2_iter->GetGausVariancePrior(group_key, class_key, par_id) };
        graph->SetPoint(count, x_value, y_value);
        graph->SetPointError(count, x_error, y_error);
        count++;
    }
    
}

#endif