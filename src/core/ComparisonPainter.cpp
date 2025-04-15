#include "ComparisonPainter.hpp"
#include "ModelObject.hpp"
#include "DataObjectBase.hpp"
#include "PotentialEntryIterator.hpp"
#include "FilePathHelper.hpp"
#include "ArrayStats.hpp"
#include "AtomClassifier.hpp"

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
#endif

#include <vector>
#include <tuple>

using ElementKeyType = GroupKeyMapping<ElementGroupClassifierTag>::type;
using ResidueKeyType = GroupKeyMapping<ResidueGroupClassifierTag>::type;

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
    auto sim_buried_charge_model_object_list{ m_ref_model_object_list_map.at("buried_charge")};

    if (m_model_object_list.size() != 0 && m_ref_model_object_list_map.size() == 2 && sim_buried_charge_model_object_list.size() > 1)
    {
        PaintSimulationGaus("figure_3_a.pdf");
        PaintGausEstimateComparison("figure_3_b.pdf");
    }
    
    if (sim_buried_charge_model_object_list.size() > 1)
    {
        PaintSimulationGausRatio("figure_2_c.pdf", sim_buried_charge_model_object_list);
    }

    if (m_model_object_list.size() == 1 && sim_buried_charge_model_object_list.size() == 1)
    {
        PainMapValueComparison("figure_2_a.pdf", m_model_object_list.at(0), sim_buried_charge_model_object_list.at(0));
    }
}

void ComparisonPainter::PaintSimulationGaus(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PaintSimulationGaus"<< std::endl;

    auto sim_buried_charge_model_object_list{ m_ref_model_object_list_map.at("buried_charge")};
    auto sim_no_charge_model_object_list{ m_ref_model_object_list_map.at("no_charge")};

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int column_size{ 4 };
    const int row_size{ 2 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 950, 600) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), column_size, row_size, 0.1, 0.01, 0.11, 0.07, 0.01, 0.005);
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
    int color_element[primary_element_size] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };

    std::unique_ptr<TGraphErrors> sim_buried_charge_graph[primary_element_size][row_size];
    std::unique_ptr<TGraphErrors> sim_no_charge_graph[primary_element_size][row_size];
    std::vector<double> x_array, y_array[row_size];
    for (int i = 0; i < primary_element_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        for (int j = 0; j < row_size; j++)
        {
            sim_buried_charge_graph[i][j] = ROOTHelper::CreateGraphErrors();
            sim_no_charge_graph[i][j] = ROOTHelper::CreateGraphErrors();

            BuildGausEstimateToBlurringWidthGraph(group_key, sim_buried_charge_graph[i][j].get(), sim_buried_charge_model_object_list, row_size-j-1);
            BuildGausEstimateToBlurringWidthGraph(group_key, sim_no_charge_graph[i][j].get(), sim_no_charge_model_object_list, row_size-j-1);

            for (int p = 0; p < sim_buried_charge_graph[i][j]->GetN(); p++)
            {
                x_array.push_back(sim_buried_charge_graph[i][j]->GetPointX(p));
                y_array[j].push_back(sim_buried_charge_graph[i][j]->GetPointY(p));
            }
            for (int p = 0; p < sim_no_charge_graph[i][j]->GetN(); p++)
            {
                x_array.push_back(sim_no_charge_graph[i][j]->GetPointX(p));
                y_array[j].push_back(sim_no_charge_graph[i][j]->GetPointY(p));
            }

            ROOTHelper::SetMarkerAttribute(sim_buried_charge_graph[i][j].get(), 20, 1.2, color_element[i]);
            ROOTHelper::SetMarkerAttribute(sim_no_charge_graph[i][j].get(), 53, 1.2, color_element[i]);
            ROOTHelper::SetLineAttribute(sim_buried_charge_graph[i][j].get(), 1, 2, color_element[i]);
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

    std::unique_ptr<TH2> frame[column_size][row_size];
    std::unique_ptr<TPaveText> title_text[column_size];
    std::unique_ptr<TLegend> legend;
    for (int i = 0; i < column_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min, x_max, 500, y_min[j], y_max[j]);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 30, 1.0, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 30, 0.005, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), y_factor*0.05/x_factor, 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 35, 1.3, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35, 0.01, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), x_factor*0.05/y_factor, 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("Blurring Width");
            frame[i][j]->GetYaxis()->SetTitle(y_axis_title[j].data());
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            sim_buried_charge_graph[i][j]->Draw("PL X0");
            sim_no_charge_graph[i][j]->Draw("PL X0");
            if (i == column_size - 1 && j == row_size - 1)
            {
                legend = ROOTHelper::CreateLegend(0.1, 0.7, 0.9, 0.99, true);
                ROOTHelper::SetLegendDefaultStyle(legend.get());
                ROOTHelper::SetTextAttribute(legend.get(), 20, 133, 12);
                ROOTHelper::SetFillAttribute(legend.get(), 1001, kWhite, 0.5);
                legend->AddEntry(sim_buried_charge_graph[i][1].get(), "Partial Charge", "pl");
                legend->AddEntry(sim_no_charge_graph[i][1].get(), "No Charge", "pl");
                legend->Draw();
            }
        }
        title_text[i] = ROOTHelper::CreatePaveText(0.01, 1.01, 0.99, 1.16, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
        ROOTHelper::SetPaveAttribute(title_text[i].get(), 0.0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[i].get(), 25, 133, 22);
        ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, color_element[i], 0.5);
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
    const int column_size{ 2 };
    const int row_size{ 1 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 600) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), column_size, row_size, 0.11, 0.11, 0.15, 0.05, 0.01, 0.01);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    
    const int primary_element_size{ 4 };
    const std::string y_axis_title[column_size]
    {
        "Amplitude Ratio A_{C#alpha} / A_{C}",
        "Width Ratio #tau_{C#alpha} / #tau_{C}"
    };

    std::unique_ptr<TGraphErrors> simulation_graph[primary_element_size][2];
    for (int i = 0; i < primary_element_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        for (int j = 0; j < 2; j++)
        {
            simulation_graph[i][j] = ROOTHelper::CreateGraphErrors();
            BuildGausEstimateToBlurringWidthGraph(group_key, simulation_graph[i][j].get(), model_list, j);
        }
    }

    std::unique_ptr<TGraphErrors> ratio_graph[column_size];
    std::vector<double> x_array, y_array;
    for (int i = 0; i < column_size; i++)
    {
        ratio_graph[i] = ROOTHelper::CreateGraphErrors();
        BuildRatioGraph(ratio_graph[i].get(), simulation_graph[0][i].get(), simulation_graph[1][i].get());

        for (int p = 0; p < ratio_graph[i]->GetN(); p++)
        {
            x_array.push_back(ratio_graph[i]->GetPointX(p));
            y_array.push_back(ratio_graph[i]->GetPointY(p));
        }
        ROOTHelper::SetMarkerAttribute(ratio_graph[i].get(), 20, 1.2, kRed+1);
        ROOTHelper::SetLineAttribute(ratio_graph[i].get(), 1, 2, kRed+1);
    }

    auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.01) };
    auto x_min{ std::get<0>(x_range) };
    auto x_max{ std::get<1>(x_range) };
    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.3) };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };

    std::unique_ptr<TH2> frame[column_size][row_size];
    for (int i = 0; i < column_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min, x_max, 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 40, 1.0, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40, 0.005, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), y_factor*0.05/x_factor, 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 40, 1.3, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 40, 0.01, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), x_factor*0.05/y_factor, 506);
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

void ComparisonPainter::PaintGausEstimateComparison(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    std::cout <<"- ComparisonPainter::PaintGausEstimateComparison"<< std::endl;

    auto sim_buried_charge_model_object_list{ m_ref_model_object_list_map.at("buried_charge")};
    auto sim_no_charge_model_object_list{ m_ref_model_object_list_map.at("no_charge")};

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int column_size{ 3 };
    const int row_size{ 1 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), column_size, row_size, 0.1, 0.01, 0.15, 0.1, 0.01, 0.01);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int primary_element_size{ 4 };
    const char * element_label[column_size]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen"
    };

    int color_element[primary_element_size] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    int marker_element[primary_element_size]{ 54, 53, 55, 59 };

    std::unique_ptr<TGraphErrors> data_graph[primary_element_size][2];
    std::unique_ptr<TGraphErrors> sim_buried_charge_graph[primary_element_size][2];
    std::unique_ptr<TGraphErrors> sim_no_charge_graph[primary_element_size][2];
    for (int i = 0; i < primary_element_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        data_graph[i][0] = ROOTHelper::CreateGraphErrors();
        sim_buried_charge_graph[i][0] = ROOTHelper::CreateGraphErrors();
        sim_no_charge_graph[i][0] = ROOTHelper::CreateGraphErrors();

        BuildAmplitudeRatioToWidthGraph(group_key, data_graph[i][0].get(), m_model_object_list);
        BuildAmplitudeRatioToWidthGraph(group_key, sim_buried_charge_graph[i][0].get(), sim_buried_charge_model_object_list);
        BuildAmplitudeRatioToWidthGraph(group_key, sim_no_charge_graph[i][0].get(), sim_no_charge_model_object_list);
    }

    int ref_id{ 3 };
    double x_min[column_size]{ 0.0 };
    double x_max[column_size]{ 1.0 };
    std::vector<double> y_array;
    for (int i = 0; i < column_size; i++)
    {
        data_graph[i][1] = ROOTHelper::CreateGraphErrors();
        sim_buried_charge_graph[i][1] = ROOTHelper::CreateGraphErrors();
        sim_no_charge_graph[i][1] = ROOTHelper::CreateGraphErrors();

        BuildRatioGraph(data_graph[i][1].get(), data_graph[i][0].get(), data_graph[ref_id][0].get());
        BuildRatioGraph(sim_buried_charge_graph[i][1].get(), sim_buried_charge_graph[i][0].get(), sim_buried_charge_graph[ref_id][0].get());
        BuildRatioGraph(sim_no_charge_graph[i][1].get(), sim_no_charge_graph[i][0].get(), sim_no_charge_graph[ref_id][0].get());

        std::vector<double> x_array;
        for (int p = 0; p < data_graph[i][1]->GetN(); p++)
        {
            x_array.emplace_back(data_graph[i][1]->GetPointX(p));
            y_array.emplace_back(data_graph[i][1]->GetPointY(p));
        }
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.15) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);

        ROOTHelper::SetMarkerAttribute(data_graph[i][1].get(), marker_element[i], 1.3, color_element[i]);
        ROOTHelper::SetMarkerAttribute(sim_buried_charge_graph[i][1].get(), marker_element[i], 1.3, color_element[i]);
        ROOTHelper::SetMarkerAttribute(sim_no_charge_graph[i][1].get(), marker_element[i], 1.3, color_element[i]);
        ROOTHelper::SetLineAttribute(data_graph[i][1].get(), 1, 2, color_element[i]);
        ROOTHelper::SetLineAttribute(sim_buried_charge_graph[i][1].get(), 1, 2, color_element[i]);
        ROOTHelper::SetLineAttribute(sim_no_charge_graph[i][1].get(), 2, 2, color_element[i]);
    }

    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.35) };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };

    std::unique_ptr<TH2> frame[column_size][row_size];
    std::unique_ptr<TPaveText> title_text[column_size];
    std::unique_ptr<TLegend> legend;
    for (int i = 0; i < column_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 35, 0.9, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 35, 0.005, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), y_factor*0.05/x_factor, 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 40, 1.1, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 40, 0.02, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), x_factor*0.05/y_factor, 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("Width");
            frame[i][j]->GetYaxis()->SetTitle("Amplitude Ratio #font[1]{A} / #font[1]{A}_{O}");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            sim_buried_charge_graph[i][1]->Draw("L X0");
            sim_no_charge_graph[i][1]->Draw("L X0");
            data_graph[i][1]->Draw("P");
        }
        title_text[i] = ROOTHelper::CreatePaveText(0.01, 1.01, 0.99, 1.13, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
        ROOTHelper::SetPaveAttribute(title_text[i].get(), 0.0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[i].get(), 40, 133, 22);
        ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, color_element[i], 0.5);
        title_text[i]->AddText(element_label[i]);
        title_text[i]->Draw();

        if (i == column_size - 1)
        {
            legend = ROOTHelper::CreateLegend(0.1, 0.7, 0.9, 0.99, true);
            ROOTHelper::SetLegendDefaultStyle(legend.get());
            ROOTHelper::SetTextAttribute(legend.get(), 30, 133, 12);
            ROOTHelper::SetFillAttribute(legend.get(), 1001, kWhite, 0.5);
            legend->AddEntry(sim_buried_charge_graph[i][1].get(), "Partial Charge", "l");
            legend->AddEntry(sim_no_charge_graph[i][1].get(), "No Charge", "l");
            legend->Draw();
        }
    }
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
    const int column_size{ 4 };
    const int row_size{ 1 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1400, 500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), column_size, row_size, 0.07, 0.01, 0.15, 0.11, 0.01, 0.01);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int primary_element_size{ 4 };
    const char * element_label[column_size]
    {
        "Alpha Carbon",
        "Carbonyl Carbon",
        "Peptide Nitrogen",
        "Carbonyl Oxygen"
    };

    int color_element[primary_element_size] { kRed+1, kOrange+1, kGreen+2, kAzure+2 };
    int marker_element[primary_element_size]{ 54, 53, 55, 59 };

    std::unique_ptr<TGraphErrors> scatter_graph[column_size];
    std::unique_ptr<TF1> fit_function[column_size];
    double r_square[column_size]{ 0.0 };
    double slope[column_size]{ 0.0 };
    double intercept[column_size]{ 0.0 };
    for (int i = 0; i < column_size; i++)
    {
        auto group_key{ m_atom_classifier->GetMainChainElementClassGroupKey(i) };
        scatter_graph[i] = ROOTHelper::CreateGraphErrors();
        BuildMapValueScatterGraph(group_key, scatter_graph[i].get(), model_sim, model_data, 15, 0.0, 1.5);
        r_square[i] = ROOTHelper::PerformLinearRegression(scatter_graph[i].get(), slope[i], intercept[i]);
        ROOTHelper::SetMarkerAttribute(scatter_graph[i].get(), marker_element[i], 1.0, color_element[i]);
        ROOTHelper::SetLineAttribute(scatter_graph[i].get(), 1, 2, color_element[i]);

        fit_function[i] = ROOTHelper::CreateFunction1D(Form("fit_%d", i), "x*[1]+[0]");
        fit_function[i]->SetParameters(intercept[i], slope[i]);
        ROOTHelper::SetLineAttribute(fit_function[i].get(), 1, 2, kRed);
    }

    double x_min[column_size]{ 0.0 };
    double x_max[column_size]{ 1.0 };
    std::vector<double> y_array;
    for (int i = 0; i < column_size; i++)
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

    std::unique_ptr<TH2> frame[column_size][row_size];
    std::unique_ptr<TPaveText> title_text[column_size];
    std::unique_ptr<TPaveText> r_square_text[column_size];
    std::unique_ptr<TPaveText> fit_info_text[column_size];
    for (int i = 0; i < column_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 35, 0.9, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 35, 0.005, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), y_factor*0.05/x_factor, 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 40, 1.1, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 40, 0.02, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), x_factor*0.05/y_factor, 506);
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
        ROOTHelper::SetPaveAttribute(title_text[i].get(), 0.0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[i].get(), 40, 133, 22);
        ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, color_element[i], 0.5);
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
    const TGraphErrors * target_graph, const TGraphErrors * reference_graph)
{
    if (target_graph->GetN() != reference_graph->GetN())
    {
        std::cerr << "Error: target and reference graph have different number of points." << std::endl;
        return;
    }
    for (int i = 0; i < target_graph->GetN(); i++)
    {
        double x_value_target, y_value_target, x_value_reference, y_value_reference;
        double x_error_target, y_error_target, x_error_reference, y_error_reference;
        target_graph->GetPoint(i, x_value_target, y_value_target);
        x_error_target = target_graph->GetErrorX(i);
        y_error_target = target_graph->GetErrorY(i);
        reference_graph->GetPoint(i, x_value_reference, y_value_reference);
        x_error_reference = reference_graph->GetErrorX(i);
        y_error_reference = reference_graph->GetErrorY(i);
        auto ratio{ y_value_target / y_value_reference };
        auto error{ CalculateErrorPropagation(y_value_target, y_value_reference, y_error_target, y_error_reference) };
        ratio_graph->SetPoint(i, x_value_target, ratio);
        ratio_graph->SetPointError(i, x_error_target, error);
    }
}

void ComparisonPainter::BuildGausEstimateToBlurringWidthGraph(
    ElementKeyType & group_key, TGraphErrors * graph,
    const std::vector<ModelObject *> & model_list, int par_id)
{
    auto count{ 0 };
    for (auto model_object : model_list)
    {
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
        if (entry_iter->IsAvailableGroupKey(group_key) == false) continue;
        auto gaus_estimate{ entry_iter->GetGausEstimatePrior(group_key, par_id) };
        auto gaus_variance{ entry_iter->GetGausVariancePrior(group_key, par_id) };
        auto x_value{ count * 0.1 + 0.05 };
        graph->SetPoint(count, x_value, gaus_estimate);
        graph->SetPointError(count, 0.0, gaus_variance);
        count++;
    }
}

void ComparisonPainter::BuildAmplitudeRatioToWidthGraph(
    ElementKeyType & group_key, TGraphErrors * graph, const std::vector<ModelObject *> & model_list)
{
    auto count{ 0 };
    for (auto model_object : model_list)
    {
        auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
        if (entry_iter->IsAvailableGroupKey(group_key) == false) continue;
        auto gaus_estimate{ entry_iter->GetGausEstimatePrior(group_key) };
        auto gaus_variance{ entry_iter->GetGausVariancePrior(group_key) };
        graph->SetPoint(count, std::get<1>(gaus_estimate), std::get<0>(gaus_estimate));
        graph->SetPointError(count, std::get<1>(gaus_variance), std::get<0>(gaus_variance));
        count++;
    }
}

void ComparisonPainter::BuildMapValueScatterGraph(
    ElementKeyType & group_key, TGraphErrors * graph, ModelObject * model1, ModelObject * model2,
    int bin_size, double x_min, double x_max)
{
    auto entry1_iter{ std::make_unique<PotentialEntryIterator>(model1) };
    auto entry2_iter{ std::make_unique<PotentialEntryIterator>(model2) };
    if (entry1_iter->IsAvailableGroupKey(group_key) == false) return;
    if (entry2_iter->IsAvailableGroupKey(group_key) == false) return;
    auto model1_atom_map{ entry1_iter->GetAtomObjectMap(group_key) };
    auto model2_atom_map{ entry2_iter->GetAtomObjectMap(group_key) };
    auto count{ 0 };
    for (auto & [atom_id, atom_object1] : model1_atom_map)
    {
        if (model2_atom_map.find(atom_id) == model2_atom_map.end()) continue;
        auto atom_object2{ model2_atom_map.at(atom_id) };
        auto atom1_iter{ std::make_unique<PotentialEntryIterator>(atom_object1) };
        auto atom2_iter{ std::make_unique<PotentialEntryIterator>(atom_object2) };
        auto data1_array{ atom1_iter->GetBinnedDistanceAndMapValueList(bin_size, x_min, x_max) };
        auto data2_array{ atom2_iter->GetBinnedDistanceAndMapValueList(bin_size, x_min, x_max) };
        for (int i = 0; i < bin_size; i++)
        {
            auto x_value{ std::get<1>(data1_array.at(i)) };
            auto y_value{ std::get<1>(data2_array.at(i)) };
            graph->SetPoint(count, x_value, y_value);
            count++;
        }
    }
}

#endif