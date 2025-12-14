#include "LocalPainter.hpp"
#include "ArrayStats.hpp"
#include "Logger.hpp"

#ifdef HAVE_ROOT
#include "ROOTHelper.hpp"
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

LocalPainter::LocalPainter(void)
{
}

LocalPainter::~LocalPainter(void)
{
}

void LocalPainter::PaintTemplate1(
    const Eigen::MatrixXd & data_matrix,
    const std::vector<double> & alpha_list,
    const std::string & x_axis_title,
    const std::string & file_path)
{
    std::vector<std::string> title_y_list{
        "Amplitude #font[2]{A}", "Width #tau"
    };

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(2);
    gStyle->SetGridColor(kGray);
    const int col_size{ 1 };
    const int row_size{ 2 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 750) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.22f, 0.15f, 0.12f, 0.02f, 0.01f, 0.01f
    );
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TGraphErrors> graph[col_size][row_size];
    std::vector<double> y_array[row_size];
    for (size_t i = 0; i < col_size; i++)
    {
        for (size_t j = 0; j < row_size; j++)
        {
            graph[i][j] = ROOTHelper::CreateGraphErrors();
            for (int p = 0; p < static_cast<int>(alpha_list.size()); p++)
            {
                auto x_value{ alpha_list.at(static_cast<size_t>(p)) };
                auto y_value{ data_matrix.col(p)(static_cast<int>(j)) };
                graph[i][j]->SetPoint(p, x_value, y_value);
                y_array[j].emplace_back(y_value);
            }
        }
    }

    double x_min[col_size]{ 0.0 };
    double x_max[col_size]{ 0.0 };
    double y_min[row_size]{ 0.0 };
    double y_max[row_size]{ 0.0 };
    for (size_t i = 0; i < col_size; i++)
    {
        x_min[i] = -0.1;
        x_max[i] = 1.1;
    }
    for (size_t j = 0; j < row_size; j++)
    {
        auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array[j], 0.38, 0.005, 0.995) };
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> resolution_text[col_size];
    std::unique_ptr<TPaveText> title_y_text[row_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto par_id{ row_size - j - 1 };
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min[par_id], y_max[par_id]);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 50.0f, 0.8f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.08/x_factor), 505);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 50.0f, 1.2f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 45.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle(x_axis_title.data());
            frame[i][j]->GetYaxis()->SetTitle("");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw("Y+");

            short color{ kRed };
            ROOTHelper::SetMarkerAttribute(graph[i][par_id].get(), 20, 1.5f, color);
            ROOTHelper::SetLineAttribute(graph[i][par_id].get(), 2, 2, color);
            ROOTHelper::SetFillAttribute(graph[i][par_id].get(), 1001, color, 0.2f);
            graph[i][par_id]->Draw("PL");

            if (i == 0)
            {
                title_y_text[j] = ROOTHelper::CreatePaveText(-0.34, 0.30, 0.00, 0.70, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(title_y_text[j].get());
                ROOTHelper::SetPaveAttribute(title_y_text[j].get(), 0, 0.1);
                ROOTHelper::SetLineAttribute(title_y_text[j].get(), 1, 0);
                ROOTHelper::SetTextAttribute(title_y_text[j].get(), 40.0f, 133, 22);
                ROOTHelper::SetFillAttribute(title_y_text[j].get(), 1001, kAzure-7, 0.5f);
                title_y_text[j]->AddText(title_y_list[static_cast<size_t>(par_id)].data());
                title_y_text[j]->Draw();
            }
        }
    }

    canvas->cd();
    auto pad_extra2{ ROOTHelper::CreatePad("pad_extra2","", 0.94, 0.10, 1.00, 0.96) };
    pad_extra2->Draw();
    pad_extra2->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra2.get());
    ROOTHelper::SetFillAttribute(pad_extra2.get(), 4000);
    auto right_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(right_title_text.get());
    ROOTHelper::SetFillAttribute(right_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(right_title_text.get(), 50.0f, 133, 22);
    right_title_text->AddText("Normalized Bias Mean #font[2]{#bar{b}}");
    auto text{ right_title_text->GetLineWith("Bias") };
    text->SetTextAngle(90.0f);
    right_title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void LocalPainter::PaintHistogram2D(
    const std::vector<double> & data_array,
    int x_bin_size, double x_min, double x_max,
    int y_bin_size, double y_min, double y_max,
    const std::string & x_axis_title,
    const std::string & y_axis_title,
    const std::string & z_axis_title,
    const std::string & file_name)
{
    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(2);
    gStyle->SetGridColor(kGray);
    gStyle->SetPalette(kLightTemperature);
    gStyle->SetNumberContours(50);
    auto canvas{ ROOTHelper::CreateCanvas("canvas","", 1200, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.15, 0.20, 0.15, 0.10);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_name);

    std::unique_ptr<TH2> hist_data{
        ROOTHelper::CreateHist2D("data","", x_bin_size, x_min, x_max, y_bin_size, y_min, y_max)
    };

    for (int i = 0; i < static_cast<int>(data_array.size()); i++)
    {
        auto x_index{ i % x_bin_size + 1 };
        auto y_index{ i / x_bin_size + 1 };
        hist_data->SetBinContent(x_index, y_index, data_array[static_cast<size_t>(i)]);
    }
    auto z_min{ ArrayStats<double>::ComputeMin(data_array.data(), data_array.size()) };
    auto z_max{ ArrayStats<double>::ComputeMax(data_array.data(), data_array.size()) };

    canvas->cd();
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    ROOTHelper::SetAxisTitleAttribute(hist_data->GetXaxis(), 60.0f, 1.0f, 133);
    ROOTHelper::SetAxisLabelAttribute(hist_data->GetXaxis(), 60.0f, 0.01f, 133);
    ROOTHelper::SetAxisTickAttribute(hist_data->GetXaxis(), 0.04f, 505);
    ROOTHelper::SetAxisTitleAttribute(hist_data->GetYaxis(), 60.0f, 1.5f, 133);
    ROOTHelper::SetAxisLabelAttribute(hist_data->GetYaxis(), 60.0f, 0.01f, 133);
    ROOTHelper::SetAxisTickAttribute(hist_data->GetYaxis(), 0.04f, 505);
    ROOTHelper::SetAxisTitleAttribute(hist_data->GetZaxis(), 60.0f, 1.3f, 133);
    ROOTHelper::SetAxisLabelAttribute(hist_data->GetZaxis(), 60.0f, 0.005f, 133);
    ROOTHelper::SetAxisTickAttribute(hist_data->GetZaxis(), 0.02f, 505);
    ROOTHelper::SetLineAttribute(hist_data.get(), 1, 0);
    hist_data->GetXaxis()->SetTitle(x_axis_title.data());
    hist_data->GetYaxis()->SetTitle(y_axis_title.data());
    hist_data->GetZaxis()->SetTitle(z_axis_title.data());
    hist_data->GetXaxis()->CenterTitle();
    hist_data->GetYaxis()->CenterTitle();
    hist_data->GetZaxis()->CenterTitle();
    hist_data->GetZaxis()->SetRangeUser(z_min, z_max);
    hist_data->SetStats(0);
    hist_data->Draw("COLZ");

    ROOTHelper::PrintCanvasPad(canvas.get(), file_name);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_name);
    Logger::Log(LogLevel::Info, " Output file: " + file_name);
    #endif
}
