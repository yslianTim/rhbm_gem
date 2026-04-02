#include <rhbm_gem/core/painter/AtomPainter.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "data/detail/ModelAnalysisAccess.hpp"
#include "PotentialPlotBuilder.hpp"
#include "detail/PainterModelAccess.hpp"
#include "detail/PainterSupport.hpp"
#include <detail/PotentialSeriesOps.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TColor.h>
#include <TMarker.h>
#include <TAxis.h>
#include <TH2.h>
#include <TF1.h>
#endif

#include <vector>
#include <tuple>

namespace rhbm_gem {

AtomPainter::AtomPainter() = default;

AtomPainter::~AtomPainter()
{

}

void AtomPainter::AddModel(ModelObject & data_object)
{
    painter_internal::PainterModelIngress::AddModel(
        *this,
        painter_internal::RequireLocalAnalyzedModel(data_object, "AtomPainter"));
}

void AtomPainter::AppendAtomObject(AtomObject & data_object)
{
    if (ModelAnalysisAccess::FindLocalEntry(data_object) == nullptr) return;
    m_atom_object_list.push_back(&data_object);
}

void AtomPainter::Painting()
{
    Logger::Log(LogLevel::Info, "AtomPainter::Painting() called.");
    Logger::Log(LogLevel::Info, "Folder path: " + m_folder_path);
    Logger::Log(LogLevel::Info, "Number of atom objects to be painted: "
                + std::to_string(m_atom_object_list.size()));

    auto label{ m_output_label.empty() ? std::string("atom") : m_output_label };
    label += ".pdf";
    
    //PaintDemoPlot("demo_plot" + label);
    PaintAtomSamplingDataSummary("atom_sampling_data_summary" + label);
}

void AtomPainter::PaintDemoPlot(const std::string & name)
{
    auto file_path{ m_folder_path + name };

    #ifdef HAVE_ROOT
    gStyle->SetGridColor(kGray);
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1200, 1500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    auto atom_object{ m_atom_object_list.at(0) };

    auto pad_main{ ROOTHelper::CreatePad("pad","", 0.00, 0.00, 1.00, 1.00) };
    pad_main->Draw();

    const auto & atom_entry{ ModelAnalysisAccess::RequireLocalEntry(*atom_object) };
    auto atom_plot_builder{ std::make_unique<PotentialPlotBuilder>(atom_object) };
    auto map_value_range{ series_ops::ComputeMapValueRange(atom_entry, 0.3) };
    auto distance_range{ series_ops::ComputeDistanceRange(atom_entry, 0.0) };

    auto map_value_graph{ atom_plot_builder->CreateDistanceToMapValueGraph() };
    auto map_value_ref_graph{ atom_plot_builder->CreateDistanceToMapValueGraph() };

    auto map_value_hist{ ROOTHelper::CreateHist2D("hist","", 15, std::get<0>(distance_range), std::get<1>(distance_range), 1000, std::get<0>(map_value_range), std::get<1>(map_value_range)) };
    for (int p = 0; p < map_value_graph->GetN(); p++)
    {
        map_value_hist->Fill(map_value_graph->GetPointX(p), map_value_graph->GetPointY(p));
    }

    pad_main->cd();
    auto frame_hist{ ROOTHelper::CreateHist2D("frame","", 100, 0.0, 1.5, 100, std::get<0>(map_value_range), std::get<1>(map_value_range)) };
    ROOTHelper::SetPadDefaultStyle(gPad);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, kWhite, 1, 0);
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.18, 0.02, 0.13, 0.02);
    ROOTHelper::SetAxisTitleAttribute(frame_hist->GetXaxis(), 60.0, 1.1f);
    ROOTHelper::SetAxisTitleAttribute(frame_hist->GetYaxis(), 60.0, 1.3f);
    ROOTHelper::SetAxisLabelAttribute(frame_hist->GetXaxis(), 60.0, 0.01f);
    ROOTHelper::SetAxisLabelAttribute(frame_hist->GetYaxis(), 60.0, 0.02f);
    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.04, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.05, 1) };
    ROOTHelper::SetAxisTickAttribute(frame_hist->GetXaxis(), static_cast<float>(x_tick_length));
    ROOTHelper::SetAxisTickAttribute(frame_hist->GetYaxis(), static_cast<float>(y_tick_length));
    frame_hist->SetStats(0);
    frame_hist->GetXaxis()->CenterTitle();
    frame_hist->GetYaxis()->CenterTitle();
    frame_hist->GetXaxis()->SetTitle("Radial Distance #font[12]{r_{j}} #[]{#AA}");
    frame_hist->GetYaxis()->SetTitle("Sampled Map Value around #font[102]{C_{#alpha}}(ALA-1) #font[12]{y_{1, j}}");
    frame_hist->Draw("");

    auto ref_line{ ROOTHelper::CreateLine(0.0, 0.0, 1.5, 0.0) };
    ROOTHelper::SetLineAttribute(ref_line.get(), 2, 1);
    ref_line->Draw();

    auto marker_size{ 1.5f };
    ROOTHelper::SetMarkerAttribute(map_value_graph.get(), 20, marker_size, kSpring-1);
    ROOTHelper::SetLineAttribute(map_value_graph.get(), 1, 1, kBlack);
    map_value_graph->Draw("P SAME");

    ROOTHelper::SetMarkerAttribute(map_value_ref_graph.get(), 24, marker_size, kBlack, 0.5);
    map_value_ref_graph->Draw("P SAME");

    map_value_hist->SetStats(0);
    map_value_hist->SetBarWidth(1.0);
    ROOTHelper::SetFillAttribute(map_value_hist.get(), 1001, kAzure-7, 0.5);
    ROOTHelper::SetLineAttribute(map_value_hist.get(), 1, 2, kAzure-5);
    map_value_hist->Draw("CANDLE2 SAME");

    const auto & estimate{ atom_entry.GetEstimateMDPDE() };
    auto amplitude{ estimate.amplitude };
    auto width{ estimate.width };
    auto gaus_func{ ROOTHelper::CreateGaus3DFunctionIn1D("gaus", amplitude, width) };
    ROOTHelper::SetLineAttribute(gaus_func.get(), 9, 4, kRed+1);
    gaus_func->Draw("SAME");

    

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);

    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    #endif
}

void AtomPainter::PaintAtomSamplingDataSummary(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, "AtomPainter::PaintAtomSamplingDataSummary");

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(2.0);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 750) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size{ 5 };

    std::unique_ptr<TPad> pad[pad_size];
    pad[0] = ROOTHelper::CreatePad("pad0","", 0.00, 0.80, 1.00, 1.00); // The title pad
    pad[1] = ROOTHelper::CreatePad("pad1","", 0.00, 0.00, 0.55, 0.80); // The left pad
    pad[2] = ROOTHelper::CreatePad("pad2","", 0.55, 0.00, 0.90, 0.60); // The bottom-middle pad
    pad[3] = ROOTHelper::CreatePad("pad3","", 0.55, 0.60, 0.90, 0.80); // The top-middle pad
    pad[4] = ROOTHelper::CreatePad("pad4","", 0.90, 0.00, 1.00, 0.60); // The bottom-right pad
    
    for (auto atom_object : m_atom_object_list)
    {
        const auto & entry_view{ ModelAnalysisAccess::RequireLocalEntry(*atom_object) };
        auto plot_builder{ std::make_unique<PotentialPlotBuilder>(atom_object) };
        auto data_graph{ plot_builder->CreateDistanceToMapValueGraph() };
        auto data_hist{ plot_builder->CreateDistanceToMapValueHistogram(15) };
        auto gaus_function_mdpde{ plot_builder->CreateAtomLocalGausFunctionMDPDE() };
        auto gaus_function_ols{ plot_builder->CreateAtomLocalGausFunctionOLS() };
        auto linear_model_mdpde{ plot_builder->CreateAtomLocalLinearModelFunctionMDPDE() };
        auto linear_model_ols{ plot_builder->CreateAtomLocalLinearModelFunctionOLS() };
        auto x_hist{ plot_builder->CreateLinearModelDataHistogram(0) };
        auto y_hist{ plot_builder->CreateLinearModelDataHistogram(1) };

        canvas->cd();
        for (int i = 0; i < pad_size; i++)
        {
            ROOTHelper::SetPadDefaultStyle(pad[i].get());
            pad[i]->Draw();
        }

        pad[0]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.01, 0.01, 0.01, 0.01);

        auto component_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
        auto component_id{ atom_object->GetComponentID() };
        auto atom_id{ atom_object->GetAtomID() };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, component_text.get(), 0.01, 0.75, 0.02, 0.02);
        ROOTHelper::SetPaveTextDefaultStyle(component_text.get());
        ROOTHelper::SetPaveAttribute(component_text.get(), 0, 0.1);
        ROOTHelper::SetFillAttribute(component_text.get(), 1001, kAzure-7);
        ROOTHelper::SetTextAttribute(component_text.get(), 75.0f, 103, 22, 0.0, kYellow-10);
        component_text->AddText((component_id + "/" + atom_id).data());
        component_text->Draw();

        auto atom_info_text{ ROOTHelper::CreatePaveText(0.26, 0.10, 0.55, 0.90, "nbNDC", false) };
        ROOTHelper::SetPaveTextDefaultStyle(atom_info_text.get());
        ROOTHelper::SetPaveAttribute(atom_info_text.get(), 0);
        ROOTHelper::SetFillAttribute(atom_info_text.get(), 4000);
        ROOTHelper::SetTextAttribute(atom_info_text.get(), 50.0f, 133);
        atom_info_text->AddText(("Serial ID: " + std::to_string(atom_object->GetSerialID())).data());
        atom_info_text->AddText(("Sequence ID: " + std::to_string(atom_object->GetSequenceID())).data());
        atom_info_text->Draw();

        auto result_text{ ROOTHelper::CreatePaveText(0.60, 0.10, 0.90, 0.90, "nbNDC", false) };
        ROOTHelper::SetPaveTextDefaultStyle(result_text.get());
        ROOTHelper::SetTextAttribute(result_text.get(), 50.0f, 133, 12, 0.0, kRed);
        ROOTHelper::SetFillAttribute(result_text.get(), 4000);
        const auto & estimate_mdpde{ entry_view.GetEstimateMDPDE() };
        auto amplitude_prior{ estimate_mdpde.amplitude };
        auto width_prior{ estimate_mdpde.width };
        result_text->AddText(Form("#font[2]{A} = %.2f", amplitude_prior));
        result_text->AddText(Form("#tau = %.2f", width_prior));
        //result_text->Draw();

        pad[1]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.10, 0.00, 0.12, 0.10);
        auto frame{ ROOTHelper::CreateHist2D("hist_0","", 100, 0.0, 1.0, 100, 0.0, 1.0) };
        ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 50.0f, 0.75f);
        ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), 50.0f, 0.005f, 133);
        ROOTHelper::SetAxisTickAttribute(frame->GetXaxis(), 0.06f, 505);
        ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 60.0f, 1.25f);
        ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 60.0f, 0.01f, 133);
        ROOTHelper::SetAxisTickAttribute(frame->GetYaxis(), 0.03f, 506);
        ROOTHelper::SetLineAttribute(frame.get(), 1, 0);
        frame->SetStats(0);
        frame->GetYaxis()->CenterTitle();
        frame->GetXaxis()->SetTitle("Radial Distance #[]{#AA}");
        frame->GetYaxis()->SetTitle("Map Value");
        auto y_range{ series_ops::ComputeMapValueRange(entry_view, 0.1) };
        auto x_min{ 0.01 };
        auto x_max{ 1.49 };
        auto y_min{ std::get<0>(y_range) };
        auto y_max{ std::get<1>(y_range) };
        frame->GetXaxis()->SetLimits(x_min, x_max);
        frame->GetYaxis()->SetLimits(y_min, y_max);
        frame->SetStats(0);
        frame->Draw();
        ROOTHelper::SetMarkerAttribute(data_graph.get(), 20, 0.8f, kAzure-7, 0.5f);
        data_graph->Draw("P");

        data_hist->SetStats(0);
        data_hist->SetBarWidth(1.0);
        ROOTHelper::SetFillAttribute(data_hist.get(), 1001, kGray, 0.3f);
        ROOTHelper::SetLineAttribute(data_hist.get(), 1, 2, kGray+2);
        data_hist->Draw("CANDLE2 SAME");
        
        ROOTHelper::SetLineAttribute(gaus_function_mdpde.get(), 2, 3, kRed);
        ROOTHelper::SetLineAttribute(gaus_function_ols.get(), 3, 3, kBlue);
        gaus_function_mdpde->Draw("SAME");
        //gaus_function_ols->Draw("SAME");

        auto legend{ ROOTHelper::CreateLegend(0.02, 0.90, 1.00, 1.00, false) };
        ROOTHelper::SetLegendDefaultStyle(legend.get());
        ROOTHelper::SetFillAttribute(legend.get(), 4000);
        ROOTHelper::SetTextAttribute(legend.get(), 40.0f, 133, 12, 0.0);
        legend->SetMargin(0.15f);
        legend->SetNColumns(2);
        legend->AddEntry(gaus_function_mdpde.get(),
            "Gaussian Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})}", "l");
        legend->AddEntry(data_graph.get(),
            "Sampling Map Value", "p");
        legend->Draw();

        auto alpha_text{ ROOTHelper::CreatePaveText(0.20, 0.10, 0.40, 0.40, "nbNDC", false) };
        ROOTHelper::SetPaveTextDefaultStyle(alpha_text.get());
        ROOTHelper::SetTextAttribute(alpha_text.get(), 50.0f, 133, 12, 0.0, kRed);
        ROOTHelper::SetFillAttribute(alpha_text.get(), 4000);
        auto alpha_r{ entry_view.GetAlphaR() };
        alpha_text->AddText(Form("#alpha_{r} = %.1f", alpha_r));
        //alpha_text->Draw();

        pad[2]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.09, 0.01, 0.12, 0.02);
        ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        auto frame_scatter{ ROOTHelper::CreateHist2D("hist_scatter","", 100, 0.0, 1.0, 100, 0.0, 1.0) };
        ROOTHelper::SetAxisTitleAttribute(frame_scatter->GetXaxis(), 40.0f, 1.1f);
        ROOTHelper::SetAxisTitleAttribute(frame_scatter->GetYaxis(), 40.0f, 1.3f);
        ROOTHelper::SetAxisLabelAttribute(frame_scatter->GetXaxis(), 40.0f, 0.01f);
        ROOTHelper::SetAxisLabelAttribute(frame_scatter->GetYaxis(), 40.0f, 0.01f);
        auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.03, 0) };
        auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.015, 1) };
        ROOTHelper::SetAxisTickAttribute(frame_scatter->GetXaxis(), static_cast<float>(x_tick_length), 505);
        ROOTHelper::SetAxisTickAttribute(frame_scatter->GetYaxis(), static_cast<float>(y_tick_length), 505);
        frame_scatter->GetXaxis()->SetTitle("x");
        frame_scatter->GetYaxis()->SetTitle("y");
        frame_scatter->GetXaxis()->SetLimits(
            x_hist->GetXaxis()->GetXmin(),
            x_hist->GetXaxis()->GetXmax()
        );
        frame_scatter->GetYaxis()->SetLimits(
            y_hist->GetXaxis()->GetXmin(),
            y_hist->GetXaxis()->GetXmax()
        );
        frame_scatter->GetXaxis()->CenterTitle();
        frame_scatter->GetYaxis()->CenterTitle();
        frame_scatter->SetStats(0);
        frame_scatter->Draw();
        auto scatter_graph{ plot_builder->CreateLinearModelDistanceToMapValueGraph() };
        ROOTHelper::SetMarkerAttribute(scatter_graph.get(), 20, 0.8f, kAzure-7, 0.5f);
        scatter_graph->Draw("P");

        ROOTHelper::SetLineAttribute(linear_model_mdpde.get(), 2, 3, kRed);
        ROOTHelper::SetLineAttribute(linear_model_ols.get(), 3, 3, kBlue);
        linear_model_mdpde->Draw("SAME");
        //linear_model_ols->Draw("SAME");

        pad[3]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.09, 0.01, 0.02, 0.02);
        ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 1, 0, 0);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        ROOTHelper::SetAxisTitleAttribute(x_hist->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisTitleAttribute(x_hist->GetYaxis(), 40.0, 1.3f);
        ROOTHelper::SetAxisLabelAttribute(x_hist->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(x_hist->GetYaxis(), 40.0, -0.01f);
        ROOTHelper::SetAxisTickAttribute(x_hist->GetXaxis(), 0.0f, 505);
        ROOTHelper::SetAxisTickAttribute(x_hist->GetYaxis(), static_cast<float>(y_tick_length), 504);
        ROOTHelper::SetFillAttribute(x_hist.get(), 1001, kAzure-7, 0.5f);
        x_hist->GetYaxis()->SetTitle("Counts");
        x_hist->GetYaxis()->CenterTitle();
        x_hist->GetYaxis()->SetRangeUser(1.0, x_hist->GetMaximum() * 1.1);
        x_hist->SetStats(0);
        x_hist->Draw("BAR");

        pad[4]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.01, 0.01, 0.12, 0.02);
        ROOTHelper::SetPadLayout(gPad, 1, 1, 1, 0, 0, 0);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        ROOTHelper::SetAxisTitleAttribute(y_hist->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisTitleAttribute(y_hist->GetYaxis(), 40.0, 1.0f);
        ROOTHelper::SetAxisLabelAttribute(y_hist->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(y_hist->GetYaxis(), 40.0, -0.02f);
        ROOTHelper::SetAxisTickAttribute(y_hist->GetXaxis(), 0.0f, 505);
        ROOTHelper::SetAxisTickAttribute(y_hist->GetYaxis(), static_cast<float>(y_tick_length), 504);
        ROOTHelper::SetFillAttribute(y_hist.get(), 1001, kAzure-7, 0.5f);
        y_hist->GetYaxis()->SetTitle("Counts");
        y_hist->GetYaxis()->CenterTitle();
        y_hist->GetYaxis()->SetRangeUser(1.0, y_hist->GetMaximum() * 1.1);
        y_hist->SetStats(0);
        y_hist->Draw("HBAR");

        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }

    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    #endif
}

} // namespace rhbm_gem
