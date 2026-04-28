#include <rhbm_gem/core/painter/AtomPainter.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include "PotentialPlotBuilder.hpp"
#include "detail/PainterModelAccess.hpp"
#include "detail/PainterSupport.hpp"
#include <rhbm_gem/utils/math/ArrayHelper.hpp>
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
    if (!LocalPotentialView::For(data_object).IsAvailable()) return;
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
    auto canvas{ root_helper::CreateCanvas("test","", 1200, 1500) };
    root_helper::SetCanvasDefaultStyle(canvas.get());
    root_helper::PrintCanvasOpen(canvas.get(), file_path);
    auto atom_object{ m_atom_object_list.at(0) };

    auto pad_main{ root_helper::CreatePad("pad","", 0.00, 0.00, 1.00, 1.00) };
    pad_main->Draw();

    const auto atom_entry{ LocalPotentialView::RequireFor(*atom_object) };
    auto atom_plot_builder{ std::make_unique<PotentialPlotBuilder>(atom_object) };
    auto map_value_range{ atom_entry.GetResponseRange(0.3) };
    auto distance_range{ atom_entry.GetDistanceRange(0.0) };

    auto map_value_graph{ atom_plot_builder->CreateDistanceToMapValueGraph() };
    auto map_value_ref_graph{ atom_plot_builder->CreateDistanceToMapValueGraph() };

    auto map_value_hist{ root_helper::CreateHist2D("hist","", 15, std::get<0>(distance_range), std::get<1>(distance_range), 1000, std::get<0>(map_value_range), std::get<1>(map_value_range)) };
    for (int p = 0; p < map_value_graph->GetN(); p++)
    {
        map_value_hist->Fill(map_value_graph->GetPointX(p), map_value_graph->GetPointY(p));
    }

    pad_main->cd();
    auto frame_hist{ root_helper::CreateHist2D("frame","", 100, 0.0, 1.5, 100, std::get<0>(map_value_range), std::get<1>(map_value_range)) };
    root_helper::SetPadDefaultStyle(gPad);
    root_helper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
    root_helper::SetPadFrameAttribute(gPad, 0, 0, 4000, kWhite, 1, 0);
    root_helper::SetPadMarginInCanvas(gPad, 0.18, 0.02, 0.13, 0.02);
    root_helper::SetAxisTitleAttribute(frame_hist->GetXaxis(), 60.0, 1.1f);
    root_helper::SetAxisTitleAttribute(frame_hist->GetYaxis(), 60.0, 1.3f);
    root_helper::SetAxisLabelAttribute(frame_hist->GetXaxis(), 60.0, 0.01f);
    root_helper::SetAxisLabelAttribute(frame_hist->GetYaxis(), 60.0, 0.02f);
    auto x_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.04, 0) };
    auto y_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.05, 1) };
    root_helper::SetAxisTickAttribute(frame_hist->GetXaxis(), static_cast<float>(x_tick_length));
    root_helper::SetAxisTickAttribute(frame_hist->GetYaxis(), static_cast<float>(y_tick_length));
    frame_hist->SetStats(0);
    frame_hist->GetXaxis()->CenterTitle();
    frame_hist->GetYaxis()->CenterTitle();
    frame_hist->GetXaxis()->SetTitle("Radial Distance #font[12]{r_{j}} #[]{#AA}");
    frame_hist->GetYaxis()->SetTitle("Sampled Map Value around #font[102]{C_{#alpha}}(ALA-1) #font[12]{y_{1, j}}");
    frame_hist->Draw("");

    auto ref_line{ root_helper::CreateLine(0.0, 0.0, 1.5, 0.0) };
    root_helper::SetLineAttribute(ref_line.get(), 2, 1);
    ref_line->Draw();

    auto marker_size{ 1.5f };
    root_helper::SetMarkerAttribute(map_value_graph.get(), 20, marker_size, kSpring-1);
    root_helper::SetLineAttribute(map_value_graph.get(), 1, 1, kBlack);
    map_value_graph->Draw("P SAME");

    root_helper::SetMarkerAttribute(map_value_ref_graph.get(), 24, marker_size, kBlack, 0.5);
    map_value_ref_graph->Draw("P SAME");

    map_value_hist->SetStats(0);
    map_value_hist->SetBarWidth(1.0);
    root_helper::SetFillAttribute(map_value_hist.get(), 1001, kAzure-7, 0.5);
    root_helper::SetLineAttribute(map_value_hist.get(), 1, 2, kAzure-5);
    map_value_hist->Draw("CANDLE2 SAME");

    const auto & estimate{ atom_entry.GetEstimateMDPDE() };
    auto amplitude{ estimate.GetAmplitude() };
    auto width{ estimate.GetWidth() };
    auto gaus_func{ root_helper::CreateGaus3DFunctionIn1D("gaus", amplitude, width) };
    root_helper::SetLineAttribute(gaus_func.get(), 9, 4, kRed+1);
    gaus_func->Draw("SAME");

    

    root_helper::PrintCanvasPad(canvas.get(), file_path);

    root_helper::PrintCanvasClose(canvas.get(), file_path);
    #endif
}

void AtomPainter::PaintAtomSamplingDataSummary(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, "AtomPainter::PaintAtomSamplingDataSummary");

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(2.0);
    gStyle->SetGridColor(kGray);

    auto canvas{ root_helper::CreateCanvas("test","", 1500, 750) };
    root_helper::SetCanvasDefaultStyle(canvas.get());
    root_helper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size{ 5 };

    std::unique_ptr<TPad> pad[pad_size];
    pad[0] = root_helper::CreatePad("pad0","", 0.00, 0.80, 1.00, 1.00); // The title pad
    pad[1] = root_helper::CreatePad("pad1","", 0.00, 0.00, 0.55, 0.80); // The left pad
    pad[2] = root_helper::CreatePad("pad2","", 0.55, 0.00, 0.90, 0.60); // The bottom-middle pad
    pad[3] = root_helper::CreatePad("pad3","", 0.55, 0.60, 0.90, 0.80); // The top-middle pad
    pad[4] = root_helper::CreatePad("pad4","", 0.90, 0.00, 1.00, 0.60); // The bottom-right pad
    
    for (auto atom_object : m_atom_object_list)
    {
        const auto entry_view{ LocalPotentialView::RequireFor(*atom_object) };
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
            root_helper::SetPadDefaultStyle(pad[i].get());
            pad[i]->Draw();
        }

        pad[0]->cd();
        root_helper::SetPadMarginInCanvas(gPad, 0.01, 0.01, 0.01, 0.01);

        auto component_text{ root_helper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
        auto component_id{ atom_object->GetComponentID() };
        auto atom_id{ atom_object->GetAtomID() };
        root_helper::SetPaveTextMarginInCanvas(gPad, component_text.get(), 0.01, 0.75, 0.02, 0.02);
        root_helper::SetPaveTextDefaultStyle(component_text.get());
        root_helper::SetPaveAttribute(component_text.get(), 0, 0.1);
        root_helper::SetFillAttribute(component_text.get(), 1001, kAzure-7);
        root_helper::SetTextAttribute(component_text.get(), 75.0f, 103, 22, 0.0, kYellow-10);
        component_text->AddText((component_id + "/" + atom_id).data());
        component_text->Draw();

        auto atom_info_text{ root_helper::CreatePaveText(0.26, 0.10, 0.55, 0.90, "nbNDC", false) };
        root_helper::SetPaveTextDefaultStyle(atom_info_text.get());
        root_helper::SetPaveAttribute(atom_info_text.get(), 0);
        root_helper::SetFillAttribute(atom_info_text.get(), 4000);
        root_helper::SetTextAttribute(atom_info_text.get(), 50.0f, 133);
        atom_info_text->AddText(("Serial ID: " + std::to_string(atom_object->GetSerialID())).data());
        atom_info_text->AddText(("Sequence ID: " + std::to_string(atom_object->GetSequenceID())).data());
        atom_info_text->Draw();

        auto result_text{ root_helper::CreatePaveText(0.60, 0.10, 0.90, 0.90, "nbNDC", false) };
        root_helper::SetPaveTextDefaultStyle(result_text.get());
        root_helper::SetTextAttribute(result_text.get(), 50.0f, 133, 12, 0.0, kRed);
        root_helper::SetFillAttribute(result_text.get(), 4000);
        const auto & estimate_mdpde{ entry_view.GetEstimateMDPDE() };
        auto amplitude_prior{ estimate_mdpde.GetAmplitude() };
        auto width_prior{ estimate_mdpde.GetWidth() };
        result_text->AddText(Form("#font[2]{A} = %.2f", amplitude_prior));
        result_text->AddText(Form("#tau = %.2f", width_prior));
        //result_text->Draw();

        pad[1]->cd();
        root_helper::SetPadMarginInCanvas(gPad, 0.10, 0.00, 0.12, 0.10);
        auto frame{ root_helper::CreateHist2D("hist_0","", 100, 0.0, 1.0, 100, 0.0, 1.0) };
        root_helper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        root_helper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        root_helper::SetAxisTitleAttribute(frame->GetXaxis(), 50.0f, 0.75f);
        root_helper::SetAxisLabelAttribute(frame->GetXaxis(), 50.0f, 0.005f, 133);
        root_helper::SetAxisTickAttribute(frame->GetXaxis(), 0.06f, 505);
        root_helper::SetAxisTitleAttribute(frame->GetYaxis(), 60.0f, 1.25f);
        root_helper::SetAxisLabelAttribute(frame->GetYaxis(), 60.0f, 0.01f, 133);
        root_helper::SetAxisTickAttribute(frame->GetYaxis(), 0.03f, 506);
        root_helper::SetLineAttribute(frame.get(), 1, 0);
        frame->SetStats(0);
        frame->GetYaxis()->CenterTitle();
        frame->GetXaxis()->SetTitle("Radial Distance #[]{#AA}");
        frame->GetYaxis()->SetTitle("Map Value");
        auto y_range{ entry_view.GetResponseRange(0.1) };
        auto x_min{ 0.01 };
        auto x_max{ 1.49 };
        auto y_min{ std::get<0>(y_range) };
        auto y_max{ std::get<1>(y_range) };
        frame->GetXaxis()->SetLimits(x_min, x_max);
        frame->GetYaxis()->SetLimits(y_min, y_max);
        frame->SetStats(0);
        frame->Draw();
        root_helper::SetMarkerAttribute(data_graph.get(), 20, 0.8f, kAzure-7, 0.5f);
        data_graph->Draw("P");

        data_hist->SetStats(0);
        data_hist->SetBarWidth(1.0);
        root_helper::SetFillAttribute(data_hist.get(), 1001, kGray, 0.3f);
        root_helper::SetLineAttribute(data_hist.get(), 1, 2, kGray+2);
        data_hist->Draw("CANDLE2 SAME");
        
        root_helper::SetLineAttribute(gaus_function_mdpde.get(), 2, 3, kRed);
        root_helper::SetLineAttribute(gaus_function_ols.get(), 3, 3, kBlue);
        gaus_function_mdpde->Draw("SAME");
        //gaus_function_ols->Draw("SAME");

        auto legend{ root_helper::CreateLegend(0.02, 0.90, 1.00, 1.00, false) };
        root_helper::SetLegendDefaultStyle(legend.get());
        root_helper::SetFillAttribute(legend.get(), 4000);
        root_helper::SetTextAttribute(legend.get(), 40.0f, 133, 12, 0.0);
        legend->SetMargin(0.15f);
        legend->SetNColumns(2);
        legend->AddEntry(gaus_function_mdpde.get(),
            "Gaussian Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})}", "l");
        legend->AddEntry(data_graph.get(),
            "Sampling Map Value", "p");
        legend->Draw();

        auto alpha_text{ root_helper::CreatePaveText(0.20, 0.10, 0.40, 0.40, "nbNDC", false) };
        root_helper::SetPaveTextDefaultStyle(alpha_text.get());
        root_helper::SetTextAttribute(alpha_text.get(), 50.0f, 133, 12, 0.0, kRed);
        root_helper::SetFillAttribute(alpha_text.get(), 4000);
        auto alpha_r{ entry_view.GetAlphaR() };
        alpha_text->AddText(Form("#alpha_{r} = %.1f", alpha_r));
        //alpha_text->Draw();

        pad[2]->cd();
        root_helper::SetPadMarginInCanvas(gPad, 0.09, 0.01, 0.12, 0.02);
        root_helper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        root_helper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        auto frame_scatter{ root_helper::CreateHist2D("hist_scatter","", 100, 0.0, 1.0, 100, 0.0, 1.0) };
        root_helper::SetAxisTitleAttribute(frame_scatter->GetXaxis(), 40.0f, 1.1f);
        root_helper::SetAxisTitleAttribute(frame_scatter->GetYaxis(), 40.0f, 1.3f);
        root_helper::SetAxisLabelAttribute(frame_scatter->GetXaxis(), 40.0f, 0.01f);
        root_helper::SetAxisLabelAttribute(frame_scatter->GetYaxis(), 40.0f, 0.01f);
        auto x_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.03, 0) };
        auto y_tick_length{ root_helper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.015, 1) };
        root_helper::SetAxisTickAttribute(frame_scatter->GetXaxis(), static_cast<float>(x_tick_length), 505);
        root_helper::SetAxisTickAttribute(frame_scatter->GetYaxis(), static_cast<float>(y_tick_length), 505);
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
        root_helper::SetMarkerAttribute(scatter_graph.get(), 20, 0.8f, kAzure-7, 0.5f);
        scatter_graph->Draw("P");

        root_helper::SetLineAttribute(linear_model_mdpde.get(), 2, 3, kRed);
        root_helper::SetLineAttribute(linear_model_ols.get(), 3, 3, kBlue);
        linear_model_mdpde->Draw("SAME");
        //linear_model_ols->Draw("SAME");

        pad[3]->cd();
        root_helper::SetPadMarginInCanvas(gPad, 0.09, 0.01, 0.02, 0.02);
        root_helper::SetPadLayout(gPad, 1, 1, 0, 1, 0, 0);
        root_helper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        root_helper::SetAxisTitleAttribute(x_hist->GetXaxis(), 0.0f);
        root_helper::SetAxisTitleAttribute(x_hist->GetYaxis(), 40.0, 1.3f);
        root_helper::SetAxisLabelAttribute(x_hist->GetXaxis(), 0.0f);
        root_helper::SetAxisLabelAttribute(x_hist->GetYaxis(), 40.0, -0.01f);
        root_helper::SetAxisTickAttribute(x_hist->GetXaxis(), 0.0f, 505);
        root_helper::SetAxisTickAttribute(x_hist->GetYaxis(), static_cast<float>(y_tick_length), 504);
        root_helper::SetFillAttribute(x_hist.get(), 1001, kAzure-7, 0.5f);
        x_hist->GetYaxis()->SetTitle("Counts");
        x_hist->GetYaxis()->CenterTitle();
        x_hist->GetYaxis()->SetRangeUser(1.0, x_hist->GetMaximum() * 1.1);
        x_hist->SetStats(0);
        x_hist->Draw("BAR");

        pad[4]->cd();
        root_helper::SetPadMarginInCanvas(gPad, 0.01, 0.01, 0.12, 0.02);
        root_helper::SetPadLayout(gPad, 1, 1, 1, 0, 0, 0);
        root_helper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        root_helper::SetAxisTitleAttribute(y_hist->GetXaxis(), 0.0f);
        root_helper::SetAxisTitleAttribute(y_hist->GetYaxis(), 40.0, 1.0f);
        root_helper::SetAxisLabelAttribute(y_hist->GetXaxis(), 0.0f);
        root_helper::SetAxisLabelAttribute(y_hist->GetYaxis(), 40.0, -0.02f);
        root_helper::SetAxisTickAttribute(y_hist->GetXaxis(), 0.0f, 505);
        root_helper::SetAxisTickAttribute(y_hist->GetYaxis(), static_cast<float>(y_tick_length), 504);
        root_helper::SetFillAttribute(y_hist.get(), 1001, kAzure-7, 0.5f);
        y_hist->GetYaxis()->SetTitle("Counts");
        y_hist->GetYaxis()->CenterTitle();
        y_hist->GetYaxis()->SetRangeUser(1.0, y_hist->GetMaximum() * 1.1);
        y_hist->SetStats(0);
        y_hist->Draw("HBAR");

        root_helper::PrintCanvasPad(canvas.get(), file_path);
    }

    root_helper::PrintCanvasClose(canvas.get(), file_path);
    #endif
}

} // namespace rhbm_gem
