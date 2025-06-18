#include "AtomPainter.hpp"
#include "AtomObject.hpp"
#include "DataObjectBase.hpp"
#include "PotentialEntryIterator.hpp"
#include "ChargeEntryIterator.hpp"
#include "FilePathHelper.hpp"
#include "AtomicInfoHelper.hpp"
#include "ArrayStats.hpp"

#ifdef HAVE_ROOT
#include "ROOTHelper.hpp"
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


AtomPainter::AtomPainter(void) :
    m_folder_path{ "./" }
{

}

AtomPainter::~AtomPainter()
{

}

void AtomPainter::SetFolder(const std::string & folder_path)
{
    m_folder_path = FilePathHelper::EnsureTrailingSlash(folder_path);
}

void AtomPainter::AddDataObject(DataObjectBase * data_object)
{
    m_atom_object_list.push_back(dynamic_cast<AtomObject *>(data_object));
}

void AtomPainter::AddReferenceDataObject(DataObjectBase * data_object, const std::string & label)
{
    m_ref_atom_object_map[label] = dynamic_cast<AtomObject *>(data_object);
}

void AtomPainter::Painting(void)
{
    std::cout <<"- AtomPainter::Painting"<<std::endl;
    std::cout <<"  Folder path: "<< m_folder_path << std::endl;
    std::cout <<"  Number of atom objects to be painted: "<< m_atom_object_list.size() << std::endl;
    
    //PaintDemoPlot("demo_plot.pdf");
    PaintRegressionCheckPlot("regression_check_plot.pdf");
}

void AtomPainter::PaintDemoPlot(const std::string & name)
{
    auto file_path{ m_folder_path + name };

    #ifdef HAVE_ROOT
    gStyle->SetGridColor(kGray);
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1200) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    auto atom_object{ m_atom_object_list.at(0) };

    auto pad_main{ ROOTHelper::CreatePad("pad","", 0.00, 0.00, 1.00, 1.00) };
    pad_main->Draw();

    auto atom_iter{ std::make_unique<PotentialEntryIterator>(atom_object) };
    auto data_array{ atom_iter->GetDistanceAndMapValueList() };
    auto map_value_range{ atom_iter->GetMapValueRange(0.3) };
    auto distance_range{ atom_iter->GetDistanceRange(0.0) };

    auto map_value_graph{ atom_iter->CreateDistanceToMapValueGraph() };
    auto map_value_ref_graph{ atom_iter->CreateDistanceToMapValueGraph() };

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

    auto amplitude{ atom_iter->GetAmplitudeEstimateMDPDE() };
    auto width{ atom_iter->GetWidthEstimateMDPDE() };
    auto gaus_func{ ROOTHelper::CreateGausFunction1D("gaus", amplitude, width) };
    ROOTHelper::SetLineAttribute(gaus_func.get(), 9, 4, kRed+1);
    gaus_func->Draw("SAME");

    

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);

    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    #endif
}

void AtomPainter::PaintRegressionCheckPlot(const std::string & name)
{
    auto file_path{ m_folder_path + name };

    const int col_size{ 2 };
    const std::string col_label[col_size]
    {
        "#font[2]{x}_{1}",
        "#font[2]{x}_{2}"
    };

    #ifdef HAVE_ROOT
    gStyle->SetGridColor(kGray);
    auto canvas{ ROOTHelper::CreateCanvas("test","", 700, 400) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, 1, 0.12f, 0.02f, 0.20f, 0.02f, 0.02f, 0.02f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TH2> frame[col_size];
    for (int i = 0; i < col_size; i++)
    {
        ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, 0);
        ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
        auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
        frame[i] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, 0),"", 100, 0.0, 1.0, 100, 0.0, 1.0);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetXaxis(), 35.0f, 0.9f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetXaxis(), 35.0f, 0.005f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 505);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetYaxis(), 40.0f, 1.0f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetYaxis(), 35.0f, 0.01f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
        ROOTHelper::SetLineAttribute(frame[i].get(), 1, 0);
        frame[i]->SetStats(0);
        frame[i]->GetXaxis()->SetTitle(col_label[i].data());
        frame[i]->GetXaxis()->CenterTitle();
        frame[i]->GetYaxis()->CenterTitle();
    }

    std::vector<std::unique_ptr<TGraphErrors>> graph_list;
    std::vector<std::unique_ptr<TF1>> function_list;
    for (auto & atom : m_atom_object_list)
    {
        auto element_label{ AtomicInfoHelper::GetLabel(atom->GetElement()) };
        auto remoteness_label{ AtomicInfoHelper::GetLabel(atom->GetRemoteness()) };
        auto label{ element_label + "_{" + remoteness_label + "}" };
        frame[0]->GetYaxis()->SetTitle(Form("#font[2]{y} (#font[102]{%s})", label.data()));
        for (int i = 0; i < col_size; i++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, 0);
            auto atom_iter{ std::make_unique<ChargeEntryIterator>(atom) };
            auto graph{ atom_iter->CreateRegressionDataGraph(i) };
            double x_min, x_max, y_min, y_max;
            graph->ComputeRange(x_min, y_min, x_max, y_max);
            auto x_range{ x_max - x_min };

            auto par1{ 0.0 };
            if      (i == 0) par1 = atom_iter->GetModelEstimateMDPDE(1);
            else if (i == 1) par1 = atom_iter->GetModelEstimateMDPDE(2);
            //else if (i == 2) par1 = atom_iter->GetModelEstimateMDPDE(2);
            auto function{ ROOTHelper::CreateFunction1D(Form("func_%d_%d", i, atom->GetSerialID()), "[0] + [1]*x") };
            function->SetRange(x_min, x_max);
            function->SetParameter(0, atom_iter->GetModelEstimateMDPDE(0)); // Intercept
            function->SetParameter(1, par1); // Slope

            frame[i]->GetXaxis()->SetLimits(x_min - 0.1*x_range, x_max + 0.1*x_range);
            frame[i]->GetYaxis()->SetLimits(y_min, y_max);
            frame[i]->Draw();

            ROOTHelper::SetLineAttribute(function.get(), 1, 2, kRed);
            //function->Draw("SAME");
            function_list.emplace_back(std::move(function));

            ROOTHelper::SetMarkerAttribute(graph.get(), 20, 1.0f, kAzure-7);
            graph->Draw("P");
            graph_list.emplace_back(std::move(graph));
        }

        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    #endif
}
