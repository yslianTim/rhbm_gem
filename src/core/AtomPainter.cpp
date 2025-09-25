#include "AtomPainter.hpp"
#include "AtomObject.hpp"
#include "DataObjectBase.hpp"
#include "PotentialEntryIterator.hpp"
#include "FilePathHelper.hpp"
#include "ArrayStats.hpp"
#include "Logger.hpp"

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
    auto atom_obj{ dynamic_cast<AtomObject *>(data_object) };
    if (atom_obj == nullptr)
    {
        throw std::runtime_error(
            "AtomPainter::AddDataObject(): invalid data_object type");
    }
    m_atom_object_list.push_back(atom_obj);
}

void AtomPainter::AddReferenceDataObject(DataObjectBase * data_object, const std::string & label)
{
    auto atom_obj{ dynamic_cast<AtomObject *>(data_object) };
    if (atom_obj == nullptr)
    {
        throw std::runtime_error(
            "AtomPainter::AddReferenceDataObject(): invalid data_object type");
    }
    m_ref_atom_object_map[label] = atom_obj;
}

void AtomPainter::Painting(void)
{
    Logger::Log(LogLevel::Info, "AtomPainter::Painting() called.");
    Logger::Log(LogLevel::Info, "Folder path: " + m_folder_path);
    Logger::Log(LogLevel::Info, "Number of atom objects to be painted: "
                + std::to_string(m_atom_object_list.size()));
    
    PaintDemoPlot("demo_plot.pdf");
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
    auto gaus_func{ ROOTHelper::CreateGaus3DFunctionIn1D("gaus", amplitude, width) };
    ROOTHelper::SetLineAttribute(gaus_func.get(), 9, 4, kRed+1);
    gaus_func->Draw("SAME");

    

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);

    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    #endif
}
