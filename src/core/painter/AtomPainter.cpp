#include <rhbm_gem/core/PainterFunctions.hpp>

#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/LocalPotentialSeries.hpp>

#include "detail/PotentialPlotBuilder.hpp"
#include "detail/PainterModelValidation.hpp"

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TColor.h>
#include <TH2.h>
#include <TF1.h>
#endif

#include <vector>

namespace rhbm_gem {

namespace {

class AtomPainter
{
    std::vector<AtomObject *> m_atom_object_list;
    std::string m_output_label;
    std::string m_folder_path{ "./" };

public:
    AtomPainter(
        const core::ModelObjectList & model_objects,
        const std::string & output_folder);
    void Run();

private:
    void AddModel(ModelObject & data_object);
    void AppendAtomObject(AtomObject & data_object);
    void PaintDemoPlot(const std::string & name);
    void PaintAtomSamplingDataSummary(const std::string & name);
    void PaintAtomSamplingDataComparison(const std::string & name);
};

AtomPainter::AtomPainter(
    const core::ModelObjectList & model_objects,
    const std::string & output_folder) :
    m_folder_path{ path_helper::EnsureTrailingSlash(output_folder) }
{
    for (auto * model_object : model_objects)
    {
        AddModel(*model_object);
    }
}

void AtomPainter::AddModel(ModelObject & data_object)
{
    painter_internal::RequireLocalAnalyzedModel(data_object, "AtomPainter");
    if (m_output_label.empty())
    {
        m_output_label = path_helper::EnsureSanitizedTag(data_object.GetKeyTag());
    }
    for (auto * atom : data_object.GetSelectedAtoms())
    {
        AppendAtomObject(*atom);
    }
}

void AtomPainter::AppendAtomObject(AtomObject & data_object)
{
    if (!AtomLocalPotentialView::For(data_object).IsAvailable()) return;
    m_atom_object_list.push_back(&data_object);
}

void AtomPainter::Run()
{
    Logger::Log(LogLevel::Info, "AtomPainter::Run() called.");
    Logger::Log(LogLevel::Info, "Folder path: " + m_folder_path);
    Logger::Log(LogLevel::Info, "Number of atom objects to be painted: "
                + std::to_string(m_atom_object_list.size()));

    auto label{ m_output_label.empty() ? std::string("atom") : m_output_label };
    label += ".pdf";
    
    //PaintDemoPlot("demo_plot" + label);
    PaintAtomSamplingDataSummary("atom_sampling_data_summary" + label);
    PaintAtomSamplingDataComparison("atom_sampling_data_comparison" + label);
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

    const auto atom_entry{ AtomLocalPotentialView::For(*atom_object) };
    auto atom_plot_builder{ std::make_unique<PotentialPlotBuilder>(atom_object) };
    auto map_value_range{
        local_potential_series::ComputeResponseRange(
            atom_entry.GetRawSamplingEntries(), 0.3)
    };
    auto distance_range{
        local_potential_series::ComputeDistanceRange(
            atom_entry.GetRawSamplingEntries(), 0.0)
    };

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

    const auto & estimate{
        atom_entry.GetJointStageModel()
    };
    auto amplitude{ estimate.GetAmplitude() };
    auto width{ estimate.GetWidth() };
    auto offset{ estimate.GetOffset() };
    auto gaus_func{ root_helper::CreateGaus3DFunctionIn1D("gaus", amplitude, width, offset) };
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
    for (auto atom_object : m_atom_object_list)
    {
        const auto entry_view{ AtomLocalPotentialView::For(*atom_object) };
        if (!entry_view.IsAvailable()) continue;
        const auto & stage = entry_view.GetJointStageEstimate();
        if (stage.source.method != EstimateMethod::JointComponents || stage.source.role != FittingRole::Target) continue;
        canvas->Clear(); canvas->cd();
        auto builder=std::make_unique<PotentialPlotBuilder>(atom_object);
        auto raw=builder->CreateDistanceToMapValueGraph(false,false);
        auto peeled=builder->CreateDistanceToMapValueGraph(false,true);
        raw->SetTitle(("Joint components: atom " + std::to_string(atom_object->GetSerialID()) + ";Radial distance (Angstrom);Map response").c_str());
        root_helper::SetMarkerAttribute(raw.get(),20,0.7f,kGray+1);
        if (raw->GetN()>0) raw->Draw("AP");
        else canvas->DrawFrame(0,-1,2,1,"Joint samples unavailable;Radial distance (Angstrom);Map response");
        root_helper::SetMarkerAttribute(peeled.get(),20,0.7f,kBlue);
        peeled->Draw("P SAME");
        auto curve=builder->CreateAtomLocalGausFunctionFinal();
        auto legend=root_helper::CreateLegend(0.55,0.68,0.95,0.94,false);
        if (raw->GetN()>0) legend->AddEntry(raw.get(),"Raw samples","p");
        if(entry_view.GetPostFitPeeling()) legend->AddEntry(peeled.get(),"Joint peeling (includes target-own charge)","p");
        if(curve) {root_helper::SetLineAttribute(curve.get(),1,2,kRed); curve->Draw("SAME"); legend->AddEntry(curve.get(),"Joint Second: Gaussian + charge C","l");}
        else legend->AddEntry(static_cast<TObject *>(nullptr),("Second unavailable: " + entry_view.GetJointStageEstimate().reason).c_str(), "");
        legend->Draw();
        root_helper::PrintCanvasPad(canvas.get(),file_path);
        canvas->Clear();
    }

    root_helper::PrintCanvasClose(canvas.get(), file_path);
    #endif
}

void AtomPainter::PaintAtomSamplingDataComparison(const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, "AtomPainter::PaintAtomSamplingDataComparison");

}

} // namespace

namespace core {

void PaintAtom(const ModelObjectList & model_objects, const std::string & output_folder)
{
    AtomPainter painter{ model_objects, output_folder };
    painter.Run();
}

} // namespace core

} // namespace rhbm_gem
