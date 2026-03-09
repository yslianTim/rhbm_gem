#include <rhbm_gem/core/painter/ModelPainter.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/PotentialEntryIterator.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/ComponentHelper.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/BondClassifier.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>
#include <rhbm_gem/utils/domain/GlobalEnumClass.hpp>
#include <rhbm_gem/utils/domain/AtomKeySystem.hpp>
#include <rhbm_gem/utils/domain/StringHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include "internal/PainterIngestionInternal.hpp"

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TColor.h>
#include <TMarker.h>
#include <TAxis.h>
#include <TH2.h>
#include <TF1.h>
#include <TLine.h>
#include <TH3F.h>
#include <TH1.h>
#endif

#include <vector>
#include <tuple>

namespace rhbm_gem {
void ModelPainter::PaintGroupWidthScatterPlot(
    ModelObject * model_object, const std::string & name, int par_id, bool draw_box_plot)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintGroupWidthScatterPlot");
    auto residue_class{ ChemicalDataHelper::GetComponentAtomClassKey() };

    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    const int col_size{ 4 };
    const int row_size{ 1 };

    auto canvas{ ROOTHelper::CreateCanvas("test","", 2000, 500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.15f, 0.02f, 0.12f, 0.02f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    
    std::vector<std::unique_ptr<TGraphErrors>> graph_list[col_size][row_size];
    std::vector<double> x_array[col_size][row_size];
    std::vector<double> y_array[col_size][row_size];
    std::vector<double> global_x_array[col_size];
    std::vector<double> global_y_array;
    for (size_t i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            for (auto residue : ChemicalDataHelper::GetStandardAminoAcidList())
            {
                auto group_key{ m_atom_classifier->GetMainChainComponentAtomClassGroupKey(i, residue) };
                if (entry_iter->IsAvailableAtomGroupKey(group_key, residue_class) == false) continue;
                auto gaus_graph
                {
                    (par_id == 0) ?
                    entry_iter->CreateCOMDistanceToGausEstimateGraph(group_key, residue_class, 1) :
                    entry_iter->CreateInRangeAtomsToGausEstimateGraph(group_key, residue_class, 5.0, 1)
                };
                for (int p = 0; p < gaus_graph->GetN(); p++)
                {
                    x_array[i][j].emplace_back(gaus_graph->GetPointX(p));
                    y_array[i][j].emplace_back(gaus_graph->GetPointY(p));
                    global_x_array[i].emplace_back(gaus_graph->GetPointX(p));
                    global_y_array.emplace_back(gaus_graph->GetPointY(p));
                }
                graph_list[i][j].emplace_back(std::move(gaus_graph));
            }
        }
    }

    double x_min[col_size];
    double x_max[col_size];
    for (int i = 0; i < col_size; i++)
    {
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(global_x_array[i], 0.05) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);
    }

    auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(global_y_array, 0.1) };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };

    std::unique_ptr<TH2D> summary_hist[col_size][row_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto x_range_tmp{ ArrayStats<double>::ComputeScalingRangeTuple(x_array[i][j], 0.0) };
            auto y_range_tmp{ ArrayStats<double>::ComputeScalingRangeTuple(y_array[i][j], 0.0) };
            summary_hist[i][j] = ROOTHelper::CreateHist2D(
                Form("hist_%d_%d", i, j), "",
                5, std::get<0>(x_range_tmp), std::get<1>(x_range_tmp),
                100, std::get<0>(y_range_tmp), std::get<1>(y_range_tmp));
            for (auto & graph : graph_list[i][j])
            {
                for (int p = 0; p < graph->GetN(); p++)
                {
                    summary_hist[i][j]->Fill(graph->GetPointX(p), graph->GetPointY(p));
                }
            }
        }
    }

    std::unique_ptr<TH2> frame[col_size][row_size];
    for (size_t i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto element_color{ AtomClassifier::GetMainChainElementColor(i) };
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), static_cast<int>(i), j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", static_cast<int>(i), j),"", 500, x_min[i], x_max[i], 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 45.0f, 1.0f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 45.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 55.0f, 1.4f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 45.0f, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            auto title{ (par_id == 0)? "Distance to C.M. #[]{#AA}" : "In Range Atoms" };
            frame[i][j]->GetXaxis()->SetTitle(title);
            frame[i][j]->GetYaxis()->SetTitle("Width");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            if (draw_box_plot == true)
            {
                ROOTHelper::SetLineAttribute(summary_hist[i][j].get(), 1, 1, element_color);
                ROOTHelper::SetFillAttribute(summary_hist[i][j].get(), 1001, element_color, 0.3f);
                summary_hist[i][j]->Draw("CANDLE3 SAME");
            }
            else
            {
                for (auto & graph : graph_list[i][j])
                {
                    ROOTHelper::SetMarkerAttribute(graph.get(), 24, 1.0f, element_color);
                    graph->Draw("P");
                }
            }
        }
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintAtomXYPosition(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomXYPosition");
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    auto normalized_z_position{ 0.5 };
    auto z_ratio_window{ 0.1 };
    auto graph{ entry_iter->CreateAtomXYPositionTomographyGraph(
        normalized_z_position, z_ratio_window, true)
    };
    auto z_position{ model_object->GetModelPosition(2, normalized_z_position) };
    auto z_window{ model_object->GetModelLength(2) * z_ratio_window };
    auto com_pos{ model_object->GetCenterOfMassPosition() };
    std::vector<double> x_array, y_array;
    x_array.reserve(static_cast<size_t>(graph->GetN()));
    y_array.reserve(static_cast<size_t>(graph->GetN()));
    for (int p = 0; p < graph->GetN(); p++)
    {
        x_array.emplace_back(graph->GetPointX(p));
        y_array.emplace_back(graph->GetPointY(p));
    }

    auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.1) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.1) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };

    canvas->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.15, 0.05, 0.12, 0.12);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);

    auto frame{ ROOTHelper::CreateHist2D("frame","", 500, x_min, x_max, 500, y_min, y_max) };
    ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 50.0f, 1.1f);
    ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 50.0f, 1.4f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), 45.0f, 0.01f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 45.0f, 0.01f);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("Atom Position #font[2]{x} #[]{#AA}");
    frame->GetYaxis()->SetTitle("Atom Position #font[2]{y} #[]{#AA}");
    frame->Draw();
    ROOTHelper::SetMarkerAttribute(graph.get(), 6, 1.0f, kAzure-7);
    graph->Draw("P X0");

    auto pos_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextMarginInCanvas(gPad, pos_text.get(), 0.45, 0.06, 0.83, 0.12);
    ROOTHelper::SetPaveTextDefaultStyle(pos_text.get());
    ROOTHelper::SetFillAttribute(pos_text.get(), 4000);
    ROOTHelper::SetTextAttribute(pos_text.get(), 35, 133, 32);
    pos_text->AddText(
        Form("Atom Position #font[2]{z} = %.1f #pm %.1f #[]{#AA}",
        z_position - com_pos.at(2), z_window*0.5));
    pos_text->Draw();

    auto title_text{ ROOTHelper::CreatePaveText(0.02, 0.89, 0.67, 0.98, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(title_text.get());
    ROOTHelper::SetPaveAttribute(title_text.get(), 0, 0.2);
    ROOTHelper::SetLineAttribute(title_text.get(), 1, 0);
    ROOTHelper::SetFillAttribute(title_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(title_text.get(), 47, 23, 22, 0.0, kYellow-10);
    title_text->AddText("Tomography of Atomic Model");
    title_text->Draw();

    auto model_text{ ROOTHelper::CreatePaveText(0.69, 0.89, 0.95, 0.98, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(model_text.get());
    ROOTHelper::SetPaveAttribute(model_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(model_text.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetTextAttribute(model_text.get(), 50, 103, 22);
    model_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
    model_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintAtomGausScatterPlot(
    ModelObject * model_object, const std::string & name, bool do_normalize)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomGausScatterPlot");
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    auto amplitude_min{ entry_iter->GetAtomGausEstimateMinimum(0, Element::OXYGEN) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1300, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TH2> frame;
    std::unordered_map<Element, std::unique_ptr<TGraphErrors>> graph_map;
    std::vector<double> x_array, y_array;
    x_array.reserve(model_object->GetNumberOfSelectedAtom());
    y_array.reserve(model_object->GetNumberOfSelectedAtom());
    for (auto & [element_type, element_name] : ChemicalDataHelper::GetElementLabelMap())
    {
        auto graph
        {
            (do_normalize == true) ?
            entry_iter->CreateNormalizedAtomGausEstimateScatterGraph(element_type, amplitude_min) :
            entry_iter->CreateAtomGausEstimateScatterGraph(element_type)
        };
        auto atomic_number{ ChemicalDataHelper::GetAtomicNumber(element_type) };
        auto marker_size{ (atomic_number <= 8) ? 1.2f : 2.0f };
        short marker_color{ (ChemicalDataHelper::IsStandardElement(element_type)) ?
            ChemicalDataHelper::GetDisplayColor(element_type) : static_cast<short>(kRed)
        };
        auto transparency{ (atomic_number <= 8) ? 0.2f : 1.0f };
        ROOTHelper::SetMarkerAttribute(graph.get(),
            ChemicalDataHelper::GetDisplayMarker(element_type), marker_size,
            marker_color, transparency);
        for (int p = 0; p < graph->GetN(); p++)
        {
            auto x{ graph->GetPointX(p) };
            auto y{ graph->GetPointY(p) };
            x_array.emplace_back(x);
            y_array.emplace_back(y);
        }
        if (graph->GetN() == 0) continue;
        graph_map.emplace(element_type, std::move(graph));
    }

    auto x_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(x_array, 0.2, 0.005, 0.995) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array, 0.2, 0.005, 0.995) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };

    canvas->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.12, 0.20, 0.12, 0.10);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);

    frame = ROOTHelper::CreateHist2D("frame","", 500, x_min, x_max, 500, y_min, y_max);
    ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 50.0f, 1.1f);
    ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 50.0f, 1.4f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), 45.0f, 0.01f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 45.0f, 0.01f);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("Amplitude #font[2]{A}");
    frame->GetYaxis()->SetTitle("Width #font[2]{#tau}");
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->Draw();
    auto legend{ ROOTHelper::CreateLegend(0.81, 0.12, 0.98, 0.90, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 40.0f, 133, 12);
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    for (int i = 1; i <= 90; i++)
    {
        auto element{ ChemicalDataHelper::GetElementFromAtomicNumber(i) };
        if (graph_map.find(element) == graph_map.end()) continue;
        auto & graph{ graph_map.at(element) };
        graph->Draw("P X0");
        legend->AddEntry(graph.get(),
                         Form("#font[102]{%s} (%d)",
                            ChemicalDataHelper::GetElementLabelMap().at(element).data(),
                            graph->GetN()), "p");
    }
    legend->Draw();

    auto title_text{ ROOTHelper::CreatePaveText(0.12, 0.91, 0.98, 0.99, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(title_text.get());
    ROOTHelper::SetPaveAttribute(title_text.get(), 0, 0.2);
    ROOTHelper::SetLineAttribute(title_text.get(), 1, 0);
    ROOTHelper::SetFillAttribute(title_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(title_text.get(), 50.0f, 23, 22, 0.0, kYellow-10);
    title_text->AddText("Atom-Based Main/Side-Chain #font[2]{A}#minus#font[2]{#tau} Scatter Plot");
    title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintBondGausScatterPlot(
    ModelObject * model_object, const std::string & name, bool do_normalize)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintBondGausScatterPlot");
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    auto amplitude_min{ entry_iter->GetBondGausEstimateMinimum(0) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1300, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TH2> frame;
    std::unordered_map<Element, std::unique_ptr<TGraphErrors>> graph_map;
    std::vector<double> x_array, y_array;
    x_array.reserve(model_object->GetNumberOfSelectedBond());
    y_array.reserve(model_object->GetNumberOfSelectedBond());
    for (auto & [element_type, element_name] : ChemicalDataHelper::GetElementLabelMap())
    {
        auto graph
        {
            (do_normalize == true) ?
            entry_iter->CreateNormalizedBondGausEstimateScatterGraph(element_type, amplitude_min) :
            entry_iter->CreateBondGausEstimateScatterGraph(element_type)
        };
        auto marker_size{ 2.0f };
        auto marker_color{ (ChemicalDataHelper::IsStandardElement(element_type)) ?
            ChemicalDataHelper::GetDisplayColor(element_type) : static_cast<short>(kRed)
        };
        ROOTHelper::SetMarkerAttribute(graph.get(),
            ChemicalDataHelper::GetDisplayMarker(element_type), marker_size, marker_color);
        for (int p = 0; p < graph->GetN(); p++)
        {
            auto x{ graph->GetPointX(p) };
            auto y{ graph->GetPointY(p) };
            x_array.emplace_back(x);
            y_array.emplace_back(y);
        }
        if (graph->GetN() == 0) continue;
        graph_map.emplace(element_type, std::move(graph));
    }

    auto x_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(x_array, 0.2, 0.005, 0.995) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array, 0.2, 0.005, 0.995) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };

    canvas->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.12, 0.20, 0.12, 0.10);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0);

    frame = ROOTHelper::CreateHist2D("frame","", 500, x_min, x_max, 500, y_min, y_max);
    ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 50.0f, 1.1f);
    ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 50.0f, 1.4f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), 45.0f, 0.01f);
    ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 45.0f, 0.01f);
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle("Amplitude #font[2]{A}");
    frame->GetYaxis()->SetTitle("Width #font[2]{#tau}");
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->Draw();
    auto legend{ ROOTHelper::CreateLegend(0.81, 0.12, 0.98, 0.90, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 40.0f, 133, 12);
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    for (int i = 1; i <= 90; i++)
    {
        auto element{ ChemicalDataHelper::GetElementFromAtomicNumber(i) };
        if (graph_map.find(element) == graph_map.end()) continue;
        auto & graph{ graph_map.at(element) };
        graph->Draw("P X0");
        legend->AddEntry(graph.get(),
                         Form("#font[102]{%s} (%d)",
                            ChemicalDataHelper::GetElementLabelMap().at(element).data(),
                            graph->GetN()), "p");
    }
    legend->Draw();

    auto title_text{ ROOTHelper::CreatePaveText(0.12, 0.91, 0.98, 0.99, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(title_text.get());
    ROOTHelper::SetPaveAttribute(title_text.get(), 0, 0.2);
    ROOTHelper::SetLineAttribute(title_text.get(), 1, 0);
    ROOTHelper::SetFillAttribute(title_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(title_text.get(), 50.0f, 23, 22, 0.0, kYellow-10);
    title_text->AddText("Bond-Based Main/Side-Chain #font[2]{A}#minus#font[2]{#tau} Scatter Plot");
    title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintAtomGausMainChain(ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomGausMainChain");
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 1 };
    const int row_size{ 3 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 900) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.08f, 0.02f, 0.10f, 0.10f, 0.02f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int main_chain_element_count{ 4 };
    const std::string y_title[row_size]
    {
        "Intensity",
        "Width",
        "Amplitude"
    };
    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> gaus_graph_map[row_size][main_chain_element_count];
    for (size_t k = 0; k < main_chain_element_count; k++)
    {
        gaus_graph_map[2][k] = entry_iter->CreateAtomGausEstimateToSequenceIDGraphMap(k, 0);
        gaus_graph_map[1][k] = entry_iter->CreateAtomGausEstimateToSequenceIDGraphMap(k, 1);
        gaus_graph_map[0][k] = entry_iter->CreateAtomGausEstimateToSequenceIDGraphMap(k, 2);
    }

    for (auto & [chain_id, gaus_graph] : gaus_graph_map[0][0])
    {
        std::vector<double> x_array;
        std::vector<double> y_array[row_size];
        double y_min[row_size], y_max[row_size];
        for (int j = 0; j < row_size; j++)
        {
            y_array[j].reserve(static_cast<size_t>(gaus_graph->GetN() * 4));
            for (size_t k = 0; k < main_chain_element_count; k++)
            {
                if (gaus_graph_map[j][k].find(chain_id) == gaus_graph_map[j][k].end()) continue;
                for (int p = 0; p < gaus_graph_map[j][k].at(chain_id)->GetN(); p++)
                {
                    if (j == 0) x_array.emplace_back(gaus_graph_map[j][k].at(chain_id)->GetPointX(p));
                    y_array[j].emplace_back(gaus_graph_map[j][k].at(chain_id)->GetPointY(p));
                }
            }
            auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array[j], 0.2) };
            y_min[j] = std::get<0>(y_range);
            y_max[j] = std::get<1>(y_range);
        }
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.05) };
        auto x_min{ std::get<0>(x_range) };
        auto x_max{ std::get<1>(x_range) };

        std::unique_ptr<TPaveText> subtitle1_text;
        std::unique_ptr<TPaveText> subtitle2_text;
        std::unique_ptr<TPaveText> subtitle3_text;
        std::unique_ptr<TLegend> legend;
        for (int i = 0; i < col_size; i++)
        {
            for (int j = 0; j < row_size; j++)
            {
                ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
                ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
                ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0, 0);
                auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
                auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
                if (frame[i][j] == nullptr)
                {
                    frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 500, 0.0, 1.0, 500, 0.0, 1.0);
                    ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 45.0f, 0.9f, 133);
                    ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.01f, 133);
                    ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 510);
                    ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 45.0f, 1.3f, 133);
                    ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 40.0f, 0.005f, 133);
                    ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.02/y_factor), 506);
                    ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
                    frame[i][j]->GetXaxis()->CenterTitle();
                    frame[i][j]->GetYaxis()->CenterTitle();
                    frame[i][j]->SetStats(0);
                }
                frame[i][j]->GetXaxis()->SetTitle(Form("Residue ID #[]{Chain %s}", chain_id.data()));
                frame[i][j]->GetYaxis()->SetTitle(y_title[j].data());
                frame[i][j]->GetXaxis()->SetLimits(x_min, x_max);
                frame[i][j]->GetYaxis()->SetLimits(y_min[j], y_max[j]);
                frame[i][j]->Draw();
                for (size_t k = 0; k < main_chain_element_count; k++)
                {
                    if (gaus_graph_map[j][k].find(chain_id) == gaus_graph_map[j][k].end()) continue;
                    auto element_color{ AtomClassifier::GetMainChainElementColor(static_cast<size_t>(k)) };
                    ROOTHelper::SetMarkerAttribute(gaus_graph_map[j][k].at(chain_id).get(), 20, 0.7f, element_color);
                    ROOTHelper::SetLineAttribute(gaus_graph_map[j][k].at(chain_id).get(), 1, 1, element_color);
                    gaus_graph_map[j][k].at(chain_id)->Draw("PL X0");
                }

                if (i == 0 && j == row_size - 1)
                {
                    subtitle1_text = ROOTHelper::CreatePaveText(0.00, 1.02, 0.15, 1.37, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle1_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle1_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle1_text.get(), 1001, kAzure-7);
                    ROOTHelper::SetTextAttribute(subtitle1_text.get(), 60.0f, 133, 22, 0.0, kYellow-10);
                    subtitle1_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
                    subtitle1_text->Draw();

                    subtitle2_text = ROOTHelper::CreatePaveText(0.16, 1.02, 0.35, 1.37, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle2_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle2_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle2_text.get(), 1001, kAzure-7, 0.5f);
                    ROOTHelper::SetTextAttribute(subtitle2_text.get(), 50.0f, 103, 22);
                    subtitle2_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
                    subtitle2_text->Draw();

                    subtitle3_text = ROOTHelper::CreatePaveText(0.36, 1.02, 0.57, 1.37, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle3_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle3_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle3_text.get(), 1001, kAzure-7, 0.5f);
                    ROOTHelper::SetTextAttribute(subtitle3_text.get(), 50.0f, 103, 22);
                    subtitle3_text->AddText(model_object->GetEmdID().data());
                    subtitle3_text->Draw();

                    legend = ROOTHelper::CreateLegend(0.58, 1.02, 0.98, 1.37, true);
                    ROOTHelper::SetLegendDefaultStyle(legend.get());
                    ROOTHelper::SetTextAttribute(legend.get(), 70.0f, 133, 12);
                    ROOTHelper::SetFillAttribute(legend.get(), 4000);
                    for (size_t k = 0; k < main_chain_element_count; k++)
                    {
                        if (gaus_graph_map[j][k].find(chain_id) == gaus_graph_map[j][k].end()) continue;
                        auto label{ "#font[102]{" + AtomClassifier::GetMainChainElementLabel(k) +"}" };
                        auto color{ AtomClassifier::GetMainChainElementColor(k) };
                        auto result{ Form("#color[%d]{%s}", color, label.data()) };
                        legend->AddEntry(gaus_graph_map[j][k].at(chain_id).get(), result, "pl");
                    }
                    legend->SetNColumns(4);
                    legend->Draw();
                }
            }
        }
        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintBondGausMainChain(ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintBondGausMainChain");
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 1 };
    const int row_size{ 3 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 900) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.08f, 0.02f, 0.10f, 0.10f, 0.02f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int main_chain_element_count{ 4 };
    const std::string y_title[row_size]
    {
        "Intensity",
        "Width",
        "Amplitude"
    };
    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> gaus_graph_map[row_size][main_chain_element_count];
    for (size_t k = 0; k < main_chain_element_count; k++)
    {
        gaus_graph_map[2][k] = entry_iter->CreateBondGausEstimateToSequenceIDGraphMap(k, 0);
        gaus_graph_map[1][k] = entry_iter->CreateBondGausEstimateToSequenceIDGraphMap(k, 1);
        gaus_graph_map[0][k] = entry_iter->CreateBondGausEstimateToSequenceIDGraphMap(k, 2);
    }

    for (auto & [chain_id, gaus_graph] : gaus_graph_map[0][0])
    {
        std::vector<double> x_array;
        std::vector<double> y_array[row_size];
        double y_min[row_size], y_max[row_size];
        for (int j = 0; j < row_size; j++)
        {
            y_array[j].reserve(static_cast<size_t>(gaus_graph->GetN() * 4));
            for (size_t k = 0; k < main_chain_element_count; k++)
            {
                if (gaus_graph_map[j][k].find(chain_id) == gaus_graph_map[j][k].end()) continue;
                for (int p = 0; p < gaus_graph_map[j][k].at(chain_id)->GetN(); p++)
                {
                    if (j == 0) x_array.emplace_back(gaus_graph_map[j][k].at(chain_id)->GetPointX(p));
                    y_array[j].emplace_back(gaus_graph_map[j][k].at(chain_id)->GetPointY(p));
                }
            }
            auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array[j], 0.2) };
            y_min[j] = std::get<0>(y_range);
            y_max[j] = std::get<1>(y_range);
        }
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.05) };
        auto x_min{ std::get<0>(x_range) };
        auto x_max{ std::get<1>(x_range) };

        std::unique_ptr<TPaveText> subtitle1_text;
        std::unique_ptr<TPaveText> subtitle2_text;
        std::unique_ptr<TPaveText> subtitle3_text;
        std::unique_ptr<TLegend> legend;
        for (int i = 0; i < col_size; i++)
        {
            for (int j = 0; j < row_size; j++)
            {
                ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
                ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
                ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0, 0);
                auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
                auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
                if (frame[i][j] == nullptr)
                {
                    frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 500, 0.0, 1.0, 500, 0.0, 1.0);
                    ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 45.0f, 0.9f, 133);
                    ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.01f, 133);
                    ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 510);
                    ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 45.0f, 1.3f, 133);
                    ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 40.0f, 0.005f, 133);
                    ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.02/y_factor), 506);
                    ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
                    frame[i][j]->GetXaxis()->CenterTitle();
                    frame[i][j]->GetYaxis()->CenterTitle();
                    frame[i][j]->SetStats(0);
                }
                frame[i][j]->GetXaxis()->SetTitle(Form("Residue ID #[]{Chain %s}", chain_id.data()));
                frame[i][j]->GetYaxis()->SetTitle(y_title[j].data());
                frame[i][j]->GetXaxis()->SetLimits(x_min, x_max);
                frame[i][j]->GetYaxis()->SetLimits(y_min[j], y_max[j]);
                frame[i][j]->Draw();
                for (size_t k = 0; k < main_chain_element_count; k++)
                {
                    if (gaus_graph_map[j][k].find(chain_id) == gaus_graph_map[j][k].end()) continue;
                    auto element_color{ BondClassifier::GetMainChainMemberColor(static_cast<size_t>(k)) };
                    ROOTHelper::SetMarkerAttribute(gaus_graph_map[j][k].at(chain_id).get(), 20, 0.7f, element_color);
                    ROOTHelper::SetLineAttribute(gaus_graph_map[j][k].at(chain_id).get(), 1, 1, element_color);
                    gaus_graph_map[j][k].at(chain_id)->Draw("PL X0");
                }

                if (i == 0 && j == row_size - 1)
                {
                    subtitle1_text = ROOTHelper::CreatePaveText(0.00, 1.02, 0.15, 1.37, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle1_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle1_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle1_text.get(), 1001, kAzure-7);
                    ROOTHelper::SetTextAttribute(subtitle1_text.get(), 60.0f, 133, 22, 0.0, kYellow-10);
                    subtitle1_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
                    subtitle1_text->Draw();

                    subtitle2_text = ROOTHelper::CreatePaveText(0.16, 1.02, 0.35, 1.37, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle2_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle2_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle2_text.get(), 1001, kAzure-7, 0.5f);
                    ROOTHelper::SetTextAttribute(subtitle2_text.get(), 50.0f, 103, 22);
                    subtitle2_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
                    subtitle2_text->Draw();

                    subtitle3_text = ROOTHelper::CreatePaveText(0.36, 1.02, 0.57, 1.37, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle3_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle3_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle3_text.get(), 1001, kAzure-7, 0.5f);
                    ROOTHelper::SetTextAttribute(subtitle3_text.get(), 50.0f, 103, 22);
                    subtitle3_text->AddText(model_object->GetEmdID().data());
                    subtitle3_text->Draw();

                    legend = ROOTHelper::CreateLegend(0.58, 1.02, 0.98, 1.37, true);
                    ROOTHelper::SetLegendDefaultStyle(legend.get());
                    ROOTHelper::SetTextAttribute(legend.get(), 70.0f, 133, 12);
                    ROOTHelper::SetFillAttribute(legend.get(), 4000);
                    for (size_t k = 0; k < main_chain_element_count; k++)
                    {
                        if (gaus_graph_map[j][k].find(chain_id) == gaus_graph_map[j][k].end()) continue;
                        auto label{ "#font[102]{" + BondClassifier::GetMainChainMemberLabel(k) +"}" };
                        auto color{ BondClassifier::GetMainChainMemberColor(k) };
                        auto result{ Form("#color[%d]{%s}", color, label.data()) };
                        legend->AddEntry(gaus_graph_map[j][k].at(chain_id).get(), result, "pl");
                    }
                    legend->SetNColumns(4);
                    legend->Draw();
                }
            }
        }
        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintAtomRankMainChain(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomRankMainChain");
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 4 };
    const int row_size{ 3 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 800) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.08f, 0.05f, 0.10f, 0.10f, 0.02f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const std::string y_title[row_size]
    {
        "Intensity",
        "Width",
        "Amplitude"
    };

    int chain_size;
    std::unique_ptr<TH2> frame[col_size][row_size];
    std::vector<std::unique_ptr<TH1D>> rank_hist_list[row_size][3];
    std::unique_ptr<TH1D> rank_reference[row_size];
    std::vector<Residue> veto_residues_list{  };
    rank_hist_list[2][0] = entry_iter->CreateMainChainAtomGausRankHistogram(0, chain_size, Residue::UNK, 2, veto_residues_list);
    rank_hist_list[1][0] = entry_iter->CreateMainChainAtomGausRankHistogram(1, chain_size, Residue::UNK, 1, veto_residues_list);
    rank_hist_list[0][0] = entry_iter->CreateMainChainAtomGausRankHistogram(2, chain_size, Residue::UNK, 0, veto_residues_list);
    rank_hist_list[2][1] = entry_iter->CreateMainChainAtomGausRankHistogram(0, chain_size, Residue::GLY);
    rank_hist_list[1][1] = entry_iter->CreateMainChainAtomGausRankHistogram(1, chain_size, Residue::GLY);
    rank_hist_list[0][1] = entry_iter->CreateMainChainAtomGausRankHistogram(2, chain_size, Residue::GLY);
    rank_hist_list[2][2] = entry_iter->CreateMainChainAtomGausRankHistogram(0, chain_size, Residue::PRO);
    rank_hist_list[1][2] = entry_iter->CreateMainChainAtomGausRankHistogram(1, chain_size, Residue::PRO);
    rank_hist_list[0][2] = entry_iter->CreateMainChainAtomGausRankHistogram(2, chain_size, Residue::PRO);
    
    for (int j = 0; j < row_size; j++)
    {
        for (auto & hist : rank_hist_list[j][0]) hist->Scale(100.0 / hist->GetEntries());
        for (auto & hist : rank_hist_list[j][1]) hist->Scale(100.0 / hist->GetEntries());
        for (auto & hist : rank_hist_list[j][2]) hist->Scale(100.0 / hist->GetEntries());
        rank_reference[j] = ROOTHelper::CreateHist1D(Form("rank_reference_%d", j),"", 4, 0.5, 4.5);
        for (int p = 1; p <= rank_reference[j]->GetNbinsX(); p++)
        {
            rank_reference[j]->SetBinContent(p, 150.0);
        }
    }

    gStyle->SetTextFont(22);
    gStyle->SetPaintTextFormat(".1f");
    std::unique_ptr<TPaveText> element_text[col_size];
    std::unique_ptr<TPaveText> par_text[row_size];
    std::unique_ptr<TLegend> legend;
    for (int i = 0; i < col_size; i++)
    {
        auto element_color{ AtomClassifier::GetMainChainElementColor(static_cast<size_t>(i)) };
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 0, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 100, 0.1, 4.9, 100, 0.01, 110.0);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 45.0f, 0.8f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.01f, 103, kCyan+3);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 5);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 28.0f, 1.3f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 25.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("Rank");
            frame[i][j]->GetYaxis()->SetTitle("Count Ratio (%)");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();

            ROOTHelper::SetFillAttribute(rank_reference[j].get(), 1001, kGray, 0.2f);
            rank_reference[j]->SetBarWidth(0.9f);
            rank_reference[j]->SetBarOffset(0.05f);
            rank_reference[j]->Draw("BAR SAME");

            auto text_size{ 0.2f * 4.0f/static_cast<float>(gPad->GetAbsHNDC()) };
            ROOTHelper::SetFillAttribute(rank_hist_list[j][0].at(static_cast<size_t>(i)).get(), 1001, kAzure-7, 0.5f);
            ROOTHelper::SetLineAttribute(rank_hist_list[j][0].at(static_cast<size_t>(i)).get(), 1, 0);
            ROOTHelper::SetMarkerAttribute(rank_hist_list[j][0].at(static_cast<size_t>(i)).get(), 1, text_size, kAzure-7);
            rank_hist_list[j][0].at(static_cast<size_t>(i))->SetBarWidth(0.25f);
            rank_hist_list[j][0].at(static_cast<size_t>(i))->SetBarOffset(0.125f);
            rank_hist_list[j][0].at(static_cast<size_t>(i))->Draw("BAR TEXT0 SAME");

            ROOTHelper::SetFillAttribute(rank_hist_list[j][1].at(static_cast<size_t>(i)).get(), 1001, kOrange-6, 0.5f);
            ROOTHelper::SetLineAttribute(rank_hist_list[j][1].at(static_cast<size_t>(i)).get(), 1, 0);
            ROOTHelper::SetMarkerAttribute(rank_hist_list[j][1].at(static_cast<size_t>(i)).get(), 1, text_size, kOrange-6);
            rank_hist_list[j][1].at(static_cast<size_t>(i))->SetBarWidth(0.25f);
            rank_hist_list[j][1].at(static_cast<size_t>(i))->SetBarOffset(0.375f);
            rank_hist_list[j][1].at(static_cast<size_t>(i))->Draw("BAR TEXT0 SAME");

            ROOTHelper::SetFillAttribute(rank_hist_list[j][2].at(static_cast<size_t>(i)).get(), 1001, kMagenta-2, 0.5f);
            ROOTHelper::SetLineAttribute(rank_hist_list[j][2].at(static_cast<size_t>(i)).get(), 1, 0);
            ROOTHelper::SetMarkerAttribute(rank_hist_list[j][2].at(static_cast<size_t>(i)).get(), 1, text_size, kMagenta-2);
            rank_hist_list[j][2].at(static_cast<size_t>(i))->SetBarWidth(0.25f);
            rank_hist_list[j][2].at(static_cast<size_t>(i))->SetBarOffset(0.625f);
            rank_hist_list[j][2].at(static_cast<size_t>(i))->Draw("BAR TEXT0 SAME");

            if (j == row_size - 1)
            {
                element_text[i] = ROOTHelper::CreatePaveText(0.01, 1.01, 0.99, 1.20, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(element_text[i].get());
                ROOTHelper::SetPaveAttribute(element_text[i].get(), 0, 0.2);
                ROOTHelper::SetTextAttribute(element_text[i].get(), 30.0f, 133, 22);
                ROOTHelper::SetFillAttribute(element_text[i].get(), 1001, element_color, 0.5f);
                element_text[i]->AddText(AtomClassifier::GetMainChainElementLabel(static_cast<size_t>(i)).data());
                element_text[i]->Draw();
            }
            if (i == col_size - 1)
            {
                par_text[j] = ROOTHelper::CreatePaveText(1.01, 0.01, 1.22, 0.99, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(par_text[j].get());
                ROOTHelper::SetPaveAttribute(par_text[j].get(), 0, 0.2);
                ROOTHelper::SetTextAttribute(par_text[j].get(), 30.0f, 133, 22, 0.0f, kYellow-10);
                ROOTHelper::SetLineAttribute(par_text[j].get(), 1, 0);
                ROOTHelper::SetFillAttribute(par_text[j].get(), 1001, kAzure-7);
                par_text[j]->AddText(y_title[j].data());
                auto text{ par_text[j]->GetLineWith(y_title[j].data()) };
                text->SetTextAngle(90.0f);
                par_text[j]->Draw();
            }
        }
    }
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}


} // namespace rhbm_gem
