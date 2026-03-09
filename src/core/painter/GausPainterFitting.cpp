#include <rhbm_gem/core/painter/GausPainter.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/PotentialEntryIterator.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
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

namespace rhbm_gem {
void GausPainter::PaintAtomGroupGausAminoAcidMainChainComponentSimple(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, "GausPainter::PaintAtomGroupGausAminoAcidMainChainComponentSimple");

    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    const auto & chemical_component_map{ model_object->GetChemicalComponentEntryMap() };

    const std::vector<Spot> spot_list{ Spot::CA, Spot::C, Spot::N };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(2.0);
    gStyle->SetGridColor(kGray);
    gStyle->SetEndErrorSize(5.0f);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 1000) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size{ 4 };
    
    std::unique_ptr<TPad> pad[pad_size];
    std::unique_ptr<TH2> frame[pad_size];

    
    pad[0] = ROOTHelper::CreatePad("pad0","", 0.00, 0.00, 1.00, 0.45);
    pad[1] = ROOTHelper::CreatePad("pad1","", 0.00, 0.45, 1.00, 0.80);
    pad[2] = ROOTHelper::CreatePad("pad2","", 0.00, 0.80, 1.00, 0.95);
    pad[3] = ROOTHelper::CreatePad("pad3","", 0.00, 0.95, 1.00, 1.00);

    frame[0] = ROOTHelper::CreateHist2D("hist_0","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[1] = ROOTHelper::CreateHist2D("hist_1","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[2] = ROOTHelper::CreateHist2D("hist_2","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[3] = ROOTHelper::CreateHist2D("hist_3","", 100, 0.0, 1.0, 100, 0.0, 1.0);

    auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    std::map<Spot, std::vector<GroupKey>> group_key_list_map;
    for (size_t i = 0; i < spot_list.size(); i++)
    {
        auto spot{ spot_list.at(i) };
        group_key_list_map.emplace(
            spot, m_atom_classifier->GetMainChainComponentAtomClassGroupKeyList(i));
        auto & group_key_list{ group_key_list_map.at(spot) };
        for (auto it = group_key_list.begin(); it != group_key_list.end(); )
        {
            if (entry_iter->IsAvailableAtomGroupKey(*it, class_key) == false)
            {
                it = group_key_list.erase(it);
            }
            else
            {
                ++it;
            }
        }
    }

    std::vector<std::string> component_id_list;
    component_id_list.reserve(ChemicalDataHelper::GetStandardAminoAcidCount());
    for (auto & group_key : group_key_list_map.at(Spot::CA))
    {
        auto residue{ entry_iter->GetResidueFromAtomGroupKey(group_key, class_key) };
        auto component_key{ static_cast<ComponentKey>(residue) };
        auto component_id{ chemical_component_map.at(component_key)->GetComponentId() };
        component_id_list.emplace_back(component_id);
    }

    std::map<Spot, std::unique_ptr<TGraphErrors>> amplitude_graph_map;
    std::map<Spot, std::unique_ptr<TGraphErrors>> width_graph_map;
    std::vector<double> amplitude_array, width_array;
    amplitude_array.reserve(group_key_list_map.at(Spot::CA).size()*4);
    width_array.reserve(group_key_list_map.at(Spot::CA).size()*4);
    std::map<Spot, short> color_list{ { Spot::CA,kRed }, { Spot::C, kBlue }, { Spot::N, kGreen+2 } };
    std::map<Spot, short> marker_list{ { Spot::CA, 72 }, { Spot::C, 71 }, { Spot::N, 73 } };
    std::map<Spot, short> line_list{ { Spot::CA, 1 }, { Spot::C, 2 }, { Spot::N, 3 } };
    
    for (auto & [spot, group_key_list] : group_key_list_map)
    {
        auto amplitude_graph{
            entry_iter->CreateAtomGausEstimateToResidueGraph(
                group_key_list, class_key, 0)
        };
        auto width_graph{
            entry_iter->CreateAtomGausEstimateToResidueGraph(
                group_key_list, class_key, 1)
        };
        for (int p = 0; p < amplitude_graph->GetN(); p++)
        {
            amplitude_array.push_back(amplitude_graph->GetPointY(p));
            width_array.push_back(width_graph->GetPointY(p));
        }
        ROOTHelper::SetMarkerAttribute(amplitude_graph.get(), marker_list[spot], 1.5f, color_list[spot]);
        ROOTHelper::SetMarkerAttribute(width_graph.get(), marker_list[spot], 1.5f, color_list[spot]);
        ROOTHelper::SetLineAttribute(amplitude_graph.get(), line_list[spot], 2, color_list[spot]);
        ROOTHelper::SetLineAttribute(width_graph.get(), line_list[spot], 2, color_list[spot]);

        amplitude_graph_map[spot] = std::move(amplitude_graph);
        width_graph_map[spot] = std::move(width_graph);
    }

    auto scaling{ 0.3 };
    auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, scaling) };
    auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, scaling) };

    auto count_hist{
        entry_iter->CreateComponentCountHistogram(group_key_list_map.at(Spot::CA), class_key)
    };

    canvas->cd();
    for (int i = 0; i < pad_size; i++)
    {
        ROOTHelper::SetPadDefaultStyle(pad[i].get());
        pad[i]->Draw();
    }

    frame[0]->GetYaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));
    frame[1]->GetYaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));

    pad[0]->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.15, 0.02, 0.11, 0.02);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    RemodelFrameInPad(frame[0].get(), pad[0].get(), 0.0, 0.015);
    RemodelAxisLabels(frame[0]->GetXaxis(), component_id_list, 90.0, 12);
    ROOTHelper::SetAxisTitleAttribute(frame[0]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTitleAttribute(frame[0]->GetYaxis(), 50.0f, 1.5f);
    ROOTHelper::SetAxisLabelAttribute(frame[0]->GetXaxis(), 50.0f, 0.13f, 103, kCyan+3);
    ROOTHelper::SetAxisLabelAttribute(frame[0]->GetYaxis(), 45.0f, 0.01f);
    frame[0]->GetYaxis()->SetTitle("Width");
    frame[0]->Draw();
    for (auto & [element, graph] : width_graph_map) graph->Draw("PL X0");

    pad[1]->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.15, 0.02, 0.02, 0.01);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    RemodelFrameInPad(frame[1].get(), pad[1].get(), 0.0, 0.015);
    RemodelAxisLabels(frame[1]->GetXaxis(), component_id_list, 90.0, 12);
    ROOTHelper::SetAxisTitleAttribute(frame[1]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTitleAttribute(frame[1]->GetYaxis(), 50.0f, 1.5f);
    ROOTHelper::SetAxisLabelAttribute(frame[1]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(frame[1]->GetYaxis(), 45.0f, 0.01f);
    frame[1]->GetYaxis()->SetTitle("Amplitude");
    frame[1]->Draw();
    for (auto & [element, graph] : amplitude_graph_map) graph->Draw("PL X0");

    pad[2]->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.15, 0.02, 0.2, 0.003);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    RemodelFrameInPad(frame[2].get(), pad[2].get(), 0.0, 0.015);
    RemodelAxisLabels(frame[2]->GetXaxis(), component_id_list, 90.0, 12);
    ROOTHelper::SetAxisTitleAttribute(frame[2]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTitleAttribute(frame[2]->GetYaxis(), 40.0f, 1.2f);
    ROOTHelper::SetAxisLabelAttribute(frame[2]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(frame[2]->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisTickAttribute(frame[2]->GetYaxis(), 0.0f, 504);
    frame[2]->GetYaxis()->SetTitle("#splitline{Member}{Counts}");
    frame[2]->GetYaxis()->SetLimits(0.5, count_hist->GetMaximum()*1.1);
    frame[2]->GetYaxis()->CenterTitle();
    frame[2]->Draw();
    gStyle->SetTextFont(132);
    ROOTHelper::SetFillAttribute(count_hist.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetLineAttribute(count_hist.get(), 0, 0);
    ROOTHelper::SetMarkerAttribute(count_hist.get(), 20, 8.0f, kAzure);
    count_hist->Draw("BAR TEXT0 SAME");

    pad[3]->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.0, 0.0, 0.0, 0.0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    auto legend{
        ROOTHelper::CreateLegend(0.05, 0.00, 1.00, 1.00, false)
    };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 35.0f, 133, 12);
    legend->AddEntry(amplitude_graph_map.at(Spot::CA).get(), "Alpha Carbon", "lp");
    legend->AddEntry(amplitude_graph_map.at(Spot::C).get(), "Carbon", "lp");
    legend->AddEntry(amplitude_graph_map.at(Spot::N).get(), "Nitrogen", "lp");
    legend->SetNColumns(3);
    legend->SetMargin(0.50f);
    legend->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);

    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void GausPainter::PaintAtomGroupGausAminoAcidMainChainStructure(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, "GausPainter::PaintAtomGroupGausAminoAcidMainChainStructure");

    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };
    const auto & chemical_component_map{ model_object->GetChemicalComponentEntryMap() };

    const std::vector<Spot> spot_list{ Spot::CA, Spot::C, Spot::N, Spot::O };
    const std::vector<Structure> structure_list{ Structure::FREE, Structure::HELX_P, Structure::SHEET };
    
    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(2.0);
    gStyle->SetGridColor(kGray);
    gStyle->SetEndErrorSize(5.0f);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 750) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size{ 8 };
    
    std::unique_ptr<TPad> pad[pad_size];
    std::unique_ptr<TH2> frame[pad_size];

    pad[0] = ROOTHelper::CreatePad("pad0","", 0.55, 0.80, 1.00, 1.00); // The title pad
    pad[1] = ROOTHelper::CreatePad("pad1","", 0.00, 0.00, 0.55, 0.45); // The bottom-left pad
    pad[2] = ROOTHelper::CreatePad("pad2","", 0.00, 0.45, 0.55, 0.80); // The top-left pad
    pad[3] = ROOTHelper::CreatePad("pad3","", 0.55, 0.00, 0.90, 0.60); // The bottom-middle pad
    pad[4] = ROOTHelper::CreatePad("pad4","", 0.55, 0.60, 0.90, 0.80); // The top-middle pad
    pad[5] = ROOTHelper::CreatePad("pad5","", 0.90, 0.00, 1.00, 0.60); // The bottom-right pad
    pad[6] = ROOTHelper::CreatePad("pad6","", 0.90, 0.60, 1.00, 0.80); // The top-right pad
    pad[7] = ROOTHelper::CreatePad("pad7","", 0.00, 0.80, 0.55, 1.00); // The left-title pad

    frame[0] = ROOTHelper::CreateHist2D("hist_0","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[1] = ROOTHelper::CreateHist2D("hist_1","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[2] = ROOTHelper::CreateHist2D("hist_2","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[3] = ROOTHelper::CreateHist2D("hist_3","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[4] = ROOTHelper::CreateHist2D("hist_4","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[5] = ROOTHelper::CreateHist2D("hist_5","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[6] = ROOTHelper::CreateHist2D("hist_6","", 100, 0.0, 1.0, 100, 0.0, 1.0);

    auto class_key{ ChemicalDataHelper::GetStructureAtomClassKey() };
    
    for (size_t i = 0; i < spot_list.size(); i++)
    {
        auto spot{ spot_list.at(i) };
        std::map<Structure, std::vector<GroupKey>> group_key_list_map;
        
        for (size_t j = 0; j < structure_list.size(); j++)
        {
            auto structure{ structure_list.at(j) };
            group_key_list_map.emplace(
                structure, m_atom_classifier->GetMainChainStructureAtomClassGroupKeyList(i, structure));
            auto & group_key_list{ group_key_list_map.at(structure) };
            for (auto it = group_key_list.begin(); it != group_key_list.end(); )
            {
                if (entry_iter->IsAvailableAtomGroupKey(*it, class_key) == false)
                {
                    it = group_key_list.erase(it);
                }
                else
                {
                    ++it;
                }
            }
        }

        std::vector<std::string> component_id_list;
        component_id_list.reserve(ChemicalDataHelper::GetStandardAminoAcidCount());
        for (auto & residue : ChemicalDataHelper::GetStandardAminoAcidList())
        {
            auto component_key{ static_cast<ComponentKey>(residue) };
            auto component_id{ chemical_component_map.at(component_key)->GetComponentId() };
            component_id_list.emplace_back(component_id);
        }

        std::map<Structure, std::unique_ptr<TGraphErrors>> amplitude_graph_map;
        std::map<Structure, std::unique_ptr<TGraphErrors>> width_graph_map;
        std::map<Structure, std::unique_ptr<TGraphErrors>> correlation_graph_map;
        std::map<Structure, std::unique_ptr<TH2D>> amplitude_hist_map;
        std::map<Structure, std::unique_ptr<TH2D>> width_hist_map;
        std::vector<double> amplitude_array, width_array;
        amplitude_array.reserve(80);
        width_array.reserve(80);
        
        for (auto & [structure, group_key_list] : group_key_list_map)
        {
            auto amplitude_graph{
                entry_iter->CreateAtomGausEstimateToResidueGraph(
                    group_key_list, class_key, 0)
            };
            auto width_graph{
                entry_iter->CreateAtomGausEstimateToResidueGraph(
                    group_key_list, class_key, 1)
            };
            auto correlation_graph{
                entry_iter->CreateAtomGausEstimateScatterGraph(group_key_list, class_key, 0, 1)
            };
            for (int p = 0; p < amplitude_graph->GetN(); p++)
            {
                amplitude_array.push_back(amplitude_graph->GetPointY(p));
                width_array.push_back(width_graph->GetPointY(p));
            }
            auto structure_color{ AtomClassifier::GetMainChainStructureColor(structure) };
            auto structure_marker{ AtomClassifier::GetMainChainStructureMarker(structure) };
            auto structure_line{ AtomClassifier::GetMainChainStructureLineStyle(structure) };
            ROOTHelper::SetMarkerAttribute(amplitude_graph.get(), structure_marker, 1.5f, structure_color);
            ROOTHelper::SetMarkerAttribute(width_graph.get(), structure_marker, 1.5f, structure_color);
            ROOTHelper::SetMarkerAttribute(correlation_graph.get(), structure_marker, 1.5f, structure_color);
            ROOTHelper::SetLineAttribute(amplitude_graph.get(), structure_line, 1, structure_color);
            ROOTHelper::SetLineAttribute(width_graph.get(), structure_line, 1, structure_color);
            ROOTHelper::SetLineAttribute(correlation_graph.get(), structure_line, 1, structure_color);

            amplitude_graph_map[structure] = std::move(amplitude_graph);
            width_graph_map[structure] = std::move(width_graph);
            correlation_graph_map[structure] = std::move(correlation_graph);
        }

        auto scaling{ 0.3 };
        auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, scaling) };
        auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, scaling) };
        auto structure_count{ structure_list.size() };
        std::vector<std::string> structure_label_list;
        structure_label_list.reserve(structure_count);
        for (size_t j = 0; j < structure_count; j++)
        {
            auto structure{ structure_list.at(j) };
            auto structure_label{ AtomClassifier::GetMainChainStructureLabel(structure) };
            auto x_value{ static_cast<double>(j) };
            structure_label_list.emplace_back(structure_label);
            std::string name_amplitude{ "amplitude_hist_"+ structure_label };
            std::string name_width{ "width_hist_"+ structure_label };
            auto amplitude_hist{
                ROOTHelper::CreateHist2D(
                    name_amplitude.data(),"",
                    100, std::get<0>(amplitude_range), std::get<1>(amplitude_range),
                    static_cast<int>(structure_count), -0.5, static_cast<int>(structure_count)-0.5)
            };
            auto width_hist{
                ROOTHelper::CreateHist2D(
                    name_width.data(),"",
                    static_cast<int>(structure_count), -0.5, static_cast<int>(structure_count)-0.5,
                    100, std::get<0>(width_range), std::get<1>(width_range))
            };
            for (int p = 0; p < amplitude_graph_map.at(structure)->GetN(); p++)
            {
                amplitude_hist->Fill(amplitude_graph_map.at(structure)->GetPointY(p), x_value);
            }
            for (int p = 0; p < width_graph_map.at(structure)->GetN(); p++)
            {
                width_hist->Fill(x_value, width_graph_map.at(structure)->GetPointY(p));
            }
            auto structure_color{ AtomClassifier::GetMainChainStructureColor(structure) };
            ROOTHelper::SetLineAttribute(amplitude_hist.get(), 1, 1, structure_color);
            ROOTHelper::SetLineAttribute(width_hist.get(), 1, 1, structure_color);
            ROOTHelper::SetFillAttribute(amplitude_hist.get(), 1001, structure_color, 0.3f);
            ROOTHelper::SetFillAttribute(width_hist.get(), 1001, structure_color, 0.3f);
            amplitude_hist->SetBarWidth(0.6f);
            width_hist->SetBarWidth(0.6f);
            amplitude_hist_map[structure] = std::move(amplitude_hist);
            width_hist_map[structure] = std::move(width_hist);
        }

        auto count_hist{ entry_iter->CreateAtomResidueCountHistogram2D(class_key) };

        canvas->cd();
        for (int j = 0; j < pad_size; j++)
        {
            ROOTHelper::SetPadDefaultStyle(pad[j].get());
            pad[j]->Draw();
        }

        pad[0]->cd();
        auto info_text{ CreateDataInfoPaveText(model_object) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, info_text.get(), 0.04, 0.18, 0.02, 0.02);
        info_text->Draw();

        auto resolution_text{ CreateResolutionPaveText(model_object) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, resolution_text.get(), 0.28, 0.01, 0.02, 0.02);
        resolution_text->Draw();

        frame[0]->GetYaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));
        frame[1]->GetYaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        frame[2]->GetYaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));
        frame[3]->GetXaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        frame[4]->GetXaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));
        frame[4]->GetYaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));

        pad[1]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.10, 0.00, 0.11, 0.02);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        RemodelFrameInPad(frame[0].get(), pad[1].get(), 0.0, 0.015);
        RemodelAxisLabels(frame[0]->GetXaxis(), component_id_list, 90.0, 12);
        ROOTHelper::SetAxisTitleAttribute(frame[0]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisTitleAttribute(frame[0]->GetYaxis(), 50.0f, 1.5f);
        ROOTHelper::SetAxisLabelAttribute(frame[0]->GetXaxis(), 40.0f, 0.13f, 103, kCyan+3);
        ROOTHelper::SetAxisLabelAttribute(frame[0]->GetYaxis(), 45.0f, 0.01f);
        frame[0]->GetYaxis()->SetTitle("Width");
        frame[0]->Draw();
        for (auto & [element, graph] : width_graph_map) graph->Draw("PL X0");

        auto spot_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, spot_text.get(), 0.005, 0.50, 0.01, 0.35);
        ROOTHelper::SetPaveTextDefaultStyle(spot_text.get());
        ROOTHelper::SetPaveAttribute(spot_text.get(), 0, 0.1);
        ROOTHelper::SetFillAttribute(spot_text.get(), 1001, kAzure-7);
        ROOTHelper::SetTextAttribute(spot_text.get(), 60.0f, 103, 22, 0.0, kYellow-10);
        auto spot_label{ AtomClassifier::GetMainChainSpotLabel(spot) };
        spot_text->AddText(spot_label.data());
        spot_text->Draw();

        pad[2]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.10, 0.00, 0.02, 0.01);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        RemodelFrameInPad(frame[1].get(), pad[2].get(), 0.0, 0.015);
        RemodelAxisLabels(frame[1]->GetXaxis(), component_id_list, 90.0, 12);
        ROOTHelper::SetAxisTitleAttribute(frame[1]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisTitleAttribute(frame[1]->GetYaxis(), 50.0f, 1.5f);
        ROOTHelper::SetAxisLabelAttribute(frame[1]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(frame[1]->GetYaxis(), 45.0f, 0.01f);
        frame[1]->GetYaxis()->SetTitle("Amplitude");
        frame[1]->Draw();
        for (auto & [element, graph] : amplitude_graph_map) graph->Draw("PL X0");

        pad[3]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.09, 0.01, 0.12, 0.02);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        RemodelFrameInPad(frame[4].get(), pad[3].get(), 0.03, 0.015);
        ROOTHelper::SetAxisTitleAttribute(frame[4]->GetXaxis(), 40.0f, 1.1f);
        ROOTHelper::SetAxisTitleAttribute(frame[4]->GetYaxis(), 40.0f, 1.2f);
        ROOTHelper::SetAxisLabelAttribute(frame[4]->GetXaxis(), 40.0f, 0.01f);
        ROOTHelper::SetAxisLabelAttribute(frame[4]->GetYaxis(), 40.0f, 0.01f);
        auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.03, 0) };
        auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, 0.015, 1) };
        ROOTHelper::SetAxisTickAttribute(frame[4]->GetXaxis(), static_cast<float>(x_tick_length), 505);
        ROOTHelper::SetAxisTickAttribute(frame[4]->GetYaxis(), static_cast<float>(y_tick_length), 505);
        frame[4]->GetXaxis()->SetTitle("Group Amplitude");
        frame[4]->GetYaxis()->SetTitle("Group Width");
        frame[4]->Draw("");
        for (auto & [element, graph] : correlation_graph_map) graph->Draw("P X0");

        pad[4]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.09, 0.01, 0.02, 0.02);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        RemodelFrameInPad(frame[3].get(), pad[4].get(), 0.0, 0.0);
        RemodelAxisLabels(frame[3]->GetYaxis(), structure_label_list, 0.0, 22);
        ROOTHelper::SetAxisTitleAttribute(frame[3]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisTitleAttribute(frame[3]->GetYaxis(), 40.0f, 1.2f);
        ROOTHelper::SetAxisLabelAttribute(frame[3]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(frame[3]->GetYaxis(), 35.0f, 0.04f);
        frame[3]->GetYaxis()->SetTitle("Element");
        frame[3]->Draw();
        for (auto & [element, hist] : amplitude_hist_map) hist->Draw("CANDLEY3 SAME");

        pad[5]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.01, 0.01, 0.12, 0.02);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        RemodelFrameInPad(frame[2].get(), pad[5].get(), 0.0, 0.0);
        RemodelAxisLabels(frame[2]->GetXaxis(), structure_label_list, 0.0, -1);
        ROOTHelper::SetAxisTitleAttribute(frame[2]->GetXaxis(), 40.0f, 1.0f);
        ROOTHelper::SetAxisTitleAttribute(frame[2]->GetYaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(frame[2]->GetXaxis(), 35.0f, 0.005f);
        ROOTHelper::SetAxisLabelAttribute(frame[2]->GetYaxis(), 0.0f);
        frame[2]->GetXaxis()->SetTitle("Element");
        frame[2]->Draw();
        for (auto & [element, hist] : width_hist_map) hist->Draw("CANDLE3 SAME");

        pad[6]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.02, 0.02, 0.20, 0.20);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        auto legend{
            ROOTHelper::CreateLegend(0.00, 0.00, 1.00, 1.00, false)
        };
        ROOTHelper::SetLegendDefaultStyle(legend.get());
        ROOTHelper::SetTextAttribute(legend.get(), 35.0f, 133, 12);
        for (auto & structure : structure_list)
        {
            auto structure_label{ AtomClassifier::GetMainChainStructureLabel(structure) };
            legend->AddEntry(
                amplitude_graph_map.at(structure).get(), structure_label.data(), "lp");
        }
        legend->SetNColumns(1);
        legend->SetMargin(0.50f);
        legend->Draw();

        pad[7]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.10, 0.00, 0.02, 0.01);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        RemodelFrameInPad(frame[6].get(), pad[7].get(), 0.0, 0.015);
        RemodelAxisLabels(frame[6]->GetXaxis(), component_id_list, 90.0, 12);
        RemodelAxisLabels(frame[6]->GetYaxis(), structure_label_list, 0.0, 22);
        ROOTHelper::SetAxisTitleAttribute(frame[6]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisTitleAttribute(frame[6]->GetYaxis(), 40.0f, 1.7f);
        ROOTHelper::SetAxisLabelAttribute(frame[6]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(frame[6]->GetYaxis(), 40.0f, 0.03f);
        ROOTHelper::SetAxisTickAttribute(frame[6]->GetYaxis(), 0.0f, 504);
        frame[6]->GetYaxis()->SetTitle("#splitline{Member}{Counts}");
        frame[6]->GetYaxis()->CenterTitle();
        frame[6]->Draw();
        gStyle->SetTextFont(132);
        gStyle->SetPalette(kLightTemperature);
        ROOTHelper::SetFillAttribute(count_hist.get(), 1001, kAzure-7, 0.5f);
        ROOTHelper::SetLineAttribute(count_hist.get(), 0, 0);
        ROOTHelper::SetMarkerAttribute(count_hist.get(), 20, 7.0f, kBlack);
        count_hist->Draw("COL TEXT SAME");
        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void GausPainter::PaintAtomLocalGausToSequenceAminoAcidMainChain(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomLocalGausToSequenceAminoAcidMainChain");
    auto entry_iter{ std::make_unique<PotentialEntryIterator>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 1 };
    const int row_size{ 2 };
    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 900) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.08f, 0.02f, 0.10f, 0.10f, 0.02f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const int main_chain_element_count{ 4 };
    const std::string y_title[row_size]
    {
        "Width",
        "Amplitude"
    };
    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> gaus_graph_map[row_size][main_chain_element_count];
    for (size_t k = 0; k < main_chain_element_count; k++)
    {
        gaus_graph_map[1][k] = entry_iter->CreateAtomGausEstimateToSequenceIDGraphMap(k, 0);
        gaus_graph_map[0][k] = entry_iter->CreateAtomGausEstimateToSequenceIDGraphMap(k, 1);
        //gaus_graph_map[0][k] = entry_iter->CreateAtomGausEstimateToSequenceIDGraphMap(k, 2);
        //gaus_graph_map[0][k] = entry_iter->CreateAtomMapValueToSequenceIDGraphMap(k);
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
            auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array[j], 0.4) };
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
                    auto spot{ AtomClassifier::GetMainChainSpot(k) };
                    auto element_color{ AtomClassifier::GetMainChainSpotColor(spot) };
                    auto element_marker{ AtomClassifier::GetMainChainSpotOpenMarker(spot) };
                    //auto element_line{ AtomClassifier::GetMainChainSpotLineStyle(spot) };
                    ROOTHelper::SetMarkerAttribute(gaus_graph_map[j][k].at(chain_id).get(), element_marker, 0.7f, element_color);
                    ROOTHelper::SetLineAttribute(gaus_graph_map[j][k].at(chain_id).get(), 1, 1, element_color);
                    gaus_graph_map[j][k].at(chain_id)->Draw("PL X0");
                    //gaus_graph_map[j][k].at(chain_id)->Draw("PL");
                }

                if (i == 0 && j == row_size - 1)
                {
                    subtitle1_text = ROOTHelper::CreatePaveText(0.00, 1.02, 0.15, 1.24, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle1_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle1_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle1_text.get(), 1001, kAzure-7);
                    ROOTHelper::SetTextAttribute(subtitle1_text.get(), 60.0f, 133, 22, 0.0, kYellow-10);
                    subtitle1_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
                    subtitle1_text->Draw();

                    subtitle2_text = ROOTHelper::CreatePaveText(0.16, 1.02, 0.35, 1.24, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle2_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle2_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle2_text.get(), 1001, kAzure-7, 0.5f);
                    ROOTHelper::SetTextAttribute(subtitle2_text.get(), 50.0f, 103, 22);
                    subtitle2_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
                    subtitle2_text->Draw();

                    subtitle3_text = ROOTHelper::CreatePaveText(0.36, 1.02, 0.57, 1.24, "nbNDC ARC", true);
                    ROOTHelper::SetPaveTextDefaultStyle(subtitle3_text.get());
                    ROOTHelper::SetPaveAttribute(subtitle3_text.get(), 0, 0.2);
                    ROOTHelper::SetFillAttribute(subtitle3_text.get(), 1001, kAzure-7, 0.5f);
                    ROOTHelper::SetTextAttribute(subtitle3_text.get(), 50.0f, 103, 22);
                    subtitle3_text->AddText(model_object->GetEmdID().data());
                    subtitle3_text->Draw();

                    legend = ROOTHelper::CreateLegend(0.58, 1.02, 0.98, 1.24, true);
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

#ifdef HAVE_ROOT
void GausPainter::RemodelFrameInPad(
    TH2 * frame, TPad * pad, double x_tick, double y_tick)
{
    pad->cd();
    ROOTHelper::SetPadLayout(pad, 1, 1, 0, 0, 0, 0);
    auto x_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, x_tick, 0) };
    auto y_tick_length{ ROOTHelper::ConvertGlobalTickLengthToPadTickLength(gPad, y_tick, 1) };
    ROOTHelper::SetAxisTickAttribute(frame->GetXaxis(), static_cast<float>(x_tick_length), 505);
    ROOTHelper::SetAxisTickAttribute(frame->GetYaxis(), static_cast<float>(y_tick_length), 505);
    frame->SetStats(0);
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
}

void GausPainter::RemodelAxisLabels(
    TAxis * axis, const std::vector<std::string> & label_list, double angle, int align)
{
    auto label_size{ static_cast<int>(label_list.size()) };
    ROOTHelper::SetAxisTickAttribute(axis, 0.0, label_size + 1);
    axis->SetLimits(-0.75, static_cast<float>(label_size) - 0.5f);
    //axis->ChangeLabel(1, -1.0, 0.0);
    for (size_t i = 0; i < label_list.size(); i++)
    {
        auto label_index{ static_cast<int>(i) + 1 };
        axis->ChangeLabel(label_index, angle, -1, align, -1, -1, label_list.at(i).data());
    }
}

std::unique_ptr<TPaveText> GausPainter::CreateDataInfoPaveText(
    ModelObject * model_object) const
{
    auto pavetext{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(pavetext.get());
    ROOTHelper::SetPaveAttribute(pavetext.get(), 0, 0.1);
    ROOTHelper::SetFillAttribute(pavetext.get(), 1001, kAzure-7, 0.5);
    ROOTHelper::SetTextAttribute(pavetext.get(), 55, 133, 12);
    pavetext->AddText(("#font[102]{PDB-" + model_object->GetPdbID() +"}").data());
    pavetext->AddText(("#font[102]{"+ model_object->GetEmdID() +"}").data());
    return pavetext;
}

std::unique_ptr<TPaveText> GausPainter::CreateResolutionPaveText(
    ModelObject * model_object) const
{
    auto pavetext{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(pavetext.get());
    ROOTHelper::SetPaveAttribute(pavetext.get(), 0, 0.1);
    ROOTHelper::SetFillAttribute(pavetext.get(), 1001, kAzure-7);
    ROOTHelper::SetLineAttribute(pavetext.get(), 1, 0);
    ROOTHelper::SetTextAttribute(pavetext.get(), 75.0f, 133, 22, 0.0, kYellow-10);
    if (model_object->GetEmdID().find("Simulation") != std::string::npos)
    {
        pavetext->AddText(Form("#sigma_{B}=%.2f", model_object->GetResolution()));
    }
    else
    {
        pavetext->AddText(Form("%.2f #AA", model_object->GetResolution()));
    }
    return pavetext;
}
#endif


} // namespace rhbm_gem
