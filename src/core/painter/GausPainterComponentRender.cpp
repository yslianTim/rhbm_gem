#include <rhbm_gem/core/painter/GausPainter.hpp>
#include <rhbm_gem/core/painter/PotentialPlotBuilder.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/ChemicalComponentEntry.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/PotentialEntryQuery.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLegend.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#endif

#include <map>
#include <vector>

namespace rhbm_gem {

void GausPainter::PaintAtomGroupGausAminoAcidMainChainComponent(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, "GausPainter::PaintAtomGroupGausAminoAcidMainChainComponent");

    auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    const auto & chemical_component_map{ model_object->GetChemicalComponentEntryMap() };
    const std::vector<Spot> spot_list{ Spot::CA, Spot::C, Spot::N };
    (void)chemical_component_map;

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
    std::map<Spot, std::unique_ptr<TGraphErrors>> correlation_graph_map;
    std::map<Spot, std::unique_ptr<TH2D>> amplitude_hist_map;
    std::map<Spot, std::unique_ptr<TH2D>> width_hist_map;
    std::vector<double> amplitude_array, width_array;
    amplitude_array.reserve(group_key_list_map.at(Spot::CA).size()*4);
    width_array.reserve(group_key_list_map.at(Spot::CA).size()*4);
    
    for (auto & [spot, group_key_list] : group_key_list_map)
    {
        auto amplitude_graph{
            plot_builder->CreateAtomGausEstimateToResidueGraph(
                group_key_list, class_key, 0)
        };
        auto width_graph{
            plot_builder->CreateAtomGausEstimateToResidueGraph(
                group_key_list, class_key, 1)
        };
        auto correlation_graph{
            plot_builder->CreateAtomGausEstimateScatterGraph(group_key_list, class_key, 1, 0)
        };
        for (int p = 0; p < amplitude_graph->GetN(); p++)
        {
            amplitude_array.push_back(amplitude_graph->GetPointY(p));
            width_array.push_back(width_graph->GetPointY(p));
        }
        auto spot_color{ AtomClassifier::GetMainChainSpotColor(spot) };
        auto spot_marker{ AtomClassifier::GetMainChainSpotOpenMarker(spot) };
        auto spot_line{ AtomClassifier::GetMainChainSpotLineStyle(spot) };
        ROOTHelper::SetMarkerAttribute(amplitude_graph.get(), spot_marker, 1.5f, spot_color);
        ROOTHelper::SetMarkerAttribute(width_graph.get(), spot_marker, 1.5f, spot_color);
        ROOTHelper::SetMarkerAttribute(correlation_graph.get(), spot_marker, 1.5f, spot_color);
        ROOTHelper::SetLineAttribute(amplitude_graph.get(), spot_line, 1, spot_color);
        ROOTHelper::SetLineAttribute(width_graph.get(), spot_line, 1, spot_color);
        ROOTHelper::SetLineAttribute(correlation_graph.get(), spot_line, 1, spot_color);

        amplitude_graph_map[spot] = std::move(amplitude_graph);
        width_graph_map[spot] = std::move(width_graph);
        correlation_graph_map[spot] = std::move(correlation_graph);
    }

    auto scaling{ 0.3 };
    auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, scaling) };
    auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, scaling) };
    auto spot_count{ spot_list.size() };
    std::vector<std::string> spot_label_list;
    spot_label_list.reserve(spot_count);
    for (size_t i = 0; i < spot_count; i++)
    {
        auto spot{ spot_list.at(i) };
        auto spot_label{ AtomClassifier::GetMainChainSpotLabel(spot) };
        auto x_value{ static_cast<double>(i) };
        spot_label_list.emplace_back(spot_label);
        std::string name_amplitude{ "amplitude_hist_"+ spot_label };
        std::string name_width{ "width_hist_"+ spot_label };
        auto amplitude_hist{
            ROOTHelper::CreateHist2D(
                name_amplitude.data(),"",
                static_cast<int>(spot_count), -0.5, static_cast<int>(spot_count)-0.5,
                100, std::get<0>(amplitude_range), std::get<1>(amplitude_range))
        };
        auto width_hist{
            ROOTHelper::CreateHist2D(
                name_width.data(),"",
                100, std::get<0>(width_range), std::get<1>(width_range),
                static_cast<int>(spot_count), -0.5, static_cast<int>(spot_count)-0.5)
        };
        for (int p = 0; p < amplitude_graph_map.at(spot)->GetN(); p++)
        {
            amplitude_hist->Fill(x_value, amplitude_graph_map.at(spot)->GetPointY(p));
        }
        for (int p = 0; p < width_graph_map.at(spot)->GetN(); p++)
        {
            width_hist->Fill(width_graph_map.at(spot)->GetPointY(p), x_value);
        }
        auto spot_color{ AtomClassifier::GetMainChainSpotColor(spot) };
        ROOTHelper::SetLineAttribute(amplitude_hist.get(), 1, 1, spot_color);
        ROOTHelper::SetLineAttribute(width_hist.get(), 1, 1, spot_color);
        ROOTHelper::SetFillAttribute(amplitude_hist.get(), 1001, spot_color, 0.3f);
        ROOTHelper::SetFillAttribute(width_hist.get(), 1001, spot_color, 0.3f);
        amplitude_hist->SetBarWidth(0.6f);
        width_hist->SetBarWidth(0.6f);
        amplitude_hist_map[spot] = std::move(amplitude_hist);
        width_hist_map[spot] = std::move(width_hist);
    }

    auto count_hist{
        plot_builder->CreateComponentCountHistogram(group_key_list_map.at(Spot::CA), class_key)
    };

    canvas->cd();
    for (int i = 0; i < pad_size; i++)
    {
        ROOTHelper::SetPadDefaultStyle(pad[i].get());
        pad[i]->Draw();
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
    frame[2]->GetYaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    frame[3]->GetXaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));
    frame[4]->GetYaxis()->SetLimits(std::get<0>(amplitude_range), std::get<1>(amplitude_range));
    frame[4]->GetXaxis()->SetLimits(std::get<0>(width_range), std::get<1>(width_range));

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
    frame[4]->GetYaxis()->SetTitle("Amplitude #font[2]{A}");
    frame[4]->GetXaxis()->SetTitle("Width #font[2]{#tau}");
    frame[4]->Draw("");
    for (auto & [element, graph] : correlation_graph_map) graph->Draw("P X0");

    pad[4]->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.09, 0.01, 0.02, 0.02);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    RemodelFrameInPad(frame[3].get(), pad[4].get(), 0.0, 0.0);
    RemodelAxisLabels(frame[3]->GetYaxis(), spot_label_list, 0.0, 22);
    ROOTHelper::SetAxisTitleAttribute(frame[3]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTitleAttribute(frame[3]->GetYaxis(), 40.0f, 1.2f);
    ROOTHelper::SetAxisLabelAttribute(frame[3]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(frame[3]->GetYaxis(), 35.0f, 0.04f, 103);
    frame[3]->GetYaxis()->SetTitle("Element");
    frame[3]->Draw();
    for (auto & [element, hist] : width_hist_map) hist->Draw("CANDLEY3 SAME");

    pad[5]->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.01, 0.01, 0.12, 0.02);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    RemodelFrameInPad(frame[2].get(), pad[5].get(), 0.0, 0.0);
    RemodelAxisLabels(frame[2]->GetXaxis(), spot_label_list, 0.0, -1);
    ROOTHelper::SetAxisTitleAttribute(frame[2]->GetXaxis(), 40.0f, 1.0f);
    ROOTHelper::SetAxisTitleAttribute(frame[2]->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(frame[2]->GetXaxis(), 35.0f, 0.005f, 103);
    ROOTHelper::SetAxisLabelAttribute(frame[2]->GetYaxis(), 0.0f);
    frame[2]->GetXaxis()->SetTitle("Element");
    frame[2]->Draw();
    for (auto & [element, hist] : amplitude_hist_map) hist->Draw("CANDLE3 SAME");

    pad[6]->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.02, 0.02, 0.20, 0.20);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    auto legend{
        ROOTHelper::CreateLegend(0.00, 0.00, 1.00, 1.00, false)
    };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetTextAttribute(legend.get(), 35.0f, 103, 12);
    for (auto & spot : spot_list)
    {
        auto spot_label{ AtomClassifier::GetMainChainSpotLabel(spot) };
        legend->AddEntry(
            amplitude_graph_map.at(spot).get(), spot_label.data(), "lp");
    }
    legend->SetNColumns(2);
    legend->SetMargin(0.50f);
    legend->Draw();

    pad[7]->cd();
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.10, 0.00, 0.02, 0.01);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    RemodelFrameInPad(frame[6].get(), pad[7].get(), 0.0, 0.015);
    RemodelAxisLabels(frame[6]->GetXaxis(), component_id_list, 90.0, 12);
    ROOTHelper::SetAxisTitleAttribute(frame[6]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisTitleAttribute(frame[6]->GetYaxis(), 40.0f, 1.2f);
    ROOTHelper::SetAxisLabelAttribute(frame[6]->GetXaxis(), 0.0f);
    ROOTHelper::SetAxisLabelAttribute(frame[6]->GetYaxis(), 0.0f);
    ROOTHelper::SetAxisTickAttribute(frame[6]->GetYaxis(), 0.0f, 504);
    frame[6]->GetYaxis()->SetTitle("#splitline{Member}{Counts}");
    frame[6]->GetYaxis()->SetLimits(0.5, count_hist->GetMaximum()*1.1);
    frame[6]->GetYaxis()->CenterTitle();
    frame[6]->Draw();
    gStyle->SetTextFont(132);
    ROOTHelper::SetFillAttribute(count_hist.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetLineAttribute(count_hist.get(), 0, 0);
    ROOTHelper::SetMarkerAttribute(count_hist.get(), 20, 7.0f, kAzure);
    count_hist->Draw("BAR TEXT0 SAME");
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);

    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

} // namespace rhbm_gem
