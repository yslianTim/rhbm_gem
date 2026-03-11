#include <rhbm_gem/core/painter/GausPainter.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/BondObject.hpp>
#include <rhbm_gem/data/object/DataObjectBase.hpp>
#include <rhbm_gem/data/object/PotentialEntryQuery.hpp>
#include <rhbm_gem/core/painter/PotentialPlotBuilder.hpp>
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

GausPainter::GausPainter() :
    m_folder_path{ "./" },
    m_atom_classifier{ std::make_unique<AtomClassifier>() },
    m_bond_classifier{ std::make_unique<BondClassifier>() }
{

}

GausPainter::~GausPainter()
{

}

void GausPainter::SetFolder(const std::string & folder_path)
{
    m_folder_path = FilePathHelper::EnsureTrailingSlash(folder_path);
}

void GausPainter::AddDataObject(DataObjectBase * data_object)
{
    painter_internal::AddDataObject<ModelObject>(
        data_object,
        m_ingest_mode,
        IngestMode::Data,
        IngestMode::Reference,
        m_ingest_label,
        "GausPainter",
        [this](ModelObject & typed_data_object) { IngestModelObject(typed_data_object); });
}

void GausPainter::AddReferenceDataObject(DataObjectBase * data_object, const std::string & label)
{
    painter_internal::AddReferenceDataObject<ModelObject>(
        data_object,
        label,
        m_ingest_mode,
        IngestMode::Reference,
        m_ingest_label,
        "GausPainter",
        [this](ModelObject & typed_data_object) { IngestModelObject(typed_data_object); });
}

void GausPainter::IngestModelObject(ModelObject & data_object)
{
    if (m_ingest_mode == IngestMode::Reference)
    {
        m_ref_model_object_list_map[m_ingest_label].push_back(&data_object);
        return;
    }
    m_model_object_list.push_back(&data_object);
}

void GausPainter::Painting()
{
    Logger::Log(LogLevel::Info, "GausPainter::Painting() called.");
    Logger::Log(LogLevel::Info, "Folder path: " + m_folder_path);
    Logger::Log(LogLevel::Info, "Number of atom objects to be painted: "
                + std::to_string(m_model_object_list.size()));

    for (auto model_object : m_model_object_list)
    {
        auto is_simulation{ model_object->GetEmdID().find("Simulation") != std::string::npos };
        auto label{ model_object->GetPdbID() };
        if (is_simulation == true)
        {
            label += "_"+ model_object->GetEmdID() +
                     "_bw"+ StringHelper::ToStringWithPrecision(model_object->GetResolution(), 2);
        }
        label += ".pdf";
        
        PaintAtomGroupGausAminoAcidMainChainComponentSimple(model_object, "atom_group_gaus_amino_acid_main_chain_component_simple_"+ label);
        PaintAtomLocalGausSummary(model_object, "atom_local_gaus_summary_"+ label);
        PaintAtomGroupGausSummary(model_object, "atom_group_gaus_summary_"+ label);
        PaintAtomGroupMapValueAminoAcidMainChainComponent(model_object, "atom_group_map_value_amino_acid_main_chain_component_"+ label);
        PaintAtomGroupGausAminoAcidMainChainComponent(model_object, "atom_group_gaus_amino_acid_main_chain_component_"+ label);
        PaintAtomGroupGausAminoAcidMainChainStructure(model_object, "atom_group_gaus_amino_acid_main_chain_structure_"+ label);
        PaintAtomLocalGausToSequenceAminoAcidMainChain(model_object, "atom_local_gaus_to_sequence_amino_acid_main_chain_"+ label);
    }
}

void GausPainter::PaintAtomLocalGausSummary(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, "GausPainter::PaintAtomLocalGausSummary");

    auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    const auto & chemical_component_map{ model_object->GetChemicalComponentEntryMap() };
    auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    auto show_outlier{ false };
    (void)chemical_component_map;
    (void)show_outlier;

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

    for (auto & [component_key, component_entry] : chemical_component_map)
    {
        const auto & component_atom_map{ component_entry->GetComponentAtomEntryMap() };
        auto component_id{ component_entry->GetComponentId() };
        for (auto & [atom_key, atom_entry] : component_atom_map)
        {
            if (atom_entry.element_type == Element::HYDROGEN) continue;
            if (atom_entry.atom_id == "OXT") continue;
            auto group_key{ AtomClassifier::GetGroupKeyInClass(component_key, atom_key) };
            if (entry_iter->IsAvailableAtomGroupKey(group_key, class_key) == false) continue;
            if (entry_iter->GetAtomObjectMap(group_key, class_key).size() < 1) continue;

            auto amplitude_hist{ plot_builder->CreateAtomGausEstimateHistogram(group_key, class_key, 0) };
            auto width_hist{ plot_builder->CreateAtomGausEstimateHistogram(group_key, class_key, 1) };

            canvas->cd();
            for (int i = 0; i < pad_size; i++)
            {
                ROOTHelper::SetPadDefaultStyle(pad[i].get());
                pad[i]->Draw();
            }

            pad[0]->cd();
            ROOTHelper::SetPadMarginInCanvas(gPad, 0.01, 0.01, 0.01, 0.01);
            auto resolution_text{ CreateResolutionPaveText(model_object) };
            ROOTHelper::SetPaveTextMarginInCanvas(gPad, resolution_text.get(), 0.83, 0.01, 0.02, 0.02);
            resolution_text->Draw();

            auto info_text{ CreateDataInfoPaveText(model_object) };
            ROOTHelper::SetPaveTextMarginInCanvas(gPad, info_text.get(), 0.59, 0.18, 0.02, 0.02);
            info_text->Draw();

            auto component_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
            ROOTHelper::SetPaveTextMarginInCanvas(gPad, component_text.get(), 0.01, 0.75, 0.02, 0.02);
            ROOTHelper::SetPaveTextDefaultStyle(component_text.get());
            ROOTHelper::SetPaveAttribute(component_text.get(), 0, 0.1);
            ROOTHelper::SetFillAttribute(component_text.get(), 1001, kAzure-7);
            ROOTHelper::SetTextAttribute(component_text.get(), 75.0f, 103, 22, 0.0, kYellow-10);
            component_text->AddText((component_id + "/" + atom_entry.atom_id).data());
            component_text->Draw();

            auto result_text{ ROOTHelper::CreatePaveText(0.26, 0.10, 0.55, 0.90, "nbNDC", false) };
            ROOTHelper::SetPaveTextDefaultStyle(result_text.get());
            ROOTHelper::SetTextAttribute(result_text.get(), 50.0f, 133, 12, 0.0, kRed);
            ROOTHelper::SetFillAttribute(result_text.get(), 4000);
            auto amplitude_prior{ entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 0) };
            auto amplitude_error{ entry_iter->GetAtomGausVariancePrior(group_key, class_key, 0) };
            auto width_prior{ entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 1) };
            auto width_error{ entry_iter->GetAtomGausVariancePrior(group_key, class_key, 1) };
            result_text->AddText(Form("#font[2]{#hat{A}} = %.2f #pm %.2f", amplitude_prior, amplitude_error));
            result_text->AddText(Form("#hat{#tau} = %.2f #pm %.2f", width_prior, width_error));
            result_text->Draw();

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
            frame_scatter->GetXaxis()->SetTitle("Local Amplitude");
            frame_scatter->GetYaxis()->SetTitle("Local Width");
            frame_scatter->GetXaxis()->SetLimits(
                amplitude_hist->GetXaxis()->GetXmin(),
                amplitude_hist->GetXaxis()->GetXmax()
            );
            frame_scatter->GetYaxis()->SetLimits(
                width_hist->GetXaxis()->GetXmin(),
                width_hist->GetXaxis()->GetXmax()
            );
            frame_scatter->GetXaxis()->CenterTitle();
            frame_scatter->GetYaxis()->CenterTitle();
            frame_scatter->SetStats(0);
            frame_scatter->Draw();
            auto scatter_graph{ plot_builder->CreateAtomGausEstimateScatterGraph(group_key, class_key, 0, 1) };
            auto outlier_scatter_graph{ plot_builder->CreateAtomGausEstimateScatterGraph(group_key, class_key, 0, 1, true) };
            ROOTHelper::SetMarkerAttribute(scatter_graph.get(), 24, 1.5, kAzure-7, 1.0f);
            ROOTHelper::SetMarkerAttribute(outlier_scatter_graph.get(), 5, 1.5, kRed+1, 1.0f);
            scatter_graph->Draw("P");
            if (show_outlier == true) outlier_scatter_graph->Draw("P");

            pad[3]->cd();
            ROOTHelper::SetPadMarginInCanvas(gPad, 0.09, 0.01, 0.02, 0.02);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            ROOTHelper::SetAxisTitleAttribute(amplitude_hist->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisTitleAttribute(amplitude_hist->GetYaxis(), 40.0, 1.3f);
            ROOTHelper::SetAxisLabelAttribute(amplitude_hist->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(amplitude_hist->GetYaxis(), 40.0, 0.01f);
            ROOTHelper::SetAxisTickAttribute(amplitude_hist->GetXaxis(), 0.0f, 505);
            ROOTHelper::SetAxisTickAttribute(amplitude_hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
            ROOTHelper::SetFillAttribute(amplitude_hist.get(), 1001, kAzure-7, 0.5f);
            amplitude_hist->GetYaxis()->SetTitle("Counts");
            amplitude_hist->GetYaxis()->CenterTitle();
            amplitude_hist->SetStats(0);
            amplitude_hist->Draw("BAR");

            pad[4]->cd();
            ROOTHelper::SetPadMarginInCanvas(gPad, 0.01, 0.01, 0.12, 0.02);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            ROOTHelper::SetAxisTitleAttribute(width_hist->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisTitleAttribute(width_hist->GetYaxis(), 40.0, 1.0f);
            ROOTHelper::SetAxisLabelAttribute(width_hist->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(width_hist->GetYaxis(), 40.0, 0.01f);
            ROOTHelper::SetAxisTickAttribute(width_hist->GetXaxis(), 0.0f, 505);
            ROOTHelper::SetAxisTickAttribute(width_hist->GetYaxis(), static_cast<float>(y_tick_length), 505);
            ROOTHelper::SetFillAttribute(width_hist.get(), 1001, kAzure-7, 0.5f);
            width_hist->GetYaxis()->SetTitle("Counts");
            width_hist->GetYaxis()->CenterTitle();
            width_hist->SetStats(0);
            width_hist->Draw("HBAR");

            pad[1]->cd();
            ROOTHelper::SetPadMarginInCanvas(gPad, 0.10, 0.00, 0.12, 0.10);
            auto member_size{ entry_iter->GetAtomObjectList(group_key, class_key).size() };
            auto outlier_size{ entry_iter->GetOutlierAtomObjectList(group_key, class_key).size() };
            std::vector<std::unique_ptr<TGraphErrors>> map_value_graph_list;
            std::vector<double> y_array;
            map_value_graph_list.reserve(member_size);
            y_array.reserve(member_size * 2);
            for (auto atom : entry_iter->GetAtomObjectList(group_key, class_key))
            {
                auto atom_iter{ std::make_unique<PotentialEntryQuery>(atom) };
                auto atom_plot_builder{ std::make_unique<PotentialPlotBuilder>(atom) };
                auto graph{ atom_plot_builder->CreateBinnedDistanceToMapValueGraph() };
                auto is_outlier{ atom_iter->IsOutlierAtom(class_key) };
                auto line_color{ kAzure-7 };
                if (show_outlier == true && is_outlier == true) line_color = kRed+1;
                ROOTHelper::SetLineAttribute(graph.get(), 1, 3, static_cast<short>(line_color), 0.3f);
                auto range_in_graph{ ROOTHelper::GetRangeInGraph(graph.get()) };
                y_array.emplace_back(std::get<0>(range_in_graph));
                y_array.emplace_back(std::get<1>(range_in_graph));
                map_value_graph_list.emplace_back(std::move(graph));
            }

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
            auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.1) };
            auto x_min_tmp{ 0.01 };
            auto x_max_tmp{ 1.49 };
            auto y_min_tmp{ std::get<0>(y_range) };
            auto y_max_tmp{ std::get<1>(y_range) };
            frame->GetXaxis()->SetLimits(x_min_tmp, x_max_tmp);
            frame->GetYaxis()->SetLimits(y_min_tmp, y_max_tmp);
            frame->SetStats(0);
            frame->Draw();
            for (auto & graph : map_value_graph_list)
            {
                graph->Draw("L X0");
            }

            auto gaus_mean{ plot_builder->CreateAtomGroupGausFunctionMean(group_key, class_key) };
            auto gaus_prior{ plot_builder->CreateAtomGroupGausFunctionPrior(group_key, class_key) };
            ROOTHelper::SetLineAttribute(gaus_prior.get(), 2, 3, kRed);
            ROOTHelper::SetLineAttribute(gaus_mean.get(), 3, 3, kBlue);
            gaus_prior->Draw("SAME");
            gaus_mean->Draw("SAME");

            auto legend{ ROOTHelper::CreateLegend(0.02, 0.90, 1.00, 1.00, false) };
            ROOTHelper::SetLegendDefaultStyle(legend.get());
            ROOTHelper::SetFillAttribute(legend.get(), 4000);
            ROOTHelper::SetTextAttribute(legend.get(), 40.0f, 133, 12, 0.0);
            legend->SetMargin(0.15f);
            legend->SetNColumns(3);
            //legend->AddEntry(gaus_prior.get(),
            //    "Gaussian Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})}", "l");
            //legend->AddEntry(map_value_graph_list.at(0).get(),
            //    "Members of Value", "l");
            legend->AddEntry(gaus_prior.get(),
                Form("#alpha_{r} = %.1f, #alpha_{g} = %.1f",
                    entry_iter->GetAtomAlphaR(group_key, class_key),
                    entry_iter->GetAtomAlphaG(group_key, class_key)), "l");
            legend->AddEntry(gaus_mean.get(),
                "#alpha_{r} = #alpha_{g} = 0", "l");
            legend->AddEntry(map_value_graph_list.at(0).get(),
                "Map Value", "l");
            legend->Draw();

            auto count_text{ ROOTHelper::CreatePaveText(0.45, 0.75, 1.00, 0.90, "nbNDC", false) };
            ROOTHelper::SetPaveTextDefaultStyle(count_text.get());
            ROOTHelper::SetTextAttribute(count_text.get(), 40.0f, 133, 32, 0.0, kGray+1);
            ROOTHelper::SetFillAttribute(count_text.get(), 4000);
            count_text->AddText(Form("Number of Members: %ld ", member_size));
            if (show_outlier == true)
            {
                count_text->AddText(Form("Group Outliers: #color[633]{%ld} ", outlier_size));
            }
            count_text->Draw();

            auto alpha_text{ ROOTHelper::CreatePaveText(0.20, 0.10, 0.40, 0.40, "nbNDC", false) };
            ROOTHelper::SetPaveTextDefaultStyle(alpha_text.get());
            ROOTHelper::SetTextAttribute(alpha_text.get(), 50.0f, 133, 12, 0.0, kRed);
            ROOTHelper::SetFillAttribute(alpha_text.get(), 4000);
            auto alpha_g{ entry_iter->GetAtomAlphaG(group_key, class_key) };
            alpha_text->AddText(Form("#alpha_{g} = %.1f", alpha_g));
            //alpha_text->Draw();
            
            ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
        }
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void GausPainter::PaintAtomGroupGausSummary(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, "GausPainter::PaintAtomGroupGausSummary");

    auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    const auto & chemical_component_map{ model_object->GetChemicalComponentEntryMap() };
    (void)chemical_component_map;
    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(2.0);
    gStyle->SetGridColor(kGray);
    gStyle->SetEndErrorSize(5.0f);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 750) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size{ 7 };
    
    std::unique_ptr<TPad> pad[pad_size];
    std::unique_ptr<TH2> frame[pad_size];

    pad[0] = ROOTHelper::CreatePad("pad0","", 0.00, 0.80, 1.00, 1.00); // The title pad
    pad[1] = ROOTHelper::CreatePad("pad1","", 0.00, 0.00, 0.55, 0.45); // The bottom-left pad
    pad[2] = ROOTHelper::CreatePad("pad2","", 0.00, 0.45, 0.55, 0.80); // The top-left pad
    pad[3] = ROOTHelper::CreatePad("pad3","", 0.55, 0.00, 0.90, 0.60); // The bottom-middle pad
    pad[4] = ROOTHelper::CreatePad("pad4","", 0.55, 0.60, 0.90, 0.80); // The top-middle pad
    pad[5] = ROOTHelper::CreatePad("pad5","", 0.90, 0.00, 1.00, 0.60); // The bottom-right pad
    pad[6] = ROOTHelper::CreatePad("pad6","", 0.90, 0.60, 1.00, 0.80); // The top-right pad

    frame[0] = ROOTHelper::CreateHist2D("hist_0","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[1] = ROOTHelper::CreateHist2D("hist_1","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[2] = ROOTHelper::CreateHist2D("hist_2","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[3] = ROOTHelper::CreateHist2D("hist_3","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[4] = ROOTHelper::CreateHist2D("hist_4","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    frame[5] = ROOTHelper::CreateHist2D("hist_5","", 100, 0.0, 1.0, 100, 0.0, 1.0);

    auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    auto struct_class_key{ ChemicalDataHelper::GetStructureAtomClassKey() };
    for (auto & [component_key, component_entry] : chemical_component_map)
    {
        if (component_entry->GetComponentType() == "non-polymer") continue;
        const auto & component_atom_map{ component_entry->GetComponentAtomEntryMap() };
        auto component_id{ component_entry->GetComponentId() };
        std::map<Element, std::map<std::string, GroupKey>> group_key_list_map;
        std::map<Element, std::map<std::string, GroupKey>> helix_group_key_list_map;
        std::map<Element, std::map<std::string, GroupKey>> sheet_group_key_list_map;
        std::vector<std::string> atom_id_list;
        atom_id_list.reserve(component_atom_map.size());
        GroupKey group_key_tmp{ 0 };
        GroupKey helix_group_key_tmp{ 0 };
        GroupKey sheet_group_key_tmp{ 0 };
        bool has_helix_component{ false };
        bool has_sheet_component{ false };
        for (auto & [atom_key, atom_entry] : component_atom_map)
        {
            if (atom_entry.element_type == Element::HYDROGEN) continue;
            if (atom_entry.atom_id == "OXT") continue;
            if (atom_entry.atom_id == "OP3") continue;
            auto group_key{ AtomClassifier::GetGroupKeyInClass(component_key, atom_key) };
            group_key_list_map[atom_entry.element_type].emplace(atom_entry.atom_id, group_key);
            if (group_key_tmp == 0) group_key_tmp = group_key;
            auto helix_group_key{
                AtomClassifier::GetGroupKeyInClass(Structure::HELX_P, component_key, atom_key)
            };
            helix_group_key_list_map[atom_entry.element_type].emplace(atom_entry.atom_id, helix_group_key);
            if (helix_group_key_tmp == 0) helix_group_key_tmp = helix_group_key;
            auto sheet_group_key{
                AtomClassifier::GetGroupKeyInClass(Structure::SHEET, component_key, atom_key)
            };
            sheet_group_key_list_map[atom_entry.element_type].emplace(atom_entry.atom_id, sheet_group_key);
            if (sheet_group_key_tmp == 0) sheet_group_key_tmp = sheet_group_key;
            atom_id_list.emplace_back(atom_entry.atom_id);
        }
        if (entry_iter->IsAvailableAtomGroupKey(group_key_tmp, class_key) == false) continue;
        if (entry_iter->IsAvailableAtomGroupKey(helix_group_key_tmp, struct_class_key) == true)
        {
            has_helix_component = true;
        }
        if (entry_iter->IsAvailableAtomGroupKey(sheet_group_key_tmp, struct_class_key) == true)
        {
            has_sheet_component = true;
        }

        std::map<Element, std::unique_ptr<TGraphErrors>> amplitude_graph_map;
        std::map<Element, std::unique_ptr<TGraphErrors>> width_graph_map;
        std::map<Element, std::unique_ptr<TGraphErrors>> correlation_graph_map;
        std::map<Element, std::unique_ptr<TGraphErrors>> helix_amplitude_graph_map;
        std::map<Element, std::unique_ptr<TGraphErrors>> sheet_amplitude_graph_map;
        std::map<Element, std::unique_ptr<TGraphErrors>> helix_width_graph_map;
        std::map<Element, std::unique_ptr<TGraphErrors>> sheet_width_graph_map;
        std::map<Element, std::unique_ptr<TH2D>> amplitude_hist_map;
        std::map<Element, std::unique_ptr<TH2D>> width_hist_map;
        std::vector<Element> element_list;
        std::vector<double> amplitude_array, width_array;
        amplitude_array.reserve(component_atom_map.size());
        width_array.reserve(component_atom_map.size());
        
        auto component_size{ entry_iter->GetAtomObjectList(group_key_tmp, class_key).size() };
        auto helix_component_size{ (has_helix_component == true) ?
            entry_iter->GetAtomObjectList(helix_group_key_tmp, struct_class_key).size() : 0
        };
        auto sheet_component_size{ (has_sheet_component == true) ?
            entry_iter->GetAtomObjectList(sheet_group_key_tmp, struct_class_key).size() : 0
        };
        for (auto & [element, group_key_list] : group_key_list_map)
        {
            auto amplitude_graph{
                plot_builder->CreateAtomGausEstimateToAtomIdGraph(
                    group_key_list, atom_id_list, class_key, 0)
            };
            auto width_graph{
                plot_builder->CreateAtomGausEstimateToAtomIdGraph(
                    group_key_list, atom_id_list, class_key, 1)
            };
            std::vector<GroupKey> group_key_list_vector;
            group_key_list_vector.reserve(group_key_list.size());
            for (auto & [atom_id, group_key] : group_key_list)
            {
                group_key_list_vector.emplace_back(group_key);
            }
            auto correlation_graph{
                plot_builder->CreateAtomGausEstimateScatterGraph(group_key_list_vector, class_key, 0, 1)
            };
            for (int p = 0; p < amplitude_graph->GetN(); p++)
            {
                amplitude_array.push_back(amplitude_graph->GetPointY(p));
                width_array.push_back(width_graph->GetPointY(p));
            }
            auto component_color{ ChemicalDataHelper::GetDisplayColor(element) };
            short component_marker{ 53 };
            ROOTHelper::SetMarkerAttribute(amplitude_graph.get(), component_marker, 1.5f, component_color);
            ROOTHelper::SetMarkerAttribute(width_graph.get(), component_marker, 1.5f, component_color);
            ROOTHelper::SetMarkerAttribute(correlation_graph.get(), component_marker, 1.5f, component_color);
            ROOTHelper::SetLineAttribute(amplitude_graph.get(), 1, 1, component_color);
            ROOTHelper::SetLineAttribute(width_graph.get(), 1, 1, component_color);
            ROOTHelper::SetLineAttribute(correlation_graph.get(), 1, 1, component_color);

            amplitude_graph_map[element] = std::move(amplitude_graph);
            width_graph_map[element] = std::move(width_graph);
            correlation_graph_map[element] = std::move(correlation_graph);
            element_list.emplace_back(element);

            if (has_helix_component == true)
            {
                auto helix_amplitude_graph{
                    plot_builder->CreateAtomGausEstimateToAtomIdGraph(
                        helix_group_key_list_map[element], atom_id_list, struct_class_key, 0)
                };
                auto helix_width_graph{
                    plot_builder->CreateAtomGausEstimateToAtomIdGraph(
                        helix_group_key_list_map[element], atom_id_list, struct_class_key, 1)
                };
                
                short helix_marker{ 20 };
                ROOTHelper::SetMarkerAttribute(helix_amplitude_graph.get(), helix_marker, 1.5f, component_color, 0.5f);
                ROOTHelper::SetMarkerAttribute(helix_width_graph.get(), helix_marker, 1.5f, component_color, 0.5f);
                ROOTHelper::SetLineAttribute(helix_amplitude_graph.get(), 1, 1, component_color, 0.5f);
                ROOTHelper::SetLineAttribute(helix_width_graph.get(), 1, 1, component_color, 0.5f);

                helix_amplitude_graph_map[element] = std::move(helix_amplitude_graph);
                helix_width_graph_map[element] = std::move(helix_width_graph);
            }

            if (has_sheet_component == true)
            {
                auto sheet_amplitude_graph{
                    plot_builder->CreateAtomGausEstimateToAtomIdGraph(
                        sheet_group_key_list_map[element], atom_id_list, struct_class_key, 0)
                };
                auto sheet_width_graph{
                    plot_builder->CreateAtomGausEstimateToAtomIdGraph(
                        sheet_group_key_list_map[element], atom_id_list, struct_class_key, 1)
                };
                
                short sheet_marker{ 52 };
                ROOTHelper::SetMarkerAttribute(sheet_amplitude_graph.get(), sheet_marker, 1.5f, component_color);
                ROOTHelper::SetMarkerAttribute(sheet_width_graph.get(), sheet_marker, 1.5f, component_color);
                ROOTHelper::SetLineAttribute(sheet_amplitude_graph.get(), 1, 1, component_color);
                ROOTHelper::SetLineAttribute(sheet_width_graph.get(), 1, 1, component_color);

                sheet_amplitude_graph_map[element] = std::move(sheet_amplitude_graph);
                sheet_width_graph_map[element] = std::move(sheet_width_graph);
            }
        }

        auto scaling{ 0.3 };
        auto amplitude_range{ ArrayStats<double>::ComputeScalingRangeTuple(amplitude_array, scaling) };
        auto width_range{ ArrayStats<double>::ComputeScalingRangeTuple(width_array, scaling) };
        auto element_count{ element_list.size() };
        std::vector<std::string> element_label_list;
        element_label_list.reserve(element_count);
        for (size_t i = 0; i < element_count; i++)
        {
            auto element{ element_list.at(i) };
            auto element_label{ ChemicalDataHelper::GetLabel(element) };
            auto x_value{ static_cast<double>(i) };
            element_label_list.emplace_back(element_label);
            std::string name_amplitude{ "amplitude_hist_"+ element_label };
            std::string name_width{ "width_hist_"+ element_label };
            auto amplitude_hist{
                ROOTHelper::CreateHist2D(
                    name_amplitude.data(),"",
                    100, std::get<0>(amplitude_range), std::get<1>(amplitude_range),
                    static_cast<int>(element_count), -0.5, static_cast<int>(element_count)-0.5)
            };
            auto width_hist{
                ROOTHelper::CreateHist2D(
                    name_width.data(),"",
                    static_cast<int>(element_count), -0.5, static_cast<int>(element_count)-0.5,
                    100, std::get<0>(width_range), std::get<1>(width_range))
            };
            for (int p = 0; p < amplitude_graph_map.at(element)->GetN(); p++)
            {
                amplitude_hist->Fill(amplitude_graph_map.at(element)->GetPointY(p), x_value);
            }
            for (int p = 0; p < width_graph_map.at(element)->GetN(); p++)
            {
                width_hist->Fill(x_value, width_graph_map.at(element)->GetPointY(p));
            }
            auto element_color{ ChemicalDataHelper::GetDisplayColor(element) };
            ROOTHelper::SetLineAttribute(amplitude_hist.get(), 1, 1, element_color);
            ROOTHelper::SetLineAttribute(width_hist.get(), 1, 1, element_color);
            ROOTHelper::SetFillAttribute(amplitude_hist.get(), 1001, element_color, 0.3f);
            ROOTHelper::SetFillAttribute(width_hist.get(), 1001, element_color, 0.3f);
            amplitude_hist->SetBarWidth(0.6f);
            width_hist->SetBarWidth(0.6f);
            amplitude_hist_map[element] = std::move(amplitude_hist);
            width_hist_map[element] = std::move(width_hist);
        }

        canvas->cd();
        for (int i = 0; i < pad_size; i++)
        {
            ROOTHelper::SetPadDefaultStyle(pad[i].get());
            pad[i]->Draw();
        }

        pad[0]->cd();
        auto component_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, component_text.get(), 0.01, 0.85, 0.02, 0.02);
        ROOTHelper::SetPaveTextDefaultStyle(component_text.get());
        ROOTHelper::SetPaveAttribute(component_text.get(), 0, 0.1);
        ROOTHelper::SetFillAttribute(component_text.get(), 1001, kAzure-7);
        ROOTHelper::SetTextAttribute(component_text.get(), 75.0f, 103, 22, 0.0, kYellow-10);
        component_text->AddText(component_entry->GetComponentId().data());
        component_text->Draw();

        auto component_info_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC", false) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, component_info_text.get(), 0.15, 0.53, 0.02, 0.02);
        ROOTHelper::SetPaveTextDefaultStyle(component_info_text.get());
        ROOTHelper::SetPaveAttribute(component_info_text.get(), 0);
        ROOTHelper::SetFillAttribute(component_info_text.get(), 4000);
        ROOTHelper::SetTextAttribute(component_info_text.get(), 30.0f, 133);
        auto component_name{ component_entry->GetComponentName() };
        StringHelper::EraseCharFromString(component_name, '\"');
        component_info_text->AddText(component_name.data());
        component_info_text->AddText(("Formula: " + component_entry->GetComponentFormula()).data());
        component_info_text->AddText(Form("Number of members: %zu (#alpha-helix: %zu, #beta-sheet: %zu)", component_size, helix_component_size, sheet_component_size));
        component_info_text->Draw();

        auto info_text{ CreateDataInfoPaveText(model_object) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, info_text.get(), 0.59, 0.18, 0.02, 0.02);
        info_text->Draw();

        auto resolution_text{ CreateResolutionPaveText(model_object) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, resolution_text.get(), 0.83, 0.01, 0.02, 0.02);
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
        RemodelAxisLabels(frame[0]->GetXaxis(), atom_id_list, 90.0, 12);
        ROOTHelper::SetAxisTitleAttribute(frame[0]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisTitleAttribute(frame[0]->GetYaxis(), 50.0f, 1.5f);
        ROOTHelper::SetAxisLabelAttribute(frame[0]->GetXaxis(), 40.0f, 0.11f, 103, kCyan+3);
        ROOTHelper::SetAxisLabelAttribute(frame[0]->GetYaxis(), 45.0f, 0.01f);
        frame[0]->GetYaxis()->SetTitle("Width");
        frame[0]->Draw();
        for (auto & [element, graph] : width_graph_map) graph->Draw("P");
        for (auto & [element, graph] : helix_width_graph_map) graph->Draw("P X0");
        for (auto & [element, graph] : sheet_width_graph_map) graph->Draw("P X0");

        pad[2]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.10, 0.00, 0.02, 0.01);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        RemodelFrameInPad(frame[1].get(), pad[2].get(), 0.0, 0.015);
        RemodelAxisLabels(frame[1]->GetXaxis(), atom_id_list, 90.0, 12);
        ROOTHelper::SetAxisTitleAttribute(frame[1]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisTitleAttribute(frame[1]->GetYaxis(), 50.0f, 1.5f);
        ROOTHelper::SetAxisLabelAttribute(frame[1]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(frame[1]->GetYaxis(), 45.0f, 0.01f);
        frame[1]->GetYaxis()->SetTitle("Amplitude");
        frame[1]->Draw();
        for (auto & [element, graph] : amplitude_graph_map) graph->Draw("P");
        for (auto & [element, graph] : helix_amplitude_graph_map) graph->Draw("P X0");
        for (auto & [element, graph] : sheet_amplitude_graph_map) graph->Draw("P X0");

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
        RemodelAxisLabels(frame[3]->GetYaxis(), element_label_list, 0.0, -1);
        ROOTHelper::SetAxisTitleAttribute(frame[3]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisTitleAttribute(frame[3]->GetYaxis(), 40.0f, 1.2f);
        ROOTHelper::SetAxisLabelAttribute(frame[3]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(frame[3]->GetYaxis(), 35.0f, 0.01f, 103);
        frame[3]->GetYaxis()->SetTitle("Element");
        frame[3]->Draw();
        for (auto & [element, hist] : amplitude_hist_map) hist->Draw("CANDLEY3 SAME");

        pad[5]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.01, 0.01, 0.12, 0.02);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        RemodelFrameInPad(frame[2].get(), pad[5].get(), 0.0, 0.0);
        RemodelAxisLabels(frame[2]->GetXaxis(), element_label_list, 0.0, -1);
        ROOTHelper::SetAxisTitleAttribute(frame[2]->GetXaxis(), 40.0f, 1.0f);
        ROOTHelper::SetAxisTitleAttribute(frame[2]->GetYaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(frame[2]->GetXaxis(), 35.0f, 0.005f, 103);
        ROOTHelper::SetAxisLabelAttribute(frame[2]->GetYaxis(), 0.0f);
        frame[2]->GetXaxis()->SetTitle("Element");
        frame[2]->Draw();
        for (auto & [element, hist] : width_hist_map) hist->Draw("CANDLE3 SAME");

        pad[6]->cd();
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.02, 0.02, 0.20, 0.20);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        auto legend{ ROOTHelper::CreateLegend(0.00, 0.00, 1.00, 1.00, false) };
        ROOTHelper::SetLegendDefaultStyle(legend.get());
        ROOTHelper::SetTextAttribute(legend.get(), 30.0f, 133, 12);
        short general_color{ kGray+2 };
        auto merge_marker{ std::make_unique<TMarker>(0.0, 0.0, 1) };
        auto helix_marker{ std::make_unique<TMarker>(0.0, 0.0, 1) };
        auto sheet_marker{ std::make_unique<TMarker>(0.0, 0.0, 1) };
        ROOTHelper::SetMarkerAttribute(merge_marker.get(), 53, 1.5f, general_color);
        ROOTHelper::SetMarkerAttribute(helix_marker.get(), 20, 1.5f, general_color);
        ROOTHelper::SetMarkerAttribute(sheet_marker.get(), 52, 1.5f, general_color);
        legend->AddEntry(merge_marker.get(), "Merge", "pe");
        if (has_helix_component == true) legend->AddEntry(helix_marker.get(), "#alpha-helix", "p");
        if (has_sheet_component == true) legend->AddEntry(sheet_marker.get(), "#beta-sheet", "p");
        legend->SetMargin(0.30f);
        legend->Draw();
        
        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void GausPainter::PaintAtomGroupMapValueAminoAcidMainChainComponent(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, "GausPainter::PaintAtomGroupMapValueAminoAcidMainChainComponent");

    auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    const auto & chemical_component_map{ model_object->GetChemicalComponentEntryMap() };
    const auto & class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    const std::vector<Spot> spot_list{ Spot::CA, Spot::C, Spot::N, Spot::O };
    const auto & standard_residue_list{ ChemicalDataHelper::GetStandardAminoAcidList() };
    bool show_outlier{ false };
    (void)chemical_component_map;
    (void)class_key;
    (void)standard_residue_list;
    (void)show_outlier;

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    const int col_size{ 5 };
    const int row_size{ 4 };

    auto canvas{ ROOTHelper::CreateCanvas("test","", 2000, 1500) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.10f, 0.02f, 0.10f, 0.12f, 0.01f, 0.01f
    );

    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);
    for (size_t k = 0; k < spot_list.size(); k++)
    {
        auto group_key_list{ m_atom_classifier->GetMainChainComponentAtomClassGroupKeyList(k) };
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

        std::map<Residue, std::vector<std::unique_ptr<TGraphErrors>>> map_value_graph_list_map;
        std::map<Residue, std::unique_ptr<TF1>> gaus_prior_map;
        std::map<Residue, std::unique_ptr<TF1>> gaus_mean_map;
        std::map<Residue, std::string> component_id_map;
        std::vector<double> global_y_array;
        for (auto & group_key : group_key_list)
        {
            auto residue{ entry_iter->GetResidueFromAtomGroupKey(group_key, class_key) };
            auto component_key{ static_cast<ComponentKey>(residue) };
            auto component_id{ chemical_component_map.at(component_key)->GetComponentId() };
            component_id_map.emplace(residue, component_id);

            auto gaus_mean{ plot_builder->CreateAtomGroupGausFunctionMean(group_key, class_key) };
            auto gaus_prior{ plot_builder->CreateAtomGroupGausFunctionPrior(group_key, class_key) };
            gaus_mean_map.emplace(residue, std::move(gaus_mean));
            gaus_prior_map.emplace(residue, std::move(gaus_prior));

            auto member_size{ entry_iter->GetAtomObjectList(group_key, class_key).size() };
            std::vector<std::unique_ptr<TGraphErrors>> map_value_graph_list;
            std::vector<double> y_array;
            map_value_graph_list.reserve(member_size);
            y_array.reserve(member_size * 2);
            for (auto atom : entry_iter->GetAtomObjectList(group_key, class_key))
            {
                auto atom_iter{ std::make_unique<PotentialEntryQuery>(atom) };
                auto atom_plot_builder{ std::make_unique<PotentialPlotBuilder>(atom) };
                auto graph{ atom_plot_builder->CreateBinnedDistanceToMapValueGraph() };
                auto is_outlier{ atom_iter->IsOutlierAtom(class_key) };
                auto line_color{ kAzure-7 };
                if (show_outlier == true && is_outlier == true) line_color = kRed+1;
                ROOTHelper::SetLineAttribute(graph.get(), 1, 3, static_cast<short>(line_color), 0.3f);
                auto range_in_graph{ ROOTHelper::GetRangeInGraph(graph.get()) };
                y_array.emplace_back(std::get<0>(range_in_graph));
                y_array.emplace_back(std::get<1>(range_in_graph));
                map_value_graph_list.emplace_back(std::move(graph));
            }
            global_y_array.insert(global_y_array.end(), y_array.begin(), y_array.end());
            map_value_graph_list_map.emplace(residue, std::move(map_value_graph_list));
        }

        double y_min{ 0.0 };
        double y_max{ 0.0 };
        auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(global_y_array, 0.20) };
        y_min = std::get<0>(y_range);
        y_max = std::get<1>(y_range);

        std::unique_ptr<TH2> frame[col_size][row_size];
        std::unique_ptr<TPaveText> info_text[col_size][row_size];
        std::unique_ptr<TPaveText> result_text[col_size][row_size];
        for (int i = 0; i < col_size; i++)
        {
            for (int j = 0; j < row_size; j++)
            {
                auto par_id{ row_size - j - 1 };
                auto residue_id{ par_id * col_size + i };
                auto residue{ standard_residue_list.at(static_cast<size_t>(residue_id)) };

                ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
                ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
                ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
                auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
                auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
                frame[i][j] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, j),"", 500, 0.01, 1.49, 500, y_min, y_max);
                ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 60.0f, 0.01f, 133);
                ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.08/x_factor), 505);
                ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 60.0f, 0.01f, 133);
                ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
                ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
                frame[i][j]->GetXaxis()->SetTitle("");
                frame[i][j]->GetYaxis()->SetTitle("");
                frame[i][j]->SetStats(0);
                frame[i][j]->Draw("");

                const auto & map_value_graph_list{ map_value_graph_list_map.at(residue) };
                for (auto & graph : map_value_graph_list)
                {
                    graph->Draw("L X0");
                }

                info_text[i][j] = ROOTHelper::CreatePaveText(0.65, 0.80, 0.99, 0.99, "nbNDC ARC", true);
                ROOTHelper::SetPaveTextDefaultStyle(info_text[i][j].get());
                ROOTHelper::SetPaveAttribute(info_text[i][j].get(), 0, 0.2);
                ROOTHelper::SetLineAttribute(info_text[i][j].get(), 1, 0);
                ROOTHelper::SetTextAttribute(info_text[i][j].get(), 50.0f, 103, 22);
                ROOTHelper::SetFillAttribute(info_text[i][j].get(), 1001, kAzure-7, 0.5f);
                info_text[i][j]->AddText(component_id_map.at(residue).data());
                info_text[i][j]->Draw();

                auto & gaus_mean{ gaus_mean_map.at(residue) };
                auto & gaus_prior{ gaus_prior_map.at(residue) };
                ROOTHelper::SetLineAttribute(gaus_mean.get(), 3, 3, kBlue);
                ROOTHelper::SetLineAttribute(gaus_prior.get(), 2, 3, kRed);
                gaus_prior->Draw("SAME");
                gaus_mean->Draw("SAME");

                result_text[i][j] = ROOTHelper::CreatePaveText(0.11, 0.10, 0.95, 0.20, "nbNDC", true);
                ROOTHelper::SetPaveTextDefaultStyle(result_text[i][j].get());
                ROOTHelper::SetLineAttribute(result_text[i][j].get(), 1, 0);
                ROOTHelper::SetTextAttribute(result_text[i][j].get(), 30.0f, 133, 22, 0.0f, kRed);
                ROOTHelper::SetFillAttribute(result_text[i][j].get(), 4000);
                result_text[i][j]->AddText(
                    Form("A = %.2f,   #tau = %.2f", gaus_prior->GetParameter(0), gaus_prior->GetParameter(1))
                );
                result_text[i][j]->Draw();
            }
        }

        canvas->cd();
        auto pad_extra_0{ ROOTHelper::CreatePad("pad_extra_0","", 0.00, 0.88, 1.00, 1.00) };
        pad_extra_0->Draw();
        pad_extra_0->cd();
        ROOTHelper::SetPadDefaultStyle(pad_extra_0.get());
        ROOTHelper::SetFillAttribute(pad_extra_0.get(), 4000);
        ROOTHelper::SetPadMarginInCanvas(gPad, 0.01, 0.01, 0.01, 0.01);

        auto spot{ spot_list.at(k) };
        auto spot_color{ AtomClassifier::GetMainChainSpotColor(spot) };
        auto spot_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC ARC", false) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, spot_text.get(), 0.01, 0.91, 0.02, 0.01);
        ROOTHelper::SetPaveTextDefaultStyle(spot_text.get());
        ROOTHelper::SetFillAttribute(spot_text.get(), 1001, spot_color, 0.3f);
        ROOTHelper::SetLineAttribute(spot_text.get(), 1, 3, spot_color);
        ROOTHelper::SetTextAttribute(spot_text.get(), 90.0f, 133, 22);
        spot_text->AddText(AtomClassifier::GetMainChainSpotLabel(spot).data());
        spot_text->Draw();

        auto resolution_text{ CreateResolutionPaveText(model_object) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, resolution_text.get(), 0.10, 0.76, 0.02, 0.01);
        resolution_text->SetTextSize(90.0f);
        resolution_text->Draw();

        auto map_info_text{ CreateDataInfoPaveText(model_object) };
        ROOTHelper::SetPaveTextMarginInCanvas(gPad, map_info_text.get(), 0.25, 0.53, 0.02, 0.01);
        map_info_text->SetTextSize(70.0f);
        map_info_text->Draw();

        auto legend{ ROOTHelper::CreateLegend(0.50, 0.02, 1.00, 0.99, false) };
        ROOTHelper::SetLegendDefaultStyle(legend.get());
        ROOTHelper::SetFillAttribute(legend.get(), 4000);
        ROOTHelper::SetTextAttribute(legend.get(), 40.0f, 133, 12, 0.0);
        legend->SetMargin(0.25f);
        legend->AddEntry(gaus_prior_map.at(Residue::ALA).get(),
            Form("Gaussian Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})} with #alpha_{r} = %.1f, #alpha_{g} = %.1f",
            entry_iter->GetAtomAlphaR(group_key_list.at(0), class_key), entry_iter->GetAtomAlphaG(group_key_list.at(0), class_key)), "l");
        legend->AddEntry(gaus_mean_map.at(Residue::ALA).get(),
            "Gaussian Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})} with #alpha_{r} = #alpha_{g} = 0", "l");
        legend->AddEntry(map_value_graph_list_map.at(Residue::ALA).front().get(),
            "Members of Map Value", "l");
        legend->Draw();

        canvas->cd();
        auto pad_extra_1{ ROOTHelper::CreatePad("pad_extra_1","", 0.00, 0.10, 0.05, 0.88) };
        pad_extra_1->Draw();
        pad_extra_1->cd();
        ROOTHelper::SetPadDefaultStyle(pad_extra_1.get());
        ROOTHelper::SetFillAttribute(pad_extra_1.get(), 4000);
        auto y_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
        ROOTHelper::SetPaveTextDefaultStyle(y_title_text.get());
        ROOTHelper::SetFillAttribute(y_title_text.get(), 4000);
        ROOTHelper::SetTextAttribute(y_title_text.get(), 70.0f, 133, 22);
        y_title_text->AddText("Normalized Map Value");
        auto y_text{ y_title_text->GetLineWith("Normalized") };
        y_text->SetTextAngle(90.0f);
        y_title_text->Draw();

        canvas->cd();
        auto pad_extra_2{ ROOTHelper::CreatePad("pad_extra_2","", 0.10, 0.00, 0.98, 0.05) };
        pad_extra_2->Draw();
        pad_extra_2->cd();
        ROOTHelper::SetPadDefaultStyle(pad_extra_2.get());
        ROOTHelper::SetFillAttribute(pad_extra_2.get(), 4000);
        auto x_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
        ROOTHelper::SetPaveTextDefaultStyle(x_title_text.get());
        ROOTHelper::SetFillAttribute(x_title_text.get(), 4000);
        ROOTHelper::SetTextAttribute(x_title_text.get(), 70.0f, 133, 22);
        x_title_text->AddText("Radial Distance #[]{#AA}");
        x_title_text->Draw();

        ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    }

    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void GausPainter::PaintAtomGroupGausAminoAcidMainChainComponent(
    ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, "GausPainter::PaintAtomGroupGausAminoAcidMainChainComponent");

    auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    const auto & chemical_component_map{ model_object->GetChemicalComponentEntryMap() };
    //const std::vector<Spot> spot_list{ Spot::CA, Spot::C, Spot::N, Spot::O };
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
