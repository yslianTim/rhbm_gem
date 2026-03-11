#include <rhbm_gem/core/painter/ModelPainter.hpp>
#include <rhbm_gem/core/painter/PotentialPlotBuilder.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/BondClassifier.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/PotentialEntryQuery.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#endif

namespace rhbm_gem {

void ModelPainter::PaintAtomMapValueMainChain(ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintAtomMapValueMainChain");
    (void)model_object;

    #ifdef HAVE_ROOT
    const auto & class_key{ ChemicalDataHelper::GetSimpleAtomClassKey() };
    auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    const int col_size{ 4 };
    const int row_size{ 1 };
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 350) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.10f, 0.02f, 0.20f, 0.20f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const size_t main_chain_element_size{ 4 };

    std::vector<std::unique_ptr<TGraphErrors>> map_value_graph_list[main_chain_element_size];
    std::unique_ptr<TF1> gaus_function[main_chain_element_size];
    double amplitude_prior[main_chain_element_size];
    double width_prior[main_chain_element_size];
    std::vector<double> y_array;
    y_array.reserve(model_object->GetNumberOfSelectedAtom());
    AtomClassifier classifier;
    for (size_t k = 0; k < main_chain_element_size; k++)
    {
        auto group_key{ classifier.GetMainChainSimpleAtomClassGroupKey(k) };
        if (entry_iter->IsAvailableAtomGroupKey(group_key, class_key) == false) continue;
        for (auto atom : entry_iter->GetAtomObjectList(group_key, class_key))
        {
            auto atom_iter{ std::make_unique<PotentialEntryQuery>(atom) };
            auto atom_plot_builder{ std::make_unique<PotentialPlotBuilder>(atom) };
            auto graph{ atom_plot_builder->CreateBinnedDistanceToMapValueGraph() };
            ROOTHelper::SetLineAttribute(graph.get(), 1, 2, static_cast<short>(kAzure-7), 0.3f);
            map_value_graph_list[k].emplace_back(std::move(graph));
            auto map_value_range{ atom_iter->GetMapValueRange(0.0) };
            y_array.emplace_back(std::get<0>(map_value_range));
            y_array.emplace_back(std::get<1>(map_value_range));
        }
        gaus_function[k] = plot_builder->CreateAtomGroupGausFunctionPrior(group_key, class_key);
        amplitude_prior[k] = entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 0);
        width_prior[k] = entry_iter->GetAtomGausEstimatePrior(group_key, class_key, 1);
    }

    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.15) };
    auto x_min{ 0.01 };
    auto x_max{ 1.49 };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };
    auto line_ref{ ROOTHelper::CreateLine(x_min, 0.0, x_max, 0.0) };

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> element_text[col_size];
    std::unique_ptr<TPaveText> result_text[col_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 100, x_min, x_max, 100, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 35.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 505);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 40.0f, 1.2f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35.0f, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->SetStats(0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("Map Value");
            frame[i][j]->Draw();

            for (auto & graph : map_value_graph_list[i])
            {
                graph->Draw("L X0");
            }

            ROOTHelper::SetLineAttribute(line_ref.get(), 2, 1, kBlack);
            line_ref->Draw();

            ROOTHelper::SetLineAttribute(gaus_function[i].get(), 2, 2, kRed);
            gaus_function[i]->Draw("SAME");

            element_text[i] = ROOTHelper::CreatePaveText(0.70, 0.75, 0.99, 0.99, "nbNDC ARC", true);
            ROOTHelper::SetPaveTextDefaultStyle(element_text[i].get());
            ROOTHelper::SetPaveAttribute(element_text[i].get(), 0, 0.2);
            ROOTHelper::SetTextAttribute(element_text[i].get(), 40.0f, 103, 22);
            ROOTHelper::SetFillAttribute(element_text[i].get(), 1001, AtomClassifier::GetMainChainElementColor(static_cast<size_t>(i)), 0.5f);
            element_text[i]->AddText(AtomClassifier::GetMainChainElementLabel(static_cast<size_t>(i)).data());
            element_text[i]->Draw();

            result_text[i] = ROOTHelper::CreatePaveText(0.05, 0.08, 0.95, 0.25, "nbNDC", true);
            ROOTHelper::SetPaveTextDefaultStyle(result_text[i].get());
            ROOTHelper::SetTextAttribute(result_text[i].get(), 20.0f, 133, 22, 0.0, kRed);
            ROOTHelper::SetFillAttribute(result_text[i].get(), 4000);
            result_text[i]->AddText(Form("#font[2]{A} = %.2f  ;  #tau = %.2f", amplitude_prior[i], width_prior[i]));
            result_text[i]->Draw();
        }
    }

    canvas->cd();
    auto pad_title{ ROOTHelper::CreatePad("pad_title","", 0.00, 0.80, 1.00, 1.00) };
    ROOTHelper::SetPadDefaultStyle(pad_title.get());
    pad_title->Draw();
    pad_title->cd();

    auto subtitle1_text{ ROOTHelper::CreatePaveText(0.10, 0.10, 0.25, 0.90, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(subtitle1_text.get());
    ROOTHelper::SetPaveAttribute(subtitle1_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle1_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(subtitle1_text.get(), 40.0f, 133, 22, 0.0, kYellow-10);
    subtitle1_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
    subtitle1_text->Draw();

    auto subtitle2_text{ ROOTHelper::CreatePaveText(0.26, 0.10, 0.45, 0.90, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(subtitle2_text.get());
    ROOTHelper::SetPaveAttribute(subtitle2_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle2_text.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetTextAttribute(subtitle2_text.get(), 35.0f, 103, 22);
    subtitle2_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
    subtitle2_text->Draw();

    auto subtitle3_text{ ROOTHelper::CreatePaveText(0.46, 0.10, 0.68, 0.90, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(subtitle3_text.get());
    ROOTHelper::SetPaveAttribute(subtitle3_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle3_text.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetTextAttribute(subtitle3_text.get(), 35.0f, 103, 22);
    subtitle3_text->AddText(model_object->GetEmdID().data());
    subtitle3_text->Draw();

    auto legend{ ROOTHelper::CreateLegend(0.69, 0.10, 0.99, 0.90, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    ROOTHelper::SetTextAttribute(legend.get(), 25.0f, 133, 12, 0.0);
    legend->SetMargin(0.15f);
    legend->AddEntry(gaus_function[0].get(),
        "Gaussian Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})}", "l");
    legend->AddEntry(map_value_graph_list[0].at(0).get(),
        "Map Value", "l");
    legend->Draw();

    canvas->cd();
    auto pad_bottom{ ROOTHelper::CreatePad("pad_bottom","", 0.10, 0.00, 0.98, 0.11) };
    ROOTHelper::SetPadDefaultStyle(pad_bottom.get());
    pad_bottom->Draw();
    pad_bottom->cd();
    auto x_title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(x_title_text.get());
    ROOTHelper::SetFillAttribute(x_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(x_title_text.get(), 35.0f, 133, 22);
    x_title_text->AddText("Radial Distance #[]{#AA}");
    x_title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void ModelPainter::PaintBondMapValueMainChain(ModelObject * model_object, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " ModelPainter::PaintBondMapValueMainChain");
    (void)model_object;

    #ifdef HAVE_ROOT
    const auto & class_key{ ChemicalDataHelper::GetSimpleBondClassKey() };
    auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };
    const int col_size{ 4 };
    const int row_size{ 1 };
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1000, 350) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.10f, 0.02f, 0.20f, 0.20f, 0.01f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    const size_t main_chain_element_size{ 4 };

    std::vector<std::unique_ptr<TGraphErrors>> map_value_graph_list[main_chain_element_size];
    std::unique_ptr<TF1> gaus_function[main_chain_element_size];
    double amplitude_prior[main_chain_element_size];
    double width_prior[main_chain_element_size];
    int member_entries[main_chain_element_size];
    std::vector<double> y_array;
    y_array.reserve(model_object->GetNumberOfSelectedBond());
    BondClassifier classifier;
    for (size_t k = 0; k < main_chain_element_size; k++)
    {
        auto group_key{ classifier.GetMainChainSimpleBondClassGroupKey(k) };
        if (entry_iter->IsAvailableBondGroupKey(group_key, class_key, true) == false) continue;

        for (auto bond : entry_iter->GetBondObjectList(group_key, class_key))
        {
            auto bond_iter{ std::make_unique<PotentialEntryQuery>(bond) };
            auto bond_plot_builder{ std::make_unique<PotentialPlotBuilder>(bond) };
            auto graph{ bond_plot_builder->CreateBinnedDistanceToMapValueGraph() };
            ROOTHelper::SetLineAttribute(graph.get(), 1, 2, static_cast<short>(kAzure-7), 0.3f);
            map_value_graph_list[k].emplace_back(std::move(graph));
            auto map_value_range{ bond_iter->GetMapValueRange(0.0) };
            y_array.emplace_back(std::get<0>(map_value_range));
            y_array.emplace_back(std::get<1>(map_value_range));
        }
        gaus_function[k] = plot_builder->CreateBondGroupGausFunctionPrior(group_key, class_key);
        amplitude_prior[k] = entry_iter->GetBondGausEstimatePrior(group_key, class_key, 0);
        width_prior[k] = entry_iter->GetBondGausEstimatePrior(group_key, class_key, 1);
        member_entries[k] = static_cast<int>(entry_iter->GetBondObjectList(group_key, class_key).size());
    }

    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.15) };
    auto x_min{ 0.01 };
    auto x_max{ 1.49 };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };
    auto line_ref{ ROOTHelper::CreateLine(x_min, 0.0, x_max, 0.0) };

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> element_text[col_size];
    std::unique_ptr<TPaveText> result_text[col_size];
    std::unique_ptr<TPaveText> count_text[col_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", i, j),"", 100, x_min, x_max, 100, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 35.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 505);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 40.0f, 1.2f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 35.0f, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->SetStats(0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("Map Value");
            frame[i][j]->Draw();

            for (auto & graph : map_value_graph_list[i])
            {
                graph->Draw("L X0");
            }

            ROOTHelper::SetLineAttribute(line_ref.get(), 2, 1, kBlack);
            line_ref->Draw();

            ROOTHelper::SetLineAttribute(gaus_function[i].get(), 2, 2, kRed);
            gaus_function[i]->Draw("SAME");

            element_text[i] = ROOTHelper::CreatePaveText(0.55, 0.75, 0.99, 0.99, "nbNDC ARC", true);
            ROOTHelper::SetPaveTextDefaultStyle(element_text[i].get());
            ROOTHelper::SetPaveAttribute(element_text[i].get(), 0, 0.2);
            ROOTHelper::SetTextAttribute(element_text[i].get(), 40.0f, 103, 22);
            ROOTHelper::SetFillAttribute(element_text[i].get(), 1001, BondClassifier::GetMainChainMemberColor(static_cast<size_t>(i)), 0.5f);
            element_text[i]->AddText(BondClassifier::GetMainChainMemberLabel(static_cast<size_t>(i)).data());
            element_text[i]->Draw();

            count_text[i] = ROOTHelper::CreatePaveText(0.03, 0.90, 0.54, 0.99, "nbNDC", true);
            ROOTHelper::SetPaveTextDefaultStyle(count_text[i].get());
            ROOTHelper::SetTextAttribute(count_text[i].get(), 15.0f, 133, 12, 0.0, kGray+2);
            ROOTHelper::SetFillAttribute(count_text[i].get(), 4000);
            count_text[i]->AddText(Form("#Entries = %d", member_entries[i]));
            count_text[i]->Draw();

            result_text[i] = ROOTHelper::CreatePaveText(0.05, 0.08, 0.95, 0.25, "nbNDC", true);
            ROOTHelper::SetPaveTextDefaultStyle(result_text[i].get());
            ROOTHelper::SetTextAttribute(result_text[i].get(), 20.0f, 133, 22, 0.0, kRed);
            ROOTHelper::SetFillAttribute(result_text[i].get(), 4000);
            result_text[i]->AddText(Form("#font[2]{A} = %.2f  ;  #tau = %.2f", amplitude_prior[i], width_prior[i]));
            result_text[i]->Draw();
        }
    }

    canvas->cd();
    auto pad_title{ ROOTHelper::CreatePad("pad_title","", 0.00, 0.80, 1.00, 1.00) };
    ROOTHelper::SetPadDefaultStyle(pad_title.get());
    pad_title->Draw();
    pad_title->cd();

    auto subtitle1_text{ ROOTHelper::CreatePaveText(0.10, 0.10, 0.25, 0.90, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(subtitle1_text.get());
    ROOTHelper::SetPaveAttribute(subtitle1_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle1_text.get(), 1001, kAzure-7);
    ROOTHelper::SetTextAttribute(subtitle1_text.get(), 40.0f, 133, 22, 0.0, kYellow-10);
    subtitle1_text->AddText(Form("%.2f #AA", model_object->GetResolution()));
    subtitle1_text->Draw();

    auto subtitle2_text{ ROOTHelper::CreatePaveText(0.26, 0.10, 0.45, 0.90, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(subtitle2_text.get());
    ROOTHelper::SetPaveAttribute(subtitle2_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle2_text.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetTextAttribute(subtitle2_text.get(), 35.0f, 103, 22);
    subtitle2_text->AddText(Form("PDB-%s", model_object->GetPdbID().data()));
    subtitle2_text->Draw();

    auto subtitle3_text{ ROOTHelper::CreatePaveText(0.46, 0.10, 0.68, 0.90, "nbNDC ARC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(subtitle3_text.get());
    ROOTHelper::SetPaveAttribute(subtitle3_text.get(), 0, 0.2);
    ROOTHelper::SetFillAttribute(subtitle3_text.get(), 1001, kAzure-7, 0.5f);
    ROOTHelper::SetTextAttribute(subtitle3_text.get(), 35.0f, 103, 22);
    subtitle3_text->AddText(model_object->GetEmdID().data());
    subtitle3_text->Draw();

    auto legend{ ROOTHelper::CreateLegend(0.69, 0.10, 0.99, 0.90, false) };
    ROOTHelper::SetLegendDefaultStyle(legend.get());
    ROOTHelper::SetFillAttribute(legend.get(), 4000);
    ROOTHelper::SetTextAttribute(legend.get(), 25.0f, 133, 12, 0.0);
    legend->SetMargin(0.15f);
    legend->AddEntry(gaus_function[0].get(),
        "Gaussian Model #color[633]{#phi (#font[1]{A},#font[1]{#tau})}", "l");
    legend->AddEntry(map_value_graph_list[0].at(0).get(),
        "Map Value", "l");
    legend->Draw();

    canvas->cd();
    auto pad_bottom{ ROOTHelper::CreatePad("pad_bottom","", 0.10, 0.00, 0.98, 0.11) };
    ROOTHelper::SetPadDefaultStyle(pad_bottom.get());
    pad_bottom->Draw();
    pad_bottom->cd();
    auto x_title_text{ ROOTHelper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(x_title_text.get());
    ROOTHelper::SetFillAttribute(x_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(x_title_text.get(), 35.0f, 133, 22);
    x_title_text->AddText("Radial Distance #[]{#AA}");
    x_title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

} // namespace rhbm_gem
