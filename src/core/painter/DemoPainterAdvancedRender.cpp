#include <rhbm_gem/core/painter/DemoPainter.hpp>
#include <rhbm_gem/core/painter/PotentialPlotBuilder.hpp>
#include <rhbm_gem/data/object/AtomClassifier.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/data/object/PotentialEntryQuery.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#endif

#include <vector>

namespace rhbm_gem {

void DemoPainter::PaintGroupGausToFSC(
    const std::vector<ModelObject *> & model_list, const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " DemoPainter::PaintGroupGausToFSC");

    for (auto & model : model_list)
    {
        if (model == nullptr) return;
        Logger::Log(LogLevel::Debug, "Find Model : "+ model->GetPdbID());
        
    }

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 4 };
    const int row_size{ 1 };

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 300) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(canvas.get(), col_size, row_size, 0.08f, 0.01f, 0.30f, 0.02f, 0.01f, 0.005f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    std::unique_ptr<TGraphErrors> graph[col_size];
    std::unique_ptr<TF1> fit_function[col_size];
    std::vector<double> x_array, y_array;
    double r_square[col_size]{ 0.0 };
    double slope[col_size]{ 0.0 };
    double intercept[col_size]{ 0.0 };
    for (size_t i = 0; i < AtomClassifier::GetMainChainMemberCount(); i++)
    {
        if (i >= col_size) continue;
        auto group_key{ m_atom_classifier->GetMainChainSimpleAtomClassGroupKey(i) };
        graph[i] = ROOTHelper::CreateGraphErrors();
        auto count{ 0 };
        for (auto model : model_list)
        {
            auto entry_iter{ std::make_unique<PotentialEntryQuery>(model) };
            auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model) };
            auto width_value{ entry_iter->GetAtomGausEstimatePrior(group_key, ChemicalDataHelper::GetSimpleAtomClassKey(), 1) };
            auto width_error{ entry_iter->GetAtomGausVariancePrior(group_key, ChemicalDataHelper::GetSimpleAtomClassKey(), 1) };
            graph[i]->SetPoint(count, model->GetResolution(), width_value);
            graph[i]->SetPointError(count, 0.0, width_error);
            count++;
        }
        for (int p = 0; p < graph[i]->GetN(); p++)
        {
            x_array.push_back(graph[i]->GetPointX(p));
            y_array.push_back(graph[i]->GetPointY(p));
        }
        auto element_color{ AtomClassifier::GetMainChainElementColor(i) };
        ROOTHelper::SetMarkerAttribute(graph[i].get(), 20, 1.2f, element_color);
        ROOTHelper::SetLineAttribute(graph[i].get(), 1, 2, element_color);

        r_square[i] = ROOTHelper::PerformLinearRegression(graph[i].get(), slope[i], intercept[i]);
        fit_function[i] = ROOTHelper::CreateFunction1D(Form("fit_%d", static_cast<int>(i)), "x*[1]+[0]");
        fit_function[i]->SetParameters(intercept[i], slope[i]);
        ROOTHelper::SetLineAttribute(fit_function[i].get(), 2, 2, kRed);
    }

    auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array, 0.2) };
    double x_min{ std::get<0>(x_range) };
    double x_max{ std::get<1>(x_range) };

    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_array, 0.3) };
    double y_min{ std::get<0>(y_range) };
    double y_max{ std::get<1>(y_range) };
    
    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> title_text[col_size];
    std::unique_ptr<TPaveText> r_square_text[col_size];
    for (size_t i = 0; i < col_size; i++)
    {
        auto element_color{ AtomClassifier::GetMainChainElementColor(i) };
        auto element_label{ AtomClassifier::GetMainChainElementLabel(i) };
        for (int j = 0; j < row_size; j++)
        {
            ROOTHelper::FindPadInCanvasPartition(canvas.get(), static_cast<int>(i), j);
            ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = ROOTHelper::CreateHist2D(Form("hist_%d_%d", static_cast<int>(i), j),"", 500, x_min, x_max, 500, y_min, y_max);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 0.0f);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.005f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 506);
            ROOTHelper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 45.0f, 1.3f, 133);
            ROOTHelper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 45.0f, 0.01f, 133);
            ROOTHelper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.08/y_factor), 506);
            ROOTHelper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("Width #tau");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw();
            graph[i]->Draw("P X0");
            fit_function[i]->SetRange(x_min, x_max);
            fit_function[i]->Draw("SAME");
        }
        title_text[i] = ROOTHelper::CreatePaveText(0.10, 0.75, 0.40, 0.99, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[i].get());
        ROOTHelper::SetPaveAttribute(title_text[i].get(), 0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[i].get(), 60.0f, 103, 22);
        ROOTHelper::SetFillAttribute(title_text[i].get(), 1001, element_color, 0.5f);
        title_text[i]->AddText(element_label.data());
        title_text[i]->Draw();

        r_square_text[i] = ROOTHelper::CreatePaveText(0.48, 0.05, 0.99, 0.25, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(r_square_text[i].get());
        ROOTHelper::SetPaveAttribute(r_square_text[i].get(), 0, 0.5);
        ROOTHelper::SetLineAttribute(r_square_text[i].get(), 1, 0);
        ROOTHelper::SetTextAttribute(r_square_text[i].get(), 40.0f, 133, 22);
        ROOTHelper::SetFillAttribute(r_square_text[i].get(), 1001, kAzure-7, 0.20f);
        r_square_text[i]->AddText(Form("R^{2} = %.2f", r_square[i]));
        r_square_text[i]->Draw();
    }
    canvas->cd();
    auto pad_extra{ ROOTHelper::CreatePad("pad_extra","", 0.08, 0.00, 0.98, 0.19) };
    pad_extra->Draw();
    pad_extra->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra.get());
    ROOTHelper::SetFillAttribute(pad_extra.get(), 4000);
    auto bottom_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(bottom_title_text.get());
    ROOTHelper::SetFillAttribute(bottom_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(bottom_title_text.get(), 40.0f, 133, 22);
    bottom_title_text->AddText("FSC Resolution #[]{#AA}");
    bottom_title_text->Draw();
    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

void DemoPainter::PaintAtomWidthScatterPlotSingle(
    ModelObject * model_object, const std::string & name, bool draw_box_plot)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, " DemoPainter::PaintAtomWidthScatterPlotSingle");
    auto class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    (void)draw_box_plot;

    if (model_object == nullptr) return;
    auto entry_iter{ std::make_unique<PotentialEntryQuery>(model_object) };
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);

    const int pad_size{ 5 };

    auto canvas{ ROOTHelper::CreateCanvas("test","", 1500, 400) };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), pad_size, 1, 0.08f, 0.08f, 0.25f, 0.13f, 0.011f, 0.011f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), file_path);

    auto normalized_z_position{ 0.5 };
    auto z_ratio_window{ 0.1 };
    auto graph_position{
        plot_builder->CreateAtomXYPositionTomographyGraph(normalized_z_position, z_ratio_window, true)
    };
    auto z_position{ model_object->GetModelPosition(2, normalized_z_position) };
    auto z_window{ model_object->GetModelLength(2) * z_ratio_window };
    auto com_pos{ model_object->GetCenterOfMassPosition() };
    std::vector<double> x_pos_array, y_pos_array;
    x_pos_array.reserve(static_cast<size_t>(graph_position->GetN()));
    y_pos_array.reserve(static_cast<size_t>(graph_position->GetN()));
    for (int p = 0; p < graph_position->GetN(); p++)
    {
        x_pos_array.emplace_back(graph_position->GetPointX(p));
        y_pos_array.emplace_back(graph_position->GetPointY(p));
    }

    auto x_pos_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_pos_array, 0.1) };
    double x_pos_min{ std::get<0>(x_pos_range) };
    double x_pos_max{ std::get<1>(x_pos_range) };

    auto y_pos_range{ ArrayStats<double>::ComputeScalingRangeTuple(y_pos_array, 0.1) };
    double y_pos_min{ std::get<0>(y_pos_range) };
    double y_pos_max{ std::get<1>(y_pos_range) };

    std::vector<std::unique_ptr<TGraphErrors>> graph_list[pad_size - 1];
    std::vector<double> x_array[pad_size - 1];
    std::vector<double> y_array;
    for (size_t i = 0; i < pad_size - 1; i++)
    {
        for (auto residue : ChemicalDataHelper::GetStandardAminoAcidList())
        {
            auto group_key{ m_atom_classifier->GetMainChainComponentAtomClassGroupKey(i, residue) };
            if (entry_iter->IsAvailableAtomGroupKey(group_key, class_key) == false) continue;
            auto gaus_graph{ plot_builder->CreateCOMDistanceToGausEstimateGraph(group_key, class_key, 1) };
            for (int p = 0; p < gaus_graph->GetN(); p++)
            {
                x_array[i].emplace_back(gaus_graph->GetPointX(p));
                y_array.emplace_back(gaus_graph->GetPointY(p));
            }
            graph_list[i].emplace_back(std::move(gaus_graph));
        }
    }

    double x_min[pad_size - 1];
    double x_max[pad_size - 1];
    for (int i = 0; i < pad_size - 1; i++)
    {
        auto x_range{ ArrayStats<double>::ComputeScalingRangeTuple(x_array[i], 0.1) };
        x_min[i] = std::get<0>(x_range);
        x_max[i] = std::get<1>(x_range);
    }

    auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(y_array, 0.1) };
    auto y_min{ std::get<0>(y_range) };
    auto y_max{ std::get<1>(y_range) };

    std::unique_ptr<TH2D> summary_hist[pad_size - 1];
    for (int i = 0; i < pad_size - 1; i++)
    {
        summary_hist[i] = ROOTHelper::CreateHist2D(
            Form("hist_%d", i), "", 5, 35.0, 65.0, 100, y_min, y_max);
        for (auto & graph : graph_list[i])
        {
            for (int p = 0; p < graph->GetN(); p++)
            {
                summary_hist[i]->Fill(graph->GetPointX(p), graph->GetPointY(p));
            }
        }
    }

    std::unique_ptr<TH2> frame[pad_size];
    ROOTHelper::FindPadInCanvasPartition(canvas.get(), 0, 0);
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    auto x_factor_tmp{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
    auto y_factor_tmp{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
    frame[0] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", 0, 0),"", 500, x_pos_min, x_pos_max, 500, y_pos_min, y_pos_max);
    ROOTHelper::SetAxisTitleAttribute(frame[0]->GetXaxis(), 50.0f, 0.8f, 133);
    ROOTHelper::SetAxisLabelAttribute(frame[0]->GetXaxis(), 45.0f, 0.005f, 133);
    ROOTHelper::SetAxisTickAttribute(frame[0]->GetXaxis(), static_cast<float>(y_factor_tmp*0.05/x_factor_tmp), 506);
    ROOTHelper::SetAxisTitleAttribute(frame[0]->GetYaxis(), 50.0f, 1.1f, 133);
    ROOTHelper::SetAxisLabelAttribute(frame[0]->GetYaxis(), 45.0f, 0.01f, 133);
    ROOTHelper::SetAxisTickAttribute(frame[0]->GetYaxis(), static_cast<float>(x_factor_tmp*0.05/y_factor_tmp), 506);
    ROOTHelper::SetLineAttribute(frame[0].get(), 1, 0);
    frame[0]->GetXaxis()->SetTitle("#font[2]{x} #[]{#AA}");
    frame[0]->GetYaxis()->SetTitle("#font[2]{y} #[]{#AA}");
    frame[0]->GetXaxis()->CenterTitle();
    frame[0]->GetYaxis()->CenterTitle();
    frame[0]->SetStats(0);
    frame[0]->Draw();
    ROOTHelper::SetMarkerAttribute(graph_position.get(), 20, 0.2f, kAzure-7);
    graph_position->Draw("P X0");

    auto circle_inner{ new TEllipse(0.0, 0.0, 40.0) };
    ROOTHelper::SetFillAttribute(circle_inner, 0);
    ROOTHelper::SetLineAttribute(circle_inner, 1, 1, kRed);
    circle_inner->Draw();
    auto circle_outer{ new TEllipse(0.0, 0.0, 60.0) };
    ROOTHelper::SetFillAttribute(circle_outer, 0);
    ROOTHelper::SetLineAttribute(circle_outer, 1, 1, kRed);
    circle_outer->Draw();

    auto pos_text{ ROOTHelper::CreatePaveText(0.08, 1.0, 1.0, 1.2, "nbNDC", true) };
    ROOTHelper::SetPaveTextDefaultStyle(pos_text.get());
    ROOTHelper::SetFillAttribute(pos_text.get(), 4000);
    ROOTHelper::SetTextAttribute(pos_text.get(), 35, 133, 32, 0.0, kAzure-7);
    auto slicing_z_center{ z_position - com_pos.at(2) };
    pos_text->AddText(
        Form("#font[2]{z} #in [%.1f, %.1f] #AA",
        slicing_z_center - z_window*0.5, slicing_z_center + z_window*0.5));
    pos_text->Draw();

    std::unique_ptr<TPaveText> title_text[pad_size - 1];
    for (int i = 1; i < pad_size; i++)
    {
        auto element_id{ i - 1 };
        ROOTHelper::FindPadInCanvasPartition(canvas.get(), i, 0);
        ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        auto x_factor{ ROOTHelper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
        auto y_factor{ ROOTHelper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
        frame[i] = ROOTHelper::CreateHist2D(Form("frame_%d_%d", i, 0),"", 500, x_min[element_id], x_max[element_id], 500, y_min, y_max);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetXaxis(), 0.0f);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetXaxis(), 45.0f, 0.005f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetXaxis(), static_cast<float>(y_factor*0.05/x_factor), 505);
        ROOTHelper::SetAxisTitleAttribute(frame[i]->GetYaxis(), 55.0f, 1.0f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame[i]->GetYaxis(), 45.0f, 0.01f, 133);
        ROOTHelper::SetAxisTickAttribute(frame[i]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
        ROOTHelper::SetLineAttribute(frame[i].get(), 1, 0);
        frame[i]->GetXaxis()->SetTitle("");
        frame[i]->GetYaxis()->SetTitle("Width #tau");
        frame[i]->GetXaxis()->CenterTitle();
        frame[i]->GetYaxis()->CenterTitle();
        frame[i]->SetStats(0);
        frame[i]->GetXaxis()->SetLimits(35.0, 65.0);
        frame[i]->Draw("Y+");

        auto element_color{ AtomClassifier::GetMainChainElementColor(static_cast<size_t>(element_id)) };
        auto element_label{ AtomClassifier::GetMainChainElementLabel(static_cast<size_t>(element_id)) };
        if (draw_box_plot == true)
        {
            ROOTHelper::SetLineAttribute(summary_hist[element_id].get(), 1, 1, element_color);
            ROOTHelper::SetFillAttribute(summary_hist[element_id].get(), 1001, element_color, 0.3f);
            summary_hist[element_id]->Draw("CANDLE3 SAME");
        }
        else
        {
            for (auto & graph : graph_list[element_id])
            {
                ROOTHelper::SetMarkerAttribute(graph.get(), 24, 0.1f, element_color);
                graph->Draw("P");
            }
        }
        title_text[element_id] = ROOTHelper::CreatePaveText(0.01, 0.95, 0.39, 1.20, "nbNDC ARC", true);
        ROOTHelper::SetPaveTextDefaultStyle(title_text[element_id].get());
        ROOTHelper::SetPaveAttribute(title_text[element_id].get(), 0, 0.2);
        ROOTHelper::SetTextAttribute(title_text[element_id].get(), 60.0f, 103, 22);
        ROOTHelper::SetFillAttribute(title_text[element_id].get(), 1001, element_color, 0.5f);
        title_text[element_id]->AddText(element_label.data());
        title_text[element_id]->Draw();
    }

    canvas->cd();
    auto pad_extra{ ROOTHelper::CreatePad("pad_extra","", 0.25, 0.00, 0.92, 0.15) };
    pad_extra->Draw();
    pad_extra->cd();
    ROOTHelper::SetPadDefaultStyle(pad_extra.get());
    ROOTHelper::SetFillAttribute(pad_extra.get(), 4000);
    auto bottom_title_text{ ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    ROOTHelper::SetPaveTextDefaultStyle(bottom_title_text.get());
    ROOTHelper::SetFillAttribute(bottom_title_text.get(), 4000);
    ROOTHelper::SetTextAttribute(bottom_title_text.get(), 50.0f, 133, 22);
    bottom_title_text->AddText("Distance to Center of Mass #font[2]{d} #[]{#AA}");
    bottom_title_text->Draw();

    ROOTHelper::PrintCanvasPad(canvas.get(), file_path);
    ROOTHelper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #endif
}

} // namespace rhbm_gem
