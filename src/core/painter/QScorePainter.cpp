#include <rhbm_gem/core/PainterFunctions.hpp>

#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/FilePathHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include "detail/PainterModelValidation.hpp"
#include "detail/PainterStyle.hpp"
#include "detail/PotentialPlotBuilder.hpp"

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#endif

#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace rhbm_gem {

namespace {

class QScorePainter
{
    std::vector<ModelObject *> m_model_object_list;
    std::string m_folder_path{ "./" };

public:
    QScorePainter(const core::ModelObjectList & model_objects, const std::string & output_folder);
    void Run();

private:
    void AddModel(ModelObject & data_object);
    void PaintAverageQScoreToSequenceSummary(ModelObject * model_object, const std::string & name);
};

QScorePainter::QScorePainter(
    const core::ModelObjectList & model_objects,
    const std::string & output_folder) :
    m_folder_path{ path_helper::EnsureTrailingSlash(output_folder) }
{
    for (auto * model_object : model_objects)
    {
        AddModel(*model_object);
    }
}

void QScorePainter::AddModel(ModelObject & data_object)
{
    painter_internal::RequireGroupedAnalyzedModel(data_object, "QScorePainter");
    m_model_object_list.push_back(&data_object);
}

void QScorePainter::Run()
{
    Logger::Log(LogLevel::Info, "QScorePainter::Run() called.");
    Logger::Log(LogLevel::Info, "Folder path: " + m_folder_path);
    Logger::Log(LogLevel::Info, "Number of model objects to be painted: "
                + std::to_string(m_model_object_list.size()));

    for (auto * model_object : m_model_object_list)
    {
        auto label{ painter_internal::BuildPainterOutputLabel(*model_object) };
        PaintAverageQScoreToSequenceSummary(model_object, "average_qscore_to_sequence_summary_" + label);
    }
}

void QScorePainter::PaintAverageQScoreToSequenceSummary(
    ModelObject * model_object,
    const std::string & name)
{
    auto file_path{ m_folder_path + name };
    Logger::Log(LogLevel::Info, "QScorePainter::PaintAverageQScoreToSequenceSummary");

    #ifdef HAVE_ROOT
    auto plot_builder{ std::make_unique<PotentialPlotBuilder>(model_object) };

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    auto canvas{ root_helper::CreateCanvas("test","", 1500, 600) };
    root_helper::SetCanvasDefaultStyle(canvas.get());
    root_helper::PrintCanvasOpen(canvas.get(), file_path);

    const int method_count{ 3 };
    std::unique_ptr<TH2> frame;
    std::unordered_map<std::string, std::unique_ptr<TGraphErrors>> qscore_graph_map[method_count];
    qscore_graph_map[0] = plot_builder->CreateAverageQScoreToSequenceIDGraphMap(true, false);
    qscore_graph_map[1] = plot_builder->CreateAverageQScoreToSequenceIDGraphMap(false, false);
    qscore_graph_map[2] = plot_builder->CreateAverageQScoreToSequenceIDGraphMap(false, true);

    for (auto & [chain_id, gaus_graph] : qscore_graph_map[0])
    {
        double y_min{ 0.5 };
        double y_max{ 1.0 };

        std::vector<double> x_array;
        x_array.reserve(static_cast<size_t>(gaus_graph->GetN()));
        for (int p = 0; p < gaus_graph->GetN(); p++)
        {
            x_array.emplace_back(gaus_graph->GetPointX(p));
        }
        auto x_range{ array_helper::ComputeScalingRangeTuple(x_array, 0.05) };
        auto x_min{ std::get<0>(x_range) };
        auto x_max{ std::get<1>(x_range) };

        std::unique_ptr<TPaveText> subtitle1_text;
        std::unique_ptr<TPaveText> subtitle2_text;
        std::unique_ptr<TPaveText> subtitle3_text;
        std::unique_ptr<TLegend> legend;
        root_helper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        root_helper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0, 0);
        root_helper::SetPadMarginInCanvas(gPad, 0.09, 0.01, 0.16, 0.12);
        if (frame == nullptr)
        {
            frame = root_helper::CreateHist2D("frame","", 500, 0.0, 1.0, 500, 0.0, 1.0);
            root_helper::SetAxisTitleAttribute(frame->GetXaxis(), 45.0f, 0.9f, 133);
            root_helper::SetAxisLabelAttribute(frame->GetXaxis(), 40.0f, 0.01f, 133);
            root_helper::SetAxisTickAttribute(frame->GetXaxis(), static_cast<float>(0.05), 510);
            root_helper::SetAxisTitleAttribute(frame->GetYaxis(), 45.0f, 1.3f, 133);
            root_helper::SetAxisLabelAttribute(frame->GetYaxis(), 40.0f, 0.005f, 133);
            root_helper::SetAxisTickAttribute(frame->GetYaxis(), static_cast<float>(0.02), 506);
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->SetStats(0);
        }
        frame->GetXaxis()->SetTitle(Form("Residue ID #[]{Chain %s}", chain_id.data()));
        frame->GetYaxis()->SetTitle("Q-Score");
        frame->GetXaxis()->SetLimits(x_min, x_max);
        frame->GetYaxis()->SetLimits(y_min, y_max);
        frame->Draw();

        short marker_list[method_count] = { 20, 24, 25 };
        short color_list[method_count] = { kBlue, kGreen+1, kRed };
        for (size_t k = 0; k < method_count; k++)
        {
            root_helper::SetMarkerAttribute(
                qscore_graph_map[k].at(chain_id).get(),
                marker_list[k], 1.2f, color_list[k]);
            root_helper::SetLineAttribute(
                qscore_graph_map[k].at(chain_id).get(),
                1, 1, color_list[k]);
            qscore_graph_map[k].at(chain_id)->Draw("PL X0");
        }

        legend = root_helper::CreateLegend(0.10, 0.90, 1.00, 1.00, false);
        root_helper::SetLegendDefaultStyle(legend.get());
        root_helper::SetFillAttribute(legend.get(), 4000);
        root_helper::SetTextAttribute(legend.get(), 40.0f, 133, 12, 0.0);
        legend->SetMargin(0.15f);
        legend->SetNColumns(2);
        legend->AddEntry(qscore_graph_map[0].at(chain_id).get(), "Original", "lp");
        legend->AddEntry(qscore_graph_map[1].at(chain_id).get(), "Using Gaussian parameters from this estimation", "lp");
        legend->Draw();
        root_helper::PrintCanvasPad(canvas.get(), file_path);
    }
    root_helper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path);
    #else
    (void)file_path;
    #endif

    auto average_qscore{ model_object->GetStandardAverageQScore() };
    Logger::Log(LogLevel::Info, Form("Average Q-Score: %.3f", average_qscore));
}

} // namespace

namespace core {

void PaintQScore(
    const ModelObjectList & model_objects,
    const std::string & output_folder)
{
    QScorePainter painter{ model_objects, output_folder };
    painter.Run();
}

} // namespace core
} // namespace rhbm_gem
