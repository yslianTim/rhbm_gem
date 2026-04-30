#include "RHBMTestPlotting.hpp"

#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <tuple>
#include <utility>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLegend.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#endif

namespace rhbm_gem::command_detail::rhbm_test_plotting {
namespace {

std::string FormatSigmaToken(double value)
{
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(3) << value;
    auto token{ stream.str() };
    std::replace(token.begin(), token.end(), '.', '_');
    std::replace(token.begin(), token.end(), '-', 'm');
    return token;
}

std::string FormatDistanceLabel(double distance)
{
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(1)
           << "Distance = " << distance << " Angstrom";
    return stream.str();
}

bool ValidateLinearizedDataset(
    const RHBMMemberDataset & dataset,
    std::string_view dataset_label)
{
    if (dataset.X.rows() != dataset.y.size())
    {
        Logger::Log(
            LogLevel::Warning,
            std::string("Skip linearized benchmark plot for ") + std::string(dataset_label)
                + " because X rows do not match y size.");
        return false;
    }
    if (dataset.X.cols() < 2)
    {
        Logger::Log(
            LogLevel::Warning,
            std::string("Skip linearized benchmark plot for ") + std::string(dataset_label)
                + " because the dataset has fewer than 2 basis columns.");
        return false;
    }
    return true;
}

LineSeries BuildLinearizedDatasetSeries(
    const RHBMMemberDataset & dataset,
    std::string_view dataset_label,
    std::string series_name)
{
    std::vector<std::pair<double, double>> points;
    points.reserve(static_cast<std::size_t>(dataset.y.size()));
    for (Eigen::Index row = 0; row < dataset.y.size(); ++row)
    {
        const double x_value{ dataset.X(row, 1) };
        const double y_value{ dataset.y(row) };
        if (!std::isfinite(x_value) || !std::isfinite(y_value))
        {
            Logger::Log(
                LogLevel::Warning,
                std::string("Drop non-finite linearized benchmark point for ")
                    + std::string(dataset_label) + ".");
            continue;
        }
        points.emplace_back(x_value, y_value);
    }
    std::sort(points.begin(), points.end());

    LineSeries series;
    series.name = std::move(series_name);
    series.x_values.reserve(points.size());
    series.y_values.reserve(points.size());
    for (const auto & point : points)
    {
        series.x_values.emplace_back(point.first);
        series.y_values.emplace_back(point.second);
    }
    return series;
}

#ifdef HAVE_ROOT
double ScaleBiasPlotX(double x, BiasXAxisMode x_axis_mode)
{
    if (x_axis_mode == BiasXAxisMode::NeighborDistance)
    {
        return -1.0 * x;
    }
    if (x_axis_mode == BiasXAxisMode::NeighborType)
    {
        return x;
    }
    return 100.0 * x;
}

short GetBiasCurveColor(BiasPlotFlavor flavor, BiasCurveKind kind)
{
    if (flavor == BiasPlotFlavor::MemberOutlier)
    {
        if (kind == BiasCurveKind::Median)
        {
            return kAzure;
        }
        if (kind == BiasCurveKind::Mdpde)
        {
            return kRed;
        }
        return kGreen + 1;
    }
    if (flavor == BiasPlotFlavor::ModelAlphaMember)
    {
        return (kind == BiasCurveKind::RequestedAlpha) ? kAzure : kRed;
    }
    if (flavor == BiasPlotFlavor::ModelAlphaData)
    {
        return (kind == BiasCurveKind::RequestedAlpha) ? kAzure : kGreen + 2;
    }
    if (kind == BiasCurveKind::Ols)
    {
        return kAzure;
    }
    if (kind == BiasCurveKind::Mdpde)
    {
        return kGreen + 2;
    }
    return kRed;
}

short GetBiasCurveMarker(BiasPlotFlavor flavor, BiasCurveKind kind)
{
    if (flavor == BiasPlotFlavor::MemberOutlier)
    {
        if (kind == BiasCurveKind::Median)
        {
            return 24;
        }
        if (kind == BiasCurveKind::Mdpde)
        {
            return 20;
        }
        return 25;
    }
    if (flavor == BiasPlotFlavor::ModelAlphaMember)
    {
        return (kind == BiasCurveKind::RequestedAlpha) ? 24 : 20;
    }
    if (kind == BiasCurveKind::Ols || kind == BiasCurveKind::RequestedAlpha)
    {
        return 24;
    }
    if (kind == BiasCurveKind::Mdpde || kind == BiasCurveKind::TrainedAlpha)
    {
        return 25;
    }
    return 20;
}

short GetBiasCurveLineStyle(BiasPlotFlavor flavor, BiasCurveKind kind)
{
    if (kind == BiasCurveKind::Ols ||
        kind == BiasCurveKind::Median ||
        kind == BiasCurveKind::RequestedAlpha)
    {
        return 2;
    }
    if (kind == BiasCurveKind::TrainedMdpde &&
        flavor != BiasPlotFlavor::ModelAlphaMember)
    {
        return 1;
    }
    return 3;
}

void ApplyBiasCurveStyle(
    TGraphErrors * graph,
    BiasPlotFlavor flavor,
    BiasCurveKind kind)
{
    const auto color{ GetBiasCurveColor(flavor, kind) };
    root_helper::SetMarkerAttribute(
        graph,
        GetBiasCurveMarker(flavor, kind),
        1.5f,
        color);
    root_helper::SetLineAttribute(
        graph,
        GetBiasCurveLineStyle(flavor, kind),
        2,
        color);
    root_helper::SetFillAttribute(graph, 1001, color, 0.2f);
}

std::vector<BiasCurveKind> GetBiasLegendOrder(BiasPlotFlavor flavor)
{
    if (flavor == BiasPlotFlavor::NeighborType)
    {
        return { BiasCurveKind::Ols, BiasCurveKind::Mdpde, BiasCurveKind::TrainedMdpde };
    }
    if (flavor == BiasPlotFlavor::ModelAlphaData ||
        flavor == BiasPlotFlavor::ModelAlphaMember)
    {
        return { BiasCurveKind::RequestedAlpha, BiasCurveKind::TrainedAlpha };
    }
    if (flavor == BiasPlotFlavor::MemberOutlier)
    {
        return { BiasCurveKind::TrainedMdpde, BiasCurveKind::Mdpde, BiasCurveKind::Median };
    }
    return { BiasCurveKind::TrainedMdpde, BiasCurveKind::Mdpde, BiasCurveKind::Ols };
}

std::string GetBiasLegendLabel(BiasPlotFlavor flavor, BiasCurveKind kind)
{
    if (flavor == BiasPlotFlavor::NeighborType)
    {
        if (kind == BiasCurveKind::Ols)
        {
            return "OLS";
        }
        if (kind == BiasCurveKind::Mdpde)
        {
            return "MDPDE";
        }
        return "MDPDE (Sampling Scheme)";
    }
    if (flavor == BiasPlotFlavor::ModelAlphaData)
    {
        return (kind == BiasCurveKind::RequestedAlpha) ?
            "MDPDE (#alpha_{r} = 0.1)" :
            "MDPDE (#alpha_{r} = 0.4)";
    }
    if (flavor == BiasPlotFlavor::ModelAlphaMember)
    {
        return (kind == BiasCurveKind::RequestedAlpha) ?
            "MDPDE (#alpha_{g} = 0.2)" :
            "MDPDE (#alpha_{g} = 0.5)";
    }
    if (flavor == BiasPlotFlavor::MemberOutlier)
    {
        if (kind == BiasCurveKind::TrainedMdpde)
        {
            return "MDPDE (#alpha_{g} from Alg.5)";
        }
        if (kind == BiasCurveKind::Mdpde)
        {
            return "MDPDE (#alpha_{g} = 0.2)";
        }
        return "MDPDE (#alpha_{g} = 0)";
    }
    if (kind == BiasCurveKind::TrainedMdpde)
    {
        return "MDPDE (w/ Sampling Scheme)";
    }
    if (kind == BiasCurveKind::Mdpde)
    {
        return "MDPDE";
    }
    return "Ordinary Least Squares";
}

TGraphErrors * FindBiasGraph(
    std::vector<std::unique_ptr<TGraphErrors>> & graph_list,
    const std::vector<BiasCurveKind> & kind_list,
    BiasCurveKind kind)
{
    for (size_t i = 0; i < kind_list.size(); i++)
    {
        if (kind_list.at(i) == kind)
        {
            return graph_list.at(i).get();
        }
    }
    return nullptr;
}

void ApplyNeighborTypeAxisLabels(TAxis * axis)
{
    axis->SetNdivisions(5, false);
    const auto label_color{ kCyan+3 };
    axis->ChangeLabel(1, -1, 0, -1, -1, -1, "");
    axis->ChangeLabel(2, -90.0, -1, -1, label_color, 103, "O");
    axis->ChangeLabel(3, -90.0, -1, -1, label_color, 103, "N");
    axis->ChangeLabel(4, -90.0, -1, -1, label_color, 103, "C");
    axis->ChangeLabel(5, -90.0, -1, -1, label_color, 103, "C_{#alpha}");
    axis->ChangeLabel(-1, -1, 0, -1, -1, -1, "");
}

std::unique_ptr<TGraphErrors> CreateDistanceToResponseGraph(
    const LocalPotentialSampleList & sampling_entries)
{
    auto graph{ root_helper::CreateGraphErrors() };
    auto count{ 0 };
    for (const auto & sample : sampling_entries)
    {
        graph->SetPoint(count, sample.distance, sample.response);
        count++;
    }
    return graph;
}

std::unique_ptr<TH2D> CreateDistanceToResponseHistogram(
    const LocalPotentialSampleList & sampling_entries, int x_bin_size, int y_bin_size = 500)
{
    auto hist{
        root_helper::CreateHist2D(
            "hist_distance_response", "Distance vs Response",
            x_bin_size, 0.0, 4.0,
            y_bin_size, 0.0, 1.0)
    };
    for (const auto & sample : sampling_entries)
    {
        hist->Fill(sample.distance, sample.response);
    }
    return hist;
}
#endif

} // namespace

BiasCurve MakeBiasCurve(BiasCurveKind kind, size_t point_capacity)
{
    BiasCurve curve{ kind, {} };
    curve.points.reserve(point_capacity);
    return curve;
}

void AppendBiasCurvePoint(
    BiasCurve & curve,
    double x,
    const rhbm_tester::BiasStatistics & bias)
{
    curve.points.emplace_back(BiasCurvePoint{ x, bias });
}

std::string FormatDataBiasPanelLabel(size_t panel_index)
{
    const double error_value[3]{ 0.0, 2.5, 5.0 };
    const auto value{ error_value[panel_index] };
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(1)
           << "#sigma_{#epsilon} = " << value << "% D_{max}";
    return stream.str();
}

std::string FormatMemberBiasPanelLabel(size_t panel_index)
{
    const std::string outlier_type_list[2]{
        "#font[2]{A}", "#tau"
    };
    return "Outlier in " + outlier_type_list[panel_index];
}

bool TryAppendBenchmarkLinearizedPanel(
    std::vector<LinePlotPanel> & panels,
    double distance,
    const RHBMMemberDataset & no_cut_dataset,
    const RHBMMemberDataset & cut_dataset)
{
    const auto label{ FormatDistanceLabel(distance) };
    if (!ValidateLinearizedDataset(no_cut_dataset, label + " (No Cut)") ||
        !ValidateLinearizedDataset(cut_dataset, label + " (Cut)"))
    {
        return false;
    }

    auto no_cut_series{
        BuildLinearizedDatasetSeries(no_cut_dataset, label + " (No Cut)", "No Cut")
    };
    auto cut_series{
        BuildLinearizedDatasetSeries(cut_dataset, label + " (Cut)", "Cut")
    };
    if (no_cut_series.x_values.empty() || cut_series.x_values.empty())
    {
        Logger::Log(
            LogLevel::Warning,
            "Skip linearized benchmark panel for " + label
                + " because at least one series has no plottable points.");
        return false;
    }

    panels.emplace_back(LinePlotPanel{
        label,
        AxisSpec{ "Linearized Response" },
        { std::move(no_cut_series), std::move(cut_series) }
    });
    return true;
}

void SaveBenchmarkLinearizedDatasetReport(
    const RHBMTestRequest & request,
    double error_sigma,
    const std::vector<LinePlotPanel> & panels)
{
    if (panels.empty())
    {
        Logger::Log(
            LogLevel::Warning,
            "Skip benchmark linearized dataset report because no valid panels were collected.");
        return;
    }

    LinePlotRequest line_plot_request;
    line_plot_request.output_path = request.output_dir /
        std::string("benchmark_cut_vs_no_cut_sigma_" + FormatSigmaToken(error_sigma) + ".pdf");
    line_plot_request.title = "";
    line_plot_request.x_axis.title = "Linearized Basis";
    line_plot_request.shared_y_axis_title = "Linearized Response";
    line_plot_request.panels = panels;
    line_plot_request.canvas_width = 2000;
    line_plot_request.canvas_height_per_panel = 200;

    const auto plot_result{ local_painter::SaveLinePlot(line_plot_request) };
    if (!plot_result.Succeeded())
    {
        Logger::Log(
            LogLevel::Warning,
            "Failed to emit benchmark linearized dataset report '"
                + line_plot_request.output_path.string() + "': " + plot_result.message);
    }
}

void SaveDataOutlierBiasPlot(
    const RHBMTestRequest & request,
    const BiasPlotRequest & plot_request)
{
    Logger::Log(LogLevel::Info, " RHBMTestCommand::SaveDataOutlierBiasPlot");
    if (plot_request.panels.empty())
    {
        Logger::Log(LogLevel::Warning, "Skip data outlier result plot because no panels were collected.");
        return;
    }

    #ifdef HAVE_ROOT
    auto file_path{ request.output_dir / plot_request.output_name };
    std::vector<std::string> title_y_list{
        "Amplitude #font[2]{A}", "Width #tau"
    };

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const size_t col_count{ plot_request.panels.size() };
    const int col_size{ static_cast<int>(col_count) };
    const int row_size{ 2 };
    const int canvas_width{ (col_size == 1) ? 900 : 1500 };
    auto canvas{ root_helper::CreateCanvas("test","", canvas_width, 750) };
    root_helper::SetCanvasDefaultStyle(canvas.get());
    root_helper::SetCanvasPartition(
        canvas.get(), col_size, row_size,
        (plot_request.flavor == BiasPlotFlavor::NeighborType) ? 0.25f : 0.17f,
        (plot_request.flavor == BiasPlotFlavor::NeighborType) ? 0.15f : 0.08f,
        0.12f,
        0.14f,
        0.01f, 0.01f
    );
    root_helper::PrintCanvasOpen(canvas.get(), file_path);

    std::vector<std::vector<std::vector<std::unique_ptr<TGraphErrors>>>> graph_list(col_count);
    std::vector<std::vector<std::vector<BiasCurveKind>>> graph_kind_list(col_count);
    for (size_t i = 0; i < col_count; i++)
    {
        graph_list.at(i).resize(row_size);
        graph_kind_list.at(i).resize(row_size);
    }
    std::vector<std::vector<double>> global_y_array(row_size);
    for (size_t i = 0; i < col_count; i++)
    {
        const auto & panel{ plot_request.panels.at(i) };
        for (size_t j = 0; j < row_size; j++)
        {
            for (const auto & curve : panel.curves)
            {
                auto graph{ root_helper::CreateGraphErrors() };
                for (size_t p = 0; p < curve.points.size(); p++)
                {
                    const auto & point{ curve.points.at(p) };
                    const auto x_value{ ScaleBiasPlotX(point.x, plot_request.x_axis_mode) };
                    const auto mean{ point.bias.mean(static_cast<int>(j)) };
                    const auto sigma{ point.bias.sigma(static_cast<int>(j)) };
                    graph->SetPoint(static_cast<int>(p), x_value, mean);
                    graph->SetPointError(static_cast<int>(p), 0.0, sigma);
                    global_y_array.at(j).emplace_back(mean);
                }
                graph_kind_list.at(i).at(j).emplace_back(curve.kind);
                graph_list.at(i).at(j).emplace_back(std::move(graph));
            }
        }
    }

    std::vector<double> x_min(col_count, 0.0);
    std::vector<double> x_max(col_count, 0.0);
    std::vector<double> y_min(row_size, 0.0);
    std::vector<double> y_max(row_size, 0.0);
    for (size_t i = 0; i < col_count; i++)
    {
        x_min.at(i) = (plot_request.flavor == BiasPlotFlavor::ModelAlphaData) ? -2.0 : -0.7;
        x_max.at(i) = (plot_request.flavor == BiasPlotFlavor::ModelAlphaData) ? 47.0 : 22.0;
        if (plot_request.x_axis_mode == BiasXAxisMode::NeighborType)
        {
            x_min.at(i) = 0.0;
            x_max.at(i) = 5.0;
        }
        if (plot_request.x_axis_mode == BiasXAxisMode::NeighborDistance)
        {
            x_min.at(i) = -2.6;
            x_max.at(i) = -0.8;
        }
    }
    for (size_t j = 0; j < row_size; j++)
    {
        auto y_range{ array_helper::ComputeScalingPercentileRangeTuple(global_y_array.at(j), 0.2, 0.005, 0.995) };
        y_min.at(j) = std::get<0>(y_range);
        y_max.at(j) = std::get<1>(y_range);
    }

    std::vector<std::vector<std::unique_ptr<TH2>>> frame(col_count);
    for (auto & frame_column : frame)
    {
        frame_column.resize(row_size);
    }
    std::vector<std::unique_ptr<TPaveText>> title_x_text(col_count);
    std::vector<std::unique_ptr<TPaveText>> title_y_text(row_size);
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto par_id{ row_size - j - 1 };
            root_helper::FindPadInCanvasPartition(canvas.get(), i, j);
            root_helper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            root_helper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ root_helper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ root_helper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            auto & frame_hist{ frame.at(static_cast<size_t>(i)).at(static_cast<size_t>(j)) };
            frame_hist = root_helper::CreateHist2D(
                Form("frame_%d_%d", i, j), "",
                500, x_min.at(static_cast<size_t>(i)), x_max.at(static_cast<size_t>(i)),
                500, y_min.at(static_cast<size_t>(par_id)), y_max.at(static_cast<size_t>(par_id)));
            root_helper::SetAxisTitleAttribute(frame_hist->GetXaxis(), 50.0f, 1.0f, 133);
            root_helper::SetAxisLabelAttribute(frame_hist->GetXaxis(), 40.0f, 0.005f, 133);
            root_helper::SetAxisTickAttribute(frame_hist->GetXaxis(), static_cast<float>(y_factor*0.08/x_factor), 505);
            root_helper::SetAxisTitleAttribute(frame_hist->GetYaxis(), 50.0f, 1.2f, 133);
            root_helper::SetAxisLabelAttribute(frame_hist->GetYaxis(), 45.0f, 0.02f, 133);
            root_helper::SetAxisTickAttribute(frame_hist->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
            if (plot_request.x_axis_mode == BiasXAxisMode::NeighborType)
            {
                ApplyNeighborTypeAxisLabels(frame_hist->GetXaxis());
            }
            root_helper::SetLineAttribute(frame_hist.get(), 1, 0);
            frame_hist->GetXaxis()->SetTitle("");
            frame_hist->GetYaxis()->SetTitle("");
            frame_hist->GetXaxis()->CenterTitle();
            frame_hist->GetYaxis()->CenterTitle();
            frame_hist->SetStats(0);
            frame_hist->Draw("Y+");

            auto & pad_graph_list{ graph_list.at(static_cast<size_t>(i)).at(static_cast<size_t>(par_id)) };
            auto & pad_graph_kind_list{ graph_kind_list.at(static_cast<size_t>(i)).at(static_cast<size_t>(par_id)) };
            for (size_t graph_index = 0; graph_index < pad_graph_list.size(); graph_index++)
            {
                const auto kind{ pad_graph_kind_list.at(graph_index) };
                auto & graph{ pad_graph_list.at(graph_index) };
                ApplyBiasCurveStyle(graph.get(), plot_request.flavor, kind);
                graph->Draw("PL3");
            }

            if (i == 0)
            {
                auto & text{ title_y_text.at(static_cast<size_t>(j)) };
                text = (plot_request.flavor == BiasPlotFlavor::NeighborType) ?
                    root_helper::CreatePaveText(-0.41, 0.30, 0.00, 0.70, "nbNDC ARC", true) :
                    root_helper::CreatePaveText(-0.68, 0.30, 0.00, 0.70, "nbNDC ARC", true);
                root_helper::SetPaveTextDefaultStyle(text.get());
                root_helper::SetPaveAttribute(text.get(), 0, 0.1);
                root_helper::SetLineAttribute(text.get(), 1, 0);
                root_helper::SetTextAttribute(text.get(), 45.0f, 133, 22);
                root_helper::SetFillAttribute(text.get(), 1001, kAzure-7, 0.5f);
                text->AddText(title_y_list.at(static_cast<size_t>(par_id)).data());
                text->Draw();
            }

            if (j == row_size - 1)
            {
                auto & text{ title_x_text.at(static_cast<size_t>(i)) };
                text = root_helper::CreatePaveText(0.10, 0.95, 0.90, 1.15, "nbNDC ARC", true);
                root_helper::SetPaveTextDefaultStyle(text.get());
                root_helper::SetPaveAttribute(text.get(), 0, 0.2);
                root_helper::SetTextAttribute(text.get(), 45.0f, 133, 22);
                root_helper::SetFillAttribute(text.get(), 1001, kRed+1, 0.5f);
                text->AddText(plot_request.panels.at(static_cast<size_t>(i)).label.data());
                if (plot_request.flavor != BiasPlotFlavor::NeighborType) text->Draw();
            }
        }
    }

    canvas->cd();
    auto pad_extra0{ root_helper::CreatePad("pad_extra0","", 0.02, 0.92, 0.98, 1.00) };
    pad_extra0->Draw();
    pad_extra0->cd();
    root_helper::SetPadDefaultStyle(pad_extra0.get());
    root_helper::SetFillAttribute(pad_extra0.get(), 4000);
    auto legend{ root_helper::CreateLegend(0.0, 0.0, 1.0, 1.0, false) };
    root_helper::SetLegendDefaultStyle(legend.get());
    root_helper::SetFillAttribute(legend.get(), 4000);
    root_helper::SetTextAttribute(legend.get(), 40.0f, 133, 12, 0.0);
    legend->SetMargin(0.25f);
    legend->SetNColumns(3);
    for (const auto kind : GetBiasLegendOrder(plot_request.flavor))
    {
        auto * graph{ FindBiasGraph(graph_list.at(0).at(0), graph_kind_list.at(0).at(0), kind) };
        if (graph != nullptr)
        {
            legend->AddEntry(graph, GetBiasLegendLabel(plot_request.flavor, kind).data(), "plf");
        }
    }
    legend->Draw();

    canvas->cd();
    auto pad_extra1{ root_helper::CreatePad("pad_extra1","", 0.17, 0.00, 0.92, 0.06) };
    pad_extra1->Draw();
    pad_extra1->cd();
    root_helper::SetPadDefaultStyle(pad_extra1.get());
    root_helper::SetFillAttribute(pad_extra1.get(), 4000);
    auto bottom_title_text{ root_helper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    root_helper::SetPaveTextDefaultStyle(bottom_title_text.get());
    root_helper::SetFillAttribute(bottom_title_text.get(), 4000);
    root_helper::SetTextAttribute(bottom_title_text.get(), 45.0f, 133, 22);
    if (plot_request.x_axis_mode == BiasXAxisMode::NeighborType)
    {
        bottom_title_text->AddText("Atom Type");
    }
    else if (plot_request.x_axis_mode == BiasXAxisMode::NeighborDistance)
    {
        bottom_title_text->AddText("Distance to Neighbor Atom (Angstrom)");
    }
    else
    {
        bottom_title_text->AddText("Data Contamination Ratio (%)");
    }
    bottom_title_text->Draw();

    canvas->cd();
    auto pad_extra2{ root_helper::CreatePad("pad_extra2","", 0.96, 0.10, 1.00, 0.86) };
    pad_extra2->Draw();
    pad_extra2->cd();
    root_helper::SetPadDefaultStyle(pad_extra2.get());
    root_helper::SetFillAttribute(pad_extra2.get(), 4000);
    auto right_title_text{ root_helper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    root_helper::SetPaveTextDefaultStyle(right_title_text.get());
    root_helper::SetFillAttribute(right_title_text.get(), 4000);
    root_helper::SetTextAttribute(right_title_text.get(), 50.0f, 133, 22);
    right_title_text->AddText("Normalized Bias");
    auto text{ right_title_text->GetLineWith("Bias") };
    text->SetTextAngle(90.0f);
    right_title_text->Draw();

    root_helper::PrintCanvasPad(canvas.get(), file_path);
    root_helper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path.string());
    #else
    (void)request;
    #endif
}

void SaveMemberOutlierBiasPlot(
    const RHBMTestRequest & request,
    const BiasPlotRequest & plot_request)
{
    Logger::Log(LogLevel::Info, " RHBMTestCommand::SaveMemberOutlierBiasPlot");

    #ifdef HAVE_ROOT
    auto file_path{ request.output_dir / plot_request.output_name };
    std::vector<std::string> title_y_list{
        "Amplitude #font[2]{A}", "Width #tau"
    };

    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 2 };
    const int row_size{ 2 };
    auto canvas{ root_helper::CreateCanvas("test","", 1500, 750) };
    root_helper::SetCanvasDefaultStyle(canvas.get());
    root_helper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.19f, 0.10f, 0.12f, 0.14f, 0.01f, 0.01f
    );
    root_helper::PrintCanvasOpen(canvas.get(), file_path);

    std::vector<std::unique_ptr<TGraphErrors>> graph_list[col_size][row_size];
    std::vector<BiasCurveKind> graph_kind_list[col_size][row_size];
    std::vector<double> global_y_array;
    for (size_t i = 0; i < col_size; i++)
    {
        const auto & panel{ plot_request.panels.at(i) };
        for (size_t j = 0; j < row_size; j++)
        {
            for (const auto & curve : panel.curves)
            {
                auto graph{ root_helper::CreateGraphErrors() };
                for (size_t p = 0; p < curve.points.size(); p++)
                {
                    const auto & point{ curve.points.at(p) };
                    const auto x_value{ ScaleBiasPlotX(point.x, plot_request.x_axis_mode) };
                    const auto mean{ point.bias.mean(static_cast<int>(j)) };
                    const auto sigma{ point.bias.sigma(static_cast<int>(j)) };
                    graph->SetPoint(static_cast<int>(p), x_value, mean);
                    graph->SetPointError(static_cast<int>(p), 0.0, sigma);
                    global_y_array.emplace_back(mean);
                }
                graph_kind_list[i][j].emplace_back(curve.kind);
                graph_list[i][j].emplace_back(std::move(graph));
            }
        }
    }

    double x_min[col_size]{ 0.0 };
    double x_max[col_size]{ 0.0 };
    double y_min[row_size]{ 0.0 };
    double y_max[row_size]{ 0.0 };
    for (size_t i = 0; i < col_size; i++)
    {
        x_min[i] = (plot_request.flavor == BiasPlotFlavor::ModelAlphaMember) ? -2.0 : -0.7;
        x_max[i] = (plot_request.flavor == BiasPlotFlavor::ModelAlphaMember) ? 47.0 : 22.0;
    }
    auto y_range{ array_helper::ComputeScalingRangeTuple(global_y_array, 0.3) };
    for (size_t j = 0; j < row_size; j++)
    {
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> title_x_text[col_size];
    std::unique_ptr<TPaveText> title_y_text[row_size];
    for (int i = 0; i < col_size; i++)
    {
        for (int j = 0; j < row_size; j++)
        {
            auto par_id{ row_size - j - 1 };
            root_helper::FindPadInCanvasPartition(canvas.get(), i, j);
            root_helper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
            root_helper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
            auto x_factor{ root_helper::GetPadXfactorInCanvasPartition(canvas.get(), gPad) };
            auto y_factor{ root_helper::GetPadYfactorInCanvasPartition(canvas.get(), gPad) };
            frame[i][j] = root_helper::CreateHist2D(Form("frame_%d_%d", i, j),"", 500, x_min[i], x_max[i], 500, y_min[par_id], y_max[par_id]);
            root_helper::SetAxisTitleAttribute(frame[i][j]->GetXaxis(), 50.0f, 1.0f, 133);
            root_helper::SetAxisLabelAttribute(frame[i][j]->GetXaxis(), 40.0f, 0.005f, 133);
            root_helper::SetAxisTickAttribute(frame[i][j]->GetXaxis(), static_cast<float>(y_factor*0.08/x_factor), 505);
            root_helper::SetAxisTitleAttribute(frame[i][j]->GetYaxis(), 50.0f, 1.2f, 133);
            root_helper::SetAxisLabelAttribute(frame[i][j]->GetYaxis(), 45.0f, 0.02f, 133);
            root_helper::SetAxisTickAttribute(frame[i][j]->GetYaxis(), static_cast<float>(x_factor*0.05/y_factor), 505);
            root_helper::SetLineAttribute(frame[i][j].get(), 1, 0);
            frame[i][j]->GetXaxis()->SetTitle("");
            frame[i][j]->GetYaxis()->SetTitle("");
            frame[i][j]->GetXaxis()->CenterTitle();
            frame[i][j]->GetYaxis()->CenterTitle();
            frame[i][j]->SetStats(0);
            frame[i][j]->Draw("Y+");

            for (size_t graph_index = 0; graph_index < graph_list[i][par_id].size(); graph_index++)
            {
                const auto kind{ graph_kind_list[i][par_id].at(graph_index) };
                auto & graph{ graph_list[i][par_id].at(graph_index) };
                ApplyBiasCurveStyle(graph.get(), plot_request.flavor, kind);
                graph->Draw("PL3");
            }

            if (i == 0)
            {
                title_y_text[j] = root_helper::CreatePaveText(-0.52, 0.30, 0.00, 0.70, "nbNDC ARC", true);
                root_helper::SetPaveTextDefaultStyle(title_y_text[j].get());
                root_helper::SetPaveAttribute(title_y_text[j].get(), 0, 0.1);
                root_helper::SetLineAttribute(title_y_text[j].get(), 1, 0);
                root_helper::SetTextAttribute(title_y_text[j].get(), 45.0f, 133, 22);
                root_helper::SetFillAttribute(title_y_text[j].get(), 1001, kAzure-7, 0.5f);
                title_y_text[j]->AddText(title_y_list[static_cast<size_t>(par_id)].data());
                title_y_text[j]->Draw();
            }

            if (j == row_size - 1)
            {
                title_x_text[i] = root_helper::CreatePaveText(0.10, 0.95, 0.90, 1.15, "nbNDC ARC", true);
                root_helper::SetPaveTextDefaultStyle(title_x_text[i].get());
                root_helper::SetPaveAttribute(title_x_text[i].get(), 0, 0.2);
                root_helper::SetTextAttribute(title_x_text[i].get(), 45.0f, 133, 22);
                root_helper::SetFillAttribute(title_x_text[i].get(), 1001, kRed+1, 0.5f);
                title_x_text[i]->AddText(plot_request.panels.at(static_cast<size_t>(i)).label.data());
                title_x_text[i]->Draw();
            }
        }
    }

    canvas->cd();
    auto pad_extra0{ root_helper::CreatePad("pad_extra0","", 0.02, 0.92, 0.98, 1.00) };
    pad_extra0->Draw();
    pad_extra0->cd();
    root_helper::SetPadDefaultStyle(pad_extra0.get());
    root_helper::SetFillAttribute(pad_extra0.get(), 4000);
    auto legend{ root_helper::CreateLegend(0.0, 0.0, 1.0, 1.0, false) };
    root_helper::SetLegendDefaultStyle(legend.get());
    root_helper::SetFillAttribute(legend.get(), 4000);
    root_helper::SetTextAttribute(legend.get(), 40.0f, 133, 12, 0.0);
    legend->SetMargin(0.25f);
    legend->SetNColumns(3);
    for (const auto kind : GetBiasLegendOrder(plot_request.flavor))
    {
        auto * graph{ FindBiasGraph(graph_list[0][0], graph_kind_list[0][0], kind) };
        if (graph != nullptr)
        {
            legend->AddEntry(graph, GetBiasLegendLabel(plot_request.flavor, kind).data(), "plf");
        }
    }
    legend->Draw();

    canvas->cd();
    auto pad_extra1{ root_helper::CreatePad("pad_extra1","", 0.19, 0.00, 0.90, 0.06) };
    pad_extra1->Draw();
    pad_extra1->cd();
    root_helper::SetPadDefaultStyle(pad_extra1.get());
    root_helper::SetFillAttribute(pad_extra1.get(), 4000);
    auto bottom_title_text{ root_helper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    root_helper::SetPaveTextDefaultStyle(bottom_title_text.get());
    root_helper::SetFillAttribute(bottom_title_text.get(), 4000);
    root_helper::SetTextAttribute(bottom_title_text.get(), 45.0f, 133, 22);
    bottom_title_text->AddText("Member Contamination Ratio (%)");
    bottom_title_text->Draw();

    canvas->cd();
    auto pad_extra2{ root_helper::CreatePad("pad_extra2","", 0.96, 0.10, 1.00, 0.86) };
    pad_extra2->Draw();
    pad_extra2->cd();
    root_helper::SetPadDefaultStyle(pad_extra2.get());
    root_helper::SetFillAttribute(pad_extra2.get(), 4000);
    auto right_title_text{ root_helper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false) };
    root_helper::SetPaveTextDefaultStyle(right_title_text.get());
    root_helper::SetFillAttribute(right_title_text.get(), 4000);
    root_helper::SetTextAttribute(right_title_text.get(), 50.0f, 133, 22);
    right_title_text->AddText("Normalized Bias");
    auto text{ right_title_text->GetLineWith("Bias") };
    text->SetTextAngle(90.0f);
    right_title_text->Draw();

    root_helper::PrintCanvasPad(canvas.get(), file_path);
    root_helper::PrintCanvasClose(canvas.get(), file_path);
    Logger::Log(LogLevel::Info, " Output file: " + file_path.string());
    #else
    (void)request;
    (void)plot_request;
    #endif
}

void SaveAtomSamplingSummary(
    const RHBMTestRequest & request,
    const std::string & name,
    const std::vector<LocalPotentialSampleList> & sampling_entries_list,
    const std::vector<double> & distance_list)
{
    #ifdef HAVE_ROOT
    auto file_path{ request.output_dir / name };
    Logger::Log(LogLevel::Info, " RHBMTestCommand::SaveAtomSamplingSummary");

    gStyle->SetLineScalePS(2.0);
    gStyle->SetGridColor(kGray);

    auto canvas{ root_helper::CreateCanvas("test","", 1000, 750) };
    root_helper::SetCanvasDefaultStyle(canvas.get());
    root_helper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size{ 1 };

    std::unique_ptr<TPad> pad[pad_size];
    pad[0] = root_helper::CreatePad("pad0","", 0.00, 0.00, 1.00, 1.00);

    size_t count{ 0 };
    for (const auto & sampling_entries : sampling_entries_list)
    {
        auto neighbor_distance{ distance_list.at(count) };
        auto data_graph{ CreateDistanceToResponseGraph(sampling_entries) };
        auto data_hist{ CreateDistanceToResponseHistogram(sampling_entries, 40) };

        canvas->cd();
        for (int i = 0; i < pad_size; i++)
        {
            root_helper::SetPadDefaultStyle(pad[i].get());
            pad[i]->Draw();
        }

        pad[0]->cd();
        root_helper::SetPadMarginInCanvas(gPad, 0.16, 0.05, 0.15, 0.05);
        auto frame{ root_helper::CreateHist2D("hist_0","", 100, 0.0, 1.0, 100, 0.0, 1.0) };
        root_helper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        root_helper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
        root_helper::SetAxisTitleAttribute(frame->GetXaxis(), 50.0f, 0.9f);
        root_helper::SetAxisLabelAttribute(frame->GetXaxis(), 50.0f, 0.005f, 133);
        root_helper::SetAxisTickAttribute(frame->GetXaxis(), 0.06f, 505);
        root_helper::SetAxisTitleAttribute(frame->GetYaxis(), 60.0f, 1.25f);
        root_helper::SetAxisLabelAttribute(frame->GetYaxis(), 50.0f, 0.01f, 133);
        root_helper::SetAxisTickAttribute(frame->GetYaxis(), 0.03f, 506);
        root_helper::SetLineAttribute(frame.get(), 1, 0);
        frame->SetStats(0);
        frame->GetYaxis()->CenterTitle();
        frame->GetXaxis()->SetTitle("Radial Distance #[]{#AA}");
        frame->GetYaxis()->SetTitle("Response");
        auto x_min{ 0.00 };
        auto x_max{ 4.00 };
        auto y_min{ 0.00 };
        auto y_max{ 1.00 };
        frame->GetXaxis()->SetLimits(x_min, x_max);
        frame->GetYaxis()->SetLimits(y_min, y_max);
        frame->SetStats(0);
        frame->Draw();
        root_helper::SetMarkerAttribute(data_graph.get(), 20, 0.8f, kAzure-7, 0.5f);
        data_graph->Draw("P");

        data_hist->SetStats(0);
        data_hist->SetBarWidth(1.0);
        root_helper::SetFillAttribute(data_hist.get(), 1001, kGray, 0.3f);
        root_helper::SetLineAttribute(data_hist.get(), 1, 2, kGray+2);
        data_hist->Draw("CANDLE2 SAME");

        auto component_text{ root_helper::CreatePaveText(0.00, 0.00, 1.00, 1.00, "nbNDC ARC", false) };
        root_helper::SetPaveTextMarginInCanvas(gPad, component_text.get(), 0.25, 0.02, 0.85, 0.02);
        root_helper::SetPaveTextDefaultStyle(component_text.get());
        root_helper::SetPaveAttribute(component_text.get(), 0, 0.1);
        root_helper::SetFillAttribute(component_text.get(), 1001, kAzure-7);
        root_helper::SetTextAttribute(component_text.get(), 75.0f, 103, 22, 0.0, kYellow-10);
        component_text->AddText(Form("Distance = %.1f #AA", neighbor_distance));
        component_text->Draw();

        count++;
        root_helper::PrintCanvasPad(canvas.get(), file_path);
    }

    root_helper::PrintCanvasClose(canvas.get(), file_path);
    #else
    (void)request;
    (void)name;
    (void)sampling_entries_list;
    (void)distance_list;
    Logger::Log(LogLevel::Info, " RHBMTestCommand::SaveAtomSamplingSummary");
    #endif
}

} // namespace rhbm_gem::command_detail::rhbm_test_plotting
