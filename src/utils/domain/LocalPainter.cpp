#include <rhbm_gem/utils/domain/LocalPainter.hpp>

#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <system_error>
#include <tuple>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>

#include <TCanvas.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLegend.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#endif

namespace rhbm_gem {
namespace {

constexpr double kMinimumAxisSpan{ 1.0 };

PlotResult BuildPlotResult(
    PlotStatus status,
    const std::filesystem::path & output_path,
    const std::string & message)
{
    return PlotResult{ status, output_path, message };
}

PlotResult FailPlot(
    PlotStatus status,
    const std::filesystem::path & output_path,
    LogLevel level,
    const std::string & message)
{
    Logger::Log(level, message);
    return BuildPlotResult(status, output_path, message);
}

bool IsFiniteRange(const std::pair<double, double> & range)
{
    return std::isfinite(range.first) &&
        std::isfinite(range.second) &&
        range.first < range.second;
}

double ExpandSpanWithFallback(double min_value, double max_value, double padding_ratio)
{
    const double span{ max_value - min_value };
    if (span > 0.0)
    {
        return span * padding_ratio;
    }
    const double anchor{
        std::max({ std::abs(min_value), std::abs(max_value), kMinimumAxisSpan })
    };
    return anchor * std::max(padding_ratio, 0.05);
}

std::pair<double, double> ExpandDegenerateRange(double min_value, double max_value)
{
    if (min_value < max_value)
    {
        return { min_value, max_value };
    }
    const double padding{ ExpandSpanWithFallback(min_value, max_value, 0.05) };
    return { min_value - padding, max_value + padding };
}

std::optional<std::pair<double, double>> ResolveAutoRange(
    const std::vector<double> & values,
    double padding_ratio)
{
    if (values.empty())
    {
        return std::nullopt;
    }

    const auto range_tuple{ ArrayStats<double>::ComputeRangeTuple(values) };
    const double min_value{ std::get<0>(range_tuple) };
    const double max_value{ std::get<1>(range_tuple) };
    if (!std::isfinite(min_value) || !std::isfinite(max_value))
    {
        return std::nullopt;
    }

    const double padding{ ExpandSpanWithFallback(min_value, max_value, padding_ratio) };
    return std::make_pair(min_value - padding, max_value + padding);
}

std::vector<double> CollectPanelYValues(const LinePlotPanel & panel)
{
    std::vector<double> values;
    for (const auto & series : panel.series)
    {
        values.insert(values.end(), series.y_values.begin(), series.y_values.end());
    }
    return values;
}

std::vector<double> CollectAllXValues(const LinePlotRequest & request)
{
    std::vector<double> values;
    for (const auto & panel : request.panels)
    {
        for (const auto & series : panel.series)
        {
            values.insert(values.end(), series.x_values.begin(), series.x_values.end());
        }
    }
    return values;
}

std::optional<std::pair<double, double>> ResolveLineAxisRange(
    const AxisSpec & axis,
    const std::vector<double> & values,
    double padding_ratio)
{
    if (axis.range.has_value())
    {
        if (!IsFiniteRange(*axis.range))
        {
            return std::nullopt;
        }
        return axis.range;
    }
    if (!axis.auto_range)
    {
        return std::nullopt;
    }
    return ResolveAutoRange(values, padding_ratio);
}

PlotResult EnsureParentDirectory(const std::filesystem::path & output_path)
{
    if (output_path.empty())
    {
        return FailPlot(
            PlotStatus::INVALID_INPUT,
            output_path,
            LogLevel::Error,
            "local_painter requires a non-empty output path.");
    }

    const auto parent_path{ output_path.parent_path() };
    if (parent_path.empty())
    {
        return BuildPlotResult(PlotStatus::SUCCESS, output_path, "");
    }

    std::error_code ec;
    std::filesystem::create_directories(parent_path, ec);
    if (ec)
    {
        return FailPlot(
            PlotStatus::OUTPUT_ERROR,
            output_path,
            LogLevel::Error,
            "local_painter failed to create output directory '" +
                parent_path.string() + "': " + ec.message());
    }
    return BuildPlotResult(PlotStatus::SUCCESS, output_path, "");
}

PlotResult ValidateLinePlotRequest(const LinePlotRequest & request)
{
    const auto directory_result{ EnsureParentDirectory(request.output_path) };
    if (!directory_result.Succeeded())
    {
        return directory_result;
    }

    if (request.panels.empty())
    {
        return FailPlot(
            PlotStatus::INVALID_INPUT,
            request.output_path,
            LogLevel::Error,
            "local_painter::SaveLinePlot() requires at least one plot panel.");
    }
    if (request.canvas_width <= 0 || request.canvas_height_per_panel <= 0)
    {
        return FailPlot(
            PlotStatus::INVALID_INPUT,
            request.output_path,
            LogLevel::Error,
            "local_painter::SaveLinePlot() requires positive canvas dimensions.");
    }

    if (!request.x_axis.auto_range && !request.x_axis.range.has_value())
    {
        return FailPlot(
            PlotStatus::INVALID_INPUT,
            request.output_path,
            LogLevel::Error,
            "local_painter::SaveLinePlot() requires an explicit x-axis range when auto_range is false.");
    }

    for (std::size_t panel_index = 0; panel_index < request.panels.size(); ++panel_index)
    {
        const auto & panel{ request.panels[panel_index] };
        if (panel.series.empty())
        {
            return FailPlot(
                PlotStatus::INVALID_INPUT,
                request.output_path,
                LogLevel::Error,
                "local_painter::SaveLinePlot() panel " + std::to_string(panel_index) +
                    " does not contain any series.");
        }
        if (!panel.y_axis.auto_range && !panel.y_axis.range.has_value())
        {
            return FailPlot(
                PlotStatus::INVALID_INPUT,
                request.output_path,
                LogLevel::Error,
                "local_painter::SaveLinePlot() panel " + std::to_string(panel_index) +
                    " requires an explicit y-axis range when auto_range is false.");
        }

        for (std::size_t series_index = 0; series_index < panel.series.size(); ++series_index)
        {
            const auto & series{ panel.series[series_index] };
            if (series.x_values.empty() || series.y_values.empty())
            {
                return FailPlot(
                    PlotStatus::INVALID_INPUT,
                    request.output_path,
                    LogLevel::Error,
                    "local_painter::SaveLinePlot() panel " + std::to_string(panel_index) +
                        " series " + std::to_string(series_index) + " is empty.");
            }
            if (series.x_values.size() != series.y_values.size())
            {
                return FailPlot(
                    PlotStatus::INVALID_INPUT,
                    request.output_path,
                    LogLevel::Error,
                    "local_painter::SaveLinePlot() panel " + std::to_string(panel_index) +
                        " series " + std::to_string(series_index) +
                        " has mismatched x/y lengths.");
            }
        }
    }

    return BuildPlotResult(PlotStatus::SUCCESS, request.output_path, "");
}

PlotResult ValidateHeatmapRequest(const HeatmapRequest & request)
{
    const auto directory_result{ EnsureParentDirectory(request.output_path) };
    if (!directory_result.Succeeded())
    {
        return directory_result;
    }

    if (request.values.rows() <= 0 || request.values.cols() <= 0)
    {
        return FailPlot(
            PlotStatus::INVALID_INPUT,
            request.output_path,
            LogLevel::Error,
            "local_painter::SaveHeatmap() requires a non-empty matrix.");
    }
    if (request.canvas_width <= 0 || request.canvas_height <= 0)
    {
        return FailPlot(
            PlotStatus::INVALID_INPUT,
            request.output_path,
            LogLevel::Error,
            "local_painter::SaveHeatmap() requires positive canvas dimensions.");
    }

    for (const auto & item : {
             std::make_pair(&request.x_axis, std::string("x-axis")),
             std::make_pair(&request.y_axis, std::string("y-axis")),
             std::make_pair(&request.z_axis, std::string("z-axis")) })
    {
        if (item.first->range.has_value() && !IsFiniteRange(*item.first->range))
        {
            return FailPlot(
                PlotStatus::INVALID_INPUT,
                request.output_path,
                LogLevel::Error,
                "local_painter::SaveHeatmap() received an invalid " + item.second + " range.");
        }
    }
    if (!request.z_axis.auto_range && !request.z_axis.range.has_value())
    {
        return FailPlot(
            PlotStatus::INVALID_INPUT,
            request.output_path,
            LogLevel::Error,
            "local_painter::SaveHeatmap() requires an explicit z-axis range when auto_range is false.");
    }

    return BuildPlotResult(PlotStatus::SUCCESS, request.output_path, "");
}

#ifdef HAVE_ROOT

short ResolveSeriesColor(const LineSeries & series, std::size_t series_index)
{
    static constexpr short kDefaultColors[]{
        static_cast<short>(kRed + 1),
        static_cast<short>(kAzure + 2),
        static_cast<short>(kGreen + 2),
        static_cast<short>(kOrange + 7),
        static_cast<short>(kMagenta + 1),
        static_cast<short>(kCyan + 2)
    };

    if (series.color_index.has_value())
    {
        return static_cast<short>(*series.color_index);
    }
    return kDefaultColors[series_index % (sizeof(kDefaultColors) / sizeof(kDefaultColors[0]))];
}

PlotResult RenderLinePlotWithRoot(const LinePlotRequest & request)
{
    gStyle->SetLineScalePS(2);
    gStyle->SetGridColor(kGray);

    const int panel_count{ static_cast<int>(request.panels.size()) };
    const int canvas_height{ std::max(1, panel_count) * request.canvas_height_per_panel };
    auto canvas{
        ROOTHelper::CreateCanvas("local_line_plot", "", request.canvas_width, canvas_height)
    };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::SetCanvasPartition(
        canvas.get(), 1, panel_count, 0.18f, 0.10f, 0.14f, 0.06f, 0.02f, 0.01f);
    ROOTHelper::PrintCanvasOpen(canvas.get(), request.output_path.string());

    const auto x_range{
        ResolveLineAxisRange(
            request.x_axis,
            CollectAllXValues(request),
            request.x_auto_padding_ratio)
    };
    if (!x_range.has_value())
    {
        return FailPlot(
            PlotStatus::INVALID_INPUT,
            request.output_path,
            LogLevel::Error,
            "local_painter::SaveLinePlot() could not determine a valid x-axis range.");
    }

    std::vector<std::unique_ptr<TH2>> frames(static_cast<std::size_t>(panel_count));
    std::vector<std::vector<std::unique_ptr<TGraphErrors>>> panel_graphs(
        static_cast<std::size_t>(panel_count));
    std::vector<std::unique_ptr<TLegend>> legends(static_cast<std::size_t>(panel_count));
    std::vector<std::unique_ptr<TPaveText>> panel_labels(static_cast<std::size_t>(panel_count));
    std::unique_ptr<TPad> shared_axis_pad;
    std::unique_ptr<TPaveText> shared_axis_title;

    for (int panel_index = 0; panel_index < panel_count; ++panel_index)
    {
        const auto & panel{ request.panels[static_cast<std::size_t>(panel_index)] };
        const auto y_range{
            ResolveLineAxisRange(
                panel.y_axis,
                CollectPanelYValues(panel),
                request.y_auto_padding_ratio)
        };
        if (!y_range.has_value())
        {
            return FailPlot(
                PlotStatus::INVALID_INPUT,
                request.output_path,
                LogLevel::Error,
                "local_painter::SaveLinePlot() could not determine a valid y-axis range for panel "
                    + std::to_string(panel_index) + ".");
        }

        const int canvas_row{ panel_count - panel_index - 1 };
        ROOTHelper::FindPadInCanvasPartition(canvas.get(), 0, canvas_row);
        ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
        ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);

        const std::string frame_name{ "local_line_frame_" + std::to_string(panel_index) };
        frames[static_cast<std::size_t>(panel_index)] = ROOTHelper::CreateHist2D(
            frame_name,
            request.title,
            500, x_range->first, x_range->second,
            500, y_range->first, y_range->second);

        auto * frame{ frames[static_cast<std::size_t>(panel_index)].get() };
        ROOTHelper::SetAxisTitleAttribute(frame->GetXaxis(), 45.0f, 1.0f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame->GetXaxis(), 38.0f, 0.01f, 133);
        ROOTHelper::SetAxisTickAttribute(frame->GetXaxis(), 0.03f, 505);
        ROOTHelper::SetAxisTitleAttribute(frame->GetYaxis(), 45.0f, 1.3f, 133);
        ROOTHelper::SetAxisLabelAttribute(frame->GetYaxis(), 38.0f, 0.01f, 133);
        ROOTHelper::SetAxisTickAttribute(frame->GetYaxis(), 0.03f, 505);
        ROOTHelper::SetLineAttribute(frame, 1, 0);
        frame->GetXaxis()->SetTitle(
            panel_index == panel_count - 1 ? request.x_axis.title.c_str() : "");
        frame->GetYaxis()->SetTitle(panel.y_axis.title.c_str());
        frame->GetXaxis()->CenterTitle();
        frame->GetYaxis()->CenterTitle();
        frame->SetStats(0);
        frame->Draw();

        auto & graphs{ panel_graphs[static_cast<std::size_t>(panel_index)] };
        graphs.reserve(panel.series.size());
        bool draw_legend{ panel.series.size() > 1 };
        for (std::size_t series_index = 0; series_index < panel.series.size(); ++series_index)
        {
            const auto & series{ panel.series[series_index] };
            auto graph{ ROOTHelper::CreateGraphErrors() };
            for (int point_index = 0; point_index < static_cast<int>(series.x_values.size()); ++point_index)
            {
                graph->SetPoint(
                    point_index,
                    series.x_values[static_cast<std::size_t>(point_index)],
                    series.y_values[static_cast<std::size_t>(point_index)]);
            }
            const short color{ ResolveSeriesColor(series, series_index) };
            ROOTHelper::SetMarkerAttribute(graph.get(), 20, 1.3f, color);
            ROOTHelper::SetLineAttribute(graph.get(), 1, 2, color);
            graph->Draw("PL SAME");
            if (series.name.empty())
            {
                draw_legend = false;
            }
            graphs.emplace_back(std::move(graph));
        }

        if (draw_legend)
        {
            auto legend{ ROOTHelper::CreateLegend(0.62, 0.72, 0.92, 0.90, false) };
            ROOTHelper::SetLegendDefaultStyle(legend.get());
            ROOTHelper::SetFillAttribute(legend.get(), 1001, 0, 0.0f);
            for (std::size_t series_index = 0; series_index < panel.series.size(); ++series_index)
            {
                legend->AddEntry(
                    graphs[series_index].get(),
                    panel.series[series_index].name.c_str(),
                    "lp");
            }
            legend->Draw();
            legends[static_cast<std::size_t>(panel_index)] = std::move(legend);
        }

        if (!panel.label.empty())
        {
            auto label{ ROOTHelper::CreatePaveText(0.02, 0.78, 0.34, 0.94, "nbNDC ARC", false) };
            ROOTHelper::SetPaveTextDefaultStyle(label.get());
            ROOTHelper::SetPaveAttribute(label.get(), 0, 0.1);
            ROOTHelper::SetLineAttribute(label.get(), 1, 0);
            ROOTHelper::SetTextAttribute(label.get(), 30.0f, 133, 12);
            ROOTHelper::SetFillAttribute(label.get(), 1001, kAzure - 7, 0.45f);
            label->AddText(panel.label.c_str());
            label->Draw();
            panel_labels[static_cast<std::size_t>(panel_index)] = std::move(label);
        }
    }

    if (!request.shared_y_axis_title.empty())
    {
        canvas->cd();
        shared_axis_pad = ROOTHelper::CreatePad(
            "local_line_plot_axis_title", "", 0.94, 0.10, 1.00, 0.96);
        shared_axis_pad->Draw();
        shared_axis_pad->cd();
        ROOTHelper::SetPadDefaultStyle(shared_axis_pad.get());
        ROOTHelper::SetFillAttribute(shared_axis_pad.get(), 4000);
        shared_axis_title = ROOTHelper::CreatePaveText(0.0, 0.0, 1.0, 1.0, "nbNDC", false);
        ROOTHelper::SetPaveTextDefaultStyle(shared_axis_title.get());
        ROOTHelper::SetFillAttribute(shared_axis_title.get(), 4000);
        ROOTHelper::SetTextAttribute(shared_axis_title.get(), 42.0f, 133, 22);
        shared_axis_title->AddText(request.shared_y_axis_title.c_str());
        if (auto * text{ shared_axis_title->GetLineWith(request.shared_y_axis_title.c_str()) })
        {
            text->SetTextAngle(90.0f);
        }
        shared_axis_title->Draw();
    }

    ROOTHelper::PrintCanvasPad(canvas.get(), request.output_path.string());
    ROOTHelper::PrintCanvasClose(canvas.get(), request.output_path.string());

    if (!std::filesystem::exists(request.output_path))
    {
        return FailPlot(
            PlotStatus::OUTPUT_ERROR,
            request.output_path,
            LogLevel::Error,
            "local_painter::SaveLinePlot() completed without producing '" +
                request.output_path.string() + "'.");
    }

    Logger::Log(LogLevel::Info, " Output file: " + request.output_path.string());
    return BuildPlotResult(PlotStatus::SUCCESS, request.output_path, "");
}

PlotResult RenderHeatmapWithRoot(const HeatmapRequest & request)
{
    gStyle->SetLineScalePS(2);
    gStyle->SetGridColor(kGray);
    gStyle->SetPalette(kLightTemperature);
    gStyle->SetNumberContours(50);

    auto canvas{
        ROOTHelper::CreateCanvas("local_heatmap", "", request.canvas_width, request.canvas_height)
    };
    ROOTHelper::SetCanvasDefaultStyle(canvas.get());
    ROOTHelper::PrintCanvasOpen(canvas.get(), request.output_path.string());

    const std::pair<double, double> x_range{
        request.x_axis.range.value_or(std::make_pair(0.0, static_cast<double>(request.values.cols())))
    };
    const std::pair<double, double> y_range{
        request.y_axis.range.value_or(std::make_pair(0.0, static_cast<double>(request.values.rows())))
    };

    std::vector<double> flat_values;
    flat_values.reserve(static_cast<std::size_t>(request.values.size()));
    for (Eigen::Index row = 0; row < request.values.rows(); ++row)
    {
        for (Eigen::Index col = 0; col < request.values.cols(); ++col)
        {
            flat_values.emplace_back(request.values(row, col));
        }
    }

    const auto z_range_spec{
        request.z_axis.range.value_or(
            ResolveAutoRange(flat_values, 0.0).value_or(std::make_pair(0.0, 1.0)))
    };
    const auto z_range{ ExpandDegenerateRange(z_range_spec.first, z_range_spec.second) };

    auto hist_data{ ROOTHelper::CreateHist2D(
        "local_heatmap_data",
        request.title,
        static_cast<int>(request.values.cols()), x_range.first, x_range.second,
        static_cast<int>(request.values.rows()), y_range.first, y_range.second) };

    for (Eigen::Index row = 0; row < request.values.rows(); ++row)
    {
        for (Eigen::Index col = 0; col < request.values.cols(); ++col)
        {
            hist_data->SetBinContent(
                static_cast<int>(col) + 1,
                static_cast<int>(row) + 1,
                request.values(row, col));
        }
    }

    canvas->cd();
    ROOTHelper::SetPadLayout(gPad, 1, 1, 0, 0, 0, 0);
    ROOTHelper::SetPadFrameAttribute(gPad, 0, 0, 4000, 0, 0, 0);
    ROOTHelper::SetPadMarginInCanvas(gPad, 0.15, 0.20, 0.15, 0.10);
    ROOTHelper::SetAxisTitleAttribute(hist_data->GetXaxis(), 60.0f, 1.0f, 133);
    ROOTHelper::SetAxisLabelAttribute(hist_data->GetXaxis(), 60.0f, 0.01f, 133);
    ROOTHelper::SetAxisTickAttribute(hist_data->GetXaxis(), 0.04f, 505);
    ROOTHelper::SetAxisTitleAttribute(hist_data->GetYaxis(), 60.0f, 1.5f, 133);
    ROOTHelper::SetAxisLabelAttribute(hist_data->GetYaxis(), 60.0f, 0.01f, 133);
    ROOTHelper::SetAxisTickAttribute(hist_data->GetYaxis(), 0.04f, 505);
    ROOTHelper::SetAxisTitleAttribute(hist_data->GetZaxis(), 60.0f, 1.3f, 133);
    ROOTHelper::SetAxisLabelAttribute(hist_data->GetZaxis(), 60.0f, 0.005f, 133);
    ROOTHelper::SetAxisTickAttribute(hist_data->GetZaxis(), 0.02f, 505);
    ROOTHelper::SetLineAttribute(hist_data.get(), 1, 0);
    hist_data->GetXaxis()->SetTitle(request.x_axis.title.c_str());
    hist_data->GetYaxis()->SetTitle(request.y_axis.title.c_str());
    hist_data->GetZaxis()->SetTitle(request.z_axis.title.c_str());
    hist_data->GetXaxis()->CenterTitle();
    hist_data->GetYaxis()->CenterTitle();
    hist_data->GetZaxis()->CenterTitle();
    hist_data->GetZaxis()->SetRangeUser(z_range.first, z_range.second);
    hist_data->SetStats(0);
    hist_data->Draw("COLZ");

    ROOTHelper::PrintCanvasPad(canvas.get(), request.output_path.string());
    ROOTHelper::PrintCanvasClose(canvas.get(), request.output_path.string());

    if (!std::filesystem::exists(request.output_path))
    {
        return FailPlot(
            PlotStatus::OUTPUT_ERROR,
            request.output_path,
            LogLevel::Error,
            "local_painter::SaveHeatmap() completed without producing '" +
                request.output_path.string() + "'.");
    }

    Logger::Log(LogLevel::Info, " Output file: " + request.output_path.string());
    return BuildPlotResult(PlotStatus::SUCCESS, request.output_path, "");
}

#endif

} // namespace
namespace local_painter {

PlotResult SaveLinePlot(const LinePlotRequest & request)
{
    const auto validation_result{ ValidateLinePlotRequest(request) };
    if (!validation_result.Succeeded())
    {
        return validation_result;
    }

#ifdef HAVE_ROOT
    try
    {
        return RenderLinePlotWithRoot(request);
    }
    catch (const std::exception & ex)
    {
        return FailPlot(
            PlotStatus::RENDER_ERROR,
            request.output_path,
            LogLevel::Error,
            "local_painter::SaveLinePlot() failed: " + std::string(ex.what()));
    }
#else
    return FailPlot(
        PlotStatus::ROOT_UNAVAILABLE,
        request.output_path,
        LogLevel::Warning,
        "local_painter::SaveLinePlot() requires ROOT support, but this build was compiled without HAVE_ROOT.");
#endif
}

PlotResult SaveHeatmap(const HeatmapRequest & request)
{
    const auto validation_result{ ValidateHeatmapRequest(request) };
    if (!validation_result.Succeeded())
    {
        return validation_result;
    }

#ifdef HAVE_ROOT
    try
    {
        return RenderHeatmapWithRoot(request);
    }
    catch (const std::exception & ex)
    {
        return FailPlot(
            PlotStatus::RENDER_ERROR,
            request.output_path,
            LogLevel::Error,
            "local_painter::SaveHeatmap() failed: " + std::string(ex.what()));
    }
#else
    return FailPlot(
        PlotStatus::ROOT_UNAVAILABLE,
        request.output_path,
        LogLevel::Warning,
        "local_painter::SaveHeatmap() requires ROOT support, but this build was compiled without HAVE_ROOT.");
#endif
}

} // namespace local_painter
} // namespace rhbm_gem
