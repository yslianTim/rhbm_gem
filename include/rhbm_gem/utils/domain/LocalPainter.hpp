#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem {

enum class PlotStatus
{
    SUCCESS,
    INVALID_INPUT,
    ROOT_UNAVAILABLE,
    OUTPUT_ERROR,
    RENDER_ERROR
};

struct AxisSpec
{
    std::string title{};
    std::optional<std::pair<double, double>> range{};
    bool auto_range{ true };
};

struct LineSeries
{
    std::string name{};
    std::vector<double> x_values;
    std::vector<double> y_values;
    std::optional<int> color_index{};
};

struct LinePlotPanel
{
    std::string label{};
    AxisSpec y_axis{};
    std::vector<LineSeries> series;
};

struct LinePlotRequest
{
    std::filesystem::path output_path{};
    AxisSpec x_axis{};
    std::vector<LinePlotPanel> panels;
    std::string title{};
    std::string shared_y_axis_title{};
    int canvas_width{ 1000 };
    int canvas_height_per_panel{ 375 };
    double x_auto_padding_ratio{ 0.05 };
    double y_auto_padding_ratio{ 0.10 };
};

struct HeatmapRequest
{
    std::filesystem::path output_path{};
    Eigen::MatrixXd values;
    AxisSpec x_axis{};
    AxisSpec y_axis{};
    AxisSpec z_axis{};
    std::string title{};
    int canvas_width{ 1200 };
    int canvas_height{ 1000 };
};

struct PlotResult
{
    PlotStatus status{ PlotStatus::SUCCESS };
    std::filesystem::path output_path{};
    std::string message{};

    bool Succeeded() const { return status == PlotStatus::SUCCESS; }
    explicit operator bool() const { return Succeeded(); }
};

// Quick plotting facade for lightweight local reports and debug visualizations.
namespace local_painter {

PlotResult SaveLinePlot(const LinePlotRequest & request);

PlotResult SaveHeatmap(const HeatmapRequest & request);

} // namespace local_painter

} // namespace rhbm_gem
