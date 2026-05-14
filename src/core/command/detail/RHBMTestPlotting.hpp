#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include <rhbm_gem/core/command/CommandSystem.hpp>
#include <rhbm_gem/utils/domain/LocalPainter.hpp>
#include <rhbm_gem/utils/hrl/RHBMTester.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>

namespace rhbm_gem::command_detail::rhbm_test_plotting {

enum class BiasCurveKind
{
    Ols,
    Median,
    Mdpde,
    TrainedMdpde,
    RequestedAlpha,
    TrainedAlpha
};

enum class BiasPlotFlavor
{
    DataOutlier,
    MemberOutlier,
    ModelAlphaData,
    ModelAlphaMember,
    NeighborType
};

enum class BiasXAxisMode
{
    ContaminationRatio,
    NeighborType
};

struct BiasCurvePoint
{
    double x{ 0.0 };
    rhbm_tester::BiasStatistics bias;
};

struct BiasCurve
{
    BiasCurveKind kind;
    std::vector<BiasCurvePoint> points;
};

struct BiasPlotPanel
{
    std::string label;
    std::vector<BiasCurve> curves;
};

struct BiasPlotRequest
{
    std::string output_name;
    BiasPlotFlavor flavor{ BiasPlotFlavor::DataOutlier };
    BiasXAxisMode x_axis_mode{ BiasXAxisMode::ContaminationRatio };
    std::vector<BiasPlotPanel> panels;
};

BiasCurve MakeBiasCurve(BiasCurveKind kind, size_t point_capacity);

void AppendBiasCurvePoint(
    BiasCurve & curve,
    double x,
    const rhbm_tester::BiasStatistics & bias);

std::string FormatDataBiasPanelLabel(size_t panel_index);

std::string FormatMemberBiasPanelLabel(size_t panel_index);

bool TryAppendBenchmarkLinearizedPanel(
    std::vector<LinePlotPanel> & panels,
    double distance,
    const RHBMMemberDataset & no_cut_dataset,
    const RHBMMemberDataset & cut_dataset);

void SaveBenchmarkLinearizedDatasetReport(
    const RHBMTestRequest & request,
    double error_sigma,
    const std::vector<LinePlotPanel> & panels);

void SaveDataOutlierBiasPlot(
    const RHBMTestRequest & request,
    const BiasPlotRequest & plot_request);

void SaveMemberOutlierBiasPlot(
    const RHBMTestRequest & request,
    const BiasPlotRequest & plot_request);

} // namespace rhbm_gem::command_detail::rhbm_test_plotting
