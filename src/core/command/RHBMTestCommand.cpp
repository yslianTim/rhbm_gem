#include "RHBMTestCommand.hpp"
#include <rhbm_gem/utils/domain/LocalPainter.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/TestDataFactory.hpp>
#include <rhbm_gem/utils/hrl/RHBMTester.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#ifdef HAVE_ROOT
#include <rhbm_gem/utils/domain/ROOTHelper.hpp>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#endif

namespace rhbm_gem {

namespace {
constexpr std::string_view kTesterOption{ "--tester" };
constexpr std::string_view kFitMinOption{ "--fit-min" };
constexpr std::string_view kFitMaxOption{ "--fit-max" };
constexpr std::string_view kFitRangeIssue{ "--fit-range" };
constexpr std::string_view kAlphaROption{ "--alpha-r" };
constexpr std::string_view kAlphaGOption{ "--alpha-g" };
constexpr int kGausParSize{ GaussianModel3D::kParameterSize };
struct RHBMTestExecutionContext
{
    const RHBMTestRequest & options;
    int thread_size;
    const std::filesystem::path & output_folder;
};

struct BetaScenarioConfig
{
    int replica_size;
    int sampling_entry_size;
    std::vector<double> alpha_r_list;
};

struct MuScenarioConfig
{
    int replica_size;
    int member_size;
    std::vector<double> alpha_g_list;
    std::vector<Eigen::VectorXd> outlier_prior_list;
};

struct NeighborDistanceScenarioConfig
{
    int replica_size;
    int sampling_entry_size;
    size_t neighbor_count;
    double rejected_angle;
    std::vector<double> error_list;
    std::vector<double> distance_list;
};

enum class ResidualCurveKind
{
    Ols,
    Median,
    Mdpde,
    TrainedMdpde,
    RequestedAlpha,
    TrainedAlpha
};

enum class ResidualPlotFlavor
{
    DataOutlier,
    MemberOutlier,
    ModelAlphaData,
    ModelAlphaMember,
    NeighborDistance
};

enum class ResidualXAxisMode
{
    ContaminationRatio,
    NeighborDistance
};

struct ResidualCurvePoint
{
    double x{ 0.0 };
    rhbm_tester::ResidualStatistics residual;
};

struct ResidualCurve
{
    ResidualCurveKind kind;
    std::vector<ResidualCurvePoint> points;
};

struct ResidualPlotPanel
{
    std::string label;
    std::vector<ResidualCurve> curves;
};

struct ResidualPlotRequest
{
    std::string output_name;
    ResidualPlotFlavor flavor{ ResidualPlotFlavor::DataOutlier };
    ResidualXAxisMode x_axis_mode{ ResidualXAxisMode::ContaminationRatio };
    std::vector<ResidualPlotPanel> panels;
};

Eigen::VectorXd MakeDefaultModelPrior()
{
    Eigen::VectorXd model_par_prior{ Eigen::VectorXd::Zero(kGausParSize) };
    model_par_prior(0) = 1.0;
    model_par_prior(1) = 0.5;
    model_par_prior(2) = 0.0;
    return model_par_prior;
}

Eigen::VectorXd MakeDefaultModelSigma()
{
    Eigen::VectorXd model_par_sigma{ Eigen::VectorXd::Zero(kGausParSize) };
    model_par_sigma(0) = 0.050;
    model_par_sigma(1) = 0.025;
    model_par_sigma(2) = 0.010;
    return model_par_sigma;
}

std::vector<double> BuildLinearSweep(int count, double step, double start = 0.0)
{
    std::vector<double> values(static_cast<size_t>(count));
    for (int i = 0; i < count; i++)
    {
        values[static_cast<size_t>(i)] = start + static_cast<double>(i) * step;
    }
    return values;
}

std::vector<double> BuildDescendingSweep(int count, double start, double step)
{
    std::vector<double> values(static_cast<size_t>(count));
    for (int i = 0; i < count; i++)
    {
        values[static_cast<size_t>(i)] = start - static_cast<double>(i) * step;
    }
    return values;
}

ResidualCurve MakeResidualCurve(ResidualCurveKind kind, size_t point_capacity)
{
    ResidualCurve curve{ kind, {} };
    curve.points.reserve(point_capacity);
    return curve;
}

void AppendResidualCurvePoint(
    ResidualCurve & curve,
    double x,
    const rhbm_tester::ResidualStatistics & residual)
{
    curve.points.emplace_back(ResidualCurvePoint{ x, residual });
}

std::string FormatDataResidualPanelLabel(size_t panel_index)
{
    const double error_value[3]{ 0.0, 2.5, 5.0 };
    const auto value{ error_value[panel_index] };
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(1)
           << "#sigma_{#epsilon} = " << value << "% D_{max}";
    return stream.str();
}

std::string FormatMemberResidualPanelLabel(size_t panel_index)
{
    const std::string outlier_type_list[2]{
        "#font[2]{A}", "#tau"
    };
    return "Outlier in " + outlier_type_list[panel_index];
}

test_data_factory::TestDataFactory BuildDataFactory(const RHBMTestExecutionContext & options)
{
    test_data_factory::TestDataFactory factory(
        linearization_service::LinearizationSpec::DefaultDataset());
    factory.SetFittingRange(options.options.fit_range_min, options.options.fit_range_max);
    return factory;
}

test_data_factory::TestDataFactory::NeighborhoodScenario BuildNeighborhoodScenario(
    const Eigen::VectorXd & model_par_prior,
    const NeighborDistanceScenarioConfig & scenario,
    double error_sigma,
    double neighbor_distance,
    bool include_sampling_summary = false)
{
    return test_data_factory::TestDataFactory::NeighborhoodScenario{
        model_par_prior,
        scenario.sampling_entry_size,
        error_sigma,
        0.0,
        1.0,
        neighbor_distance,
        scenario.neighbor_count,
        scenario.rejected_angle,
        include_sampling_summary,
        0.0,
        4.0,
        scenario.replica_size
    };
}

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
    const RHBMTestExecutionContext & options,
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

    LinePlotRequest request;
    request.output_path = options.output_folder /
        std::string("benchmark_cut_vs_no_cut_sigma_" + FormatSigmaToken(error_sigma) + ".pdf");
    request.title = "";
    request.x_axis.title = "Linearized Basis";
    request.shared_y_axis_title = "Linearized Response";
    request.panels = panels;
    request.canvas_width = 2000;
    request.canvas_height_per_panel = 200;

    const auto plot_result{ local_painter::SaveLinePlot(request) };
    if (!plot_result.Succeeded())
    {
        Logger::Log(
            LogLevel::Warning,
            "Failed to emit benchmark linearized dataset report '"
                + request.output_path.string() + "': " + plot_result.message);
    }
}

#ifdef HAVE_ROOT
double ScaleResidualPlotX(double x, ResidualXAxisMode x_axis_mode)
{
    if (x_axis_mode == ResidualXAxisMode::NeighborDistance)
    {
        return -1.0 * x;
    }
    return 100.0 * x;
}

bool ShouldDrawResidualCurve(ResidualPlotFlavor, ResidualCurveKind)
{
    return true;
}

short GetResidualCurveColor(ResidualPlotFlavor flavor, ResidualCurveKind kind)
{
    if (flavor == ResidualPlotFlavor::MemberOutlier)
    {
        if (kind == ResidualCurveKind::Median)
        {
            return kAzure;
        }
        if (kind == ResidualCurveKind::Mdpde)
        {
            return kRed;
        }
        return kGreen + 1;
    }
    if (flavor == ResidualPlotFlavor::ModelAlphaMember)
    {
        return (kind == ResidualCurveKind::RequestedAlpha) ? kAzure : kRed;
    }
    if (flavor == ResidualPlotFlavor::ModelAlphaData)
    {
        return (kind == ResidualCurveKind::RequestedAlpha) ? kAzure : kGreen + 2;
    }
    if (kind == ResidualCurveKind::Ols)
    {
        return kAzure;
    }
    if (kind == ResidualCurveKind::Mdpde)
    {
        return kGreen + 2;
    }
    return kRed;
}

short GetResidualCurveMarker(ResidualPlotFlavor flavor, ResidualCurveKind kind)
{
    if (flavor == ResidualPlotFlavor::MemberOutlier)
    {
        if (kind == ResidualCurveKind::Median)
        {
            return 24;
        }
        if (kind == ResidualCurveKind::Mdpde)
        {
            return 20;
        }
        return 25;
    }
    if (flavor == ResidualPlotFlavor::ModelAlphaMember)
    {
        return (kind == ResidualCurveKind::RequestedAlpha) ? 24 : 20;
    }
    if (kind == ResidualCurveKind::Ols || kind == ResidualCurveKind::RequestedAlpha)
    {
        return 24;
    }
    if (kind == ResidualCurveKind::Mdpde || kind == ResidualCurveKind::TrainedAlpha)
    {
        return 25;
    }
    return 20;
}

short GetResidualCurveLineStyle(ResidualPlotFlavor flavor, ResidualCurveKind kind)
{
    if (kind == ResidualCurveKind::Ols ||
        kind == ResidualCurveKind::Median ||
        kind == ResidualCurveKind::RequestedAlpha)
    {
        return 2;
    }
    if (kind == ResidualCurveKind::TrainedMdpde &&
        flavor != ResidualPlotFlavor::ModelAlphaMember)
    {
        return 3;
    }
    return 1;
}

void ApplyResidualCurveStyle(
    TGraphErrors * graph,
    ResidualPlotFlavor flavor,
    ResidualCurveKind kind)
{
    const auto color{ GetResidualCurveColor(flavor, kind) };
    root_helper::SetMarkerAttribute(
        graph,
        GetResidualCurveMarker(flavor, kind),
        1.5f,
        color);
    root_helper::SetLineAttribute(
        graph,
        GetResidualCurveLineStyle(flavor, kind),
        2,
        color);
    root_helper::SetFillAttribute(graph, 1001, color, 0.2f);
}

std::vector<ResidualCurveKind> GetResidualLegendOrder(ResidualPlotFlavor flavor)
{
    if (flavor == ResidualPlotFlavor::ModelAlphaData ||
        flavor == ResidualPlotFlavor::ModelAlphaMember)
    {
        return { ResidualCurveKind::RequestedAlpha, ResidualCurveKind::TrainedAlpha };
    }
    if (flavor == ResidualPlotFlavor::MemberOutlier)
    {
        return { ResidualCurveKind::TrainedMdpde, ResidualCurveKind::Mdpde, ResidualCurveKind::Median };
    }
    return { ResidualCurveKind::TrainedMdpde, ResidualCurveKind::Mdpde, ResidualCurveKind::Ols };
}

std::string GetResidualLegendLabel(ResidualPlotFlavor flavor, ResidualCurveKind kind)
{
    if (flavor == ResidualPlotFlavor::ModelAlphaData)
    {
        return (kind == ResidualCurveKind::RequestedAlpha) ?
            "MDPDE (#alpha_{r} = 0.1)" :
            "MDPDE (#alpha_{r} = 0.4)";
    }
    if (flavor == ResidualPlotFlavor::ModelAlphaMember)
    {
        return (kind == ResidualCurveKind::RequestedAlpha) ?
            "MDPDE (#alpha_{g} = 0.2)" :
            "MDPDE (#alpha_{g} = 0.5)";
    }
    if (flavor == ResidualPlotFlavor::MemberOutlier)
    {
        if (kind == ResidualCurveKind::TrainedMdpde)
        {
            return "MDPDE (#alpha_{g} from Alg.5)";
        }
        if (kind == ResidualCurveKind::Mdpde)
        {
            return "MDPDE (#alpha_{g} = 0.2)";
        }
        return "MDPDE (#alpha_{g} = 0)";
    }
    if (kind == ResidualCurveKind::TrainedMdpde)
    {
        return "MDPDE (w/ Sampling Scheme)";
    }
    if (kind == ResidualCurveKind::Mdpde)
    {
        return "MDPDE";
    }
    return "Ordinary Least Squares";
}

TGraphErrors * FindResidualGraph(
    std::vector<std::unique_ptr<TGraphErrors>> & graph_list,
    const std::vector<ResidualCurveKind> & kind_list,
    ResidualCurveKind kind)
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
#endif

} // namespace

void PrintDataOutlierResult(
    const RHBMTestExecutionContext & options,
    const ResidualPlotRequest & request
);

void PrintMemberOutlierResult(
    const RHBMTestExecutionContext & options,
    const ResidualPlotRequest & request
);

void PrintAtomSamplingDataSummary(
    const RHBMTestExecutionContext & options,
    const std::string & name,
    const std::vector<LocalPotentialSampleList> & sampling_entries_list,
    const std::vector<double> & distance_list
);

#ifdef HAVE_ROOT
std::unique_ptr<TGraphErrors> CreateDistanceToResponseGraph(
    const LocalPotentialSampleList & sampling_entries);

std::unique_ptr<TH2D> CreateDistanceToResponseHistogram(
    const LocalPotentialSampleList & sampling_entries,
    int x_bin_size, int y_bin_size = 500);
#endif

void RunSimulationTestOnBenchMark(const RHBMTestExecutionContext & options)
{
    ScopeTimer timer("RHBMTestCommand::RunSimulationTestOnBenchMark");

    const auto scenario{ NeighborDistanceScenarioConfig{
        1,
        50,
        3,
        45.0,
        std::vector<double>{ 0.0 },
        BuildDescendingSweep(16, 2.5, 0.1)
    } };
    const auto model_par_prior{ MakeDefaultModelPrior() };
    auto data_factory{ BuildDataFactory(options) };
    for (auto error_sigma : scenario.error_list)
    {
        std::vector<LinePlotPanel> linearized_panels;
        linearized_panels.reserve(scenario.distance_list.size());
        for (double distance : scenario.distance_list)
        {
            const auto test_input{
                data_factory.BuildNeighborhoodTestInput(
                    BuildNeighborhoodScenario(
                        model_par_prior,
                        scenario,
                        error_sigma,
                        distance))
            };
            rhbm_tester::BetaReplicaResidual residual_result;
            rhbm_tester::RunSingleBetaMDPDETest(
                residual_result,
                test_input.cut_datasets.front(),
                test_input.gaus_true,
                options.options.alpha_r,
                options.thread_size);
            TryAppendBenchmarkLinearizedPanel(
                linearized_panels,
                distance,
                test_input.no_cut_datasets.front(),
                test_input.cut_datasets.front());

            std::ostringstream stream;
            stream << "Distance: " << distance
                   << " , OLS: "
                   << residual_result.ols_residual(0) << " , "
                   << residual_result.ols_residual(1) << " , "
                   << residual_result.ols_residual(2)
                   << " , MDPDE: "
                   << residual_result.mdpde_residual(0) << " , "
                   << residual_result.mdpde_residual(1) << " , "
                   << residual_result.mdpde_residual(2)
                   << " (Alpha-R = " << options.options.alpha_r << ")";
            Logger::Log(LogLevel::Info, stream.str());
        }
        SaveBenchmarkLinearizedDatasetReport(options, error_sigma, linearized_panels);
    }
}

void RunSimulationTestOnDataOutlier(const RHBMTestExecutionContext & options)
{
    ScopeTimer timer("RHBMTestCommand::RunSimulationTestOnDataOutlier");

    const auto scenario{ BetaScenarioConfig{
        100,
        1000,
        std::vector<double>{ options.options.alpha_r }
    } };
    const auto model_par_prior{ MakeDefaultModelPrior() };
    auto data_factory{ BuildDataFactory(options) };
    std::vector<double> error_list{ 0.1, 0.2, 0.3 };
    const auto outlier_list{ BuildLinearSweep(9, 0.025) };
    ResidualPlotRequest plot_request;
    plot_request.output_name = "bias_outlier_in_data.pdf";
    plot_request.flavor = ResidualPlotFlavor::DataOutlier;
    plot_request.x_axis_mode = ResidualXAxisMode::ContaminationRatio;
    plot_request.panels.reserve(error_list.size());

    for (size_t panel_index = 0; panel_index < error_list.size(); panel_index++)
    {
        const auto error_sigma{ error_list.at(panel_index) };
        ResidualPlotPanel panel;
        panel.label = FormatDataResidualPanelLabel(panel_index);
        panel.curves.emplace_back(MakeResidualCurve(ResidualCurveKind::Ols, outlier_list.size()));
        panel.curves.emplace_back(MakeResidualCurve(ResidualCurveKind::Mdpde, outlier_list.size()));
        panel.curves.emplace_back(MakeResidualCurve(ResidualCurveKind::TrainedMdpde, outlier_list.size()));
        for (size_t i = 0; i < outlier_list.size(); i++)
        {
            rhbm_tester::BetaMDPDETestResidual residual;
            const auto test_input{
                data_factory.BuildBetaTestInput(test_data_factory::TestDataFactory::BetaScenario{
                    model_par_prior,
                    scenario.sampling_entry_size,
                    error_sigma,
                    outlier_list.at(i),
                    scenario.replica_size
                })
            };
            rhbm_tester::RunBetaMDPDETest(
                residual,
                scenario.alpha_r_list,
                test_input,
                options.thread_size
            );

            AppendResidualCurvePoint(panel.curves.at(0), outlier_list.at(i), residual.ols);
            AppendResidualCurvePoint(
                panel.curves.at(1),
                outlier_list.at(i),
                residual.mdpde.requested_alpha.front());
            AppendResidualCurvePoint(
                panel.curves.at(2),
                outlier_list.at(i),
                residual.mdpde.trained_alpha);
        }
        plot_request.panels.emplace_back(std::move(panel));
    }

    PrintDataOutlierResult(options, plot_request);
}

void RunSimulationTestOnMemberOutlier(const RHBMTestExecutionContext & options)
{
    ScopeTimer timer("RHBMTestCommand::RunSimulationTestOnMemberOutlier");

    const auto scenario{ MuScenarioConfig{
        100,
        100,
        std::vector<double>{ options.options.alpha_g },
        std::vector<Eigen::VectorXd>{
            Eigen::VectorXd::Zero(kGausParSize),
            Eigen::VectorXd::Zero(kGausParSize)
        }
    } };
    auto outlier_prior_list{ scenario.outlier_prior_list };
    outlier_prior_list.at(0)(0) = 1.50;
    outlier_prior_list.at(0)(1) = 0.50;
    outlier_prior_list.at(0)(2) = 0.10;
    outlier_prior_list.at(1)(0) = 1.00;
    outlier_prior_list.at(1)(1) = 1.00;
    outlier_prior_list.at(1)(2) = 0.10;

    const auto model_par_prior{ MakeDefaultModelPrior() };
    const auto model_par_sigma{ MakeDefaultModelSigma() };
    auto data_factory{ BuildDataFactory(options) };
    const auto outlier_list{ BuildLinearSweep(9, 0.025) };
    ResidualPlotRequest plot_request;
    plot_request.output_name = "bias_outlier_in_member.pdf";
    plot_request.flavor = ResidualPlotFlavor::MemberOutlier;
    plot_request.x_axis_mode = ResidualXAxisMode::ContaminationRatio;
    plot_request.panels.reserve(outlier_prior_list.size());

    for (size_t panel_index = 0; panel_index < outlier_prior_list.size(); panel_index++)
    {
        const auto & outlier_prior{ outlier_prior_list.at(panel_index) };
        ResidualPlotPanel panel;
        panel.label = FormatMemberResidualPanelLabel(panel_index);
        panel.curves.emplace_back(MakeResidualCurve(ResidualCurveKind::Median, outlier_list.size()));
        panel.curves.emplace_back(MakeResidualCurve(ResidualCurveKind::Mdpde, outlier_list.size()));
        panel.curves.emplace_back(MakeResidualCurve(ResidualCurveKind::TrainedMdpde, outlier_list.size()));
        for (size_t i = 0; i < outlier_list.size(); i++)
        {
            rhbm_tester::MuMDPDETestResidual residual;
            const auto test_input{
                data_factory.BuildMuTestInput(test_data_factory::TestDataFactory::MuScenario{
                    scenario.member_size,
                    model_par_prior,
                    model_par_sigma,
                    outlier_prior,
                    model_par_sigma,
                    outlier_list.at(i),
                    scenario.replica_size
                })
            };
            rhbm_tester::RunMuMDPDETest(
                residual,
                scenario.alpha_g_list,
                test_input,
                options.thread_size
            );

            AppendResidualCurvePoint(panel.curves.at(0), outlier_list.at(i), residual.median);
            AppendResidualCurvePoint(
                panel.curves.at(1),
                outlier_list.at(i),
                residual.mdpde.requested_alpha.front());
            AppendResidualCurvePoint(
                panel.curves.at(2),
                outlier_list.at(i),
                residual.mdpde.trained_alpha);
        }
        plot_request.panels.emplace_back(std::move(panel));
    }

    PrintMemberOutlierResult(options, plot_request);
}

void RunSimulationTestOnModelAlphaData(const RHBMTestExecutionContext & options)
{
    ScopeTimer timer("RHBMTestCommand::RunSimulationTestOnModelAlphaData");

    const auto scenario{ BetaScenarioConfig{
        100,
        1000,
        std::vector<double>{ options.options.alpha_r }
    } };
    const auto model_par_prior{ MakeDefaultModelPrior() };
    auto data_factory{ BuildDataFactory(options) };
    std::vector<double> error_list{ 0.1, 0.2, 0.3 };
    const auto outlier_list{ BuildLinearSweep(10, 0.05) };
    ResidualPlotRequest plot_request;
    plot_request.output_name = "bias_outlier_with_alpha_in_data.pdf";
    plot_request.flavor = ResidualPlotFlavor::ModelAlphaData;
    plot_request.x_axis_mode = ResidualXAxisMode::ContaminationRatio;
    plot_request.panels.reserve(error_list.size());

    for (size_t panel_index = 0; panel_index < error_list.size(); panel_index++)
    {
        const auto error_sigma{ error_list.at(panel_index) };
        ResidualPlotPanel panel;
        panel.label = FormatDataResidualPanelLabel(panel_index);
        panel.curves.emplace_back(MakeResidualCurve(ResidualCurveKind::RequestedAlpha, outlier_list.size()));
        panel.curves.emplace_back(MakeResidualCurve(ResidualCurveKind::TrainedAlpha, outlier_list.size()));
        for (size_t i = 0; i < outlier_list.size(); i++)
        {
            rhbm_tester::BetaMDPDETestResidual residual;
            const auto test_input{
                data_factory.BuildBetaTestInput(test_data_factory::TestDataFactory::BetaScenario{
                    model_par_prior,
                    scenario.sampling_entry_size,
                    error_sigma,
                    outlier_list.at(i),
                    scenario.replica_size
                })
            };
            rhbm_tester::RunBetaMDPDETest(
                residual,
                scenario.alpha_r_list,
                test_input,
                options.thread_size
            );

            AppendResidualCurvePoint(
                panel.curves.at(0),
                outlier_list.at(i),
                residual.mdpde.requested_alpha.front());
            AppendResidualCurvePoint(
                panel.curves.at(1),
                outlier_list.at(i),
                residual.mdpde.trained_alpha);
        }
        plot_request.panels.emplace_back(std::move(panel));
    }

    PrintDataOutlierResult(options, plot_request);
}

void RunSimulationTestOnModelAlphaMember(const RHBMTestExecutionContext & options)
{
    ScopeTimer timer("RHBMTestCommand::RunSimulationTestOnModelAlphaMember");

    const auto scenario{ MuScenarioConfig{
        100,
        100,
        std::vector<double>{ options.options.alpha_g },
        std::vector<Eigen::VectorXd>{
            Eigen::VectorXd::Zero(kGausParSize),
            Eigen::VectorXd::Zero(kGausParSize)
        }
    } };
    auto outlier_prior_list{ scenario.outlier_prior_list };
    outlier_prior_list.at(0)(0) = 1.50;
    outlier_prior_list.at(0)(1) = 0.50;
    outlier_prior_list.at(0)(2) = 0.10;
    outlier_prior_list.at(1)(0) = 1.00;
    outlier_prior_list.at(1)(1) = 1.00;
    outlier_prior_list.at(1)(2) = 0.10;

    const auto model_par_prior{ MakeDefaultModelPrior() };
    const auto model_par_sigma{ MakeDefaultModelSigma() };
    auto data_factory{ BuildDataFactory(options) };
    const auto outlier_list{ BuildLinearSweep(10, 0.05) };
    ResidualPlotRequest plot_request;
    plot_request.output_name = "bias_outlier_with_alpha_in_member.pdf";
    plot_request.flavor = ResidualPlotFlavor::ModelAlphaMember;
    plot_request.x_axis_mode = ResidualXAxisMode::ContaminationRatio;
    plot_request.panels.reserve(outlier_prior_list.size());

    for (size_t panel_index = 0; panel_index < outlier_prior_list.size(); panel_index++)
    {
        const auto & outlier_prior{ outlier_prior_list.at(panel_index) };
        ResidualPlotPanel panel;
        panel.label = FormatMemberResidualPanelLabel(panel_index);
        panel.curves.emplace_back(MakeResidualCurve(ResidualCurveKind::RequestedAlpha, outlier_list.size()));
        panel.curves.emplace_back(MakeResidualCurve(ResidualCurveKind::TrainedAlpha, outlier_list.size()));
        for (size_t i = 0; i < outlier_list.size(); i++)
        {
            rhbm_tester::MuMDPDETestResidual residual;
            const auto test_input{
                data_factory.BuildMuTestInput(test_data_factory::TestDataFactory::MuScenario{
                    scenario.member_size,
                    model_par_prior,
                    model_par_sigma,
                    outlier_prior,
                    model_par_sigma,
                    outlier_list.at(i),
                    scenario.replica_size
                })
            };
            rhbm_tester::RunMuMDPDETest(
                residual,
                scenario.alpha_g_list,
                test_input,
                options.thread_size
            );

            AppendResidualCurvePoint(
                panel.curves.at(0),
                outlier_list.at(i),
                residual.mdpde.requested_alpha.front());
            AppendResidualCurvePoint(
                panel.curves.at(1),
                outlier_list.at(i),
                residual.mdpde.trained_alpha);
        }
        plot_request.panels.emplace_back(std::move(panel));
    }

    PrintMemberOutlierResult(options, plot_request);
}

void RunSimulationTestOnNeighborDistance(const RHBMTestExecutionContext & options)
{
    ScopeTimer timer("RHBMTestCommand::RunSimulationTestOnNeighborDistance");

    const auto scenario{ NeighborDistanceScenarioConfig{
        10,
        50,
        2,
        45.0,
        std::vector<double>{ 0.0, 0.025, 0.05 },
        BuildDescendingSweep(16, 2.5, 0.1)
    } };
    const auto model_par_prior{ MakeDefaultModelPrior() };
    auto data_factory{ BuildDataFactory(options) };
    ResidualPlotRequest plot_request;
    plot_request.output_name = "bias_from_neighbor_atom.pdf";
    plot_request.flavor = ResidualPlotFlavor::NeighborDistance;
    plot_request.x_axis_mode = ResidualXAxisMode::NeighborDistance;
    plot_request.panels.reserve(scenario.error_list.size());

    bool is_print_sampling_summary{ false };
    for (size_t panel_index = 0; panel_index < scenario.error_list.size(); panel_index++)
    {
        const auto error_sigma{ scenario.error_list.at(panel_index) };
        ResidualPlotPanel panel;
        panel.label = FormatDataResidualPanelLabel(panel_index);
        panel.curves.emplace_back(MakeResidualCurve(ResidualCurveKind::Ols, scenario.distance_list.size()));
        panel.curves.emplace_back(MakeResidualCurve(ResidualCurveKind::Mdpde, scenario.distance_list.size()));
        panel.curves.emplace_back(MakeResidualCurve(ResidualCurveKind::TrainedMdpde, scenario.distance_list.size()));
        std::vector<LocalPotentialSampleList> sampling_entries_list;
        sampling_entries_list.reserve(scenario.distance_list.size());
        for (size_t i = 0; i < scenario.distance_list.size(); i++)
        {
            rhbm_tester::NeighborhoodMDPDETestResidual residual;
            const auto test_input{
                data_factory.BuildNeighborhoodTestInput(
                    BuildNeighborhoodScenario(
                        model_par_prior,
                        scenario,
                        error_sigma,
                        scenario.distance_list.at(i),
                        true))
            };
            rhbm_tester::RunBetaMDPDEWithNeighborhoodTest(
                residual,
                test_input,
                options.thread_size,
                scenario.rejected_angle
            );

            sampling_entries_list.emplace_back(test_input.sampling_summaries.front());

            AppendResidualCurvePoint(
                panel.curves.at(0),
                scenario.distance_list.at(i),
                residual.no_cut_ols);
            AppendResidualCurvePoint(
                panel.curves.at(1),
                scenario.distance_list.at(i),
                residual.no_cut_mdpde);
            AppendResidualCurvePoint(
                panel.curves.at(2),
                scenario.distance_list.at(i),
                residual.cut_mdpde);

            Logger::Log(LogLevel::Info,
                std::string("Distance: ") + std::to_string(scenario.distance_list.at(i))
                + " , OLS: " + std::to_string(residual.no_cut_ols.mean(0)) + " +- " + std::to_string(residual.no_cut_ols.sigma(0))
                + " , MDPDE: " + std::to_string(residual.no_cut_mdpde.mean(0)) + " +- " + std::to_string(residual.no_cut_mdpde.sigma(0))
                + " , Train: " + std::to_string(residual.cut_mdpde.mean(0)) + " +- " + std::to_string(residual.cut_mdpde.sigma(0))
                + " (Alpha-R = " + std::to_string(residual.trained_alpha_r_average) + ")"
            );
        }
        plot_request.panels.emplace_back(std::move(panel));

        if (!is_print_sampling_summary)
        {
            PrintAtomSamplingDataSummary(
                options,
                "neighbor_atom_sampling_entries.pdf",
                sampling_entries_list, scenario.distance_list
            );
            is_print_sampling_summary = true;
        }
    }

    PrintDataOutlierResult(options, plot_request);
}

void PrintDataOutlierResult(
    const RHBMTestExecutionContext & options,
    const ResidualPlotRequest & request)
{
    auto file_path{ options.output_folder / request.output_name };
    Logger::Log(LogLevel::Info, " RHBMTestCommand::PrintDataOutlierResult");

    std::vector<std::string> title_y_list{
        "Amplitude #font[2]{A}", "Width #tau"
    };

    #ifdef HAVE_ROOT
    gStyle->SetLineScalePS(1.5);
    gStyle->SetGridColor(kGray);
    const int col_size{ 3 };
    const int row_size{ 2 };
    auto canvas{ root_helper::CreateCanvas("test","", 1500, 750) };
    root_helper::SetCanvasDefaultStyle(canvas.get());
    root_helper::SetCanvasPartition(
        canvas.get(), col_size, row_size, 0.17f, 0.08f, 0.12f, 0.14f, 0.01f, 0.01f
    );
    root_helper::PrintCanvasOpen(canvas.get(), file_path);

    std::vector<std::unique_ptr<TGraphErrors>> graph_list[col_size][row_size];
    std::vector<ResidualCurveKind> graph_kind_list[col_size][row_size];
    std::vector<double> global_y_array[row_size];
    for (size_t i = 0; i < col_size; i++)
    {
        const auto & panel{ request.panels.at(i) };
        for (size_t j = 0; j < row_size; j++)
        {
            for (const auto & curve : panel.curves)
            {
                auto graph{ root_helper::CreateGraphErrors() };
                for (size_t p = 0; p < curve.points.size(); p++)
                {
                    const auto & point{ curve.points.at(p) };
                    const auto x_value{ ScaleResidualPlotX(point.x, request.x_axis_mode) };
                    const auto mean{ point.residual.mean(static_cast<int>(j)) };
                    const auto sigma{ point.residual.sigma(static_cast<int>(j)) };
                    graph->SetPoint(static_cast<int>(p), x_value, mean);
                    graph->SetPointError(static_cast<int>(p), 0.0, sigma);
                    global_y_array[j].emplace_back(mean);
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
        x_min[i] = (request.flavor == ResidualPlotFlavor::ModelAlphaData) ? -2.0 : -0.7;
        x_max[i] = (request.flavor == ResidualPlotFlavor::ModelAlphaData) ? 47.0 : 22.0;
        if (request.x_axis_mode == ResidualXAxisMode::NeighborDistance)
        {
            x_min[i] = -2.6;
            x_max[i] = -0.8;
        }
    }
    for (size_t j = 0; j < row_size; j++)
    {
        auto y_range{ array_helper::ComputeScalingPercentileRangeTuple(global_y_array[j], 0.38, 0.005, 0.995) };
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> resolution_text[col_size];
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
                ApplyResidualCurveStyle(graph.get(), request.flavor, kind);
                if (ShouldDrawResidualCurve(request.flavor, kind))
                {
                    graph->Draw("PL3");
                }
            }

            if (i == 0)
            {
                title_y_text[j] = root_helper::CreatePaveText(-0.68, 0.30, 0.00, 0.70, "nbNDC ARC", true);
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
                title_x_text[i]->AddText(request.panels.at(static_cast<size_t>(i)).label.data());
                title_x_text[i]->Draw();
            }
        }
    }

    canvas->cd();
    auto pad_extra0{ root_helper::CreatePad("pad_extra0","", 0.02, 0.92, 0.98, 1.00) };
    //auto pad_extra0{ root_helper::CreatePad("pad_extra0","", 0.20, 0.92, 0.98, 1.00) };
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
    for (const auto kind : GetResidualLegendOrder(request.flavor))
    {
        auto * graph{ FindResidualGraph(graph_list[0][0], graph_kind_list[0][0], kind) };
        if (graph != nullptr)
        {
            legend->AddEntry(graph, GetResidualLegendLabel(request.flavor, kind).data(), "plf");
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
    if (request.x_axis_mode == ResidualXAxisMode::NeighborDistance)
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
    #endif
}

void PrintMemberOutlierResult(
    const RHBMTestExecutionContext & options,
    const ResidualPlotRequest & request)
{
    auto file_path{ options.output_folder / request.output_name };
    Logger::Log(LogLevel::Info, " RHBMTestCommand::PrintMemberOutlierResult");

    std::vector<std::string> title_y_list{
        "Amplitude #font[2]{A}", "Width #tau"
    };

    #ifdef HAVE_ROOT
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
    std::vector<ResidualCurveKind> graph_kind_list[col_size][row_size];
    std::vector<double> global_y_array;
    for (size_t i = 0; i < col_size; i++)
    {
        const auto & panel{ request.panels.at(i) };
        for (size_t j = 0; j < row_size; j++)
        {
            for (const auto & curve : panel.curves)
            {
                auto graph{ root_helper::CreateGraphErrors() };
                for (size_t p = 0; p < curve.points.size(); p++)
                {
                    const auto & point{ curve.points.at(p) };
                    const auto x_value{ ScaleResidualPlotX(point.x, request.x_axis_mode) };
                    const auto mean{ point.residual.mean(static_cast<int>(j)) };
                    const auto sigma{ point.residual.sigma(static_cast<int>(j)) };
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
        x_min[i] = (request.flavor == ResidualPlotFlavor::ModelAlphaMember) ? -2.0 : -0.7;
        x_max[i] = (request.flavor == ResidualPlotFlavor::ModelAlphaMember) ? 47.0 : 22.0;
    }
    auto y_range{ array_helper::ComputeScalingRangeTuple(global_y_array, 0.3) };
    for (size_t j = 0; j < row_size; j++)
    {
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> resolution_text[col_size];
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
                ApplyResidualCurveStyle(graph.get(), request.flavor, kind);
                if (ShouldDrawResidualCurve(request.flavor, kind))
                {
                    graph->Draw("PL3");
                }
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
                title_x_text[i]->AddText(request.panels.at(static_cast<size_t>(i)).label.data());
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
    for (const auto kind : GetResidualLegendOrder(request.flavor))
    {
        auto * graph{ FindResidualGraph(graph_list[0][0], graph_kind_list[0][0], kind) };
        if (graph != nullptr)
        {
            legend->AddEntry(graph, GetResidualLegendLabel(request.flavor, kind).data(), "plf");
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
    #endif
}

void PrintAtomSamplingDataSummary(
    const RHBMTestExecutionContext & options,
    const std::string & name,
    const std::vector<LocalPotentialSampleList> & sampling_entries_list,
    const std::vector<double> & distance_list
)
{
    auto file_path{ options.output_folder / name };
    Logger::Log(LogLevel::Info, " RHBMTestCommand::PrintAtomSamplingDataSummary");

    #ifdef HAVE_ROOT

    gStyle->SetLineScalePS(2.0);
    gStyle->SetGridColor(kGray);

    auto canvas{ root_helper::CreateCanvas("test","", 1000, 750) };
    root_helper::SetCanvasDefaultStyle(canvas.get());
    root_helper::PrintCanvasOpen(canvas.get(), file_path);
    const int pad_size{ 1 };

    std::unique_ptr<TPad> pad[pad_size];
    pad[0] = root_helper::CreatePad("pad0","", 0.00, 0.00, 1.00, 1.00); // The left pad
    
    size_t count{ 0 };
    for (auto sampling_entries : sampling_entries_list)
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
    #endif
}

#ifdef HAVE_ROOT
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
    const LocalPotentialSampleList & sampling_entries, int x_bin_size, int y_bin_size)
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

bool RunRHBMTestWorkflow(const RHBMTestExecutionContext & options)
{
    switch (options.options.tester_choice)
    {
    case TesterType::BENCHMARK:
        RunSimulationTestOnBenchMark(options);
        return true;
    case TesterType::DATA_OUTLIER:
        RunSimulationTestOnDataOutlier(options);
        return true;
    case TesterType::MEMBER_OUTLIER:
        RunSimulationTestOnMemberOutlier(options);
        return true;
    case TesterType::MODEL_ALPHA_DATA:
        RunSimulationTestOnModelAlphaData(options);
        return true;
    case TesterType::MODEL_ALPHA_MEMBER:
        RunSimulationTestOnModelAlphaMember(options);
        return true;
    case TesterType::NEIGHBOR_DISTANCE:
        RunSimulationTestOnNeighborDistance(options);
        return true;
    default:
        Logger::Log(
            LogLevel::Error,
            "Invalid tester choice reached execution path: ["
                + std::to_string(static_cast<int>(options.options.tester_choice)) + "]");
        return false;
    }
}


} // namespace rhbm_gem::detail

namespace rhbm_gem {

RHBMTestCommand::RHBMTestCommand() :
    CommandWithRequest<RHBMTestRequest>{}
{
}

void RHBMTestCommand::NormalizeRequest()
{
    auto & request{ MutableRequest() };
    CoerceEnum(
        request.tester_choice,
        kTesterOption,
        TesterType::BENCHMARK,
        "Tester choice");
    CoerceFiniteNonNegativeScalar(
        request.fit_range_min,
        kFitMinOption,
        0.0,
        LogLevel::Error,
        "Minimum fitting range");
    CoerceFiniteNonNegativeScalar(
        request.fit_range_max,
        kFitMaxOption,
        1.0,
        LogLevel::Error,
        "Maximum fitting range");
    CoerceFinitePositiveScalar(
        request.alpha_r,
        kAlphaROption,
        0.1,
        LogLevel::Error,
        "Alpha-R");
    CoerceFinitePositiveScalar(
        request.alpha_g,
        kAlphaGOption,
        0.2,
        LogLevel::Error,
        "Alpha-G");
}

void RHBMTestCommand::ValidateOptions()
{
    const auto & request{ RequestOptions() };
    RequireCondition(
        request.fit_range_min <= request.fit_range_max,
        kFitRangeIssue,
        "Expected --fit-min <= --fit-max.");
}

bool RHBMTestCommand::ExecuteImpl()
{
    return RunRHBMTestWorkflow(RHBMTestExecutionContext{
        RequestOptions(),
        ThreadSize(),
        OutputFolder(),
    });
}

} // namespace rhbm_gem
