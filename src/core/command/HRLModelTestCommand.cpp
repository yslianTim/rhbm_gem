#include "HRLModelTestCommand.hpp"
#include <rhbm_gem/utils/domain/LocalPainter.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/hrl/GaussianLinearizationService.hpp>
#include <rhbm_gem/utils/hrl/HRLModelTestDataFactory.hpp>
#include <rhbm_gem/utils/hrl/HRLModelTester.hpp>
#include <rhbm_gem/utils/math/ArrayStats.hpp>

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
} // namespace

HRLModelTestCommand::HRLModelTestCommand() :
    CommandWithRequest<HRLModelTestRequest>{}
{
}

void HRLModelTestCommand::NormalizeRequest()
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

void HRLModelTestCommand::ValidateOptions()
{
    const auto & request{ RequestOptions() };
    RequireCondition(
        request.fit_range_min <= request.fit_range_max,
        kFitRangeIssue,
        "Expected --fit-min <= --fit-max.");
}

} // namespace rhbm_gem

namespace rhbm_gem::detail {

struct HRLModelTestExecutionContext
{
    const HRLModelTestRequest & options;
    int thread_size;
    const std::filesystem::path & output_folder;
};

namespace {

constexpr int kGausParSize{ 3 };
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

std::filesystem::path BuildOutputPath(
    const HRLModelTestExecutionContext & options,
    std::string_view stem)
{
    return options.output_folder / std::string(stem);
}

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

HRLModelTestDataFactory BuildDataFactory(const HRLModelTestExecutionContext & options)
{
    HRLModelTestDataFactory factory(
        kGausParSize,
        rhbm_gem::GaussianLinearizationSpec::DefaultDataset());
    factory.SetFittingRange(options.options.fit_range_min, options.options.fit_range_max);
    return factory;
}

HRLModelTestDataFactory::NeighborhoodScenario BuildNeighborhoodScenario(
    const Eigen::VectorXd & model_par_prior,
    const NeighborDistanceScenarioConfig & scenario,
    double error_sigma,
    double neighbor_distance,
    bool include_sampling_summary = false)
{
    return HRLModelTestDataFactory::NeighborhoodScenario{
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
    const HRLMemberDataset & dataset,
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
    const HRLMemberDataset & dataset,
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
    const HRLMemberDataset & no_cut_dataset,
    const HRLMemberDataset & cut_dataset)
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
    const HRLModelTestExecutionContext & options,
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
    request.output_path = BuildOutputPath(
        options,
        "benchmark_linearized_cut_vs_no_cut_sigma_" + FormatSigmaToken(error_sigma) + ".pdf");
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

} // namespace

void PrintDataOutlierResult(
    const HRLModelTestExecutionContext & options,
    const std::string & name,
    const std::vector<double> & outlier_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_ols_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_mdpde_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_train_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_ols_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_mdpde_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_train_list,
    bool is_neighbor_distance = false
);

void PrintMemberOutlierResult(
    const HRLModelTestExecutionContext & options,
    const std::string & name,
    const std::vector<double> & outlier_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_median_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_mdpde_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_train_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_median_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_mdpde_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_train_list
);

void PrintAtomSamplingDataSummary(
    const HRLModelTestExecutionContext & options,
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

void RunSimulationTestOnBenchMark(const HRLModelTestExecutionContext & options)
{
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnBenchMark");

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
    HRLModelTester tester(kGausParSize);

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
            HRLModelTester::BetaReplicaResidual residual_result;
            tester.RunSingleBetaMDPDETest(
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

void RunSimulationTestOnDataOutlier(const HRLModelTestExecutionContext & options)
{
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnDataOutlier");

    const auto scenario{ BetaScenarioConfig{
        100,
        1000,
        std::vector<double>{ options.options.alpha_r }
    } };
    const auto model_par_prior{ MakeDefaultModelPrior() };
    auto data_factory{ BuildDataFactory(options) };
    HRLModelTester tester(kGausParSize);
    std::vector<double> error_list{ 0.1, 0.2, 0.3 };
    std::vector<Eigen::MatrixXd> mean_matrix_ols_list;
    std::vector<Eigen::MatrixXd> mean_matrix_mdpde_list;
    std::vector<Eigen::MatrixXd> mean_matrix_train_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_ols_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_mdpde_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_train_list;

    const auto outlier_list{ BuildLinearSweep(9, 0.025) };
    const auto outlier_size{ static_cast<int>(outlier_list.size()) };

    for (auto error_sigma : error_list)
    {
        Eigen::MatrixXd mean_matrix_ols{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd mean_matrix_mdpde{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd mean_matrix_train{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd sigma_matrix_ols{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd sigma_matrix_mdpde{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd sigma_matrix_train{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        for (int i = 0; i < outlier_size; i++)
        {
            std::vector<Eigen::VectorXd> residual_mean_ols_list;
            std::vector<Eigen::VectorXd> residual_mean_mdpde_list;
            std::vector<Eigen::VectorXd> residual_sigma_ols_list;
            std::vector<Eigen::VectorXd> residual_sigma_mdpde_list;
            const auto test_input{
                data_factory.BuildBetaTestInput(HRLModelTestDataFactory::BetaScenario{
                    model_par_prior,
                    scenario.sampling_entry_size,
                    error_sigma,
                    outlier_list[static_cast<size_t>(i)],
                    scenario.replica_size
                })
            };
            tester.RunBetaMDPDETest(
                scenario.alpha_r_list,
                residual_mean_ols_list, residual_mean_mdpde_list,
                residual_sigma_ols_list, residual_sigma_mdpde_list,
                test_input,
                options.thread_size
            );

            mean_matrix_ols.col(i) = residual_mean_ols_list.front();
            sigma_matrix_ols.col(i) = residual_sigma_ols_list.front();
            mean_matrix_mdpde.col(i) = residual_mean_mdpde_list.front();
            sigma_matrix_mdpde.col(i) = residual_sigma_mdpde_list.front();
            mean_matrix_train.col(i) = residual_mean_mdpde_list.back();
            sigma_matrix_train.col(i) = residual_sigma_mdpde_list.back();
        }
        mean_matrix_ols_list.emplace_back(mean_matrix_ols);
        mean_matrix_mdpde_list.emplace_back(mean_matrix_mdpde);
        mean_matrix_train_list.emplace_back(mean_matrix_train);
        sigma_matrix_ols_list.emplace_back(sigma_matrix_ols);
        sigma_matrix_mdpde_list.emplace_back(sigma_matrix_mdpde);
        sigma_matrix_train_list.emplace_back(sigma_matrix_train);
    }

    PrintDataOutlierResult(
        options,
        "bias_outlier_in_data.pdf",
        outlier_list,
        mean_matrix_ols_list, mean_matrix_mdpde_list, mean_matrix_train_list,
        sigma_matrix_ols_list, sigma_matrix_mdpde_list, sigma_matrix_train_list
    );
}

void RunSimulationTestOnMemberOutlier(const HRLModelTestExecutionContext & options)
{
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnMemberOutlier");

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
    HRLModelTester tester(kGausParSize);

    std::vector<Eigen::MatrixXd> mean_matrix_median_list;
    std::vector<Eigen::MatrixXd> mean_matrix_mdpde_list;
    std::vector<Eigen::MatrixXd> mean_matrix_train_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_median_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_mdpde_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_train_list;

    const auto outlier_list{ BuildLinearSweep(9, 0.025) };
    const auto outlier_size{ static_cast<int>(outlier_list.size()) };

    for (auto outlier_prior : outlier_prior_list)
    {
        Eigen::MatrixXd mean_matrix_median{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd mean_matrix_mdpde{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd mean_matrix_train{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd sigma_matrix_median{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd sigma_matrix_mdpde{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd sigma_matrix_train{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        for (int i = 0; i < outlier_size; i++)
        {
            std::vector<Eigen::VectorXd> residual_mean_median_list;
            std::vector<Eigen::VectorXd> residual_mean_mdpde_list;
            std::vector<Eigen::VectorXd> residual_sigma_median_list;
            std::vector<Eigen::VectorXd> residual_sigma_mdpde_list;
            const auto test_input{
                data_factory.BuildMuTestInput(HRLModelTestDataFactory::MuScenario{
                    scenario.member_size,
                    model_par_prior,
                    model_par_sigma,
                    outlier_prior,
                    model_par_sigma,
                    outlier_list[static_cast<size_t>(i)],
                    scenario.replica_size
                })
            };
            tester.RunMuMDPDETest(
                scenario.alpha_g_list,
                residual_mean_median_list, residual_mean_mdpde_list,
                residual_sigma_median_list, residual_sigma_mdpde_list,
                test_input,
                options.thread_size
            );

            mean_matrix_median.col(i) = residual_mean_median_list.front();
            mean_matrix_mdpde.col(i) = residual_mean_mdpde_list.front();
            mean_matrix_train.col(i) = residual_mean_mdpde_list.back();
            sigma_matrix_median.col(i) = residual_sigma_median_list.front();
            sigma_matrix_mdpde.col(i) = residual_sigma_mdpde_list.front();
            sigma_matrix_train.col(i) = residual_sigma_mdpde_list.back();
        }
        mean_matrix_median_list.emplace_back(mean_matrix_median);
        mean_matrix_mdpde_list.emplace_back(mean_matrix_mdpde);
        mean_matrix_train_list.emplace_back(mean_matrix_train);
        sigma_matrix_median_list.emplace_back(sigma_matrix_median);
        sigma_matrix_mdpde_list.emplace_back(sigma_matrix_mdpde);
        sigma_matrix_train_list.emplace_back(sigma_matrix_train);
    }

    PrintMemberOutlierResult(
        options,
        "bias_outlier_in_member.pdf",
        outlier_list,
        mean_matrix_median_list, mean_matrix_mdpde_list, mean_matrix_train_list,
        sigma_matrix_median_list, sigma_matrix_mdpde_list, sigma_matrix_train_list
    );
}

void RunSimulationTestOnModelAlphaData(const HRLModelTestExecutionContext & options)
{
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnModelAlphaData");

    const auto scenario{ BetaScenarioConfig{
        100,
        1000,
        std::vector<double>{ options.options.alpha_r }
    } };
    const auto model_par_prior{ MakeDefaultModelPrior() };
    auto data_factory{ BuildDataFactory(options) };
    HRLModelTester tester(kGausParSize);
    std::vector<double> error_list{ 0.1, 0.2, 0.3 };
    std::vector<Eigen::MatrixXd> mean_matrix_alpha1_list;
    std::vector<Eigen::MatrixXd> mean_matrix_alpha2_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_alpha1_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_alpha2_list;

    const auto outlier_list{ BuildLinearSweep(10, 0.05) };
    const auto outlier_size{ static_cast<int>(outlier_list.size()) };

    for (auto error_sigma : error_list)
    {
        Eigen::MatrixXd mean_matrix_alpha1{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd mean_matrix_alpha2{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd sigma_matrix_alpha1{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd sigma_matrix_alpha2{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        for (int i = 0; i < outlier_size; i++)
        {
            std::vector<Eigen::VectorXd> residual_mean_ols_list;
            std::vector<Eigen::VectorXd> residual_mean_mdpde_list;
            std::vector<Eigen::VectorXd> residual_sigma_ols_list;
            std::vector<Eigen::VectorXd> residual_sigma_mdpde_list;
            const auto test_input{
                data_factory.BuildBetaTestInput(HRLModelTestDataFactory::BetaScenario{
                    model_par_prior,
                    scenario.sampling_entry_size,
                    error_sigma,
                    outlier_list[static_cast<size_t>(i)],
                    scenario.replica_size
                })
            };
            tester.RunBetaMDPDETest(
                scenario.alpha_r_list,
                residual_mean_ols_list, residual_mean_mdpde_list,
                residual_sigma_ols_list, residual_sigma_mdpde_list,
                test_input,
                options.thread_size
            );

            mean_matrix_alpha1.col(i) = residual_mean_mdpde_list.front();
            mean_matrix_alpha2.col(i) = residual_mean_mdpde_list.back();
            sigma_matrix_alpha1.col(i) = residual_sigma_mdpde_list.front();
            sigma_matrix_alpha2.col(i) = residual_sigma_mdpde_list.back();
        }
        mean_matrix_alpha1_list.emplace_back(mean_matrix_alpha1);
        mean_matrix_alpha2_list.emplace_back(mean_matrix_alpha2);
        sigma_matrix_alpha1_list.emplace_back(sigma_matrix_alpha1);
        sigma_matrix_alpha2_list.emplace_back(sigma_matrix_alpha2);
    }

    PrintDataOutlierResult(
        options,
        "bias_outlier_with_alpha_in_data.pdf",
        outlier_list,
        mean_matrix_alpha1_list, mean_matrix_alpha2_list, mean_matrix_alpha2_list,
        sigma_matrix_alpha1_list, sigma_matrix_alpha2_list, sigma_matrix_alpha2_list
    );
}

void RunSimulationTestOnModelAlphaMember(const HRLModelTestExecutionContext & options)
{
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnModelAlphaMember");

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
    HRLModelTester tester(kGausParSize);

    std::vector<Eigen::MatrixXd> mean_matrix_alpha1_list;
    std::vector<Eigen::MatrixXd> mean_matrix_alpha2_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_alpha1_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_alpha2_list;

    const auto outlier_list{ BuildLinearSweep(10, 0.05) };
    const auto outlier_size{ static_cast<int>(outlier_list.size()) };

    for (auto outlier_prior : outlier_prior_list)
    {
        Eigen::MatrixXd mean_matrix_alpha1{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd mean_matrix_alpha2{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd sigma_matrix_alpha1{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        Eigen::MatrixXd sigma_matrix_alpha2{ Eigen::MatrixXd::Zero(kGausParSize, outlier_size) };
        for (int i = 0; i < outlier_size; i++)
        {
            std::vector<Eigen::VectorXd> residual_mean_median_list;
            std::vector<Eigen::VectorXd> residual_mean_mdpde_list;
            std::vector<Eigen::VectorXd> residual_sigma_median_list;
            std::vector<Eigen::VectorXd> residual_sigma_mdpde_list;
            const auto test_input{
                data_factory.BuildMuTestInput(HRLModelTestDataFactory::MuScenario{
                    scenario.member_size,
                    model_par_prior,
                    model_par_sigma,
                    outlier_prior,
                    model_par_sigma,
                    outlier_list[static_cast<size_t>(i)],
                    scenario.replica_size
                })
            };
            tester.RunMuMDPDETest(
                scenario.alpha_g_list,
                residual_mean_median_list, residual_mean_mdpde_list,
                residual_sigma_median_list, residual_sigma_mdpde_list,
                test_input,
                options.thread_size
            );

            mean_matrix_alpha1.col(i) = residual_mean_mdpde_list.front();
            mean_matrix_alpha2.col(i) = residual_mean_mdpde_list.back();
            sigma_matrix_alpha1.col(i) = residual_sigma_mdpde_list.front();
            sigma_matrix_alpha2.col(i) = residual_sigma_mdpde_list.back();
        }
        mean_matrix_alpha1_list.emplace_back(mean_matrix_alpha1);
        mean_matrix_alpha2_list.emplace_back(mean_matrix_alpha2);
        sigma_matrix_alpha1_list.emplace_back(sigma_matrix_alpha1);
        sigma_matrix_alpha2_list.emplace_back(sigma_matrix_alpha2);
    }

    PrintMemberOutlierResult(
        options,
        "bias_outlier_with_alpha_in_member.pdf",
        outlier_list,
        mean_matrix_alpha1_list, mean_matrix_alpha2_list, mean_matrix_alpha2_list,
        sigma_matrix_alpha1_list, sigma_matrix_alpha2_list, sigma_matrix_alpha2_list
    );
}

void RunSimulationTestOnNeighborDistance(const HRLModelTestExecutionContext & options)
{
    ScopeTimer timer("HRLModelTestCommand::RunSimulationTestOnNeighborDistance");

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
    HRLModelTester tester(kGausParSize);

    std::vector<Eigen::MatrixXd> mean_matrix_ols_list;
    std::vector<Eigen::MatrixXd> mean_matrix_mdpde_list;
    std::vector<Eigen::MatrixXd> mean_matrix_train_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_ols_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_mdpde_list;
    std::vector<Eigen::MatrixXd> sigma_matrix_train_list;

    const auto distance_size{ static_cast<int>(scenario.distance_list.size()) };
    bool is_print_sampling_summary{ false };
    for (auto error_sigma : scenario.error_list)
    {
        Eigen::MatrixXd mean_matrix_ols{ Eigen::MatrixXd::Zero(kGausParSize, distance_size) };
        Eigen::MatrixXd mean_matrix_mdpde{ Eigen::MatrixXd::Zero(kGausParSize, distance_size) };
        Eigen::MatrixXd mean_matrix_train{ Eigen::MatrixXd::Zero(kGausParSize, distance_size) };
        Eigen::MatrixXd sigma_matrix_ols{ Eigen::MatrixXd::Zero(kGausParSize, distance_size) };
        Eigen::MatrixXd sigma_matrix_mdpde{ Eigen::MatrixXd::Zero(kGausParSize, distance_size) };
        Eigen::MatrixXd sigma_matrix_train{ Eigen::MatrixXd::Zero(kGausParSize, distance_size) };
        std::vector<LocalPotentialSampleList> sampling_entries_list;
        sampling_entries_list.reserve(static_cast<size_t>(distance_size));
        double training_alpha_r_average{ 0.0 };
        for (int i = 0; i < distance_size; i++)
        {
            std::vector<Eigen::VectorXd> residual_mean_list;
            std::vector<Eigen::VectorXd> residual_sigma_list;
            const auto test_input{
                data_factory.BuildNeighborhoodTestInput(
                    BuildNeighborhoodScenario(
                        model_par_prior,
                        scenario,
                        error_sigma,
                        scenario.distance_list[static_cast<size_t>(i)],
                        true))
            };
            tester.RunBetaMDPDEWithNeighborhoodTest(
                residual_mean_list, residual_sigma_list,
                test_input,
                training_alpha_r_average,
                options.thread_size,
                scenario.rejected_angle
            );

            sampling_entries_list.emplace_back(test_input.sampling_summaries.front());

            mean_matrix_ols.col(i) = residual_mean_list.at(0);
            sigma_matrix_ols.col(i) = residual_sigma_list.at(0);
            mean_matrix_mdpde.col(i) = residual_mean_list.at(1);
            sigma_matrix_mdpde.col(i) = residual_sigma_list.at(1);
            mean_matrix_train.col(i) = residual_mean_list.at(2);
            sigma_matrix_train.col(i) = residual_sigma_list.at(2);

            Logger::Log(LogLevel::Info,
                std::string("Distance: ") + std::to_string(scenario.distance_list[static_cast<size_t>(i)])
                + " , OLS: " + std::to_string(mean_matrix_ols.col(i)(0)) + " +- " + std::to_string(sigma_matrix_ols.col(i)(0))
                + " , MDPDE: " + std::to_string(mean_matrix_mdpde.col(i)(0)) + " +- " + std::to_string(sigma_matrix_mdpde.col(i)(0))
                + " , Train: " + std::to_string(mean_matrix_train.col(i)(0)) + " +- " + std::to_string(sigma_matrix_train.col(i)(0))
                + " (Alpha-R = " + std::to_string(training_alpha_r_average) + ")"
            );
        }
        mean_matrix_ols_list.emplace_back(mean_matrix_ols);
        mean_matrix_mdpde_list.emplace_back(mean_matrix_mdpde);
        mean_matrix_train_list.emplace_back(mean_matrix_train);
        sigma_matrix_ols_list.emplace_back(sigma_matrix_ols);
        sigma_matrix_mdpde_list.emplace_back(sigma_matrix_mdpde);
        sigma_matrix_train_list.emplace_back(sigma_matrix_train);

        if (!is_print_sampling_summary)
        {
            detail::PrintAtomSamplingDataSummary(
                options,
                "neighbor_atom_sampling_entries.pdf",
                sampling_entries_list, scenario.distance_list
            );
            is_print_sampling_summary = true;
        }
    }

    PrintDataOutlierResult(
        options,
        "bias_from_neighbor_atom.pdf",
        scenario.distance_list,
        mean_matrix_ols_list, mean_matrix_mdpde_list, mean_matrix_train_list,
        sigma_matrix_ols_list, sigma_matrix_mdpde_list, sigma_matrix_train_list, true
    );
}

void PrintDataOutlierResult(
    const HRLModelTestExecutionContext & options,
    const std::string & name,
    const std::vector<double> & outlier_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_ols_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_mdpde_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_train_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_ols_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_mdpde_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_train_list,
    bool is_neighbor_distance)
{
    auto file_path{ BuildOutputPath(options, name) };
    Logger::Log(LogLevel::Info, " HRLModelTestCommand::PrintDataOutlierResult");

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

    std::vector<std::unique_ptr<TGraphErrors>> graph_ols_list[col_size][row_size];
    std::vector<std::unique_ptr<TGraphErrors>> graph_mdpde_list[col_size][row_size];
    std::vector<std::unique_ptr<TGraphErrors>> graph_train_list[col_size][row_size];
    std::vector<double> y_array[col_size][row_size];
    std::vector<double> global_y_array[row_size];
    for (size_t i = 0; i < col_size; i++)
    {
        for (size_t j = 0; j < row_size; j++)
        {
            auto graph_ols{ root_helper::CreateGraphErrors() };
            auto graph_mdpde{ root_helper::CreateGraphErrors() };
            auto graph_train{ root_helper::CreateGraphErrors() };
            for (int p = 0; p < static_cast<int>(outlier_list.size()); p++)
            {
                auto x_value{ (is_neighbor_distance) ?
                    outlier_list.at(static_cast<size_t>(p)) * -1.0:
                    outlier_list.at(static_cast<size_t>(p)) * 100.0
                };
                auto mean_ols{ mean_matrix_ols_list.at(i).col(p)(static_cast<int>(j)) };
                auto sigma_ols{ sigma_matrix_ols_list.at(i).col(p)(static_cast<int>(j)) };
                auto mean_mdpde{ mean_matrix_mdpde_list.at(i).col(p)(static_cast<int>(j)) };
                auto sigma_mdpde{ sigma_matrix_mdpde_list.at(i).col(p)(static_cast<int>(j)) };
                auto mean_train{ mean_matrix_train_list.at(i).col(p)(static_cast<int>(j)) };
                auto sigma_train{ sigma_matrix_train_list.at(i).col(p)(static_cast<int>(j)) };
                graph_ols->SetPoint(p, x_value, mean_ols);
                graph_ols->SetPointError(p, 0.0, sigma_ols);
                graph_mdpde->SetPoint(p, x_value, mean_mdpde);
                graph_mdpde->SetPointError(p, 0.0, sigma_mdpde);
                graph_train->SetPoint(p, x_value, mean_train);
                graph_train->SetPointError(p, 0.0, sigma_train);
                y_array[i][j].emplace_back(mean_ols);
                y_array[i][j].emplace_back(mean_mdpde);
                y_array[i][j].emplace_back(mean_train);
                global_y_array[j].emplace_back(mean_ols);
                global_y_array[j].emplace_back(mean_mdpde);
                global_y_array[j].emplace_back(mean_train);
            }
            graph_ols_list[i][j].emplace_back(std::move(graph_ols));
            graph_mdpde_list[i][j].emplace_back(std::move(graph_mdpde));
            graph_train_list[i][j].emplace_back(std::move(graph_train));
        }
    }

    double x_min[col_size]{ 0.0 };
    double x_max[col_size]{ 0.0 };
    double y_min[row_size]{ 0.0 };
    double y_max[row_size]{ 0.0 };
    for (size_t i = 0; i < col_size; i++)
    {
        x_min[i] = (options.options.tester_choice == TesterType::MODEL_ALPHA_DATA) ? -2.0 : -0.7;
        x_max[i] = (options.options.tester_choice == TesterType::MODEL_ALPHA_DATA) ? 47.0 : 22.0;
        if (is_neighbor_distance)
        {
            x_min[i] = -2.6;
            x_max[i] = -0.8;
        }
    }
    for (size_t j = 0; j < row_size; j++)
    {
        auto y_range{ ArrayStats<double>::ComputeScalingPercentileRangeTuple(global_y_array[j], 0.38, 0.005, 0.995) };
        y_min[j] = std::get<0>(y_range);
        y_max[j] = std::get<1>(y_range);
    }

    std::unique_ptr<TH2> frame[col_size][row_size];
    std::unique_ptr<TPaveText> resolution_text[col_size];
    std::unique_ptr<TPaveText> title_x_text[col_size];
    std::unique_ptr<TPaveText> title_y_text[row_size];
    double error_value[3]{ 0.0, 2.5, 5.0 };
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

            short color_ols{ kAzure };
            short color_mdpde{ kGreen+2 };
            short color_train{ kRed };
            for (auto & graph : graph_ols_list[i][par_id])
            {
                root_helper::SetMarkerAttribute(graph.get(), 24, 1.5f, color_ols);
                root_helper::SetLineAttribute(graph.get(), 2, 2, color_ols);
                root_helper::SetFillAttribute(graph.get(), 1001, color_ols, 0.2f);
                graph->Draw("PL3");
            }
            for (auto & graph : graph_mdpde_list[i][par_id])
            {
                root_helper::SetMarkerAttribute(graph.get(), 25, 1.5f, color_mdpde);
                root_helper::SetLineAttribute(graph.get(), 1, 2, color_mdpde);
                root_helper::SetFillAttribute(graph.get(), 1001, color_mdpde, 0.2f);
                graph->Draw("PL3");
            }
            for (auto & graph : graph_train_list[i][par_id])
            {
                root_helper::SetMarkerAttribute(graph.get(), 20, 1.5f, color_train);
                root_helper::SetLineAttribute(graph.get(), 3, 2, color_train);
                root_helper::SetFillAttribute(graph.get(), 1001, color_train, 0.2f);
                if (options.options.tester_choice != TesterType::MODEL_ALPHA_DATA)
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
                title_x_text[i]->AddText(Form("#sigma_{#epsilon} = %.1f%% D_{max}", error_value[i]));
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
    if (options.options.tester_choice == TesterType::MODEL_ALPHA_DATA)
    {
        legend->AddEntry(graph_ols_list[0][0].front().get(),
            "MDPDE (#alpha_{r} = 0.1)", "plf");
        legend->AddEntry(graph_mdpde_list[0][0].front().get(),
            "MDPDE (#alpha_{r} = 0.4)", "plf");
    }
    else
    {
        legend->AddEntry(graph_train_list[0][0].front().get(),
            "MDPDE (w/ Sampling Scheme)", "plf");
        legend->AddEntry(graph_mdpde_list[0][0].front().get(),
            "MDPDE", "plf");
        legend->AddEntry(graph_ols_list[0][0].front().get(),
            "Ordinary Least Squares", "plf");
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
    if (is_neighbor_distance)
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
    const HRLModelTestExecutionContext & options,
    const std::string & name,
    const std::vector<double> & outlier_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_median_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_mdpde_list,
    const std::vector<Eigen::MatrixXd> & mean_matrix_train_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_median_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_mdpde_list,
    const std::vector<Eigen::MatrixXd> & sigma_matrix_train_list)
{
    auto file_path{ BuildOutputPath(options, name) };
    Logger::Log(LogLevel::Info, " HRLModelTestCommand::PrintMemberOutlierResult");
    (void)outlier_list;
    (void)mean_matrix_median_list;
    (void)mean_matrix_mdpde_list;
    (void)mean_matrix_train_list;
    (void)sigma_matrix_median_list;
    (void)sigma_matrix_mdpde_list;
    (void)sigma_matrix_train_list;

    std::vector<std::string> title_y_list{
        "Amplitude #font[2]{A}", "Width #tau"
    };
    std::vector<std::string> outlier_type_list{
        "#font[2]{A}", "#tau"
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

    std::vector<std::unique_ptr<TGraphErrors>> graph_ols_list[col_size][row_size];
    std::vector<std::unique_ptr<TGraphErrors>> graph_mdpde_list[col_size][row_size];
    std::vector<std::unique_ptr<TGraphErrors>> graph_train_list[col_size][row_size];
    std::vector<double> global_y_array;
    for (size_t i = 0; i < col_size; i++)
    {
        for (size_t j = 0; j < row_size; j++)
        {
            auto graph_ols{ root_helper::CreateGraphErrors() };
            auto graph_mdpde{ root_helper::CreateGraphErrors() };
            auto graph_train{ root_helper::CreateGraphErrors() };
            for (int p = 0; p < static_cast<int>(outlier_list.size()); p++)
            {
                auto x_value{ outlier_list.at(static_cast<size_t>(p)) * 100.0 };
                auto mean_ols{ mean_matrix_median_list.at(i).col(p)(static_cast<int>(j)) };
                auto sigma_ols{ sigma_matrix_median_list.at(i).col(p)(static_cast<int>(j)) };
                auto mean_mdpde{ mean_matrix_mdpde_list.at(i).col(p)(static_cast<int>(j)) };
                auto sigma_mdpde{ sigma_matrix_mdpde_list.at(i).col(p)(static_cast<int>(j)) };
                auto mean_train{ mean_matrix_train_list.at(i).col(p)(static_cast<int>(j)) };
                auto sigma_train{ sigma_matrix_train_list.at(i).col(p)(static_cast<int>(j)) };
                graph_ols->SetPoint(p, x_value, mean_ols);
                graph_ols->SetPointError(p, 0.0, sigma_ols);
                graph_mdpde->SetPoint(p, x_value, mean_mdpde);
                graph_mdpde->SetPointError(p, 0.0, sigma_mdpde);
                graph_train->SetPoint(p, x_value, mean_train);
                graph_train->SetPointError(p, 0.0, sigma_train);
                global_y_array.emplace_back(mean_ols);
                global_y_array.emplace_back(mean_mdpde);
                global_y_array.emplace_back(mean_train);
            }
            graph_ols_list[i][j].emplace_back(std::move(graph_ols));
            graph_mdpde_list[i][j].emplace_back(std::move(graph_mdpde));
            graph_train_list[i][j].emplace_back(std::move(graph_train));
        }
    }

    double x_min[col_size]{ 0.0 };
    double x_max[col_size]{ 0.0 };
    double y_min[row_size]{ 0.0 };
    double y_max[row_size]{ 0.0 };
    for (size_t i = 0; i < col_size; i++)
    {
        x_min[i] = (options.options.tester_choice == TesterType::MODEL_ALPHA_MEMBER) ? -2.0 : -0.7;
        x_max[i] = (options.options.tester_choice == TesterType::MODEL_ALPHA_MEMBER) ? 47.0 : 22.0;
    }
    auto y_range{ ArrayStats<double>::ComputeScalingRangeTuple(global_y_array, 0.3) };
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

            short color_ols{ kAzure };
            short color_mdpde{ kRed };
            short color_train{ kGreen+1 };
            for (auto & graph : graph_ols_list[i][par_id])
            {
                root_helper::SetMarkerAttribute(graph.get(), 24, 1.5f, color_ols);
                root_helper::SetLineAttribute(graph.get(), 2, 2, color_ols);
                root_helper::SetFillAttribute(graph.get(), 1001, color_ols, 0.2f);
                graph->Draw("PL3");
            }
            for (auto & graph : graph_mdpde_list[i][par_id])
            {
                root_helper::SetMarkerAttribute(graph.get(), 20, 1.5f, color_mdpde);
                root_helper::SetLineAttribute(graph.get(), 1, 2, color_mdpde);
                root_helper::SetFillAttribute(graph.get(), 1001, color_mdpde, 0.2f);
                graph->Draw("PL3");
            }
            for (auto & graph : graph_train_list[i][par_id])
            {
                root_helper::SetMarkerAttribute(graph.get(), 25, 1.5f, color_train);
                root_helper::SetLineAttribute(graph.get(), 3, 2, color_train);
                root_helper::SetFillAttribute(graph.get(), 1001, color_train, 0.2f);
                if (options.options.tester_choice != TesterType::MODEL_ALPHA_MEMBER)
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
                title_x_text[i]->AddText(
                    Form("Outlier in %s", outlier_type_list[static_cast<size_t>(i)].data())
                );
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
    if (options.options.tester_choice == TesterType::MODEL_ALPHA_MEMBER)
    {
        legend->AddEntry(graph_ols_list[0][0].front().get(),
            "MDPDE (#alpha_{g} = 0.2)", "plf");
        legend->AddEntry(graph_mdpde_list[0][0].front().get(),
            "MDPDE (#alpha_{g} = 0.5)", "plf");
    }
    else
    {
        legend->AddEntry(graph_train_list[0][0].front().get(),
            "MDPDE (#alpha_{g} from Alg.5)", "plf");
        legend->AddEntry(graph_mdpde_list[0][0].front().get(),
            "MDPDE (#alpha_{g} = 0.2)", "plf");
        legend->AddEntry(graph_ols_list[0][0].front().get(),
            "MDPDE (#alpha_{g} = 0)", "plf");
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
    const HRLModelTestExecutionContext & options,
    const std::string & name,
    const std::vector<LocalPotentialSampleList> & sampling_entries_list,
    const std::vector<double> & distance_list
)
{
    auto file_path{ BuildOutputPath(options, name) };
    Logger::Log(LogLevel::Info, " HRLModelTestCommand::PrintAtomSamplingDataSummary");

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
        auto data_graph{ detail::CreateDistanceToResponseGraph(sampling_entries) };
        auto data_hist{ detail::CreateDistanceToResponseHistogram(sampling_entries, 40) };

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

bool RunHRLModelTestWorkflow(const HRLModelTestExecutionContext & options)
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

bool HRLModelTestCommand::ExecuteImpl()
{
    return detail::RunHRLModelTestWorkflow(detail::HRLModelTestExecutionContext{
        RequestOptions(),
        ThreadSize(),
        OutputFolder(),
    });
}

} // namespace rhbm_gem
