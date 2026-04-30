#include "RHBMTestCommand.hpp"
#include "detail/RHBMTestPlotting.hpp"

#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/TestDataFactory.hpp>
#include <rhbm_gem/utils/hrl/RHBMTester.hpp>

#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace rhbm_gem {

namespace rhbm_test_plotting = command_detail::rhbm_test_plotting;
using rhbm_test_plotting::AppendBiasCurvePoint;
using rhbm_test_plotting::BiasCurveKind;
using rhbm_test_plotting::BiasPlotFlavor;
using rhbm_test_plotting::BiasPlotPanel;
using rhbm_test_plotting::BiasPlotRequest;
using rhbm_test_plotting::BiasXAxisMode;
using rhbm_test_plotting::FormatDataBiasPanelLabel;
using rhbm_test_plotting::FormatMemberBiasPanelLabel;
using rhbm_test_plotting::MakeBiasCurve;
using rhbm_test_plotting::TryAppendBenchmarkLinearizedPanel;

namespace {
constexpr std::string_view kTesterOption{ "--tester" };
constexpr std::string_view kFitMinOption{ "--fit-min" };
constexpr std::string_view kFitMaxOption{ "--fit-max" };
constexpr std::string_view kFitRangeIssue{ "--fit-range" };
constexpr std::string_view kAlphaROption{ "--alpha-r" };
constexpr std::string_view kAlphaGOption{ "--alpha-g" };
constexpr int kGausParSize{ GaussianModel3D::ParameterSize() };

struct NeighborDistanceSweepConfig
{
    test_data_factory::NeighborhoodScenario base_scenario;
    std::vector<double> error_list;
    std::vector<double> distance_list;
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

test_data_factory::TestDataBuildOptions BuildTestDataOptions(
    const RHBMTestRequest & request)
{
    return test_data_factory::TestDataBuildOptions{
        linearization_service::LinearizationSpec::AtomDecode(),
        linearization_service::LinearizationRange{request.fit_range_min, request.fit_range_max}
    };
}

test_data_factory::NeighborhoodScenario BuildNeighborDistanceScenario(
    const test_data_factory::NeighborhoodScenario & base_scenario,
    double error_sigma,
    double neighbor_distance,
    bool include_sampling_summary = false)
{
    auto scenario{ base_scenario };
    scenario.data_error_sigma = error_sigma;
    scenario.neighbor_distance = neighbor_distance;
    scenario.include_sampling_summary = include_sampling_summary;
    return scenario;
}

test_data_factory::BetaScenario BuildBaseBetaScenario(
    const Eigen::VectorXd & model_par_prior,
    std::vector<double> alpha_r_list)
{
    test_data_factory::BetaScenario scenario;
    scenario.gaus_true = GaussianModel3D::FromVector(model_par_prior);
    scenario.sampling_entry_size = 50;
    scenario.data_error_sigma = 1.0;
    scenario.outlier_ratio = 0.0;
    scenario.replica_size = 10;
    scenario.requested_alpha_r_list = std::move(alpha_r_list);
    scenario.alpha_training = true;
    return scenario;
}

test_data_factory::MuScenario BuildBaseMuScenario(
    const Eigen::VectorXd & model_par_prior,
    const Eigen::VectorXd & model_par_sigma,
    std::vector<double> alpha_g_list)
{
    test_data_factory::MuScenario scenario;
    scenario.member_size = 100;
    scenario.gaus_prior = model_par_prior;
    scenario.gaus_sigma = model_par_sigma;
    scenario.outlier_prior = model_par_prior;
    scenario.outlier_sigma = model_par_sigma;
    scenario.outlier_ratio = 0.0;
    scenario.replica_size = 100;
    scenario.requested_alpha_g_list = std::move(alpha_g_list);
    scenario.alpha_training = true;
    return scenario;
}

test_data_factory::NeighborhoodScenario BuildBaseNeighborDistanceScenario(
    const Eigen::VectorXd & model_par_prior,
    size_t neighbor_count)
{
    test_data_factory::NeighborhoodScenario scenario;
    scenario.gaus_true = GaussianModel3D::FromVector(model_par_prior);
    scenario.sampling_entry_size = 50;
    scenario.data_error_sigma = 1.0;
    scenario.radius_min = 0.0;
    scenario.radius_max = 1.0;
    scenario.neighbor_distance = 2.0;
    scenario.neighbor_count = neighbor_count;
    scenario.rejected_angle = 45.0;
    scenario.include_sampling_summary = false;
    scenario.summary_radius_min = 0.0;
    scenario.summary_radius_max = 4.0;
    scenario.replica_size = 10;
    return scenario;
}

} // namespace

void RunSimulationTestOnBenchMark(const RHBMTestRequest & request)
{
    ScopeTimer timer("RHBMTestCommand::RunSimulationTestOnBenchMark");

    const auto error_sigma{ 0.1 };
    const auto model_par_prior{ MakeDefaultModelPrior() };
    const auto test_data_options{ BuildTestDataOptions(request) };

    std::vector<test_data_factory::AtomNeighborType> neighbor_type_list{
        test_data_factory::AtomNeighborType::O,
        test_data_factory::AtomNeighborType::N,
        test_data_factory::AtomNeighborType::C,
        test_data_factory::AtomNeighborType::CA
    };

    BiasPlotRequest plot_request;
    plot_request.output_name = "bias_from_neighbor_atom.pdf";
    plot_request.flavor = BiasPlotFlavor::NeighborType;
    plot_request.x_axis_mode = BiasXAxisMode::NeighborType;
    plot_request.panels.reserve(1);

    BiasPlotPanel panel;
    panel.label = "Neighbor Atom Type";
    panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Ols, neighbor_type_list.size()));
    panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Mdpde, neighbor_type_list.size()));
    panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::TrainedMdpde, neighbor_type_list.size()));

    for (const auto neighbor_type : neighbor_type_list)
    {
        test_data_factory::AtomNeighborhoodScenario base_scenario;
        base_scenario.gaus_true = GaussianModel3D::FromVector(model_par_prior);
        base_scenario.sampling_entry_size = 50;
        base_scenario.data_error_sigma = error_sigma;
        base_scenario.radius_min = test_data_options.fitting_range.min;
        base_scenario.radius_max = test_data_options.fitting_range.max;
        base_scenario.neighbor_type = neighbor_type;
        base_scenario.rejected_angle = 15.0;
        base_scenario.include_sampling_summary = false;
        base_scenario.summary_radius_min = 0.0;
        base_scenario.summary_radius_max = 4.0;
        base_scenario.replica_size = 10;

        std::vector<LinePlotPanel> linearized_panels;
        linearized_panels.reserve(1);

        const auto test_input{
            test_data_factory::BuildAtomNeighborhoodTestInput(base_scenario, test_data_options)
        };
        rhbm_tester::BetaMDPDETestBias no_cut_result;
        rhbm_tester::BetaMDPDETestBias cut_result;
        rhbm_tester::RunBetaMDPDETest(no_cut_result, test_input.no_cut_input, request.job_count);
        rhbm_tester::RunBetaMDPDETest(cut_result, test_input.cut_input, request.job_count);
        TryAppendBenchmarkLinearizedPanel(
            linearized_panels,
            0.0,
            test_input.no_cut_input.replica_datasets.front(),
            test_input.cut_input.replica_datasets.front());

        std::ostringstream stream;
        stream  << " OLS   (No Cut): " << std::setprecision(3) << std::fixed
                << no_cut_result.ols.mean(0) << " , "
                << no_cut_result.ols.mean(1)
                << " , MDPDE (No Cut): "
                << no_cut_result.mdpde.trained_alpha.value().mean(0) << " , "
                << no_cut_result.mdpde.trained_alpha.value().mean(1)
                << " , MDPDE (Cut):    "
                << cut_result.mdpde.trained_alpha.value().mean(0) << " , "
                << cut_result.mdpde.trained_alpha.value().mean(1)
                << " (Alpha-R = " << cut_result.mdpde.trained_alpha_average.value() << ")";
        Logger::Log(LogLevel::Info, stream.str());

        const auto neighbor_type_value{ static_cast<double>(neighbor_type) };
        AppendBiasCurvePoint(
            panel.curves.at(0),
            neighbor_type_value,
            no_cut_result.ols);
        AppendBiasCurvePoint(
            panel.curves.at(1),
            neighbor_type_value,
            no_cut_result.mdpde.trained_alpha.value());
        AppendBiasCurvePoint(
            panel.curves.at(2),
            neighbor_type_value,
            cut_result.mdpde.trained_alpha.value());
        //rhbm_test_plotting::SaveBenchmarkLinearizedDatasetReport(request, error_sigma, linearized_panels);
    }

    plot_request.panels.emplace_back(std::move(panel));

    rhbm_test_plotting::SaveDataOutlierBiasPlot(request, plot_request);
}

void RunSimulationTestOnDataOutlier(const RHBMTestRequest & request)
{
    ScopeTimer timer("RHBMTestCommand::RunSimulationTestOnDataOutlier");

    const auto model_par_prior{ MakeDefaultModelPrior() };
    const auto base_scenario{
        BuildBaseBetaScenario(model_par_prior, std::vector<double>{ request.alpha_r })
    };
    const auto test_data_options{ BuildTestDataOptions(request) };
    std::vector<double> error_list{ 0.1, 0.2, 0.3 };
    const auto outlier_list{ BuildLinearSweep(9, 0.025) };
    BiasPlotRequest plot_request;
    plot_request.output_name = "bias_outlier_in_data.pdf";
    plot_request.flavor = BiasPlotFlavor::DataOutlier;
    plot_request.x_axis_mode = BiasXAxisMode::ContaminationRatio;
    plot_request.panels.reserve(error_list.size());

    for (size_t panel_index = 0; panel_index < error_list.size(); panel_index++)
    {
        const auto error_sigma{ error_list.at(panel_index) };
        BiasPlotPanel panel;
        panel.label = FormatDataBiasPanelLabel(panel_index);
        panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Ols, outlier_list.size()));
        panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Mdpde, outlier_list.size()));
        panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::TrainedMdpde, outlier_list.size()));
        for (size_t i = 0; i < outlier_list.size(); i++)
        {
            rhbm_tester::BetaMDPDETestBias bias;
            auto beta_scenario{ base_scenario };
            beta_scenario.data_error_sigma = error_sigma;
            beta_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{
                test_data_factory::BuildBetaTestInput(beta_scenario, test_data_options)
            };
            rhbm_tester::RunBetaMDPDETest(bias, test_input, request.job_count);

            AppendBiasCurvePoint(panel.curves.at(0), outlier_list.at(i), bias.ols);
            AppendBiasCurvePoint(
                panel.curves.at(1),
                outlier_list.at(i),
                bias.mdpde.requested_alpha.front());
            if (bias.mdpde.trained_alpha.has_value())
            {
                AppendBiasCurvePoint(
                    panel.curves.at(2),
                    outlier_list.at(i),
                    bias.mdpde.trained_alpha.value());
            }
        }
        plot_request.panels.emplace_back(std::move(panel));
    }

    rhbm_test_plotting::SaveDataOutlierBiasPlot(request, plot_request);
}

void RunSimulationTestOnMemberOutlier(const RHBMTestRequest & request)
{
    ScopeTimer timer("RHBMTestCommand::RunSimulationTestOnMemberOutlier");

    auto outlier_prior_list{ std::vector<Eigen::VectorXd>{
        Eigen::VectorXd::Zero(kGausParSize),
        Eigen::VectorXd::Zero(kGausParSize)
    } };
    outlier_prior_list.at(0)(0) = 1.50;
    outlier_prior_list.at(0)(1) = 0.50;
    outlier_prior_list.at(0)(2) = 0.10;
    outlier_prior_list.at(1)(0) = 1.00;
    outlier_prior_list.at(1)(1) = 1.00;
    outlier_prior_list.at(1)(2) = 0.10;

    const auto model_par_prior{ MakeDefaultModelPrior() };
    const auto model_par_sigma{ MakeDefaultModelSigma() };
    const auto base_scenario{
        BuildBaseMuScenario(model_par_prior, model_par_sigma, std::vector<double>{ request.alpha_g })
    };
    const auto test_data_options{ BuildTestDataOptions(request) };
    const auto outlier_list{ BuildLinearSweep(9, 0.025) };
    BiasPlotRequest plot_request;
    plot_request.output_name = "bias_outlier_in_member.pdf";
    plot_request.flavor = BiasPlotFlavor::MemberOutlier;
    plot_request.x_axis_mode = BiasXAxisMode::ContaminationRatio;
    plot_request.panels.reserve(outlier_prior_list.size());

    for (size_t panel_index = 0; panel_index < outlier_prior_list.size(); panel_index++)
    {
        const auto & outlier_prior{ outlier_prior_list.at(panel_index) };
        BiasPlotPanel panel;
        panel.label = FormatMemberBiasPanelLabel(panel_index);
        panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Median, outlier_list.size()));
        panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Mdpde, outlier_list.size()));
        panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::TrainedMdpde, outlier_list.size()));
        for (size_t i = 0; i < outlier_list.size(); i++)
        {
            rhbm_tester::MuMDPDETestBias bias;
            auto mu_scenario{ base_scenario };
            mu_scenario.outlier_prior = outlier_prior;
            mu_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{
                test_data_factory::BuildMuTestInput(mu_scenario, test_data_options)
            };
            rhbm_tester::RunMuMDPDETest(bias, test_input, request.job_count);

            AppendBiasCurvePoint(panel.curves.at(0), outlier_list.at(i), bias.median);
            AppendBiasCurvePoint(
                panel.curves.at(1),
                outlier_list.at(i),
                bias.mdpde.requested_alpha.front());
            if (bias.mdpde.trained_alpha.has_value())
            {
                AppendBiasCurvePoint(
                    panel.curves.at(2),
                    outlier_list.at(i),
                    bias.mdpde.trained_alpha.value());
            }
        }
        plot_request.panels.emplace_back(std::move(panel));
    }

    rhbm_test_plotting::SaveMemberOutlierBiasPlot(request, plot_request);
}

void RunSimulationTestOnModelAlphaData(const RHBMTestRequest & request)
{
    ScopeTimer timer("RHBMTestCommand::RunSimulationTestOnModelAlphaData");

    const auto model_par_prior{ MakeDefaultModelPrior() };
    const auto base_scenario{
        BuildBaseBetaScenario(model_par_prior, std::vector<double>{ request.alpha_r })
    };
    const auto test_data_options{ BuildTestDataOptions(request) };
    std::vector<double> error_list{ 0.1, 0.2, 0.3 };
    const auto outlier_list{ BuildLinearSweep(10, 0.05) };
    BiasPlotRequest plot_request;
    plot_request.output_name = "bias_outlier_with_alpha_in_data.pdf";
    plot_request.flavor = BiasPlotFlavor::ModelAlphaData;
    plot_request.x_axis_mode = BiasXAxisMode::ContaminationRatio;
    plot_request.panels.reserve(error_list.size());

    for (size_t panel_index = 0; panel_index < error_list.size(); panel_index++)
    {
        const auto error_sigma{ error_list.at(panel_index) };
        BiasPlotPanel panel;
        panel.label = FormatDataBiasPanelLabel(panel_index);
        panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::RequestedAlpha, outlier_list.size()));
        panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::TrainedAlpha, outlier_list.size()));
        for (size_t i = 0; i < outlier_list.size(); i++)
        {
            rhbm_tester::BetaMDPDETestBias bias;
            auto beta_scenario{ base_scenario };
            beta_scenario.data_error_sigma = error_sigma;
            beta_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{
                test_data_factory::BuildBetaTestInput(beta_scenario, test_data_options)
            };
            rhbm_tester::RunBetaMDPDETest(bias, test_input, request.job_count);

            AppendBiasCurvePoint(
                panel.curves.at(0),
                outlier_list.at(i),
                bias.mdpde.requested_alpha.front());
            if (bias.mdpde.trained_alpha.has_value())
            {
                AppendBiasCurvePoint(
                    panel.curves.at(1),
                    outlier_list.at(i),
                    bias.mdpde.trained_alpha.value());
            }
        }
        plot_request.panels.emplace_back(std::move(panel));
    }

    rhbm_test_plotting::SaveDataOutlierBiasPlot(request, plot_request);
}

void RunSimulationTestOnModelAlphaMember(const RHBMTestRequest & request)
{
    ScopeTimer timer("RHBMTestCommand::RunSimulationTestOnModelAlphaMember");

    auto outlier_prior_list{ std::vector<Eigen::VectorXd>{
        Eigen::VectorXd::Zero(kGausParSize),
        Eigen::VectorXd::Zero(kGausParSize)
    } };
    outlier_prior_list.at(0)(0) = 1.50;
    outlier_prior_list.at(0)(1) = 0.50;
    outlier_prior_list.at(0)(2) = 0.10;
    outlier_prior_list.at(1)(0) = 1.00;
    outlier_prior_list.at(1)(1) = 1.00;
    outlier_prior_list.at(1)(2) = 0.10;

    const auto model_par_prior{ MakeDefaultModelPrior() };
    const auto model_par_sigma{ MakeDefaultModelSigma() };
    const auto base_scenario{
        BuildBaseMuScenario(model_par_prior, model_par_sigma, std::vector<double>{ request.alpha_g })
    };
    const auto test_data_options{ BuildTestDataOptions(request) };
    const auto outlier_list{ BuildLinearSweep(10, 0.05) };
    BiasPlotRequest plot_request;
    plot_request.output_name = "bias_outlier_with_alpha_in_member.pdf";
    plot_request.flavor = BiasPlotFlavor::ModelAlphaMember;
    plot_request.x_axis_mode = BiasXAxisMode::ContaminationRatio;
    plot_request.panels.reserve(outlier_prior_list.size());

    for (size_t panel_index = 0; panel_index < outlier_prior_list.size(); panel_index++)
    {
        const auto & outlier_prior{ outlier_prior_list.at(panel_index) };
        BiasPlotPanel panel;
        panel.label = FormatMemberBiasPanelLabel(panel_index);
        panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::RequestedAlpha, outlier_list.size()));
        panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::TrainedAlpha, outlier_list.size()));
        for (size_t i = 0; i < outlier_list.size(); i++)
        {
            rhbm_tester::MuMDPDETestBias bias;
            auto mu_scenario{ base_scenario };
            mu_scenario.outlier_prior = outlier_prior;
            mu_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{
                test_data_factory::BuildMuTestInput(mu_scenario, test_data_options)
            };
            rhbm_tester::RunMuMDPDETest(bias, test_input, request.job_count);

            AppendBiasCurvePoint(
                panel.curves.at(0),
                outlier_list.at(i),
                bias.mdpde.requested_alpha.front());
            if (bias.mdpde.trained_alpha.has_value())
            {
                AppendBiasCurvePoint(
                    panel.curves.at(1),
                    outlier_list.at(i),
                    bias.mdpde.trained_alpha.value());
            }
        }
        plot_request.panels.emplace_back(std::move(panel));
    }

    rhbm_test_plotting::SaveMemberOutlierBiasPlot(request, plot_request);
}

void RunSimulationTestOnNeighborDistance(const RHBMTestRequest & request)
{
    ScopeTimer timer("RHBMTestCommand::RunSimulationTestOnNeighborDistance");

    const auto model_par_prior{ MakeDefaultModelPrior() };
    const auto sweep{ NeighborDistanceSweepConfig{
        BuildBaseNeighborDistanceScenario(model_par_prior, 2),
        std::vector<double>{ 0.0, 0.025, 0.05 },
        BuildDescendingSweep(16, 2.5, 0.1)
    } };
    const auto test_data_options{ BuildTestDataOptions(request) };
    BiasPlotRequest plot_request;
    plot_request.output_name = "bias_from_neighbor_atom.pdf";
    plot_request.flavor = BiasPlotFlavor::NeighborDistance;
    plot_request.x_axis_mode = BiasXAxisMode::NeighborDistance;
    plot_request.panels.reserve(sweep.error_list.size());

    bool is_print_sampling_summary{ false };
    for (size_t panel_index = 0; panel_index < sweep.error_list.size(); panel_index++)
    {
        const auto error_sigma{ sweep.error_list.at(panel_index) };
        BiasPlotPanel panel;
        panel.label = FormatDataBiasPanelLabel(panel_index);
        panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Ols, sweep.distance_list.size()));
        panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Mdpde, sweep.distance_list.size()));
        panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::TrainedMdpde, sweep.distance_list.size()));
        std::vector<LocalPotentialSampleList> sampling_entries_list;
        sampling_entries_list.reserve(sweep.distance_list.size());
        for (size_t i = 0; i < sweep.distance_list.size(); i++)
        {
            const auto test_input{
                test_data_factory::BuildNeighborhoodTestInput(
                    BuildNeighborDistanceScenario(
                        sweep.base_scenario,
                        error_sigma,
                        sweep.distance_list.at(i),
                        true),
                    test_data_options)
            };
            rhbm_tester::BetaMDPDETestBias no_cut_result;
            rhbm_tester::BetaMDPDETestBias cut_result;
            rhbm_tester::RunBetaMDPDETest(no_cut_result, test_input.no_cut_input, request.job_count);
            rhbm_tester::RunBetaMDPDETest(cut_result, test_input.cut_input, request.job_count);
            sampling_entries_list.emplace_back(test_input.sampling_summaries.front());

            AppendBiasCurvePoint(
                panel.curves.at(0),
                sweep.distance_list.at(i),
                no_cut_result.ols);
            AppendBiasCurvePoint(
                panel.curves.at(1),
                sweep.distance_list.at(i),
                no_cut_result.mdpde.trained_alpha.value());
            AppendBiasCurvePoint(
                panel.curves.at(2),
                sweep.distance_list.at(i),
                cut_result.mdpde.trained_alpha.value());

            Logger::Log(LogLevel::Info,
                std::string("Distance: ") + std::to_string(sweep.distance_list.at(i))
                + " , OLS: " + std::to_string(no_cut_result.ols.mean(0)) + " +- " + std::to_string(no_cut_result.ols.sigma(0))
                + " , MDPDE: " + std::to_string(no_cut_result.mdpde.trained_alpha.value().mean(0)) + " +- " + std::to_string(no_cut_result.mdpde.trained_alpha.value().sigma(0))
                + " , Train: " + std::to_string(cut_result.mdpde.trained_alpha.value().mean(0)) + " +- " + std::to_string(cut_result.mdpde.trained_alpha.value().sigma(0))
                + " (Alpha-R = " + std::to_string(cut_result.mdpde.trained_alpha_average.value()) + ")"
            );
        }
        plot_request.panels.emplace_back(std::move(panel));

        if (!is_print_sampling_summary)
        {
            rhbm_test_plotting::SaveAtomSamplingSummary(
                request,
                "neighbor_atom_sampling_entries.pdf",
                sampling_entries_list, sweep.distance_list
            );
            is_print_sampling_summary = true;
        }
    }

    rhbm_test_plotting::SaveDataOutlierBiasPlot(request, plot_request);
}

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
    const auto & request{ RequestOptions() };
    switch (request.tester_choice)
    {
    case TesterType::BENCHMARK:
        RunSimulationTestOnBenchMark(request);
        return true;
    case TesterType::DATA_OUTLIER:
        RunSimulationTestOnDataOutlier(request);
        return true;
    case TesterType::MEMBER_OUTLIER:
        RunSimulationTestOnMemberOutlier(request);
        return true;
    case TesterType::MODEL_ALPHA_DATA:
        RunSimulationTestOnModelAlphaData(request);
        return true;
    case TesterType::MODEL_ALPHA_MEMBER:
        RunSimulationTestOnModelAlphaMember(request);
        return true;
    case TesterType::NEIGHBOR_DISTANCE:
        RunSimulationTestOnNeighborDistance(request);
        return true;
    default:
        Logger::Log(
            LogLevel::Error,
            "Invalid tester choice reached execution path: ["
                + std::to_string(static_cast<int>(request.tester_choice)) + "]");
        return false;
    }
}

} // namespace rhbm_gem
