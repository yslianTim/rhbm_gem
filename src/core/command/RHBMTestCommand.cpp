#include "detail/CommandBase.hpp"
#include "detail/RHBMTestPlotting.hpp"

#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/SampleFilter.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/TestDataFactory.hpp>
#include <rhbm_gem/utils/hrl/RHBMTester.hpp>

#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace rhbm_gem::core {

class RHBMTestCommand final : public CommandBase<RHBMTestRequest>
{
public:
    RHBMTestCommand();

private:
    void NormalizeAndValidateRequest(RHBMTestRequest & request) override;
    void ValidatePreparedRequest(const RHBMTestRequest & request) override;
    bool ExecuteImpl(const RHBMTestRequest & request) override;
};

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

namespace {

GaussianModel3D MakeDefaultModelPrior()
{
    return GaussianModel3D{ 1.0, 0.5, 0.0 };
}

GaussianModel3DUncertainty MakeDefaultModelSigma()
{
    return GaussianModel3DUncertainty{ 0.050, 0.025, 0.010 };
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

std::vector<RHBMMemberDataset> BuildReplicaDatasets(
    const std::vector<LocalPotentialSampleList> & replica_sampling_entries,
    double fit_range_min,
    double fit_range_max,
    LocalGaussianFitModel local_fit_model)
{
    std::vector<RHBMMemberDataset> replica_datasets;
    replica_datasets.reserve(replica_sampling_entries.size());
    for (const auto & sampling_entries : replica_sampling_entries)
    {
        replica_datasets.emplace_back(
            rhbm_helper::BuildMemberDataset(
                sampling_entries,
                fit_range_min,
                fit_range_max,
                local_fit_model));
    }
    return replica_datasets;
}

test_data_factory::LocalTestData BuildSelectedLocalTestData(
    const test_data_factory::LocalTestData & input)
{
    test_data_factory::LocalTestData selected_input;
    selected_input.gaus_true = input.gaus_true;
    selected_input.replica_sampling_entries.reserve(input.replica_sampling_entries.size());
    for (const auto & sampling_entries : input.replica_sampling_entries)
    {
        LocalPotentialSampleList selected_entries;
        selected_entries.reserve(sampling_entries.size());
        for (const auto & sample : sampling_entries)
        {
            if (sample.point.is_selected) selected_entries.emplace_back(sample);
        }
        auto filtered_entries{ sample_filter::FilterLocalPotentialSampleList(std::move(selected_entries)) };
        selected_input.replica_sampling_entries.emplace_back(std::move(filtered_entries));
    }
    return selected_input;
}

const rhbm_tester::BiasStatistics & SelectBenchmarkMdpdeBias(const rhbm_tester::BetaMDPDETestBias & bias)
{
    if (bias.mdpde.trained_alpha.has_value())
    {
        return bias.mdpde.trained_alpha.value();
    }
    return bias.mdpde.requested_alpha;
}

double SelectBenchmarkAlphaR(const rhbm_tester::BetaMDPDETestBias & bias, double requested_alpha_r)
{
    if (bias.mdpde.trained_alpha_median.has_value())
    {
        return bias.mdpde.trained_alpha_median.value();
    }
    return requested_alpha_r;
}

rhbm_tester::BetaMDPDETestOptions MakeBetaTestOptions(const RHBMTestRequest & request)
{
    rhbm_tester::BetaMDPDETestOptions options;
    options.requested_alpha_r = request.alpha_r;
    options.alpha_training = true;
    options.thread_size = request.job_count;
    return options;
}

rhbm_tester::MuMDPDETestOptions MakeMuTestOptions(const RHBMTestRequest & request)
{
    rhbm_tester::MuMDPDETestOptions options;
    options.requested_alpha_g = request.alpha_g;
    options.alpha_training = true;
    options.thread_size = request.job_count;
    return options;
}

} // namespace

void RunSimulationTestOnBenchMark(const RHBMTestRequest & request)
{
    const auto error_sigma{ 0.01 };
    const auto model_prior{ MakeDefaultModelPrior() };
    const auto beta_options{ MakeBetaTestOptions(request) };

    std::vector<Spot> spot_list{ Spot::O, Spot::N, Spot::C, Spot::CA };

    BiasPlotRequest plot_request;
    plot_request.output_name = "bias_from_neighbor_atom.pdf";
    plot_request.flavor = BiasPlotFlavor::NeighborType;
    plot_request.x_axis_mode = BiasXAxisMode::NeighborType;
    plot_request.panels.reserve(1);

    BiasPlotPanel panel;
    panel.label = "Neighbor Atom Type";
    panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Ols, spot_list.size()));
    panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Mdpde, spot_list.size()));
    panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::TrainedMdpde, spot_list.size()));

    std::vector<LinePlotPanel> linearized_panels;
    linearized_panels.reserve(spot_list.size());
    for (size_t i = 0; i < spot_list.size(); i++)
    {
        test_data_factory::AtomModelScenario base_scenario;
        base_scenario.gaus_true = model_prior;
        base_scenario.data_error_sigma = error_sigma;
        base_scenario.spot = spot_list.at(i);
        base_scenario.replica_size = 10;

        const auto input{ test_data_factory::BuildLocalTestData(base_scenario) };
        const auto cut_input{ BuildSelectedLocalTestData(input) };
        const auto no_cut_result{
            rhbm_tester::RunBetaMDPDETest(input, beta_options)
        };
        const auto cut_result{
            rhbm_tester::RunBetaMDPDETest(cut_input, beta_options)
        };
        const auto datasets{
            BuildReplicaDatasets(
                input.replica_sampling_entries,
                request.fit_range_min,
                request.fit_range_max,
                LocalGaussianFitModel::LogQuadratic)
        };
        const auto cut_datasets{
            BuildReplicaDatasets(
                cut_input.replica_sampling_entries,
                request.fit_range_min,
                request.fit_range_max,
                LocalGaussianFitModel::LogQuadratic)
        };
        rhbm_test_plotting::TryAppendBenchmarkLinearizedPanel(
            linearized_panels,
            0.0,
            datasets.front(),
            cut_datasets.front());

        const auto & no_cut_mdpde_bias{ SelectBenchmarkMdpdeBias(no_cut_result) };
        const auto & cut_mdpde_bias{ SelectBenchmarkMdpdeBias(cut_result) };
        std::ostringstream stream;
        stream  << " OLS   (No Cut): " << std::setprecision(3) << std::fixed
                << no_cut_result.ols.mean(0) << " , "
                << no_cut_result.ols.mean(1)
                << " , MDPDE (No Cut): "
                << no_cut_mdpde_bias.mean(0) << " , "
                << no_cut_mdpde_bias.mean(1)
                << " , MDPDE (Cut):    "
                << cut_mdpde_bias.mean(0) << " , "
                << cut_mdpde_bias.mean(1)
                << " (Alpha-R = " << SelectBenchmarkAlphaR(cut_result, request.alpha_r) << ")";
        Logger::Log(LogLevel::Info, stream.str());

        const auto spot_axis_value{ static_cast<double>(i + 1) };
        AppendBiasCurvePoint(panel.curves.at(0), spot_axis_value, no_cut_result.ols);
        AppendBiasCurvePoint(panel.curves.at(1), spot_axis_value, no_cut_mdpde_bias);
        AppendBiasCurvePoint(panel.curves.at(2), spot_axis_value, cut_mdpde_bias);
    }
    plot_request.panels.emplace_back(std::move(panel));
    rhbm_test_plotting::SaveBenchmarkLinearizedDatasetReport(request, error_sigma, linearized_panels);
    rhbm_test_plotting::SaveDataOutlierBiasPlot(request, plot_request);
}

void RunSimulationTestOnDataOutlier(const RHBMTestRequest & request)
{
    ScopeTimer timer("RHBMTestCommand::RunSimulationTestOnDataOutlier");

    const auto model_prior{ MakeDefaultModelPrior() };
    const auto beta_options{ MakeBetaTestOptions(request) };
    test_data_factory::LocalScenario base_scenario;
    base_scenario.gaus_true = model_prior;
    base_scenario.sampling_entry_size = 50;
    base_scenario.data_error_sigma = 0.0;
    base_scenario.outlier_ratio = 0.0;
    base_scenario.replica_size = 10;
    std::vector<double> error_list{ 0.0, 0.05, 0.1 };
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
            auto beta_scenario{ base_scenario };
            beta_scenario.data_error_sigma = error_sigma;
            beta_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{ test_data_factory::BuildLocalTestData(beta_scenario) };
            const auto bias{ rhbm_tester::RunBetaMDPDETest(test_input, beta_options) };

            AppendBiasCurvePoint(panel.curves.at(0), outlier_list.at(i), bias.ols);
            AppendBiasCurvePoint(
                panel.curves.at(1),
                outlier_list.at(i),
                bias.mdpde.requested_alpha);
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

    const std::vector<GaussianModel3D> outlier_prior_list{
        GaussianModel3D{ 1.50, 0.50, 0.10 },
        GaussianModel3D{ 1.00, 1.00, 0.10 }
    };

    const auto model_prior{ MakeDefaultModelPrior() };
    const auto model_sigma{ MakeDefaultModelSigma() };
    const auto mu_options{ MakeMuTestOptions(request) };
    test_data_factory::GroupScenario base_scenario;
    base_scenario.member_size = 100;
    base_scenario.inlier_distribution = { model_prior, model_sigma };
    base_scenario.outlier_distribution = { model_prior, model_sigma };
    base_scenario.outlier_ratio = 0.0;
    base_scenario.replica_size = 100;
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
            auto mu_scenario{ base_scenario };
            mu_scenario.outlier_distribution.mean = outlier_prior;
            mu_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{
                test_data_factory::BuildGroupTestData(mu_scenario)
            };
            const auto bias{ rhbm_tester::RunMuMDPDETest(test_input, mu_options) };

            AppendBiasCurvePoint(panel.curves.at(0), outlier_list.at(i), bias.median);
            AppendBiasCurvePoint(
                panel.curves.at(1),
                outlier_list.at(i),
                bias.mdpde.requested_alpha);
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

    const auto model_prior{ MakeDefaultModelPrior() };
    const auto beta_options{ MakeBetaTestOptions(request) };
    test_data_factory::LocalScenario base_scenario;
    base_scenario.gaus_true = model_prior;
    base_scenario.sampling_entry_size = 50;
    base_scenario.data_error_sigma = 1.0;
    base_scenario.outlier_ratio = 0.0;
    base_scenario.replica_size = 10;
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
            auto beta_scenario{ base_scenario };
            beta_scenario.data_error_sigma = error_sigma;
            beta_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{ test_data_factory::BuildLocalTestData(beta_scenario) };
            const auto bias{ rhbm_tester::RunBetaMDPDETest(test_input, beta_options) };

            AppendBiasCurvePoint(
                panel.curves.at(0),
                outlier_list.at(i),
                bias.mdpde.requested_alpha);
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

    const std::vector<GaussianModel3D> outlier_prior_list{
        GaussianModel3D{ 1.50, 0.50, 0.10 },
        GaussianModel3D{ 1.00, 1.00, 0.10 }
    };

    const auto model_prior{ MakeDefaultModelPrior() };
    const auto model_sigma{ MakeDefaultModelSigma() };
    const auto mu_options{ MakeMuTestOptions(request) };
    test_data_factory::GroupScenario base_scenario;
    base_scenario.member_size = 100;
    base_scenario.inlier_distribution = { model_prior, model_sigma };
    base_scenario.outlier_distribution = { model_prior, model_sigma };
    base_scenario.outlier_ratio = 0.0;
    base_scenario.replica_size = 100;
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
            auto mu_scenario{ base_scenario };
            mu_scenario.outlier_distribution.mean = outlier_prior;
            mu_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{
                test_data_factory::BuildGroupTestData(mu_scenario)
            };
            const auto bias{ rhbm_tester::RunMuMDPDETest(test_input, mu_options) };

            AppendBiasCurvePoint(
                panel.curves.at(0),
                outlier_list.at(i),
                bias.mdpde.requested_alpha);
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

RHBMTestCommand::RHBMTestCommand() : CommandBase<RHBMTestRequest>{}
{
}

void RHBMTestCommand::NormalizeAndValidateRequest(RHBMTestRequest & request)
{
    RequireEnum(request, &RHBMTestRequest::tester_choice);
    RequireFiniteNonNegativeScalar(request, &RHBMTestRequest::fit_range_min);
    RequireFiniteNonNegativeScalar(request, &RHBMTestRequest::fit_range_max);
    RequireFinitePositiveScalar(request, &RHBMTestRequest::alpha_r);
    RequireFinitePositiveScalar(request, &RHBMTestRequest::alpha_g);
}

void RHBMTestCommand::ValidatePreparedRequest(const RHBMTestRequest & request)
{
    RequirePrepareCondition(
        request.fit_range_min <= request.fit_range_max,
        "Expected --fit-min <= --fit-max.");
}

bool RHBMTestCommand::ExecuteImpl(const RHBMTestRequest & request)
{
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
    default:
        Logger::Log(LogLevel::Error,
            "Invalid tester choice reached execution path: ["
                + std::to_string(static_cast<int>(request.tester_choice)) + "]");
        return false;
    }
}

namespace command_internal {

CommandResult ExecuteRHBMTestCommand(const RHBMTestRequest & request)
{
    RHBMTestCommand command;
    return command.ExecuteRequest(request);
}

} // namespace command_internal

} // namespace rhbm_gem::core
