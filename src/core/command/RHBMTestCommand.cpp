#include "detail/CommandBase.hpp"
#include "detail/RHBMTestPlotting.hpp"

#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/domain/ScopeTimer.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/TestDataFactory.hpp>
#include <rhbm_gem/utils/hrl/RHBMTester.hpp>

#include <iomanip>
#include <stdexcept>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

namespace rhbm_gem {

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

Eigen::VectorXd MakeDefaultModelPrior()
{
    Eigen::VectorXd model_par_prior{ Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize()) };
    model_par_prior(0) = 1.0;
    model_par_prior(1) = 0.5;
    model_par_prior(2) = 0.0;
    return model_par_prior;
}

Eigen::VectorXd MakeDefaultModelSigma()
{
    Eigen::VectorXd model_par_sigma{ Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize()) };
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

const rhbm_tester::BiasStatistics & SelectBenchmarkMdpdeBias(const rhbm_tester::BetaMDPDETestBias & bias)
{
    if (bias.mdpde.trained_alpha.has_value())
    {
        return bias.mdpde.trained_alpha.value();
    }
    if (!bias.mdpde.requested_alpha.empty())
    {
        return bias.mdpde.requested_alpha.front();
    }
    throw std::invalid_argument(
        "Benchmark MDPDE output requires trained alpha or at least one requested alpha.");
}

double SelectBenchmarkAlphaR(
    const rhbm_tester::BetaMDPDETestBias & bias,
    double requested_alpha_r)
{
    if (bias.mdpde.trained_alpha_average.has_value())
    {
        return bias.mdpde.trained_alpha_average.value();
    }
    return requested_alpha_r;
}

} // namespace

void RunSimulationTestOnBenchMark(const RHBMTestRequest & request)
{
    const auto error_sigma{ 0.01 };
    const auto model_par_prior{ MakeDefaultModelPrior() };

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
        base_scenario.gaus_true = GaussianModel3D::FromVector(model_par_prior);
        base_scenario.data_error_sigma = error_sigma;
        base_scenario.spot = spot_list.at(i);
        base_scenario.rejected_angle = 15.0;
        base_scenario.replica_size = 10;
        base_scenario.requested_alpha_r_list = { request.alpha_r };
        base_scenario.alpha_training = true;
        base_scenario.local_fit_model = LocalGaussianFitModel::LogQuadratic;

        const auto test_input{ test_data_factory::BuildAtomModelTestData(base_scenario) };
        const auto no_cut_result{
            rhbm_tester::RunBetaMDPDETest(test_input.no_cut_input, request.job_count)
        };
        const auto cut_result{
            rhbm_tester::RunBetaMDPDETest(test_input.cut_input, request.job_count)
        };
        const auto no_cut_datasets{
            BuildReplicaDatasets(
                test_input.no_cut_input.replica_sampling_entries,
                request.fit_range_min,
                request.fit_range_max,
                test_input.no_cut_input.local_fit_model)
        };
        const auto cut_datasets{
            BuildReplicaDatasets(
                test_input.cut_input.replica_sampling_entries,
                request.fit_range_min,
                request.fit_range_max,
                test_input.cut_input.local_fit_model)
        };
        rhbm_test_plotting::TryAppendBenchmarkLinearizedPanel(
            linearized_panels,
            0.0,
            no_cut_datasets.front(),
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

    const auto model_par_prior{ MakeDefaultModelPrior() };
    test_data_factory::LocalScenario base_scenario;
    base_scenario.gaus_true = GaussianModel3D::FromVector(model_par_prior);
    base_scenario.sampling_entry_size = 50;
    base_scenario.data_error_sigma = 0.0;
    base_scenario.outlier_ratio = 0.0;
    base_scenario.replica_size = 10;
    base_scenario.requested_alpha_r_list = { request.alpha_r };
    base_scenario.alpha_training = true;
    base_scenario.local_fit_model = LocalGaussianFitModel::LogQuadratic;
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
            const auto bias{ rhbm_tester::RunBetaMDPDETest(test_input, request.job_count) };

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
        Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize()),
        Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize())
    } };
    outlier_prior_list.at(0)(0) = 1.50;
    outlier_prior_list.at(0)(1) = 0.50;
    outlier_prior_list.at(0)(2) = 0.10;
    outlier_prior_list.at(1)(0) = 1.00;
    outlier_prior_list.at(1)(1) = 1.00;
    outlier_prior_list.at(1)(2) = 0.10;

    const auto model_par_prior{ MakeDefaultModelPrior() };
    const auto model_par_sigma{ MakeDefaultModelSigma() };
    test_data_factory::GroupScenario base_scenario;
    base_scenario.member_size = 100;
    base_scenario.gaus_prior = model_par_prior;
    base_scenario.gaus_sigma = model_par_sigma;
    base_scenario.outlier_prior = model_par_prior;
    base_scenario.outlier_sigma = model_par_sigma;
    base_scenario.outlier_ratio = 0.0;
    base_scenario.replica_size = 100;
    base_scenario.requested_alpha_g_list = { request.alpha_g };
    base_scenario.alpha_training = true;
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
            mu_scenario.outlier_prior = outlier_prior;
            mu_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{
                test_data_factory::BuildGroupTestData(mu_scenario)
            };
            const auto bias{ rhbm_tester::RunMuMDPDETest(test_input, request.job_count) };

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
    test_data_factory::LocalScenario base_scenario;
    base_scenario.gaus_true = GaussianModel3D::FromVector(model_par_prior);
    base_scenario.sampling_entry_size = 50;
    base_scenario.data_error_sigma = 1.0;
    base_scenario.outlier_ratio = 0.0;
    base_scenario.replica_size = 10;
    base_scenario.requested_alpha_r_list = { request.alpha_r };
    base_scenario.alpha_training = true;
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
            const auto bias{ rhbm_tester::RunBetaMDPDETest(test_input, request.job_count) };

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
        Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize()),
        Eigen::VectorXd::Zero(GaussianModel3D::ParameterSize())
    } };
    outlier_prior_list.at(0)(0) = 1.50;
    outlier_prior_list.at(0)(1) = 0.50;
    outlier_prior_list.at(0)(2) = 0.10;
    outlier_prior_list.at(1)(0) = 1.00;
    outlier_prior_list.at(1)(1) = 1.00;
    outlier_prior_list.at(1)(2) = 0.10;

    const auto model_par_prior{ MakeDefaultModelPrior() };
    const auto model_par_sigma{ MakeDefaultModelSigma() };
    test_data_factory::GroupScenario base_scenario;
    base_scenario.member_size = 100;
    base_scenario.gaus_prior = model_par_prior;
    base_scenario.gaus_sigma = model_par_sigma;
    base_scenario.outlier_prior = model_par_prior;
    base_scenario.outlier_sigma = model_par_sigma;
    base_scenario.outlier_ratio = 0.0;
    base_scenario.replica_size = 100;
    base_scenario.requested_alpha_g_list = { request.alpha_g };
    base_scenario.alpha_training = true;
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
            mu_scenario.outlier_prior = outlier_prior;
            mu_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{
                test_data_factory::BuildGroupTestData(mu_scenario)
            };
            const auto bias{ rhbm_tester::RunMuMDPDETest(test_input, request.job_count) };

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

} // namespace rhbm_gem
