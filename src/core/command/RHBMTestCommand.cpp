#include "detail/CommandBase.hpp"
#include "detail/RHBMTestPlotting.hpp"

#include <rhbm_gem/core/EstimatorTester.hpp>
#include <rhbm_gem/core/TestDataFactory.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include <array>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>
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
    return GaussianModel3D{ 1.0, 0.5, -0.1 };
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

LocalTestOptions MakeLocalTestOptions(const RHBMTestRequest & request)
{
    LocalTestOptions options;
    options.requested_alpha_r = request.alpha_r;
    options.alpha_training = true;
    options.thread_size = request.job_count;
    options.quiet_mode = true;
    return options;
}

GroupTestOptions MakeGroupTestOptions(const RHBMTestRequest & request)
{
    GroupTestOptions options;
    options.requested_alpha_g = request.alpha_g;
    options.alpha_training = true;
    options.thread_size = request.job_count;
    options.quiet_mode = true;
    return options;
}

struct AtomicModelTestCase
{
    Spot spot;
    Element element;
    double charge{ 0.0 };
};

} // namespace

void RunSimulationTestOnBenchMark(const RHBMTestRequest & request)
{
    const std::vector<double> error_list{ 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 };
    const auto width_prior{ 0.5 };
    const GaussianModel3D gaus_truth{ 8.0, width_prior, -0.1 };

    FitOptions options;
    options.distance_min = request.fit_range_min;
    options.distance_max = request.fit_range_max;
    options.thread_size = request.job_count;
    options.quiet_mode = true;

    ElectricPotential potential_model;
    potential_model.SetModelChoice(0);
    potential_model.SetBlurringWidth(width_prior);

    PotentialModelScenario base_scenario;
    base_scenario.gaus_true = gaus_truth;
    base_scenario.potential_model = potential_model;
    base_scenario.spot = Spot::UNK;
    base_scenario.element = Element::OXYGEN;
    base_scenario.charge = -0.1;
    base_scenario.replica_size = 100;

    BiasPlotRequest plot_request;
    plot_request.output_name = "bias_from_benchmark.pdf";
    plot_request.flavor = BiasPlotFlavor::DataOutlier;
    plot_request.x_axis_mode = BiasXAxisMode::ErrorSigma;
    plot_request.panels.reserve(1);

    BiasPlotPanel panel;
    panel.label = "Independent Atom";
    panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Ols, error_list.size()));
    panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Mdpde, error_list.size()));

    for (const auto error_sigma : error_list)
    {
        auto scenario{ base_scenario };
        scenario.data_error_sigma = error_sigma;

        const auto input{ BuildPotentialModelTestData(scenario) };
        const auto bias{ RunAtomicModelLocalEstimationTest(input, request.alpha_r, options) };

        std::ostringstream stream;
        stream  << " error_sigma = " << std::setprecision(2) << std::fixed << error_sigma
                << ", OLS: " << std::setprecision(2)
                << bias.ols.mean(0) << " , "
                << bias.ols.mean(1) << " , "
                << bias.ols.mean(2)
                << " , MDPDE: "
                << bias.mdpde.mean(0) << " , "
                << bias.mdpde.mean(1) << " , "
                << bias.mdpde.mean(2);
        Logger::Log(LogLevel::Info, stream.str());

        AppendBiasCurvePoint(panel.curves.at(0), error_sigma, bias.ols);
        AppendBiasCurvePoint(panel.curves.at(1), error_sigma, bias.mdpde);
    }

    plot_request.panels.emplace_back(std::move(panel));
    rhbm_test_plotting::SaveDataOutlierBiasPlot(request, plot_request);
}

void RunSimulationTestOnAtomicModel(const RHBMTestRequest & request)
{
    const auto error_sigma{ 0.01 };
    const auto width_prior{ 0.5 };
    FitOptions options;
    options.distance_min = request.fit_range_min;
    options.distance_max = request.fit_range_max;
    options.thread_size = request.job_count;
    options.quiet_mode = true;

    ElectricPotential potential_model;
    potential_model.SetModelChoice(1);
    potential_model.SetBlurringWidth(width_prior);

    const std::array<AtomicModelTestCase, 3> benchmark_cases{{
        { Spot::UNK, Element::CARBON, 0.0 },
        { Spot::UNK, Element::NITROGEN, 0.0 },
        { Spot::UNK, Element::OXYGEN, 0.0 }
    }};

    const std::array<AtomicModelTestCase, 5> test_cases{{
        { Spot::UNK, Element::CARBON, 0.0 },
        { Spot::O, Element::OXYGEN, 0.0 },
        { Spot::N, Element::NITROGEN, 0.0 },
        { Spot::C, Element::CARBON, 0.0 },
        { Spot::CA, Element::CARBON, 0.0 }
    }};

    std::unordered_map<Element, GaussianModel3D> reference_gaus_by_element;
    reference_gaus_by_element.reserve(benchmark_cases.size());
    for (const auto & benchmark_case : benchmark_cases)
    {
        PotentialModelScenario scenario;
        scenario.potential_model = potential_model;
        scenario.data_error_sigma = 0.0;
        scenario.spot = benchmark_case.spot;
        scenario.element = benchmark_case.element;
        scenario.charge = benchmark_case.charge;
        scenario.replica_size = 1;

        const auto input{ BuildPotentialModelTestData(scenario) };
        const auto reference_gaus{ EstimateAtomicModelFirstStageMean(input, options) };
        reference_gaus_by_element.emplace(benchmark_case.element, reference_gaus);

        std::ostringstream stream;
        stream  << " Benchmark: " << static_cast<int>(benchmark_case.element)
                << " one-gaussian reference: " << std::setprecision(3) << std::fixed
                << reference_gaus.GetAmplitude() << " , "
                << reference_gaus.GetWidth() << " , "
                << reference_gaus.GetOffset();
        Logger::Log(LogLevel::Info, stream.str());
    }

    BiasPlotRequest plot_request;
    plot_request.output_name = "bias_from_neighbor_atom_atomic_model.pdf";
    plot_request.flavor = BiasPlotFlavor::NeighborType;
    plot_request.x_axis_mode = BiasXAxisMode::NeighborType;
    plot_request.panels.reserve(1);

    BiasPlotPanel panel;
    panel.label = "Neighbor Atom Type";
    panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Ols, test_cases.size()));
    panel.curves.emplace_back(MakeBiasCurve(BiasCurveKind::Mdpde, test_cases.size()));

    for (size_t i = 0; i < test_cases.size(); i++)
    {
        const auto & test_case{ test_cases.at(i) };
        PotentialModelScenario scenario;
        scenario.gaus_true = reference_gaus_by_element.at(test_case.element);
        scenario.potential_model = potential_model;
        scenario.data_error_sigma = error_sigma;
        scenario.spot = test_case.spot;
        scenario.element = test_case.element;
        scenario.charge = test_case.charge;
        scenario.replica_size = 10;

        const auto input{ BuildPotentialModelTestData(scenario) };
        const auto result_1{ RunAtomicModelFirstStageEstimationTest(input, options) };
        const auto result_2{ RunAtomicModelFullEstimationTest(input, options) };

        std::ostringstream stream;
        stream  << " Method 1: " << std::setprecision(2) << std::fixed
                << result_1.mean(0) << " , "
                << result_1.mean(1) << " , "
                << result_1.mean(2)
                << " , Method 2: "
                << result_2.mean(0) << " , "
                << result_2.mean(1) << " , "
                << result_2.mean(2);
        Logger::Log(LogLevel::Info, stream.str());

        const auto spot_axis_value{ static_cast<double>(i + 1) };
        AppendBiasCurvePoint(panel.curves.at(0), spot_axis_value, result_1);
        AppendBiasCurvePoint(panel.curves.at(1), spot_axis_value, result_2);
    }

    plot_request.panels.emplace_back(std::move(panel));
    rhbm_test_plotting::SaveDataOutlierBiasPlot(request, plot_request);
}

void RunSimulationTestOnDataOutlier(const RHBMTestRequest & request)
{
    const auto model_prior{ MakeDefaultModelPrior() };
    const auto local_options{ MakeLocalTestOptions(request) };
    LocalScenario base_scenario;
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
            auto local_scenario{ base_scenario };
            local_scenario.data_error_sigma = error_sigma;
            local_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{ BuildLocalTestData(local_scenario) };
            const auto bias{ RunLocalEstimationTest(test_input, local_options) };

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
    const std::vector<GaussianModel3D> outlier_prior_list{
        GaussianModel3D{ 1.50, 0.50, 0.10 },
        GaussianModel3D{ 1.00, 1.00, 0.10 }
    };

    const auto model_prior{ MakeDefaultModelPrior() };
    const auto model_sigma{ MakeDefaultModelSigma() };
    const auto group_options{ MakeGroupTestOptions(request) };
    GroupScenario base_scenario;
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
            auto group_scenario{ base_scenario };
            group_scenario.outlier_distribution.mean = outlier_prior;
            group_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{
                BuildGroupTestData(group_scenario)
            };
            const auto bias{ RunGroupEstimationTest(test_input, group_options) };

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
    const auto model_prior{ MakeDefaultModelPrior() };
    const auto local_options{ MakeLocalTestOptions(request) };
    LocalScenario base_scenario;
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
            auto local_scenario{ base_scenario };
            local_scenario.data_error_sigma = error_sigma;
            local_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{ BuildLocalTestData(local_scenario) };
            const auto bias{ RunLocalEstimationTest(test_input, local_options) };

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
    const std::vector<GaussianModel3D> outlier_prior_list{
        GaussianModel3D{ 1.50, 0.50, 0.10 },
        GaussianModel3D{ 1.00, 1.00, 0.10 }
    };

    const auto model_prior{ MakeDefaultModelPrior() };
    const auto model_sigma{ MakeDefaultModelSigma() };
    const auto group_options{ MakeGroupTestOptions(request) };
    GroupScenario base_scenario;
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
            auto group_scenario{ base_scenario };
            group_scenario.outlier_distribution.mean = outlier_prior;
            group_scenario.outlier_ratio = outlier_list.at(i);
            const auto test_input{
                BuildGroupTestData(group_scenario)
            };
            const auto bias{ RunGroupEstimationTest(test_input, group_options) };

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
    case TesterType::ATOMIC_MODEL:
        RunSimulationTestOnAtomicModel(request);
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
