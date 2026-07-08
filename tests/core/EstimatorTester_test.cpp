#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/TestDataFactory.hpp>
#include <rhbm_gem/core/EstimatorTester.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

namespace {
namespace rt = rhbm_gem::core;
namespace tdf = rhbm_gem::core;
namespace rg = rhbm_gem;

tdf::GaussianParameterDistribution MakeDistribution(
    const rg::GaussianModel3D & mean,
    const rg::GaussianModel3DUncertainty & sigma = rg::GaussianModel3DUncertainty{ 0.05, 0.025, 0.01 })
{
    return tdf::GaussianParameterDistribution{ mean, sigma };
}

rt::LocalTestOptions MakeLocalOptions(
    double alpha_r,
    bool alpha_training)
{
    rt::LocalTestOptions options;
    options.requested_alpha_r = alpha_r;
    options.alpha_training = alpha_training;
    options.thread_size = 1;
    return options;
}

rt::GroupTestOptions MakeGroupOptions(
    double alpha_g,
    bool alpha_training)
{
    rt::GroupTestOptions options;
    options.requested_alpha_g = alpha_g;
    options.alpha_training = alpha_training;
    options.thread_size = 1;
    return options;
}

rt::FitOptions MakeSecondStageOptions()
{
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;
    return options;
}

void ExpectBiasStatisticSize(const rt::BiasStatistics & bias)
{
    EXPECT_EQ(bias.mean.size(), rg::GaussianModel3D::ParameterSize());
    EXPECT_EQ(bias.sigma.size(), rg::GaussianModel3D::ParameterSize());
}

double Distance(
    const std::array<float, 3> & lhs,
    const std::array<float, 3> & rhs)
{
    const auto dx{ static_cast<double>(lhs.at(0) - rhs.at(0)) };
    const auto dy{ static_cast<double>(lhs.at(1) - rhs.at(1)) };
    const auto dz{ static_cast<double>(lhs.at(2) - rhs.at(2)) };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double CalculateSelectedAtomResponseMeanSquaredError(const rg::ModelObject & model_object)
{
    double squared_error_sum{ 0.0 };
    std::size_t sample_count{ 0 };
    const auto & selected_atoms{ model_object.GetSelectedAtoms() };
    for (const auto * atom : selected_atoms)
    {
        const auto local_view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
        for (const auto & sample : local_view.GetSamplingEntries(false))
        {
            double fitted_response{ 0.0 };
            for (const auto * fitted_atom : selected_atoms)
            {
                const auto fitted_view{ rg::AtomLocalPotentialView::RequireFor(*fitted_atom) };
                fitted_response += fitted_view.GetEstimateMDPDE().ResponseAtDistance(
                    Distance(sample.point.position, fitted_atom->GetPosition()));
            }
            const auto residual{ static_cast<double>(sample.response) - fitted_response };
            squared_error_sum += residual * residual;
            sample_count++;
        }
    }
    return squared_error_sum / static_cast<double>(sample_count);
}

std::unique_ptr<rg::ModelObject> BuildSecondStageFallbackDiagnosticModel()
{
    ElectricPotential potential_model;
    potential_model.SetModelChoice(0);
    potential_model.SetBlurringWidth(0.5);
    auto input{
        tdf::BuildPotentialModelTestData(tdf::PotentialModelScenario{
            Spot::UNK,
            Element::OXYGEN,
            -0.1,
            rg::GaussianModel3D{ 8.0, 0.5, -0.1 },
            potential_model,
            0.0,
            1,
            42
        })
    };
    auto model{ std::move(input.replica_model_objects.front()) };

    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;
    rt::RunLocalAlphaTraining(*model, options);
    rt::RunFirstStageLocalFitting(*model, options);

    auto analysis{ model->EditAnalysis() };
    auto * atom{ model->GetSelectedAtoms().front() };
    auto sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false)
    };
    sampling_entries.front().response = std::numeric_limits<float>::quiet_NaN();
    analysis.EnsureAtomLocalPotential(*atom).SetSamplingEntries(std::move(sampling_entries));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildSecondStageScaleDiagnosticModel()
{
    ElectricPotential potential_model;
    potential_model.SetModelChoice(0);
    potential_model.SetBlurringWidth(0.5);
    auto input{
        tdf::BuildPotentialModelTestData(tdf::PotentialModelScenario{
            Spot::UNK,
            Element::OXYGEN,
            -0.1,
            rg::GaussianModel3D{ 8.0, 0.5, -0.1 },
            potential_model,
            0.0,
            1,
            42
        })
    };
    auto model{ std::move(input.replica_model_objects.front()) };
    const auto options{ MakeSecondStageOptions() };
    rt::RunLocalAlphaTraining(*model, options);
    rt::RunFirstStageLocalFitting(*model, options);
    return model;
}

void SetSelectedAtomEstimateModel(
    rg::ModelObject & model_object,
    const rg::GaussianModel3D & model)
{
    auto analysis{ model_object.EditAnalysis() };
    for (auto * atom : model_object.GetSelectedAtoms())
    {
        auto result{
            rg::AtomLocalPotentialView::RequireFor(*atom).GetGaussianResult()
        };
        result.ols = rg::GaussianModel3DWithUncertainty{
            model,
            result.ols.GetStandardDeviationModel()
        };
        result.mdpde = rg::GaussianModel3DWithUncertainty{
            model,
            result.mdpde.GetStandardDeviationModel()
        };
        analysis.EnsureAtomLocalPotential(*atom).SetGaussianResult(std::move(result));
    }
}

void RewriteSamplingResponsesFromSelectedAtomEstimates(rg::ModelObject & model_object)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto & selected_atoms{ model_object.GetSelectedAtoms() };
    for (auto * atom : selected_atoms)
    {
        auto sampling_entries{
            rg::AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false)
        };
        for (auto & sample : sampling_entries)
        {
            double response{ 0.0 };
            for (const auto * fitted_atom : selected_atoms)
            {
                const auto fitted_view{
                    rg::AtomLocalPotentialView::RequireFor(*fitted_atom)
                };
                response += fitted_view.GetEstimateMDPDE().ResponseAtDistance(
                    Distance(sample.point.position, fitted_atom->GetPosition()));
            }
            sample.response = static_cast<float>(response);
        }
        analysis.EnsureAtomLocalPotential(*atom).SetSamplingEntries(
            std::move(sampling_entries));
    }
}

std::unique_ptr<rg::AtomObject> MakeSecondStageAtom(
    int serial_id,
    Spot spot,
    Element element,
    const std::array<float, 3> & position)
{
    auto atom{ std::make_unique<rg::AtomObject>() };
    atom->SetSerialID(serial_id);
    atom->SetComponentKey(1);
    atom->SetAtomKey(static_cast<AtomKey>(spot));
    atom->SetElement(element);
    atom->SetSpot(spot);
    atom->SetPosition(position);
    return atom;
}

rg::LocalGaussianResult MakeSecondStageGaussianResult(const rg::GaussianModel3D & model)
{
    rg::LocalGaussianResult result;
    result.ols = rg::GaussianModel3DWithUncertainty{
        model,
        rg::GaussianModel3DUncertainty{}
    };
    result.mdpde = rg::GaussianModel3DWithUncertainty{
        model,
        rg::GaussianModel3DUncertainty{}
    };
    return result;
}

LocalPotentialSampleList BuildNearCollinearSamples(
    const rg::AtomObject & target_atom,
    const std::vector<rg::AtomObject *> & atom_list,
    const std::vector<rg::GaussianModel3D> & truth_model_list)
{
    const std::array<std::array<float, 3>, 6> direction_list{
        std::array<float, 3>{ 1.0F, 0.0F, 0.0F },
        std::array<float, 3>{ -1.0F, 0.0F, 0.0F },
        std::array<float, 3>{ 0.0F, 1.0F, 0.0F },
        std::array<float, 3>{ 0.0F, -1.0F, 0.0F },
        std::array<float, 3>{ 0.0F, 0.0F, 1.0F },
        std::array<float, 3>{ 0.0F, 0.0F, -1.0F }
    };
    const std::array<float, 4> radius_list{ 0.15F, 0.35F, 0.65F, 0.95F };

    LocalPotentialSampleList sample_list;
    const auto target_position{ target_atom.GetPosition() };
    sample_list.reserve(direction_list.size() * radius_list.size());
    for (const auto radius : radius_list)
    {
        for (const auto & direction : direction_list)
        {
            SamplingPoint point;
            point.distance = radius;
            for (std::size_t axis = 0; axis < point.position.size(); axis++)
            {
                point.position.at(axis) =
                    target_position.at(axis) + radius * direction.at(axis);
            }

            double response{ 0.0 };
            for (std::size_t i = 0; i < atom_list.size(); i++)
            {
                response += truth_model_list.at(i).ResponseAtDistance(
                    Distance(point.position, atom_list.at(i)->GetPosition()));
            }
            sample_list.emplace_back(LocalPotentialSample{
                static_cast<float>(response),
                point
            });
        }
    }
    return sample_list;
}

std::unique_ptr<rg::ModelObject> BuildNearCollinearSecondStageModel()
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    atom_list.emplace_back(MakeSecondStageAtom(
        1,
        Spot::C,
        Element::CARBON,
        std::array<float, 3>{ 0.0F, 0.0F, 0.0F }));
    atom_list.emplace_back(MakeSecondStageAtom(
        2,
        Spot::O,
        Element::OXYGEN,
        std::array<float, 3>{ 1.0e-4F, 0.0F, 0.0F }));

    auto model{ std::make_unique<rg::ModelObject>(std::move(atom_list)) };
    model->SelectAllAtoms();
    auto analysis{ model->EditAnalysis() };
    const std::vector<rg::GaussianModel3D> truth_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
        rg::GaussianModel3D{ 5.5, 0.55, -0.15 }
    };
    const rg::GaussianModel3D initial_model{ 5.75, 0.55, 0.0 };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    for (auto * atom : selected_atoms)
    {
        auto local_editor{ analysis.EnsureAtomLocalPotential(*atom) };
        local_editor.SetAlphaR(0.0);
        local_editor.SetGaussianResult(MakeSecondStageGaussianResult(initial_model));
        local_editor.SetSamplingEntries(
            BuildNearCollinearSamples(*atom, selected_atoms, truth_model_list));
    }
    return model;
}

std::unique_ptr<rg::ModelObject> BuildSeparatedSecondStageClusterModelWithSuspiciousLeftCluster()
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    atom_list.emplace_back(MakeSecondStageAtom(
        1,
        Spot::C,
        Element::CARBON,
        std::array<float, 3>{ 0.0F, 0.0F, 0.0F }));
    atom_list.emplace_back(MakeSecondStageAtom(
        2,
        Spot::O,
        Element::OXYGEN,
        std::array<float, 3>{ 1.0e-4F, 0.0F, 0.0F }));
    atom_list.emplace_back(MakeSecondStageAtom(
        3,
        Spot::N,
        Element::NITROGEN,
        std::array<float, 3>{ 10.0F, 0.0F, 0.0F }));
    atom_list.emplace_back(MakeSecondStageAtom(
        4,
        Spot::CA,
        Element::CARBON,
        std::array<float, 3>{ 10.0001F, 0.0F, 0.0F }));

    auto model{ std::make_unique<rg::ModelObject>(std::move(atom_list)) };
    model->SelectAllAtoms();
    auto analysis{ model->EditAnalysis() };
    const std::vector<rg::GaussianModel3D> truth_model_list{
        rg::GaussianModel3D{ 6.0, 0.55, 0.20 },
        rg::GaussianModel3D{ 5.5, 0.55, -0.15 },
        rg::GaussianModel3D{ 6.2, 0.55, 0.18 },
        rg::GaussianModel3D{ 5.7, 0.55, -0.12 }
    };
    const rg::GaussianModel3D initial_model{ 5.8, 0.55, 0.0 };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    for (auto * atom : selected_atoms)
    {
        auto local_editor{ analysis.EnsureAtomLocalPotential(*atom) };
        local_editor.SetAlphaR(0.0);
        local_editor.SetGaussianResult(MakeSecondStageGaussianResult(initial_model));
        local_editor.SetSamplingEntries(
            BuildNearCollinearSamples(*atom, selected_atoms, truth_model_list));
    }

    auto * suspicious_atom{ selected_atoms.at(0) };
    const auto suspicious_position{ suspicious_atom->GetPosition() };
    auto suspicious_sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*suspicious_atom).GetSamplingEntries(false)
    };
    suspicious_sampling_entries.resize(256);
    suspicious_sampling_entries.front().response = 0.0F;
    suspicious_sampling_entries.front().point.distance = 0.0F;
    suspicious_sampling_entries.front().point.position = suspicious_position;
    for (std::size_t i = 1; i < suspicious_sampling_entries.size(); i++)
    {
        auto & sample{ suspicious_sampling_entries.at(i) };
        const auto response_scale{
            0.5F + 0.5F * static_cast<float>(i) /
                static_cast<float>(suspicious_sampling_entries.size())
        };
        sample.response = std::numeric_limits<float>::max() * response_scale;
        sample.point.position = suspicious_position;
        sample.point.position.at(0) += 100.0F;
        sample.point.distance = 100.0F;
    }
    analysis.EnsureAtomLocalPotential(*suspicious_atom).SetSamplingEntries(
        std::move(suspicious_sampling_entries));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildSecondStageSuspiciousOffsetDiagnosticModel()
{
    ElectricPotential potential_model;
    potential_model.SetModelChoice(0);
    potential_model.SetBlurringWidth(0.5);
    auto input{
        tdf::BuildPotentialModelTestData(tdf::PotentialModelScenario{
            Spot::UNK,
            Element::OXYGEN,
            -0.1,
            rg::GaussianModel3D{ 8.0, 0.5, -0.1 },
            potential_model,
            0.0,
            1,
            42
        })
    };
    auto model{ std::move(input.replica_model_objects.front()) };

    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;
    rt::RunLocalAlphaTraining(*model, options);
    rt::RunFirstStageLocalFitting(*model, options);

    const auto & atom_list{ model->GetSelectedAtoms() };
    auto * target_atom{ atom_list.at(0) };
    auto target_position{ target_atom->GetPosition() };

    auto analysis{ model->EditAnalysis() };
    auto target_sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*target_atom).GetSamplingEntries(false)
    };
    target_sampling_entries.resize(256);
    target_sampling_entries.front().response = 0.0F;
    target_sampling_entries.front().point.distance = 0.0F;
    target_sampling_entries.front().point.position = target_position;
    for (std::size_t i = 1; i < target_sampling_entries.size(); i++)
    {
        auto & sample{ target_sampling_entries.at(i) };
        const auto response_scale{
            0.5F + 0.5F * static_cast<float>(i) /
                static_cast<float>(target_sampling_entries.size())
        };
        sample.response = std::numeric_limits<float>::max() * response_scale;
        sample.point.position = target_position;
        sample.point.position.at(0) += 100.0F;
        sample.point.distance = 100.0F;
    }
    analysis.EnsureAtomLocalPotential(*target_atom).SetSamplingEntries(
        std::move(target_sampling_entries));
    return model;
}

std::unique_ptr<rg::ModelObject> BuildSecondStageCoupledSuspiciousOffsetDiagnosticModel()
{
    auto model{ BuildNearCollinearSecondStageModel() };
    const auto & atom_list{ model->GetSelectedAtoms() };
    auto * target_atom{ atom_list.at(0) };
    const auto target_position{ target_atom->GetPosition() };

    auto analysis{ model->EditAnalysis() };
    auto target_sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*target_atom).GetSamplingEntries(false)
    };
    target_sampling_entries.resize(256);
    target_sampling_entries.front().response = 0.0F;
    target_sampling_entries.front().point.distance = 0.0F;
    target_sampling_entries.front().point.position = target_position;
    for (std::size_t i = 1; i < target_sampling_entries.size(); i++)
    {
        auto & sample{ target_sampling_entries.at(i) };
        const auto response_scale{
            0.5F + 0.5F * static_cast<float>(i) /
                static_cast<float>(target_sampling_entries.size())
        };
        sample.response = std::numeric_limits<float>::max() * response_scale;
        sample.point.position = target_position;
        sample.point.position.at(0) += 100.0F;
        sample.point.distance = 100.0F;
    }
    analysis.EnsureAtomLocalPotential(*target_atom).SetSamplingEntries(
        std::move(target_sampling_entries));
    return model;
}

void ExpectSelectedAtomEstimatesAreFinite(const rg::ModelObject & model_object)
{
    for (const auto * atom : model_object.GetSelectedAtoms())
    {
        const auto model{
            rg::AtomLocalPotentialView::RequireFor(*atom).GetEstimateMDPDE()
        };
        EXPECT_TRUE(std::isfinite(model.GetAmplitude()));
        EXPECT_TRUE(std::isfinite(model.GetWidth()));
        EXPECT_TRUE(std::isfinite(model.GetOffset()));
    }
}

void ExpectSelectedAtomUpdatedSamplesArePresent(const rg::ModelObject & model_object)
{
    for (const auto * atom : model_object.GetSelectedAtoms())
    {
        const auto local_view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
        EXPECT_FALSE(local_view.GetSamplingEntries(false, true).empty());
    }
}

} // namespace

TEST(EstimatorTesterTest, RunLocalEstimationTestPopulatesBiasOutputs)
{
    constexpr double alpha_r{ 0.5 };
    const auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };

    const auto bias{
        rt::RunLocalEstimationTest(test_input, MakeLocalOptions(alpha_r, true))
    };

    ExpectBiasStatisticSize(bias.ols);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha);
    ASSERT_TRUE(bias.mdpde.trained_alpha.has_value());
    ExpectBiasStatisticSize(bias.mdpde.trained_alpha.value());
    ASSERT_TRUE(bias.mdpde.trained_alpha_median.has_value());
    EXPECT_GE(bias.mdpde.trained_alpha_median.value(), 0.0);
}

TEST(EstimatorTesterTest, RunGroupEstimationTestPopulatesBiasOutputs)
{
    constexpr double alpha_g{ 0.2 };
    const auto test_input{
        tdf::BuildGroupTestData(tdf::GroupScenario{
            12,
            10,
            MakeDistribution(rg::GaussianModel3D{ 1.0, 0.5, 0.1 }),
            MakeDistribution(rg::GaussianModel3D{ 1.5, 0.5, 0.1 }),
            0.1,
            2,
            33
        })
    };

    const auto bias{
        rt::RunGroupEstimationTest(test_input, MakeGroupOptions(alpha_g, true))
    };

    ExpectBiasStatisticSize(bias.median);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha);
    ASSERT_TRUE(bias.mdpde.trained_alpha.has_value());
    ExpectBiasStatisticSize(bias.mdpde.trained_alpha.value());
    ASSERT_TRUE(bias.mdpde.trained_alpha_median.has_value());
    EXPECT_GE(bias.mdpde.trained_alpha_median.value(), 0.0);
}

TEST(EstimatorTesterTest, RunLocalEstimationTestSkipsTrainedAlphaWhenDisabled)
{
    constexpr double alpha_r{ 0.5 };
    const auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };

    const auto bias{
        rt::RunLocalEstimationTest(test_input, MakeLocalOptions(alpha_r, false))
    };

    ExpectBiasStatisticSize(bias.ols);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha);
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_median.has_value());
}

TEST(EstimatorTesterTest, RunGroupEstimationTestSkipsTrainedAlphaWhenDisabled)
{
    constexpr double alpha_g{ 0.2 };
    const auto test_input{
        tdf::BuildGroupTestData(tdf::GroupScenario{
            12,
            10,
            MakeDistribution(rg::GaussianModel3D{ 1.0, 0.5, 0.1 }),
            MakeDistribution(rg::GaussianModel3D{ 1.5, 0.5, 0.1 }),
            0.1,
            2,
            33
        })
    };

    const auto bias{
        rt::RunGroupEstimationTest(test_input, MakeGroupOptions(alpha_g, false))
    };

    ExpectBiasStatisticSize(bias.median);
    ExpectBiasStatisticSize(bias.mdpde.requested_alpha);
    EXPECT_FALSE(bias.mdpde.trained_alpha.has_value());
    EXPECT_FALSE(bias.mdpde.trained_alpha_median.has_value());
}

TEST(EstimatorTesterTest, RunLocalEstimationTestRejectsNonFiniteTruth)
{
    auto test_input{
        tdf::BuildLocalTestData(tdf::LocalScenario{
            rg::GaussianModel3D{ 1.0, 0.5, 0.0 },
            10,
            0.01,
            0.0,
            2,
            42
        })
    };
    test_input.gaus_true = rg::GaussianModel3D{
        std::numeric_limits<double>::quiet_NaN(),
        0.5,
        0.0
    };

    EXPECT_THROW(
        rt::RunLocalEstimationTest(test_input, MakeLocalOptions(0.5, true)),
        std::invalid_argument
    );
}

TEST(EstimatorTesterTest, RunPotentialFittingWorkflowStoresUpdatedSamplingEntries)
{
    ElectricPotential potential_model;
    potential_model.SetModelChoice(0);
    potential_model.SetBlurringWidth(0.5);
    auto input{
        tdf::BuildPotentialModelTestData(tdf::PotentialModelScenario{
            Spot::UNK,
            Element::OXYGEN,
            -0.1,
            rg::GaussianModel3D{ 8.0, 0.5, -0.1 },
            potential_model,
            0.0,
            1,
            42
        })
    };
    auto model{ std::move(input.replica_model_objects.front()) };
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;

    rt::RunPotentialFittingWorkflow(*model, options);

    ExpectSelectedAtomUpdatedSamplesArePresent(*model);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingHandlesNearCollinearAtoms)
{
    auto model{ BuildNearCollinearSecondStageModel() };
    const auto initial_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };
    const auto options{ MakeSecondStageOptions() };

    rt::RunSecondStageLocalFitting(*model, options);

    const auto fitted_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };
    const auto tolerance{ 1.0e-3 * std::max(initial_error, 1.0) };
    EXPECT_LE(fitted_error, initial_error + tolerance);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingLogsAndersonAccelerationMode)
{
    auto model{ BuildNearCollinearSecondStageModel() };
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_TRUE(
        output.find("acceleration = aa") != std::string::npos ||
        output.find("acceleration = damped-aa") != std::string::npos);
    const auto fallback_position{
        output.find("switching from Anderson acceleration to damped fixed-point fallback")
    };
    ASSERT_NE(fallback_position, std::string::npos);
    EXPECT_NE(
        output.rfind('\n', fallback_position),
        std::string::npos);
    EXPECT_NE(
        output.find("Objective backtracking rejected all attempts; backtracking retry"),
        std::string::npos);
    const auto fixed_point_position{
        output.find("acceleration = damped-fixed-point", fallback_position)
    };
    ASSERT_NE(fixed_point_position, std::string::npos);
    const auto next_anderson_position{
        output.find("acceleration = aa", fallback_position)
    };
    const auto next_damped_anderson_position{
        output.find("acceleration = damped-aa", fallback_position)
    };
    if (next_anderson_position != std::string::npos)
    {
        EXPECT_GT(next_anderson_position, fixed_point_position);
    }
    if (next_damped_anderson_position != std::string::npos)
    {
        EXPECT_GT(next_damped_anderson_position, fixed_point_position);
    }
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingKeepsAndersonAccelerationWithSeparatedRollbackCluster)
{
    auto model{ BuildSeparatedSecondStageClusterModelWithSuspiciousLeftCluster() };
    const auto initial_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string output{ testing::internal::GetCapturedStdout() };
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    const auto fitted_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };
    const auto tolerance{ 1.0e-3 * std::max(initial_error, 1.0) };
    EXPECT_LE(fitted_error, initial_error + tolerance);
    EXPECT_TRUE(
        output.find("acceleration = aa") != std::string::npos ||
        output.find("acceleration = damped-aa") != std::string::npos);
    EXPECT_NE(
        error_output.find("suspicious offset atom-events = "),
        std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingImprovesBadFiniteEntryScale)
{
    auto model{ BuildSecondStageScaleDiagnosticModel() };
    SetSelectedAtomEstimateModel(
        *model,
        rg::GaussianModel3D{ 1.0e4, 0.25, 1.0e3 });
    const auto entry_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    const auto fitted_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };
    ExpectSelectedAtomEstimatesAreFinite(*model);
    EXPECT_LT(fitted_error, entry_error);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingHandlesNearPerfectEntryScale)
{
    auto model{ BuildSecondStageScaleDiagnosticModel() };
    RewriteSamplingResponsesFromSelectedAtomEstimates(*model);
    const auto entry_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };

    rt::RunSecondStageLocalFitting(*model, MakeSecondStageOptions());

    const auto fitted_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };
    ExpectSelectedAtomEstimatesAreFinite(*model);
    EXPECT_LE(fitted_error, entry_error + 1.0e-8);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingLogsFallbackSummary)
{
    auto model{ BuildSecondStageFallbackDiagnosticModel() };
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = false;

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(
        error_output.find("Second-stage local fitting fallback summary:"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("joint offset fallback iterations = 1"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("refit fallback atom-events = 1"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("refit fallback distinct atoms = 1"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("suspicious offset atom-events = "),
        std::string::npos);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingQuietModeSuppressesFallbackSummary)
{
    auto model{ BuildSecondStageFallbackDiagnosticModel() };
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_EQ(
        error_output.find("Second-stage local fitting fallback summary:"),
        std::string::npos);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingRollsBackSuspiciousJointOffset)
{
    auto model{ BuildSecondStageSuspiciousOffsetDiagnosticModel() };
    auto * target_atom{ model->GetSelectedAtoms().front() };
    const auto previous_offset{
        rg::AtomLocalPotentialView::RequireFor(*target_atom).GetEstimateMDPDE().GetOffset()
    };

    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;

    rt::RunSecondStageLocalFitting(*model, options);

    const auto fitted_offset{
        rg::AtomLocalPotentialView::RequireFor(*target_atom).GetEstimateMDPDE().GetOffset()
    };
    EXPECT_NEAR(fitted_offset, previous_offset, 1.0e-12);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingRollsBackSuspiciousJointOffsetCluster)
{
    auto model{ BuildSecondStageCoupledSuspiciousOffsetDiagnosticModel() };
    const auto & atom_list{ model->GetSelectedAtoms() };
    const std::array<double, 2> previous_offset_list{
        rg::AtomLocalPotentialView::RequireFor(*atom_list.at(0)).GetEstimateMDPDE().GetOffset(),
        rg::AtomLocalPotentialView::RequireFor(*atom_list.at(1)).GetEstimateMDPDE().GetOffset()
    };

    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = false;

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    for (std::size_t i = 0; i < previous_offset_list.size(); i++)
    {
        const auto fitted_offset{
            rg::AtomLocalPotentialView::RequireFor(*atom_list.at(i)).GetEstimateMDPDE().GetOffset()
        };
        EXPECT_NEAR(fitted_offset, previous_offset_list.at(i), 1.0e-12);
    }
    EXPECT_NE(
        error_output.find("suspicious offset distinct atoms = 2"),
        std::string::npos);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}

TEST(EstimatorTesterTest, RunSecondStageLocalFittingLogsSuspiciousOffsetSummary)
{
    auto model{ BuildSecondStageSuspiciousOffsetDiagnosticModel() };
    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = false;

    const auto previous_log_level{ Logger::GetLogLevel() };
    Logger::SetLogLevel(LogLevel::Warning);
    testing::internal::CaptureStderr();
    rt::RunSecondStageLocalFitting(*model, options);
    const std::string error_output{ testing::internal::GetCapturedStderr() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(
        error_output.find("Second-stage local fitting fallback summary:"),
        std::string::npos);
    EXPECT_NE(
        error_output.find("suspicious offset atom-events = "),
        std::string::npos);
    EXPECT_NE(
        error_output.find("suspicious offset distinct atoms = "),
        std::string::npos);
}
