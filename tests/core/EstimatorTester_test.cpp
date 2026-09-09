#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "core/command/detail/LocalFittingFeatures.hpp"
#include "core/detail/IterationProcess.hpp"
#include "core/detail/FittingRanges.hpp"
#include "core/detail/PreparedLocalGaussianFit.hpp"
#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/TestDataFactory.hpp>
#include <rhbm_gem/core/EstimatorTester.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>

namespace {
namespace rt = rhbm_gem::core;
namespace rt_detail = rhbm_gem::core::detail;
namespace tdf = rhbm_gem::core;
namespace rg = rhbm_gem;
using rg::FittingStage;

using GaussianParameterGetter = double (rg::GaussianModel3D::*)() const;

int ComputeExpectedParameterRank(
    const rg::AtomObject & atom,
    const std::vector<rg::AtomObject *> & comparison_atoms,
    FittingStage stage,
    GaussianParameterGetter parameter_getter)
{
    const auto & current_model{
        rg::AtomLocalPotentialView::For(atom).GetEstimateMDPDE(stage)
    };
    const auto current_value{ (current_model.*parameter_getter)() };
    return 1 + static_cast<int>(std::count_if(
        comparison_atoms.begin(),
        comparison_atoms.end(),
        [stage, parameter_getter, current_value](const rg::AtomObject * comparison_atom)
        {
            const auto & comparison_model{
                rg::AtomLocalPotentialView::For(*comparison_atom).GetEstimateMDPDE(stage)
            };
            return (comparison_model.*parameter_getter)() > current_value;
        }));
}

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
    const std::array<double, 3> & lhs,
    const std::array<double, 3> & rhs)
{
    const auto dx{ lhs.at(0) - rhs.at(0) };
    const auto dy{ lhs.at(1) - rhs.at(1) };
    const auto dz{ lhs.at(2) - rhs.at(2) };
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double CalculateSelectedAtomResponseMeanSquaredError(const rg::ModelObject & model_object)
{
    double squared_error_sum{ 0.0 };
    std::size_t sample_count{ 0 };
    const auto & selected_atoms{ model_object.GetSelectedAtoms() };
    for (const auto * atom : selected_atoms)
    {
        const auto local_view{ rg::AtomLocalPotentialView::For(*atom) };
        for (const auto & sample : local_view.GetRawSamplingEntries(false))
        {
            double fitted_response{ 0.0 };
            for (const auto * fitted_atom : selected_atoms)
            {
                const auto fitted_view{ rg::AtomLocalPotentialView::For(*fitted_atom) };
                fitted_response += fitted_view.GetEstimateMDPDE(
                    FittingStage::Second).ResponseAtDistance(
                    Distance(sample.point.position, fitted_atom->GetPosition()));
            }
            const auto residual{ sample.response - fitted_response };
            squared_error_sum += residual * residual;
            sample_count++;
        }
    }
    return squared_error_sum / static_cast<double>(sample_count);
}

void CopyFirstStageLocalResultsToSecond(rg::ModelObject & model_object)
{
    model_object.EditAnalysis().CopyLocalFittingStageResult(
        FittingStage::First,
        FittingStage::Second);
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
    rt::RunLocalAlphaTraining(*model, options, FittingStage::First);
    rt::RunFixedOffsetLocalFitting(*model, options, FittingStage::First);
    CopyFirstStageLocalResultsToSecond(*model);
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
            rg::AtomLocalPotentialView::For(*atom).GetGaussianResult(
                FittingStage::Second)
        };
        result.ols = rg::GaussianModel3DWithUncertainty{
            model,
            result.ols.GetStandardDeviationModel()
        };
        result.mdpde = rg::GaussianModel3DWithUncertainty{
            model,
            result.mdpde.GetStandardDeviationModel()
        };
        analysis.SetAtomLocalGaussianResult(
            FittingStage::Second,
            *atom,
            std::move(result));
    }
}

void RewriteSamplingResponsesFromSelectedAtomEstimates(rg::ModelObject & model_object)
{
    auto analysis{ model_object.EditAnalysis() };
    const auto & selected_atoms{ model_object.GetSelectedAtoms() };
    for (auto * atom : selected_atoms)
    {
        auto raw_sampling_entries{
            rg::AtomLocalPotentialView::For(*atom)
                .GetRawSamplingEntries(false)
        };
        for (auto & sample : raw_sampling_entries)
        {
            double response{ 0.0 };
            for (const auto * fitted_atom : selected_atoms)
            {
                const auto fitted_view{
                    rg::AtomLocalPotentialView::For(*fitted_atom)
                };
                response += fitted_view.GetEstimateMDPDE(
                    FittingStage::Second).ResponseAtDistance(
                    Distance(sample.point.position, fitted_atom->GetPosition()));
            }
            sample.response = response;
        }
        analysis.SetAtomLocalRawSamplingEntries(
            *atom, std::move(raw_sampling_entries));
    }
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
    options.thread_size = 2;
    options.quiet_mode = true;
    rt::RunLocalAlphaTraining(*model, options, FittingStage::First);
    rt::RunFixedOffsetLocalFitting(*model, options, FittingStage::First);
    CopyFirstStageLocalResultsToSecond(*model);

    const auto & atom_list{ model->GetSelectedAtoms() };
    auto * target_atom{ atom_list.at(0) };
    auto target_position{ target_atom->GetPosition() };

    auto analysis{ model->EditAnalysis() };
    auto target_raw_sampling_entries{
        rg::AtomLocalPotentialView::For(*target_atom)
            .GetRawSamplingEntries(false)
    };
    target_raw_sampling_entries.resize(256);
    target_raw_sampling_entries.front().response = 0.0;
    target_raw_sampling_entries.front().point.distance = 0.0;
    target_raw_sampling_entries.front().point.position = target_position;
    for (std::size_t i = 1; i < target_raw_sampling_entries.size(); i++)
    {
        auto & sample{ target_raw_sampling_entries.at(i) };
        const auto response_scale{
            0.5 + 0.5 * static_cast<double>(i) /
                static_cast<double>(target_raw_sampling_entries.size())
        };
        sample.response = std::numeric_limits<double>::max() * response_scale;
        sample.point.position = target_position;
        sample.point.position.at(0) += 100.0;
        sample.point.distance = 100.0;
    }
    analysis.SetAtomLocalRawSamplingEntries(
        *target_atom, std::move(target_raw_sampling_entries));
    return model;
}

void ExpectSelectedAtomEstimatesAreFinite(const rg::ModelObject & model_object)
{
    for (const auto * atom : model_object.GetSelectedAtoms())
    {
        const auto model{
            rg::AtomLocalPotentialView::For(*atom).GetEstimateMDPDE(
                FittingStage::Second)
        };
        EXPECT_TRUE(std::isfinite(model.GetAmplitude()));
        EXPECT_TRUE(std::isfinite(model.GetWidth()));
        EXPECT_TRUE(std::isfinite(model.GetOffset()));
    }
}

} // namespace

TEST(EstimatorTesterTest, PreparedLocalGaussianDatasetMatchesLegacyBuilder)
{
    LocalPotentialSampleList samples{
        { 0.6, SamplingPoint{ 0.0 } },
        { 0.8, SamplingPoint{ 0.25 } },
        { -1.0, SamplingPoint{ 0.5 } },
        { 0.7, SamplingPoint{ 1.0 } },
        { 0.9, SamplingPoint{ 1.25 } }
    };
    constexpr double range_min{ 0.25 };
    constexpr double range_max{ 1.0 };
    const rg::GaussianModel3D offset_model{ 1.0, 0.5, 0.1 };
    auto adjusted_samples{ samples };
    std::vector<double> response_list;
    response_list.reserve(samples.size());
    for (std::size_t i = 0; i < samples.size(); i++)
    {
        const auto evaluation{
            offset_model.EvaluateAtDistance(samples.at(i).point.distance)
        };
        adjusted_samples.at(i).response = samples.at(i).response -
            (evaluation.response - evaluation.signal);
        response_list.emplace_back(samples.at(i).response);
    }

    const auto legacy_dataset{
        rhbm_gem::rhbm_helper::BuildMemberDataset(
            adjusted_samples,
            range_min,
            range_max)
    };
    const rt_detail::PreparedLocalGaussianDesign design{
        samples,
        range_min,
        range_max
    };
    const auto prepared_dataset{
        design.BuildDataset(response_list, offset_model)
    };

    EXPECT_TRUE(prepared_dataset.X.isApprox(legacy_dataset.X, 0.0));
    EXPECT_TRUE(prepared_dataset.y.isApprox(legacy_dataset.y, 0.0));

    const rt_detail::PreparedLocalGaussianDesign signal_design{
        samples, 0.0, rt_detail::kSignalDistanceMax
    };
    const auto signal_dataset{ signal_design.BuildDataset(response_list, offset_model) };
    ASSERT_EQ(signal_dataset.X.rows(), 3);
    EXPECT_DOUBLE_EQ(signal_dataset.X(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(signal_dataset.X(2, 1), -0.5);
    const rt::FitOptions options;
    const auto fitted{ rt::EstimateLocalGaussian(samples, 0.0, options, offset_model) };
    const auto expected{ signal_design.Estimate(response_list, 0.0, options.thread_size, offset_model) };
    EXPECT_DOUBLE_EQ(fitted.mdpde.GetModel().GetAmplitude(), expected.mdpde.GetModel().GetAmplitude());
    EXPECT_DOUBLE_EQ(fitted.mdpde.GetModel().GetWidth(), expected.mdpde.GetModel().GetWidth());


    std::fill(response_list.begin(), response_list.end(), -1.0);
    for (auto & sample : adjusted_samples) sample.response = -1.0;
    const auto legacy_fallback{
        rhbm_gem::rhbm_helper::BuildMemberDataset(
            adjusted_samples,
            range_min,
            range_max)
    };
    const auto prepared_fallback{
        design.BuildDataset(response_list, offset_model)
    };
    EXPECT_TRUE(prepared_fallback.X.isApprox(legacy_fallback.X, 0.0));
    EXPECT_TRUE(prepared_fallback.y.isApprox(legacy_fallback.y, 0.0));
}

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

TEST(EstimatorTesterTest, GroupFittingUsesSecondLocalInputsWithoutChangingLocalStages)
{
    std::vector<std::unique_ptr<rg::AtomObject>> atoms;
    for (int i = 0; i < 3; ++i)
    {
        auto atom{ std::make_unique<rg::AtomObject>() };
        atom->SetSerialID(i + 1);
        atom->SetComponentKey(1);
        atom->SetAtomKey(static_cast<AtomKey>(Spot::CA));
        atom->SetElement(Element::CARBON);
        atom->SetSpot(Spot::CA);
        atom->SetPosition({ 10.0 * static_cast<double>(i), 0.0, 0.0 });
        atoms.emplace_back(std::move(atom));
    }
    rg::ModelObject model{ std::move(atoms) };
    model.SelectAllAtoms();
    auto analysis{ model.EditAnalysis() };
    analysis.InitializeFromSelection();
    const auto view{ model.GetAnalysisView() };
    const auto group_keys{ view.CollectAtomGroupKeys() };
    ASSERT_EQ(group_keys.size(), 1u);
    const auto group_key{ group_keys.front() };
    const auto & members{ view.GetAtomObjectList(group_key) };
    constexpr std::array stages{ FittingStage::First, FittingStage::Second };
    std::vector<std::array<rg::LocalGaussianResult, 2>> local_results;
    std::vector<rg::GroupGaussianMemberInput> expected_inputs;
    for (const auto * atom : members)
    {
        const double index{ static_cast<double>(atom->GetSerialID()) };
        const rg::GaussianModel3D final_model{ 1.0 + 0.2 * index, 0.5 + 0.04 * index, 0.1 * index };
        std::array<rg::LocalGaussianResult, 2> results;
        for (std::size_t i = 0; i < stages.size(); ++i)
        {
            results[i].alpha_r = 0.1 * static_cast<double>(i + 1);
            results[i].mdpde = rg::GaussianModel3DWithUncertainty{
                i == 1 ? final_model : rg::GaussianModel3D{ 20.0 + index, 0.9, 2.0 + static_cast<double>(i) },
                rg::GaussianModel3DUncertainty{}
            };
            results[i].ols = results[i].mdpde;
            analysis.SetAtomLocalGaussianResult(stages[i], *atom, results[i]);
        }
        LocalPotentialSampleList samples;
        for (int i = 0; i < 16; ++i)
        {
            const double distance{ static_cast<double>(i) / 15.0 };
            samples.emplace_back(LocalPotentialSample{
                final_model.ResponseAtDistance(distance) + 0.005 * std::sin(4.0 * distance + index),
                SamplingPoint{ distance, { distance, 0.0, 0.0 }, i % 3 != 0 }
            });
        }
        analysis.SetAtomLocalPeelingSamplingEntries(*atom, samples);
        auto raw_samples{ samples };
        for (auto & sample : raw_samples) sample.response += 5.0;
        analysis.SetAtomLocalRawSamplingEntries(*atom, std::move(raw_samples));
        expected_inputs.emplace_back(rg::GroupGaussianMemberInput{
            samples, results[1].alpha_r, final_model
        });
        local_results.emplace_back(std::move(results));
        EXPECT_FALSE(rg::AtomLocalPotentialView::For(*atom).GetGroupMemberResult().has_value());
    }

    constexpr double alpha_g{ 0.2 };
    analysis.InitializeGroupAlpha(alpha_g);
    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = true;
    const auto expected{ rt::EstimateGroupGaussian(expected_inputs, alpha_g, options) };
    rt::RunGroupPotentialFitting(model, options);

    EXPECT_EQ(view.CollectAtomGroupKeys(), group_keys);
    EXPECT_DOUBLE_EQ(view.GetAtomAlphaG(group_key), alpha_g);
    for (int parameter = 0; parameter < 3; ++parameter)
    {
        EXPECT_NEAR(view.GetAtomGroupMean(group_key).GetModelParameter(parameter),
            expected.mean.GetModelParameter(parameter), 1e-12);
        EXPECT_NEAR(view.GetAtomGroupMDPDE(group_key).GetModelParameter(parameter),
            expected.mdpde.GetModelParameter(parameter), 1e-12);
        EXPECT_NEAR(view.GetAtomGroupPrior(group_key).GetModelParameter(parameter),
            expected.prior.GetModelParameter(parameter), 1e-12);
        EXPECT_NEAR(view.GetAtomGroupPriorWithUncertainty(group_key).GetModelStandardDeviation(parameter),
            expected.prior.GetModelStandardDeviation(parameter), 1e-12);
    }
    for (std::size_t i = 0; i < members.size(); ++i)
    {
        const auto local_view{ rg::AtomLocalPotentialView::For(*members[i]) };
        const auto & member{ local_view.GetGroupMemberResult() };
        ASSERT_TRUE(member.has_value());
        EXPECT_EQ(member->is_outlier, expected.member_results[i].is_outlier);
        EXPECT_NEAR(member->statistical_distance, expected.member_results[i].statistical_distance, 1e-12);
        for (int parameter = 0; parameter < 3; ++parameter)
        {
            EXPECT_NEAR(member->posterior.GetModelParameter(parameter),
                expected.member_results[i].posterior.GetModelParameter(parameter), 1e-12);
            EXPECT_NEAR(member->posterior.GetModelStandardDeviation(parameter),
                expected.member_results[i].posterior.GetModelStandardDeviation(parameter), 1e-12);
            for (std::size_t stage = 0; stage < stages.size(); ++stage)
            {
                const auto & local_result{ local_view.GetGaussianResult(stages[stage]) };
                EXPECT_DOUBLE_EQ(local_result.alpha_r, local_results[i][stage].alpha_r);
                EXPECT_DOUBLE_EQ(local_result.ols.GetModelParameter(parameter),
                    local_results[i][stage].ols.GetModelParameter(parameter));
                EXPECT_DOUBLE_EQ(local_result.mdpde.GetModelParameter(parameter),
                    local_results[i][stage].mdpde.GetModelParameter(parameter));
            }
        }
    }
}

TEST(
    EstimatorTesterTest,
    RunPotentialFittingWorkflowProducesSingleGroupResultAfterLocalStages)
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
    ASSERT_EQ(model->GetSelectedAtoms().size(), 1u);
    const auto initial_view{
        rg::AtomLocalPotentialView::For(
            *model->GetSelectedAtoms().front())
    };
    ASSERT_FALSE(initial_view.GetRawSamplingEntries(false).empty());
    ASSERT_TRUE(initial_view.GetPeelingSamplingEntries(false).empty());

    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = true;
    auto analysis{ model->EditAnalysis() };
    analysis.InitializeLocalFittingSeedModels();
    rt::RunLocalAlphaTraining(*model, options, FittingStage::First);
    rt::RunFixedOffsetLocalFitting(*model, options, FittingStage::First);
    analysis.CopyLocalFittingStageResult(FittingStage::First, FittingStage::Second);
    rt_detail::RunSecondStageIterations(*model, options);
    const auto expected_local{ initial_view.GetGaussianResult(FittingStage::Second) };
    const auto expected_samples{ initial_view.GetPeelingSamplingEntries(false) };

    rt::RunPotentialFittingWorkflow(*model, options);

    const auto fitted_view{
        rg::AtomLocalPotentialView::For(
            *model->GetSelectedAtoms().front())
    };
    const auto actual_local{ fitted_view.GetGaussianResult(FittingStage::Second) };
    EXPECT_DOUBLE_EQ(actual_local.alpha_r, expected_local.alpha_r);
    for (int parameter = 0; parameter < 3; ++parameter)
    {
        EXPECT_DOUBLE_EQ(actual_local.ols.GetModelParameter(parameter),
            expected_local.ols.GetModelParameter(parameter));
        EXPECT_DOUBLE_EQ(actual_local.mdpde.GetModelParameter(parameter),
            expected_local.mdpde.GetModelParameter(parameter));
    }
    const auto actual_samples{ fitted_view.GetPeelingSamplingEntries(false) };
    ASSERT_EQ(actual_samples.size(), expected_samples.size());
    for (std::size_t i = 0; i < actual_samples.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(actual_samples[i].response, expected_samples[i].response);
        EXPECT_DOUBLE_EQ(actual_samples[i].point.distance, expected_samples[i].point.distance);
    }
    EXPECT_FALSE(fitted_view.GetPeelingSamplingEntries(false).empty());
    const auto analysis_view{ model->GetAnalysisView() };
    for (const auto stage : { FittingStage::First, FittingStage::Second })
    {
        const auto & local_model{ fitted_view.GetEstimateMDPDE(stage) };
        EXPECT_TRUE(std::isfinite(local_model.GetAmplitude()));
        EXPECT_GT(local_model.GetAmplitude(), 0.0);
        EXPECT_TRUE(std::isfinite(local_model.GetWidth()));
        EXPECT_GT(local_model.GetWidth(), 0.0);
        EXPECT_TRUE(std::isfinite(local_model.GetOffset()));
    }
    const auto group_keys{
        analysis_view.CollectAtomGroupKeys()
    };
    ASSERT_EQ(group_keys.size(), 1u);
    const auto group_key{ group_keys.front() };
    EXPECT_EQ(
        analysis_view.GetAtomObjectList(group_key).size(),
        1u);
    EXPECT_TRUE(std::isfinite(
        analysis_view.GetAtomGroupPrior(group_key).GetAmplitude()));
    EXPECT_GT(
        analysis_view.GetAtomGroupPrior(group_key).GetWidth(),
        0.0);
    EXPECT_TRUE(std::isfinite(
        analysis_view.GetAtomAlphaG(group_key)));
    EXPECT_TRUE(
        fitted_view.GetGroupMemberResult().has_value());
}

TEST(EstimatorTesterTest, LocalFittingResultRanksUseThreeNearestAtomsInSecondStage)
{
    ElectricPotential potential_model;
    potential_model.SetModelChoice(0);
    potential_model.SetBlurringWidth(0.5);
    auto input{
        tdf::BuildPotentialModelTestData(tdf::PotentialModelScenario{
            Spot::CA,
            Element::CARBON,
            -0.1,
            rg::GaussianModel3D{ 8.0, 0.5, -0.1 },
            potential_model,
            0.0,
            1,
            42
        })
    };
    auto model{ std::move(input.replica_model_objects.front()) };
    const auto & selected_atoms{ model->GetSelectedAtoms() };
    ASSERT_EQ(selected_atoms.size(), 5u);

    auto options{ MakeSecondStageOptions() };
    rt::RunPotentialFittingWorkflow(*model, options);

    const auto rows{ rt_detail::BuildLocalFittingFeatureRows(*model, true) };
    constexpr std::array<GaussianParameterGetter, 3> parameter_getters{
        &rg::GaussianModel3D::GetAmplitude,
        &rg::GaussianModel3D::GetWidth,
        &rg::GaussianModel3D::GetOffset
    };
    constexpr std::array<std::size_t, 3> rank_feature_indices{ 10, 11, 12 };

    std::size_t row_count{ 0 };
    std::size_t verified_neighbor_set_count{ 0 };
    for (const auto & row : rows)
    {
        const auto serial_id{ row.serial_id };
        const auto atom_iter{ std::find_if(
            selected_atoms.begin(),
            selected_atoms.end(),
            [serial_id](const rg::AtomObject * atom)
            {
                return atom->GetSerialID() == serial_id;
            })
        };
        ASSERT_NE(atom_iter, selected_atoms.end());
        EXPECT_FALSE(row.residue.empty());
        EXPECT_EQ(row.spot, (*atom_iter)->GetAtomID());
        for (std::size_t feature = 10; feature < row.features.size(); ++feature)
        {
            const auto rank{ row.features[feature] };
            EXPECT_GE(rank, 1);
            EXPECT_LE(rank, 4);
        }

        std::vector<rg::AtomObject *> comparison_atoms;
        for (auto * comparison_atom : selected_atoms)
        {
            if (comparison_atom != *atom_iter)
            {
                comparison_atoms.emplace_back(comparison_atom);
            }
        }
        std::sort(
            comparison_atoms.begin(),
            comparison_atoms.end(),
            [atom = *atom_iter](const rg::AtomObject * lhs, const rg::AtomObject * rhs)
            {
                const auto lhs_distance{ Distance(atom->GetPosition(), lhs->GetPosition()) };
                const auto rhs_distance{ Distance(atom->GetPosition(), rhs->GetPosition()) };
                if (lhs_distance != rhs_distance)
                {
                    return lhs_distance < rhs_distance;
                }
                return lhs->GetSerialID() < rhs->GetSerialID();
            });
        const auto third_neighbor_distance{
            Distance((*atom_iter)->GetPosition(), comparison_atoms[2]->GetPosition())
        };
        const auto fourth_neighbor_distance{
            Distance((*atom_iter)->GetPosition(), comparison_atoms[3]->GetPosition())
        };
        if (std::abs(third_neighbor_distance - fourth_neighbor_distance) < 1e-12)
        {
            ++row_count;
            continue;
        }
        comparison_atoms.resize(3);
        comparison_atoms.emplace_back(*atom_iter);

        for (std::size_t parameter = 0; parameter < parameter_getters.size(); ++parameter)
        {
            const auto expected_rank{ ComputeExpectedParameterRank(
                **atom_iter,
                comparison_atoms,
                FittingStage::Second,
                parameter_getters[parameter])
            };
            const auto actual_rank{ row.features[rank_feature_indices[parameter]] };
            EXPECT_EQ(actual_rank, expected_rank);
        }
        ++verified_neighbor_set_count;
        ++row_count;
    }
    EXPECT_EQ(row_count, selected_atoms.size());
    EXPECT_GT(verified_neighbor_set_count, 0u);
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

TEST(EstimatorTesterTest, RunSecondStageIterationsImprovesBadFiniteEntryScale)
{
    auto model{ BuildSecondStageScaleDiagnosticModel() };
    SetSelectedAtomEstimateModel(
        *model,
        rg::GaussianModel3D{ 1.0e4, 0.25, 1.0e3 });
    const auto entry_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };

    rt_detail::RunSecondStageIterations(*model, MakeSecondStageOptions());

    const auto fitted_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };
    ExpectSelectedAtomEstimatesAreFinite(*model);
    EXPECT_LT(fitted_error, entry_error);
}

TEST(EstimatorTesterTest, RunSecondStageIterationsHandlesNearPerfectEntryScale)
{
    auto model{ BuildSecondStageScaleDiagnosticModel() };
    RewriteSamplingResponsesFromSelectedAtomEstimates(*model);
    const auto entry_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };

    rt_detail::RunSecondStageIterations(*model, MakeSecondStageOptions());

    const auto fitted_error{
        CalculateSelectedAtomResponseMeanSquaredError(*model)
    };
    ExpectSelectedAtomEstimatesAreFinite(*model);
    EXPECT_LE(fitted_error, entry_error + 1.0e-8);
}

TEST(EstimatorTesterTest, RunSecondStageIterationsRollsBackSuspiciousJointOffset)
{
    auto model{ BuildSecondStageSuspiciousOffsetDiagnosticModel() };
    auto * target_atom{ model->GetSelectedAtoms().front() };
    const auto previous_offset{
        rg::AtomLocalPotentialView::For(*target_atom)
            .GetEstimateMDPDE(FittingStage::Second).GetOffset()
    };

    rt::FitOptions options;
    options.thread_size = 1;
    options.quiet_mode = true;

    rt_detail::RunSecondStageIterations(*model, options);

    const auto fitted_offset{
        rg::AtomLocalPotentialView::For(*target_atom)
            .GetEstimateMDPDE(FittingStage::Second).GetOffset()
    };
    EXPECT_NEAR(fitted_offset, previous_offset, 1.0e-12);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}
