#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "support/CommandTestHelpers.hpp"
#include <rhbm_gem/core/TestDataFactory.hpp>
#include <rhbm_gem/core/EstimatorTester.hpp>
#include "core/detail/GaussianEstimatorStages.hpp"
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>

namespace {
namespace rt = rhbm_gem::core;
namespace tdf = rhbm_gem::core;
namespace rg = rhbm_gem;

std::vector<std::string> SplitCsvLine(const std::string & line)
{
    std::vector<std::string> cells;
    std::istringstream stream{ line };
    std::string cell;
    while (std::getline(stream, cell, ','))
    {
        cells.emplace_back(std::move(cell));
    }
    return cells;
}

std::string ReadTextFile(const std::filesystem::path & path)
{
    std::ifstream input{ path };
    return std::string{
        std::istreambuf_iterator<char>{ input },
        std::istreambuf_iterator<char>{}
    };
}

bool HasTwoFractionalDigits(const std::string & value)
{
    const auto decimal_position{ value.find('.') };
    return decimal_position != std::string::npos
        && decimal_position + 3 == value.size()
        && value.at(decimal_position + 1) >= '0'
        && value.at(decimal_position + 1) <= '9'
        && value.at(decimal_position + 2) >= '0'
        && value.at(decimal_position + 2) <= '9';
}

TEST(EstimatorTesterTest, CalculateLocalFittingPeelingRatioUsesSignedResponseSums)
{
    const LocalPotentialSampleList raw_sampling_entries{
        LocalPotentialSample{ 2.0F, SamplingPoint{} },
        LocalPotentialSample{ 3.0F, SamplingPoint{} }
    };
    const LocalPotentialSampleList peeling_sampling_entries{
        LocalPotentialSample{ 1.0F, SamplingPoint{} },
        LocalPotentialSample{ 1.0F, SamplingPoint{} }
    };

    const auto ratio{ rt::detail::CalculateLocalFittingPeelingRatio(
        raw_sampling_entries,
        peeling_sampling_entries,
        true) };

    ASSERT_TRUE(ratio.has_value());
    EXPECT_DOUBLE_EQ(*ratio, 0.6);

    const auto negative_ratio{ rt::detail::CalculateLocalFittingPeelingRatio(
        LocalPotentialSampleList{ LocalPotentialSample{ 1.0F, SamplingPoint{} } },
        LocalPotentialSampleList{ LocalPotentialSample{ 2.0F, SamplingPoint{} } },
        true) };
    ASSERT_TRUE(negative_ratio.has_value());
    EXPECT_DOUBLE_EQ(*negative_ratio, -1.0);

    const auto above_one_ratio{ rt::detail::CalculateLocalFittingPeelingRatio(
        LocalPotentialSampleList{ LocalPotentialSample{ 1.0F, SamplingPoint{} } },
        LocalPotentialSampleList{ LocalPotentialSample{ -1.0F, SamplingPoint{} } },
        true) };
    ASSERT_TRUE(above_one_ratio.has_value());
    EXPECT_DOUBLE_EQ(*above_one_ratio, 2.0);
}

TEST(EstimatorTesterTest, CalculateLocalFittingPeelingRatioRejectsUnavailableInputs)
{
    const LocalPotentialSampleList valid_entries{
        LocalPotentialSample{ 1.0F, SamplingPoint{} }
    };
    const LocalPotentialSampleList zero_sum_entries{
        LocalPotentialSample{ 1.0F, SamplingPoint{} },
        LocalPotentialSample{ -1.0F, SamplingPoint{} }
    };
    const LocalPotentialSampleList non_finite_entries{
        LocalPotentialSample{
            std::numeric_limits<float>::quiet_NaN(),
            SamplingPoint{}
        }
    };

    EXPECT_FALSE(rt::detail::CalculateLocalFittingPeelingRatio(
        {}, valid_entries, true).has_value());
    EXPECT_FALSE(rt::detail::CalculateLocalFittingPeelingRatio(
        valid_entries, {}, true).has_value());
    EXPECT_FALSE(rt::detail::CalculateLocalFittingPeelingRatio(
        zero_sum_entries, valid_entries, true).has_value());
    EXPECT_FALSE(rt::detail::CalculateLocalFittingPeelingRatio(
        non_finite_entries, valid_entries, true).has_value());
    EXPECT_FALSE(rt::detail::CalculateLocalFittingPeelingRatio(
        valid_entries, non_finite_entries, true).has_value());
    EXPECT_FALSE(rt::detail::CalculateLocalFittingPeelingRatio(
        valid_entries, valid_entries, false).has_value());
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
        for (const auto & sample : local_view.GetRawSamplingEntries(false))
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

void SetSelectedAtomPosteriorFromMdpde(rg::ModelObject & model_object)
{
    auto analysis{ model_object.EditAnalysis() };
    for (auto * atom : model_object.GetSelectedAtoms())
    {
        auto result{
            rg::AtomLocalPotentialView::RequireFor(*atom).GetGaussianResult()
        };
        result.posterior = result.mdpde;
        analysis.EnsureAtomLocalPotential(*atom).SetGaussianResult(
            std::move(result));
    }
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
    rt::RunLocalAlphaTraining(*model, options, false);
    rt::RunFixedOffsetLocalFitting(*model, options, false);
    SetSelectedAtomPosteriorFromMdpde(*model);
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
        auto raw_sampling_entries{
            rg::AtomLocalPotentialView::RequireFor(*atom)
                .GetRawSamplingEntries(false)
        };
        for (auto & sample : raw_sampling_entries)
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
        analysis.EnsureAtomLocalPotential(*atom).SetRawSamplingEntries(
            std::move(raw_sampling_entries));
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
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;
    rt::RunLocalAlphaTraining(*model, options, false);
    rt::RunFixedOffsetLocalFitting(*model, options, false);
    SetSelectedAtomPosteriorFromMdpde(*model);

    const auto & atom_list{ model->GetSelectedAtoms() };
    auto * target_atom{ atom_list.at(0) };
    auto target_position{ target_atom->GetPosition() };

    auto analysis{ model->EditAnalysis() };
    auto target_raw_sampling_entries{
        rg::AtomLocalPotentialView::RequireFor(*target_atom)
            .GetRawSamplingEntries(false)
    };
    target_raw_sampling_entries.resize(256);
    target_raw_sampling_entries.front().response = 0.0F;
    target_raw_sampling_entries.front().point.distance = 0.0F;
    target_raw_sampling_entries.front().point.position = target_position;
    for (std::size_t i = 1; i < target_raw_sampling_entries.size(); i++)
    {
        auto & sample{ target_raw_sampling_entries.at(i) };
        const auto response_scale{
            0.5F + 0.5F * static_cast<float>(i) /
                static_cast<float>(target_raw_sampling_entries.size())
        };
        sample.response = std::numeric_limits<float>::max() * response_scale;
        sample.point.position = target_position;
        sample.point.position.at(0) += 100.0F;
        sample.point.distance = 100.0F;
    }
    analysis.EnsureAtomLocalPotential(*target_atom).SetRawSamplingEntries(
        std::move(target_raw_sampling_entries));
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

TEST(
    EstimatorTesterTest,
    RunPotentialFittingWorkflowBootstrapsGroupPosteriorFromAdjustedSamples)
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
        rg::AtomLocalPotentialView::RequireFor(
            *model->GetSelectedAtoms().front())
    };
    ASSERT_FALSE(initial_view.GetRawSamplingEntries(false).empty());
    ASSERT_TRUE(initial_view.GetPeelingSamplingEntries(false).empty());

    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = false;
    command_test::ScopedTempDir temp_dir{ "local_fitting_result_table" };
    const auto csv_path{ temp_dir.path() / "local_fitting_result.csv" };
    options.local_fitting_result_csv_path = csv_path;
    testing::internal::CaptureStdout();
    rt::RunPotentialFittingWorkflow(*model, options);
    const std::string out{ testing::internal::GetCapturedStdout() };

    const auto fitted_view{
        rg::AtomLocalPotentialView::RequireFor(
            *model->GetSelectedAtoms().front())
    };
    EXPECT_FALSE(fitted_view.GetPeelingSamplingEntries(false).empty());
    EXPECT_NE(
        out.find(
            "Selected second-stage initial seeds = 1, sources = "
            "group-posterior:1, group-prior:0, group-median:0, "
            "global-median:0."),
        std::string::npos);
    EXPECT_EQ(out.find("stop_reason=no-valid-seed"), std::string::npos);
    const auto count_occurrences = [&](const std::string & text)
    {
        std::size_t count{ 0 };
        for (std::size_t position = 0;
            (position = out.find(text, position)) != std::string::npos;
            position += text.size())
        {
            count++;
        }
        return count;
    };
    EXPECT_EQ(count_occurrences("Run atom group fitting."), 3U);

    constexpr std::string_view csv_header{
        "serial id,residue,spot,neighbor count,peeling ratio,"
        "amplitude 1st,amplitude 2nd,amplitude 3rd,"
        "width 1st,width 2nd,width 3rd,"
        "offset 1st,offset 2nd,offset 3rd"
    };
    EXPECT_EQ(out.find(csv_header), std::string::npos);
    ASSERT_TRUE(std::filesystem::exists(csv_path));
    const auto csv_content{ ReadTextFile(csv_path) };
    ASSERT_EQ(csv_content.find(csv_header), 0U);
    const auto row_begin{ csv_header.size() + 1 };
    const auto row_end{ csv_content.find('\n', row_begin) };
    const auto row{ SplitCsvLine(csv_content.substr(row_begin, row_end - row_begin)) };
    ASSERT_EQ(row.size(), 14U);
    const auto * atom{ model->GetSelectedAtoms().front() };
    EXPECT_EQ(std::stoi(row.at(0)), atom->GetSerialID());
    EXPECT_FALSE(row.at(1).empty());
    EXPECT_EQ(row.at(2), atom->GetAtomID());
    EXPECT_EQ(std::stoul(row.at(3)), 0U);
    const auto expected_peeling_ratio{
        rt::detail::CalculateLocalFittingPeelingRatio(
            fitted_view.GetRawSamplingEntries(false),
            fitted_view.GetPeelingSamplingEntries(false),
            true)
    };
    ASSERT_TRUE(expected_peeling_ratio.has_value());
    EXPECT_NEAR(std::stod(row.at(4)), *expected_peeling_ratio, 0.0051);
    for (std::size_t column = 4; column < row.size(); column++)
    {
        EXPECT_TRUE(HasTwoFractionalDigits(row.at(column)));
        EXPECT_TRUE(std::isfinite(std::stod(row.at(column))));
    }
    const auto & final_model{ fitted_view.GetEstimateMDPDE() };
    EXPECT_NEAR(std::stod(row.at(7)), final_model.GetAmplitude(), 0.0051);
    EXPECT_NEAR(std::stod(row.at(10)), final_model.GetWidth(), 0.0051);
    EXPECT_NEAR(std::stod(row.at(13)), final_model.GetOffset(), 0.0051);

    {
        std::ofstream stale_output{ csv_path };
        stale_output << "stale content";
    }
    options.quiet_mode = true;
    testing::internal::CaptureStdout();
    rt::RunPotentialFittingWorkflow(*model, options);
    const std::string quiet_out{ testing::internal::GetCapturedStdout() };
    EXPECT_EQ(quiet_out.find(csv_header), std::string::npos);
    const auto quiet_csv_content{ ReadTextFile(csv_path) };
    EXPECT_EQ(quiet_csv_content.find(csv_header), 0U);
    EXPECT_EQ(quiet_csv_content.find("stale content"), std::string::npos);
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
