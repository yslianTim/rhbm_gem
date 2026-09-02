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

void SetSelectedAtomPosteriorFromMdpde(rg::ModelObject & model_object)
{
    auto analysis{ model_object.EditAnalysis() };
    for (auto * atom : model_object.GetSelectedAtoms())
    {
        auto result{
            rg::AtomLocalPotentialView::For(*atom).GetGaussianResult(
                FittingStage::First)
        };
        result.posterior = result.mdpde;
        analysis.SetAtomLocalGaussianResult(
            FittingStage::Second,
            *atom,
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
    rt::RunLocalAlphaTraining(*model, options, FittingStage::First);
    rt::RunFixedOffsetLocalFitting(*model, options, FittingStage::First);
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
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 2;
    options.quiet_mode = true;
    rt::RunLocalAlphaTraining(*model, options, FittingStage::First);
    rt::RunFixedOffsetLocalFitting(*model, options, FittingStage::First);
    SetSelectedAtomPosteriorFromMdpde(*model);

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
        rg::AtomLocalPotentialView::For(
            *model->GetSelectedAtoms().front())
    };
    ASSERT_FALSE(initial_view.GetRawSamplingEntries(false).empty());
    ASSERT_TRUE(initial_view.GetPeelingSamplingEntries(false).empty());

    auto options{ MakeSecondStageOptions() };
    options.quiet_mode = true;
    command_test::ScopedTempDir temp_dir{ "local_fitting_result_table" };
    const auto csv_path{ temp_dir.path() / "local_fitting_result.csv" };
    options.result_csv_path = csv_path;
    rt::RunPotentialFittingWorkflow(*model, options);

    const auto fitted_view{
        rg::AtomLocalPotentialView::For(
            *model->GetSelectedAtoms().front())
    };
    EXPECT_FALSE(fitted_view.GetPeelingSamplingEntries(false).empty());
    const auto analysis_view{ model->GetAnalysisView() };
    for (const auto stage : {
             FittingStage::First,
             FittingStage::Second,
             FittingStage::Third })
    {
        const auto group_keys{ analysis_view.CollectAtomGroupKeys(stage) };
        ASSERT_EQ(group_keys.size(), 1u);
        const auto group_key{ group_keys.front() };
        EXPECT_EQ(analysis_view.GetAtomObjectList(stage, group_key).size(), 1u);
        EXPECT_TRUE(std::isfinite(
            analysis_view.GetAtomGroupPrior(stage, group_key).GetAmplitude()));
        EXPECT_TRUE(std::isfinite(analysis_view.GetAtomAlphaG(stage, group_key)));
        EXPECT_TRUE(fitted_view.GetGaussianResult(stage).posterior.has_value());
    }

    constexpr std::string_view csv_header{
        "serial id,residue,spot,neighbor count,peeling ratio,"
        "amplitude 1st,amplitude 2nd,amplitude 3rd,"
        "width 1st,width 2nd,width 3rd,"
        "offset 1st,offset 2nd,offset 3rd,"
        "amplitude rank 1st,amplitude rank 2nd,amplitude rank 3rd,"
        "width rank 1st,width rank 2nd,width rank 3rd,"
        "offset rank 1st,offset rank 2nd,offset rank 3rd"
    };
    ASSERT_TRUE(std::filesystem::exists(csv_path));
    const auto csv_content{ ReadTextFile(csv_path) };
    ASSERT_EQ(csv_content.find(csv_header), 0U);
    const auto row_begin{ csv_header.size() + 1 };
    const auto row_end{ csv_content.find('\n', row_begin) };
    const auto row{ SplitCsvLine(csv_content.substr(row_begin, row_end - row_begin)) };
    ASSERT_EQ(row.size(), 23U);
    const auto * atom{ model->GetSelectedAtoms().front() };
    EXPECT_EQ(std::stoi(row.at(0)), atom->GetSerialID());
    EXPECT_FALSE(row.at(1).empty());
    EXPECT_EQ(row.at(2), atom->GetAtomID());
    EXPECT_EQ(
        std::stoi(row.at(3)),
        fitted_view.GetNeighborCountForPeeling());
    const auto expected_peeling_ratio{ fitted_view.GetLocalFittingPeelingRatio(true) };
    ASSERT_TRUE(expected_peeling_ratio.has_value());
    EXPECT_NEAR(std::stod(row.at(4)), *expected_peeling_ratio, 0.0051);
    for (std::size_t column = 4; column < 14; column++)
    {
        EXPECT_TRUE(HasTwoFractionalDigits(row.at(column)));
        EXPECT_TRUE(std::isfinite(std::stod(row.at(column))));
    }
    const auto & final_model{
        fitted_view.GetEstimateMDPDE(FittingStage::Third)
    };
    const auto & first_model{
        fitted_view.GetEstimateMDPDE(FittingStage::First)
    };
    const auto & second_model{
        fitted_view.GetEstimateMDPDE(FittingStage::Second)
    };
    EXPECT_NEAR(std::stod(row.at(5)), first_model.GetAmplitude(), 0.0051);
    EXPECT_NEAR(std::stod(row.at(6)), second_model.GetAmplitude(), 0.0051);
    EXPECT_NEAR(std::stod(row.at(7)), final_model.GetAmplitude(), 0.0051);
    EXPECT_NEAR(std::stod(row.at(8)), first_model.GetWidth(), 0.0051);
    EXPECT_NEAR(std::stod(row.at(9)), second_model.GetWidth(), 0.0051);
    EXPECT_NEAR(std::stod(row.at(10)), final_model.GetWidth(), 0.0051);
    EXPECT_NEAR(std::stod(row.at(11)), first_model.GetOffset(), 0.0051);
    EXPECT_NEAR(std::stod(row.at(12)), second_model.GetOffset(), 0.0051);
    EXPECT_NEAR(std::stod(row.at(13)), final_model.GetOffset(), 0.0051);
    for (std::size_t column = 14; column < row.size(); ++column)
    {
        EXPECT_EQ(std::stoi(row.at(column)), 1);
    }

    {
        std::ofstream stale_output{ csv_path };
        stale_output << "stale content";
    }
    rt::RunPotentialFittingWorkflow(*model, options);
    const auto quiet_csv_content{ ReadTextFile(csv_path) };
    EXPECT_EQ(quiet_csv_content.find(csv_header), 0U);
    EXPECT_EQ(quiet_csv_content.find("stale content"), std::string::npos);
}

TEST(EstimatorTesterTest, LocalFittingResultRanksUseThreeNearestAtomsAcrossAllStages)
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
    command_test::ScopedTempDir temp_dir{ "local_fitting_result_ranks" };
    const auto csv_path{ temp_dir.path() / "local_fitting_result.csv" };
    options.result_csv_path = csv_path;
    rt::RunPotentialFittingWorkflow(*model, options);

    std::istringstream csv{ ReadTextFile(csv_path) };
    std::string line;
    ASSERT_TRUE(std::getline(csv, line));
    constexpr std::array stages{
        FittingStage::First,
        FittingStage::Second,
        FittingStage::Third
    };
    constexpr std::array<GaussianParameterGetter, 3> parameter_getters{
        &rg::GaussianModel3D::GetAmplitude,
        &rg::GaussianModel3D::GetWidth,
        &rg::GaussianModel3D::GetOffset
    };
    constexpr std::array<std::size_t, 3> rank_column_starts{ 14, 17, 20 };

    std::size_t row_count{ 0 };
    std::size_t verified_neighbor_set_count{ 0 };
    while (std::getline(csv, line))
    {
        const auto row{ SplitCsvLine(line) };
        ASSERT_EQ(row.size(), 23u);
        const auto serial_id{ std::stoi(row.at(0)) };
        const auto atom_iter{ std::find_if(
            selected_atoms.begin(),
            selected_atoms.end(),
            [serial_id](const rg::AtomObject * atom)
            {
                return atom->GetSerialID() == serial_id;
            })
        };
        ASSERT_NE(atom_iter, selected_atoms.end());
        for (std::size_t column = 14; column < row.size(); ++column)
        {
            const auto rank{ std::stoi(row.at(column)) };
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
            for (std::size_t stage = 0; stage < stages.size(); ++stage)
            {
                const auto expected_rank{ ComputeExpectedParameterRank(
                    **atom_iter,
                    comparison_atoms,
                    stages[stage],
                    parameter_getters[parameter])
                };
                const auto actual_rank{
                    std::stoi(row.at(rank_column_starts[parameter] + stage))
                };
                EXPECT_EQ(actual_rank, expected_rank);
            }
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
        rg::AtomLocalPotentialView::For(*target_atom)
            .GetEstimateMDPDE(FittingStage::Second).GetOffset()
    };

    rt::FitOptions options;
    options.distance_min = 0.0;
    options.distance_max = 1.0;
    options.thread_size = 1;
    options.quiet_mode = true;

    rt::RunSecondStageLocalFitting(*model, options);

    const auto fitted_offset{
        rg::AtomLocalPotentialView::For(*target_atom)
            .GetEstimateMDPDE(FittingStage::Second).GetOffset()
    };
    EXPECT_NEAR(fitted_offset, previous_offset, 1.0e-12);
    ExpectSelectedAtomEstimatesAreFinite(*model);
}
