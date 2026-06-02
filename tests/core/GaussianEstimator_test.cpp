#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/data/object/AtomLocalPotentialView.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisEditor.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/ChemicalDataHelper.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>

namespace rg = rhbm_gem;
namespace ge = rhbm_gem::core;
namespace ls = rhbm_gem::linearization_service;

namespace
{

LocalPotentialSampleList MakeSampleEntries(double log_response_shift = 0.0)
{
    LocalPotentialSampleList sample_entries;
    sample_entries.reserve(6);
    for (int i = 0; i < 6; i++)
    {
        const auto distance{ static_cast<float>(0.1 * static_cast<double>(i)) };
        const auto response{
            static_cast<float>(std::exp(1.0 + log_response_shift - 0.5 * distance * distance))
        };
        sample_entries.emplace_back(LocalPotentialSample{ response, SamplingPoint{ distance } });
    }
    return sample_entries;
}

LocalPotentialSampleList MakeAlphaTrainingSampleEntries(
    std::size_t sample_size,
    double log_response_shift = 0.0)
{
    LocalPotentialSampleList sample_entries;
    sample_entries.reserve(sample_size);
    for (std::size_t i = 0; i < sample_size; i++)
    {
        const auto distance{ static_cast<float>(0.05 * static_cast<double>(i)) };
        const auto response{
            static_cast<float>(std::exp(1.0 + log_response_shift - 0.5 * distance * distance))
        };
        sample_entries.emplace_back(LocalPotentialSample{ response, SamplingPoint{ distance } });
    }
    return sample_entries;
}

LocalPotentialSampleList BuildShiftedSampleEntries(
    const LocalPotentialSampleList & sample_entries,
    double intercept)
{
    LocalPotentialSampleList shifted_sample_entries;
    shifted_sample_entries.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        shifted_sample_entries.emplace_back(
            LocalPotentialSample{
                static_cast<float>(static_cast<double>(sample.response) - intercept),
                sample.point
            }
        );
    }
    return shifted_sample_entries;
}

ge::FitOptions MakeOptions()
{
    ge::FitOptions options;
    options.thread_size = 1;
    return options;
}

rg::LocalGaussianFitModel MakeUnsupportedFitModel()
{
    return static_cast<rg::LocalGaussianFitModel>(999);
}

std::vector<LocalPotentialSampleList> MakeIdenticalSampleGroup(std::size_t member_size)
{
    std::vector<LocalPotentialSampleList> sample_group;
    sample_group.reserve(member_size);
    for (std::size_t i = 0; i < member_size; i++)
    {
        sample_group.emplace_back(MakeSampleEntries());
    }
    return sample_group;
}

std::vector<rg::LocalGaussianResult> EstimateMemberResults(
    const std::vector<LocalPotentialSampleList> & sample_group,
    const ge::FitOptions & options)
{
    std::vector<rg::LocalGaussianResult> member_results;
    member_results.reserve(sample_group.size());
    for (const auto & sample_entries : sample_group)
    {
        member_results.emplace_back(
            ge::EstimateLocalGaussianWithIntercept(sample_entries, 0.0, options));
    }
    return member_results;
}

std::vector<std::vector<rg::RHBMParameterVector>> BuildBetaGroupList(
    const std::vector<std::vector<rg::LocalGaussianResult>> & member_result_list)
{
    std::vector<std::vector<rg::RHBMParameterVector>> beta_group_list;
    beta_group_list.reserve(member_result_list.size());
    for (const auto & group_results : member_result_list)
    {
        std::vector<rg::RHBMParameterVector> beta_list;
        beta_list.reserve(group_results.size());
        for (const auto & member_result : group_results)
        {
            beta_list.emplace_back(
                ls::EncodeGaussianToParameterVector(member_result.mdpde.GetModel()));
        }
        beta_group_list.emplace_back(std::move(beta_list));
    }
    return beta_group_list;
}

std::unique_ptr<rg::ModelObject> MakeLocalFittingModel()
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    atom_list.reserve(2);
    for (int i = 0; i < 2; i++)
    {
        auto atom{ std::make_unique<rg::AtomObject>() };
        atom->SetSerialID(i + 1);
        atom->SetPosition(static_cast<float>(10 * i), 0.0f, 0.0f);
        atom_list.emplace_back(std::move(atom));
    }

    auto model{ std::make_unique<rg::ModelObject>(std::move(atom_list)) };
    model->SelectAllAtoms();

    auto analysis{ model->EditAnalysis() };
    for (std::size_t i = 0; i < model->GetSelectedAtoms().size(); i++)
    {
        auto local_editor{ analysis.EnsureAtomLocalPotential(*model->GetSelectedAtoms().at(i)) };
        local_editor.SetSamplingEntries(MakeSampleEntries(0.1 * static_cast<double>(i)));
        local_editor.SetAlphaR(0.2);
    }
    return model;
}

std::unique_ptr<rg::ModelObject> MakeLocalAlphaTrainingModel(
    const std::vector<LocalPotentialSampleList> & sample_entries_list,
    double initial_alpha_r)
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    atom_list.reserve(sample_entries_list.size());
    for (std::size_t i = 0; i < sample_entries_list.size(); i++)
    {
        auto atom{ std::make_unique<rg::AtomObject>() };
        atom->SetSerialID(static_cast<int>(i + 1));
        atom->SetComponentKey(7);
        atom->SetAtomKey(3);
        atom->SetPosition(static_cast<float>(10 * i), 0.0f, 0.0f);
        atom_list.emplace_back(std::move(atom));
    }

    auto model{ std::make_unique<rg::ModelObject>(std::move(atom_list)) };
    model->SelectAllAtoms();

    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    for (std::size_t i = 0; i < model->GetSelectedAtoms().size(); i++)
    {
        auto local_editor{ analysis.EnsureAtomLocalPotential(*model->GetSelectedAtoms().at(i)) };
        local_editor.SetSamplingEntries(sample_entries_list.at(i));
        local_editor.SetAlphaR(initial_alpha_r);
    }
    return model;
}

std::unique_ptr<rg::ModelObject> MakeGroupAlphaTrainingModel(std::size_t member_size)
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    atom_list.reserve(member_size);
    for (std::size_t i = 0; i < member_size; i++)
    {
        auto atom{ std::make_unique<rg::AtomObject>() };
        atom->SetSerialID(static_cast<int>(i + 1));
        atom->SetComponentKey(7);
        atom->SetAtomID("CA");
        atom->SetAtomKey(static_cast<AtomKey>(Spot::CA));
        atom->SetPosition(static_cast<float>(10 * i), 0.0f, 0.0f);
        atom_list.emplace_back(std::move(atom));
    }

    auto model{ std::make_unique<rg::ModelObject>(std::move(atom_list)) };
    model->SelectAllAtoms();

    auto analysis{ model->EditAnalysis() };
    analysis.RebuildAtomGroupsFromSelection();
    for (std::size_t i = 0; i < model->GetSelectedAtoms().size(); i++)
    {
        auto local_editor{ analysis.EnsureAtomLocalPotential(*model->GetSelectedAtoms().at(i)) };
        local_editor.SetSamplingEntries(MakeSampleEntries(0.02 * static_cast<double>(i)));
        local_editor.SetAlphaR(0.2);
    }
    return model;
}

std::vector<std::vector<rg::LocalGaussianResult>> CollectMainChainComponentMemberResults(
    const rg::ModelObject & model)
{
    const auto analysis_view{ model.GetAnalysisView() };
    const auto component_class_key{ ChemicalDataHelper::GetComponentAtomClassKey() };
    std::vector<std::vector<rg::LocalGaussianResult>> member_result_list;
    for (const auto group_key : analysis_view.CollectAtomGroupKeys(component_class_key))
    {
        const auto & atom_list{ analysis_view.GetAtomObjectList(group_key, component_class_key) };
        if (atom_list.size() < 10) continue;
        if (atom_list.front()->IsMainChainAtom() == false) continue;

        std::vector<rg::LocalGaussianResult> group_member_results;
        group_member_results.reserve(atom_list.size());
        for (const auto * atom : atom_list)
        {
            group_member_results.emplace_back(
                rg::AtomLocalPotentialView::RequireFor(*atom).GetGaussianResult());
        }
        member_result_list.emplace_back(std::move(group_member_results));
    }
    return member_result_list;
}

void ExpectAllAtomGroupsHaveAlphaG(const rg::ModelObject & model, double expected_alpha_g)
{
    const auto analysis_view{ model.GetAnalysisView() };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        for (const auto group_key : analysis_view.CollectAtomGroupKeys(class_key))
        {
            EXPECT_DOUBLE_EQ(expected_alpha_g, analysis_view.GetAtomAlphaG(group_key, class_key));
        }
    }
}

double GetFirstAtomGroupAlphaG(const rg::ModelObject & model)
{
    const auto analysis_view{ model.GetAnalysisView() };
    const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(0) };
    const auto group_keys{ analysis_view.CollectAtomGroupKeys(class_key) };
    if (group_keys.empty())
    {
        throw std::runtime_error("Test model has no atom groups.");
    }
    return analysis_view.GetAtomAlphaG(group_keys.front(), class_key);
}

void SetAllAtomGroupsAlphaG(rg::ModelObject & model, double alpha_g)
{
    auto analysis{ model.EditAnalysis() };
    const auto analysis_view{ model.GetAnalysisView() };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        for (const auto group_key : analysis_view.CollectAtomGroupKeys(class_key))
        {
            analysis.SetAtomGroupAlphaG(group_key, class_key, alpha_g);
        }
    }
}

void ExpectAllAtomGroupsHavePotentialFittingOutput(const rg::ModelObject & model)
{
    const auto analysis_view{ model.GetAnalysisView() };
    std::size_t group_count{ 0 };
    for (size_t i = 0; i < ChemicalDataHelper::GetGroupAtomClassCount(); i++)
    {
        const auto & class_key{ ChemicalDataHelper::GetGroupAtomClassKey(i) };
        for (const auto group_key : analysis_view.CollectAtomGroupKeys(class_key))
        {
            (void)analysis_view.GetAtomGroupMean(group_key, class_key);
            (void)analysis_view.GetAtomGroupMDPDE(group_key, class_key);
            (void)analysis_view.GetAtomGroupPriorWithUncertainty(group_key, class_key);

            const auto & atom_list{ analysis_view.GetAtomObjectList(group_key, class_key) };
            ASSERT_FALSE(atom_list.empty());
            for (const auto * atom : atom_list)
            {
                const auto local_view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
                EXPECT_TRUE(local_view.FindAnnotation(class_key).has_value());
                EXPECT_TRUE(local_view.GetGaussianResult().fit_result.has_value());
            }
            group_count++;
        }
    }
    EXPECT_GT(group_count, 0U);
}

} // namespace

TEST(GaussianEstimatorTest, SampleListAlphaRReturnsFiniteAlpha)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list{ MakeSampleEntries() };
    const auto alpha_r{
        ge::TrainAlphaR(sample_entries_list, options)
    };
    const rg::rhbm_trainer::RHBMTrainingOptions trainer_options;

    EXPECT_TRUE(std::isfinite(alpha_r));
    EXPECT_GE(alpha_r, trainer_options.alpha_min);
    EXPECT_LE(alpha_r, trainer_options.alpha_max);
}

TEST(GaussianEstimatorTest, AlphaRMatchesTrainingFunctionBestAlpha)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list{
        MakeSampleEntries(),
        MakeSampleEntries(0.2)
    };
    const std::vector<rg::RHBMMemberDataset> dataset_list{
        rg::rhbm_helper::BuildMemberDataset(
            sample_entries_list.at(0),
            options.distance_min,
            options.distance_max,
            options.local_fit_model),
        rg::rhbm_helper::BuildMemberDataset(
            sample_entries_list.at(1),
            options.distance_min,
            options.distance_max,
            options.local_fit_model)
    };
    rg::rhbm_trainer::RHBMTrainingOptions trainer_options;
    trainer_options.execution_options.quiet_mode = true;
    trainer_options.execution_options.thread_size = options.thread_size;

    const auto expected{
        rg::rhbm_trainer::CrossValidationAlphaR(
            dataset_list,
            trainer_options).best_alpha
    };
    const auto actual{
        ge::TrainAlphaR(sample_entries_list, options)
    };

    EXPECT_DOUBLE_EQ(actual, expected);
}

TEST(GaussianEstimatorTest, AlphaGMatchesTrainingFunctionBestAlpha)
{
    const auto options{ MakeOptions() };
    const auto sample_group{ MakeIdenticalSampleGroup(10) };
    const std::vector<std::vector<rg::LocalGaussianResult>> member_result_list{
        EstimateMemberResults(sample_group, options)
    };
    const auto beta_group_list{ BuildBetaGroupList(member_result_list) };
    rg::rhbm_trainer::RHBMTrainingOptions trainer_options;
    trainer_options.execution_options.quiet_mode = true;
    trainer_options.execution_options.thread_size = options.thread_size;

    const auto expected{
        rg::rhbm_trainer::CrossValidationAlphaG(
            beta_group_list,
            trainer_options).best_alpha
    };
    const auto actual{
        ge::TrainAlphaG(member_result_list, options)
    };

    EXPECT_DOUBLE_EQ(actual, expected);
}

TEST(GaussianEstimatorTest, RunLocalPotentialFittingUpdatesSelectedAtomLocalEntries)
{
    auto model{ MakeLocalFittingModel() };
    const auto options{ MakeOptions() };
    const auto expected_sample_size{ MakeSampleEntries().size() };

    ge::RunLocalPotentialFitting(*model, options);

    for (const auto * atom : model->GetSelectedAtoms())
    {
        const auto local_view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
        const auto & result{ local_view.GetGaussianResult() };

        EXPECT_DOUBLE_EQ(0.2, result.alpha_r);
        EXPECT_TRUE(result.fit_result.has_value());
        EXPECT_EQ(expected_sample_size, local_view.GetSamplingEntries(false).size());
    }
}

TEST(GaussianEstimatorTest, RunLocalPotentialFittingStopsAfterConvergence)
{
    auto model{ MakeLocalFittingModel() };
    const auto options{ MakeOptions() };
    const auto expected_sample_size{ MakeSampleEntries().size() };
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    ge::RunLocalPotentialFitting(*model, options);
    const auto output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_EQ(std::string::npos, output.find("(20/20)"));
    for (const auto * atom : model->GetSelectedAtoms())
    {
        const auto local_view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
        EXPECT_TRUE(local_view.GetGaussianResult().fit_result.has_value());
        EXPECT_EQ(expected_sample_size, local_view.GetSamplingEntries(false).size());
    }
}

TEST(GaussianEstimatorTest, RunLocalAlphaTrainingUpdatesComponentGroupAlphaR)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list{
        MakeAlphaTrainingSampleEntries(12),
        MakeAlphaTrainingSampleEntries(12, 0.2)
    };
    const auto expected_alpha_r{ ge::TrainAlphaR(sample_entries_list, options) };
    auto model{ MakeLocalAlphaTrainingModel(sample_entries_list, 0.2) };

    ge::RunLocalAlphaTraining(*model, options);

    for (const auto * atom : model->GetSelectedAtoms())
    {
        EXPECT_DOUBLE_EQ(
            expected_alpha_r,
            rg::AtomLocalPotentialView::RequireFor(*atom).GetAlphaR());
    }
}

TEST(GaussianEstimatorTest, RunLocalAlphaTrainingSkipsGroupsWithoutEnoughSamples)
{
    const auto options{ MakeOptions() };
    constexpr double initial_alpha_r{ 0.2 };
    const std::vector<LocalPotentialSampleList> sample_entries_list{
        MakeAlphaTrainingSampleEntries(6),
        MakeAlphaTrainingSampleEntries(6, 0.2)
    };
    auto model{ MakeLocalAlphaTrainingModel(sample_entries_list, initial_alpha_r) };

    ge::RunLocalAlphaTraining(*model, options);

    for (const auto * atom : model->GetSelectedAtoms())
    {
        EXPECT_DOUBLE_EQ(
            initial_alpha_r,
            rg::AtomLocalPotentialView::RequireFor(*atom).GetAlphaR());
    }
}

TEST(GaussianEstimatorTest, RunGroupAlphaTrainingUpdatesAllAtomGroupAlphaG)
{
    const auto options{ MakeOptions() };
    auto model{ MakeGroupAlphaTrainingModel(10) };

    ge::RunGroupAlphaTraining(*model, options);

    const auto member_result_list{ CollectMainChainComponentMemberResults(*model) };
    const rg::rhbm_trainer::RHBMTrainingOptions trainer_options;
    const auto alpha_g{ GetFirstAtomGroupAlphaG(*model) };
    ASSERT_FALSE(member_result_list.empty());
    EXPECT_GE(alpha_g, trainer_options.alpha_min);
    EXPECT_LE(alpha_g, trainer_options.alpha_max);
    ExpectAllAtomGroupsHaveAlphaG(*model, alpha_g);
}

TEST(GaussianEstimatorTest, RunGroupAlphaTrainingUsesFallbackAlphaGWithoutEnoughMembers)
{
    const auto options{ MakeOptions() };
    const std::vector<std::vector<rg::LocalGaussianResult>> empty_member_result_list;
    const auto expected_alpha_g{ ge::TrainAlphaG(empty_member_result_list, options) };
    auto model{ MakeGroupAlphaTrainingModel(6) };

    ge::RunGroupAlphaTraining(*model, options);

    ExpectAllAtomGroupsHaveAlphaG(*model, expected_alpha_g);
}

TEST(GaussianEstimatorTest, RunGroupPotentialFittingWritesGroupResultsAndMemberAnnotations)
{
    const auto options{ MakeOptions() };
    auto model{ MakeGroupAlphaTrainingModel(10) };
    SetAllAtomGroupsAlphaG(*model, 0.3);
    ge::RunLocalPotentialFitting(*model, options);

    ge::RunGroupPotentialFitting(*model, options);

    ExpectAllAtomGroupsHavePotentialFittingOutput(*model);
}

TEST(GaussianEstimatorTest, RunGroupPotentialFittingPreservesCorrectedLocalSamplingEntries)
{
    const auto options{ MakeOptions() };
    auto model{ MakeGroupAlphaTrainingModel(10) };
    SetAllAtomGroupsAlphaG(*model, 0.3);
    ge::RunLocalPotentialFitting(*model, options);
    const auto expected_sample_size{
        rg::AtomLocalPotentialView::RequireFor(*model->GetSelectedAtoms().front())
            .GetSamplingEntries(false)
            .size()
    };

    ge::RunGroupPotentialFitting(*model, options);

    for (const auto * atom : model->GetSelectedAtoms())
    {
        EXPECT_EQ(
            expected_sample_size,
            rg::AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false).size());
    }
    ExpectAllAtomGroupsHavePotentialFittingOutput(*model);
}

TEST(GaussianEstimatorTest, RejectsEmptyAlphaRTrainingInputs)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> empty_sample_entries_list;

    EXPECT_THROW(
        ge::TrainAlphaR(empty_sample_entries_list, options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, EmptyAlphaGTrainingInputReturnsFallbackAlpha)
{
    const auto options{ MakeOptions() };
    const std::vector<std::vector<rg::LocalGaussianResult>> empty_member_result_list;

    const auto alpha_g{
        ge::TrainAlphaG(empty_member_result_list, options)
    };
    const rg::rhbm_trainer::RHBMTrainingOptions trainer_options;

    EXPECT_DOUBLE_EQ(alpha_g, trainer_options.alpha_min);
}

TEST(GaussianEstimatorTest, EstimateLocalGaussianRejectsInvalidAlphaR)
{
    const auto options{ MakeOptions() };
    const auto sample_entries{ MakeSampleEntries() };

    EXPECT_THROW(
        ge::EstimateLocalGaussian(
            sample_entries, -std::numeric_limits<double>::min(), options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, EstimateLocalGaussianWithInterceptRejectsInvalidAlphaR)
{
    const auto options{ MakeOptions() };
    const auto sample_entries{ MakeSampleEntries() };

    EXPECT_THROW(
        ge::EstimateLocalGaussianWithIntercept(
            sample_entries, -std::numeric_limits<double>::min(), options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, EstimateGroupGaussianRejectsInvalidAlphaG)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list;
    const std::vector<rg::LocalGaussianResult> member_result_list;

    EXPECT_THROW(
        ge::EstimateGroupGaussian(
            sample_entries_list, member_result_list,
            std::numeric_limits<double>::quiet_NaN(), options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, TrainAlphaRRejectsInvalidFitRange)
{
    auto options{ MakeOptions() };
    options.distance_min = 1.0;
    options.distance_max = 0.5;
    const std::vector<LocalPotentialSampleList> sample_entries_list{ MakeSampleEntries() };

    EXPECT_THROW(
        ge::TrainAlphaR(sample_entries_list, options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, TrainAlphaRRejectsUnsupportedLocalFitModel)
{
    auto options{ MakeOptions() };
    options.local_fit_model = MakeUnsupportedFitModel();
    const std::vector<LocalPotentialSampleList> sample_entries_list{ MakeSampleEntries() };

    EXPECT_THROW(
        ge::TrainAlphaR(sample_entries_list, options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, TrainAlphaGRejectsUnsupportedLocalFitModel)
{
    auto options{ MakeOptions() };
    options.local_fit_model = MakeUnsupportedFitModel();
    const auto sample_group{ MakeIdenticalSampleGroup(10) };
    const std::vector<std::vector<rg::LocalGaussianResult>> member_result_list{
        EstimateMemberResults(sample_group, MakeOptions())
    };

    EXPECT_THROW(
        ge::TrainAlphaG(member_result_list, options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, EstimateLocalGaussianMatchesDirectHelperPath)
{
    const auto options{ MakeOptions() };
    const auto sample_entries{ MakeSampleEntries() };
    constexpr double alpha_r{ 0.2 };
    const auto actual{
        ge::EstimateLocalGaussian(sample_entries, alpha_r, options)
    };
    const auto dataset{
        rg::rhbm_helper::BuildMemberDataset(
            sample_entries,
            options.distance_min,
            options.distance_max,
            options.local_fit_model)
    };
    const auto expected_fit{
        rg::rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset)
    };

    const auto expected_ols{ ls::DecodeParameterVector(expected_fit.beta_ols) };
    const auto expected_mdpde{ ls::DecodeParameterVector(expected_fit.beta_mdpde) };
    EXPECT_NEAR(expected_ols.GetAmplitude(), actual.ols.GetModel().GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_ols.GetWidth(), actual.ols.GetModel().GetWidth(), 1e-12);
    EXPECT_DOUBLE_EQ(0.0, actual.ols.GetModel().GetIntercept());
    EXPECT_NEAR(expected_mdpde.GetAmplitude(), actual.mdpde.GetModel().GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_mdpde.GetWidth(), actual.mdpde.GetModel().GetWidth(), 1e-12);
    EXPECT_DOUBLE_EQ(0.0, actual.mdpde.GetModel().GetIntercept());
    EXPECT_DOUBLE_EQ(alpha_r, actual.alpha_r);
    ASSERT_TRUE(actual.fit_result.has_value());
    EXPECT_TRUE(actual.fit_result->beta_ols.isApprox(expected_fit.beta_ols, 1e-12));
    EXPECT_TRUE(actual.fit_result->beta_mdpde.isApprox(expected_fit.beta_mdpde, 1e-12));
}

TEST(GaussianEstimatorTest, EstimateLocalGaussianWithInterceptMatchesHelperPath)
{
    const auto options{ MakeOptions() };
    const auto sample_entries{ MakeSampleEntries() };
    constexpr double alpha_r{ 0.2 };
    const auto actual{
        ge::EstimateLocalGaussianWithIntercept(sample_entries, alpha_r, options)
    };
    const auto intercept{ actual.mdpde.GetModel().GetIntercept() };
    const auto shifted_sample_entries{ BuildShiftedSampleEntries(sample_entries, intercept) };
    const auto dataset{
        rg::rhbm_helper::BuildMemberDataset(
            shifted_sample_entries,
            options.distance_min,
            options.distance_max,
            options.local_fit_model)
    };
    const auto expected_fit{
        rg::rhbm_helper::EstimateBetaMDPDE(alpha_r, dataset)
    };

    const auto expected_ols{ ls::DecodeParameterVector(expected_fit.beta_ols).WithIntercept(intercept) };
    const auto expected_mdpde{ ls::DecodeParameterVector(expected_fit.beta_mdpde).WithIntercept(intercept) };
    EXPECT_NEAR(expected_ols.GetAmplitude(), actual.ols.GetModel().GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_ols.GetWidth(), actual.ols.GetModel().GetWidth(), 1e-12);
    EXPECT_NEAR(expected_ols.GetIntercept(), actual.ols.GetModel().GetIntercept(), 1e-12);
    EXPECT_NEAR(expected_mdpde.GetAmplitude(), actual.mdpde.GetModel().GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_mdpde.GetWidth(), actual.mdpde.GetModel().GetWidth(), 1e-12);
    EXPECT_NEAR(expected_mdpde.GetIntercept(), actual.mdpde.GetModel().GetIntercept(), 1e-12);
    EXPECT_DOUBLE_EQ(alpha_r, actual.alpha_r);
    ASSERT_TRUE(actual.fit_result.has_value());
    EXPECT_TRUE(actual.fit_result->beta_ols.isApprox(expected_fit.beta_ols, 1e-12));
    EXPECT_TRUE(actual.fit_result->beta_mdpde.isApprox(expected_fit.beta_mdpde, 1e-12));
}

TEST(GaussianEstimatorTest, EstimateLocalGaussianRejectsUnsupportedLocalFitModel)
{
    auto options{ MakeOptions() };
    options.local_fit_model = MakeUnsupportedFitModel();
    const auto sample_entries{ MakeSampleEntries() };

    EXPECT_THROW(
        ge::EstimateLocalGaussian(sample_entries, 0.0, options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, EstimateLocalGaussianWithInterceptRejectsUnsupportedLocalFitModel)
{
    auto options{ MakeOptions() };
    options.local_fit_model = MakeUnsupportedFitModel();
    const auto sample_entries{ MakeSampleEntries() };

    EXPECT_THROW(
        ge::EstimateLocalGaussianWithIntercept(sample_entries, 0.0, options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, EstimateGroupGaussianMatchesHelperPath)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list{
        MakeSampleEntries(),
        MakeSampleEntries(),
        MakeSampleEntries()
    };
    const std::vector<double> alpha_r_list{ 0.0, 0.1, 0.2 };
    std::vector<rg::LocalGaussianResult> member_result_list;
    constexpr double alpha_g{ 0.3 };

    std::vector<rg::RHBMMemberDataset> dataset_list;
    std::vector<rg::RHBMBetaEstimateResult> fit_result_list;
    dataset_list.reserve(sample_entries_list.size());
    fit_result_list.reserve(sample_entries_list.size());
    member_result_list.reserve(sample_entries_list.size());
    for (std::size_t i = 0; i < sample_entries_list.size(); i++)
    {
        auto member_result{
            ge::EstimateLocalGaussianWithIntercept(
                sample_entries_list.at(i), alpha_r_list.at(i), options)
        };
        ASSERT_TRUE(member_result.fit_result.has_value());
        const auto intercept{ member_result.mdpde.GetModel().GetIntercept() };
        const auto shifted_sample_entries{
            BuildShiftedSampleEntries(sample_entries_list.at(i), intercept)
        };
        auto dataset{
            rg::rhbm_helper::BuildMemberDataset(
                shifted_sample_entries,
                options.distance_min,
                options.distance_max,
                options.local_fit_model)
        };
        fit_result_list.emplace_back(*member_result.fit_result);
        member_result_list.emplace_back(std::move(member_result));
        dataset_list.emplace_back(std::move(dataset));
    }
    const auto expected_raw{
        rg::rhbm_helper::EstimateGroup(
            alpha_g,
            rg::rhbm_helper::BuildGroupInput(dataset_list, fit_result_list))
    };

    const auto actual{
        ge::EstimateGroupGaussian(
            sample_entries_list, member_result_list, alpha_g, options)
    };

    const auto expected_mean{ ls::DecodeParameterVector(expected_raw.mu_mean) };
    const auto expected_mdpde{ ls::DecodeParameterVector(expected_raw.mu_mdpde) };
    const auto expected_prior{
        ls::DecodeParameterVector(expected_raw.mu_prior, expected_raw.capital_lambda)
    };
    EXPECT_NEAR(expected_mean.GetAmplitude(), actual.mean.GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_mean.GetWidth(), actual.mean.GetWidth(), 1e-12);
    EXPECT_NEAR(expected_mdpde.GetAmplitude(), actual.mdpde.GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_mdpde.GetWidth(), actual.mdpde.GetWidth(), 1e-12);
    EXPECT_NEAR(expected_prior.GetModel().GetAmplitude(), actual.prior.GetModel().GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_prior.GetModel().GetWidth(), actual.prior.GetModel().GetWidth(), 1e-12);
    ASSERT_EQ(sample_entries_list.size(), actual.member_results.size());
}

TEST(GaussianEstimatorTest, EstimateGroupGaussianRejectsUnsupportedLocalFitModel)
{
    auto options{ MakeOptions() };
    options.local_fit_model = MakeUnsupportedFitModel();
    const std::vector<LocalPotentialSampleList> sample_entries_list{
        MakeSampleEntries(),
        MakeSampleEntries()
    };
    const std::vector<rg::LocalGaussianResult> member_result_list{
        ge::EstimateLocalGaussianWithIntercept(sample_entries_list.at(0), 0.0, MakeOptions()),
        ge::EstimateLocalGaussianWithIntercept(sample_entries_list.at(1), 0.0, MakeOptions())
    };

    EXPECT_THROW(
        ge::EstimateGroupGaussian(
            sample_entries_list, member_result_list, 0.0, options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, EstimateGroupGaussianRejectsUnsupportedMemberResults)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list{
        MakeSampleEntries(),
        MakeSampleEntries()
    };
    std::vector<rg::LocalGaussianResult> member_result_list{
        ge::EstimateLocalGaussianWithIntercept(sample_entries_list.at(0), 0.0, options),
        ge::EstimateLocalGaussianWithIntercept(sample_entries_list.at(1), 0.0, options)
    };
    member_result_list.front().fit_model = MakeUnsupportedFitModel();

    EXPECT_THROW(
        ge::EstimateGroupGaussian(
            sample_entries_list, member_result_list, 0.0, options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, TrainAlphaGRejectsUnsupportedMemberResults)
{
    const auto options{ MakeOptions() };
    const auto sample_group{ MakeIdenticalSampleGroup(10) };
    auto member_results{ EstimateMemberResults(sample_group, options) };
    member_results.front().fit_model = MakeUnsupportedFitModel();
    const std::vector<std::vector<rg::LocalGaussianResult>> member_result_list{
        member_results
    };

    EXPECT_THROW(
        ge::TrainAlphaG(member_result_list, options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, EstimateGroupGaussianRejectsInconsistentMemberCount)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list{
        MakeSampleEntries(),
        MakeSampleEntries()
    };
    const std::vector<rg::LocalGaussianResult> member_result_list{
        ge::EstimateLocalGaussianWithIntercept(sample_entries_list.front(), 0.0, options)
    };

    EXPECT_THROW(
        ge::EstimateGroupGaussian(
            sample_entries_list, member_result_list, 0.0, options),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, EstimateGroupGaussianRejectsMissingTransientFitState)
{
    const auto options{ MakeOptions() };
    const std::vector<LocalPotentialSampleList> sample_entries_list{
        MakeSampleEntries(),
        MakeSampleEntries()
    };
    const std::vector<rg::LocalGaussianResult> member_result_list(sample_entries_list.size());

    EXPECT_THROW(
        ge::EstimateGroupGaussian(
            sample_entries_list, member_result_list, 0.0, options),
        std::invalid_argument);
}
