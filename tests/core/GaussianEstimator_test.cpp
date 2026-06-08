#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
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
#include <rhbm_gem/utils/domain/SampleFilter.hpp>
#include <rhbm_gem/utils/hrl/LinearizationService.hpp>
#include <rhbm_gem/utils/hrl/RHBMHelper.hpp>
#include <rhbm_gem/utils/hrl/RHBMTrainer.hpp>
#include <rhbm_gem/utils/hrl/RHBMTypes.hpp>
#include <rhbm_gem/utils/domain/SamplingTypes.hpp>
#include <rhbm_gem/utils/math/ArrayHelper.hpp>

namespace rg = rhbm_gem;
namespace ge = rhbm_gem::core;
namespace ls = rhbm_gem::linearization_service;
namespace sf = rhbm_gem::sample_filter;

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

LocalPotentialSampleList MakeInterceptIterationSampleEntries(int seed)
{
    const auto phase{ 0.173 * static_cast<double>(seed) };
    const auto amplitude{
        0.5 + 1.4 * (0.5 + 0.5 * std::sin(0.31 * static_cast<double>(seed)))
    };
    const auto width{
        0.25 + 0.5 *
            (0.5 + 0.5 * std::sin(0.47 * static_cast<double>(seed) + 0.2))
    };
    const auto intercept{
        -0.08 + 0.5 *
            (0.5 + 0.5 * std::sin(0.23 * static_cast<double>(seed) + 0.7))
    };
    const rg::GaussianModel3D signal_model{ amplitude, width, 0.0 };

    LocalPotentialSampleList sample_entries;
    for (int shell = 0; shell <= 15; shell++)
    {
        const auto distance{ static_cast<float>(0.1 * static_cast<double>(shell)) };
        for (int sample_index = 0; sample_index < 9; sample_index++)
        {
            auto response{
                signal_model.SignalAtDistance(static_cast<double>(distance)) + intercept +
                0.12 * std::sin(
                    1.91 * static_cast<double>(shell) +
                    2.37 * static_cast<double>(sample_index) + phase)
            };
            if ((shell * 11 + sample_index * 7 + seed) % 13 == 0)
            {
                response += 0.75 * std::sin(
                    0.83 * static_cast<double>(shell) +
                    1.41 * static_cast<double>(sample_index) +
                    0.19 * static_cast<double>(seed));
            }
            sample_entries.emplace_back(
                LocalPotentialSample{
                    static_cast<float>(response),
                    SamplingPoint{ distance }
                });
        }
    }
    return sample_entries;
}

double EstimateInitialInterceptForTest(const LocalPotentialSampleList & sample_entries)
{
    float maximum_distance{ 0.0f };
    for (const auto & sample : sample_entries)
    {
        maximum_distance = std::max(maximum_distance, sample.point.distance);
    }

    std::vector<double> response_list;
    for (const auto & sample : sample_entries)
    {
        if (sample.point.distance == maximum_distance)
        {
            response_list.emplace_back(static_cast<double>(sample.response));
        }
    }
    return rg::array_helper::ComputeMedian(
        rg::array_helper::ComputeSmallestProportionValues(response_list, 0.10));
}

double EstimateResidualInterceptForTest(
    const LocalPotentialSampleList & sample_entries,
    const rg::RHBMBetaEstimateResult & fit_result)
{
    const auto signal_model{
        ls::DecodeParameterVector(fit_result.beta_mdpde)
    };
    const auto median_sample_entries{
        sf::BuildMedianResponseSampleEntriesByRadius(sample_entries)
    };
    std::vector<double> residual_list;
    for (const auto & sample : median_sample_entries)
    {
        if (sample.point.distance < 1.0f || sample.point.distance > 1.5f) continue;
        residual_list.emplace_back(
            static_cast<double>(static_cast<float>(
                static_cast<double>(sample.response) -
                signal_model.SignalAtDistance(
                    static_cast<double>(sample.point.distance)))));
    }
    return rg::array_helper::ComputeMedian(residual_list);
}

struct InterceptIterationSimulation
{
    bool converged{ false };
    bool cycle_detected{ false };
    bool max_iterations_reached{ false };
    std::size_t cycle_period{ 0 };
    int iteration_count{ 0 };
    double final_intercept{ 0.0 };
};

InterceptIterationSimulation SimulateInterceptIteration(
    const LocalPotentialSampleList & sample_entries,
    double alpha_r,
    const ge::FitOptions & options,
    double intercept_initial,
    double damping_weight)
{
    const rg::RHBMExecutionOptions execution_options;
    auto intercept{ intercept_initial };
    std::vector<double> intercept_history;
    intercept_history.reserve(64);
    auto best_intercept{ intercept };
    auto best_defect{ std::numeric_limits<double>::infinity() };
    for (int t = 0; t < execution_options.max_iterations; t++)
    {
        const auto result{
            ge::EstimateLocalGaussian(sample_entries, alpha_r, options, intercept)
        };
        if (!intercept_history.empty() &&
            std::abs(intercept - intercept_history.back()) < execution_options.tolerance)
        {
            return InterceptIterationSimulation{
                true, false, false, 0, t + 1, intercept
            };
        }
        for (std::size_t period = 2; period <= intercept_history.size(); period++)
        {
            const auto cycle_begin{ intercept_history.size() - period };
            if (std::abs(intercept - intercept_history[cycle_begin]) >=
                execution_options.tolerance)
            {
                continue;
            }
            double final_intercept{ 0.0 };
            for (std::size_t i = cycle_begin; i < intercept_history.size(); i++)
            {
                final_intercept += intercept_history[i];
            }
            final_intercept /= static_cast<double>(period);
            return InterceptIterationSimulation{
                false,
                true,
                false,
                period,
                t + 1,
                final_intercept
            };
        }

        const auto raw_intercept{
            EstimateResidualInterceptForTest(sample_entries, *result.fit_result)
        };
        const auto defect{ std::abs(raw_intercept - intercept) };
        if (defect < best_defect)
        {
            best_intercept = intercept;
            best_defect = defect;
        }
        if (t + 1 == execution_options.max_iterations)
        {
            return InterceptIterationSimulation{
                false,
                false,
                true,
                0,
                execution_options.max_iterations,
                best_intercept
            };
        }
        if (intercept_history.size() == 64)
        {
            intercept_history.erase(intercept_history.begin());
        }
        intercept_history.emplace_back(intercept);
        intercept += damping_weight * (raw_intercept - intercept);
    }
    throw std::logic_error("Intercept iteration simulation exited unexpectedly.");
}

void VerifyInterceptCycleRefit(
    const ge::FitOptions & options,
    int seed,
    std::size_t expected_period)
{
    constexpr double alpha_r{ 0.2 };
    const auto sample_entries{ MakeInterceptIterationSampleEntries(seed) };
    const auto intercept_initial{ EstimateInitialInterceptForTest(sample_entries) };
    const auto simulation{
        SimulateInterceptIteration(
            sample_entries, alpha_r, options, intercept_initial, 0.5)
    };
    ASSERT_TRUE(simulation.cycle_detected);
    ASSERT_EQ(expected_period, simulation.cycle_period);
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Debug);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    const auto actual{
        ge::EstimateLocalGaussianWithIntercept(
            sample_entries, alpha_r, options, intercept_initial)
    };
    const auto error_output{ testing::internal::GetCapturedStderr() };
    const auto output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(std::string::npos, output.find("Cycle detected"));
    EXPECT_NE(
        std::string::npos,
        output.find("period " + std::to_string(expected_period)));
    EXPECT_EQ(std::string::npos, error_output.find("Maximum iterations reached"));
    EXPECT_NEAR(
        simulation.final_intercept,
        actual.mdpde.GetModel().GetIntercept(),
        1e-12);

    const auto expected{
        ge::EstimateLocalGaussian(
            sample_entries, alpha_r, options, simulation.final_intercept)
    };
    EXPECT_NEAR(
        expected.mdpde.GetModel().GetAmplitude(),
        actual.mdpde.GetModel().GetAmplitude(),
        1e-12);
    EXPECT_NEAR(
        expected.mdpde.GetModel().GetWidth(),
        actual.mdpde.GetModel().GetWidth(),
        1e-12);
    ASSERT_TRUE(expected.fit_result.has_value());
    ASSERT_TRUE(actual.fit_result.has_value());
    EXPECT_TRUE(actual.fit_result->beta_mdpde.isApprox(
        expected.fit_result->beta_mdpde, 1e-12));
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

ge::FitOptions MakeOptions()
{
    ge::FitOptions options;
    options.thread_size = 1;
    return options;
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

double Distance(
    const std::array<float, 3> & lhs,
    const std::array<float, 3> & rhs)
{
    return static_cast<double>(rg::array_helper::ComputeNorm<float>(lhs, rhs));
}

LocalPotentialSampleList MakeCoupledSampleEntries(
    const rg::AtomObject & atom,
    const std::vector<rg::AtomObject *> & atom_list)
{
    const std::vector<std::array<float, 3>> direction_list{
        { 1.0f, 0.0f, 0.0f },
        { -1.0f, 0.0f, 0.0f },
        { 0.0f, 1.0f, 0.0f },
        { 0.0f, -1.0f, 0.0f },
        { 0.0f, 0.0f, 1.0f },
        { 0.0f, 0.0f, -1.0f }
    };
    const auto center{ atom.GetPosition() };
    constexpr double background{ 0.015 };
    constexpr double width{ 0.50 };
    LocalPotentialSampleList sample_entries;
    for (int shell = 0; shell <= 15; shell++)
    {
        const auto radius{ static_cast<float>(0.1 * static_cast<double>(shell)) };
        for (const auto & direction : direction_list)
        {
            SamplingPoint point{ radius };
            point.position = {
                center[0] + radius * direction[0],
                center[1] + radius * direction[1],
                center[2] + radius * direction[2]
            };

            double response{ background };
            for (const auto * source_atom : atom_list)
            {
                const rg::GaussianModel3D source_model{
                    1.0 + 0.02 * static_cast<double>(source_atom->GetSerialID()),
                    width,
                    0.0
                };
                response += source_model.SignalAtDistance(
                    Distance(point.position, source_atom->GetPosition()));
            }
            sample_entries.emplace_back(
                LocalPotentialSample{ static_cast<float>(response), point });
        }
    }
    return sample_entries;
}

std::unique_ptr<rg::ModelObject> MakeCoupledLocalFittingModel()
{
    std::vector<std::unique_ptr<rg::AtomObject>> atom_list;
    atom_list.reserve(10);
    for (int i = 0; i < 10; i++)
    {
        auto atom{ std::make_unique<rg::AtomObject>() };
        atom->SetSerialID(i + 1);
        atom->SetPosition(1.2f * static_cast<float>(i), 0.0f, 0.0f);
        atom_list.emplace_back(std::move(atom));
    }

    auto model{ std::make_unique<rg::ModelObject>(std::move(atom_list)) };
    model->SelectAllAtoms();

    auto analysis{ model->EditAnalysis() };
    const auto selected_atoms{ model->GetSelectedAtoms() };
    for (auto * atom : selected_atoms)
    {
        auto local_editor{ analysis.EnsureAtomLocalPotential(*atom) };
        local_editor.SetSamplingEntries(MakeCoupledSampleEntries(*atom, selected_atoms));
        local_editor.SetAlphaR(0.2);
    }
    return model;
}

using GaussianSnapshotForTest = std::unordered_map<const rg::AtomObject *, rg::GaussianModel3D>;

GaussianSnapshotForTest BuildGaussianSnapshotForTest(
    const std::vector<rg::AtomObject *> & atom_list)
{
    GaussianSnapshotForTest snapshot;
    snapshot.reserve(atom_list.size());
    for (const auto * atom : atom_list)
    {
        const auto local_view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
        snapshot.emplace(atom, local_view.GetGaussianResult().mdpde.GetModel());
    }
    return snapshot;
}

LocalPotentialSampleList UpdateSampleListWithFittedGaussianForTest(
    const rg::AtomObject & atom,
    const LocalPotentialSampleList & sample_entries,
    const GaussianSnapshotForTest & snapshot)
{
    LocalPotentialSampleList updated_list;
    updated_list.reserve(sample_entries.size());
    for (const auto & sample : sample_entries)
    {
        auto response_value{ sample.response };
        for (const auto * neighbor_atom : atom.FindNeighborAtoms())
        {
            const auto gaussian_iter{ snapshot.find(neighbor_atom) };
            if (gaussian_iter == snapshot.end()) continue;

            response_value -= static_cast<float>(
                gaussian_iter->second.SignalAtDistance(
                    Distance(sample.point.position, neighbor_atom->GetPosition())));
        }
        updated_list.emplace_back(LocalPotentialSample{ response_value, sample.point });
    }
    return updated_list;
}

double EstimateNextLocalFittingMaxStep(
    const rg::ModelObject & model,
    const std::vector<LocalPotentialSampleList> & original_sample_entries_list,
    const ge::FitOptions & options)
{
    const auto atom_list{ model.GetSelectedAtoms() };
    const auto snapshot{ BuildGaussianSnapshotForTest(atom_list) };
    double max_step{ 0.0 };
    for (std::size_t i = 0; i < atom_list.size(); i++)
    {
        const auto & atom{ *atom_list.at(i) };
        const auto local_view{ rg::AtomLocalPotentialView::RequireFor(atom) };
        const auto sample_entries{
            UpdateSampleListWithFittedGaussianForTest(
                atom,
                original_sample_entries_list.at(i),
                snapshot)
        };
        const auto current_model{ local_view.GetGaussianResult().mdpde.GetModel() };
        const auto next_result{
            ge::EstimateLocalGaussianWithIntercept(
                sample_entries,
                local_view.GetAlphaR(),
                options,
                current_model.GetIntercept())
        };
        const auto step{
            (next_result.mdpde.GetModel().ToVector() - current_model.ToVector()).norm()
        };
        if (step > max_step)
        {
            max_step = step;
        }
    }
    return max_step;
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
            std::vector<double> member_intercepts;
            member_intercepts.reserve(atom_list.size());
            for (const auto * atom : atom_list)
            {
                const auto local_view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
                const auto local_intercept{
                    local_view.GetGaussianResult().mdpde.GetModel().GetIntercept()
                };
                member_intercepts.emplace_back(local_intercept);
                const auto annotation{ local_view.FindAnnotation(class_key) };
                ASSERT_TRUE(annotation.has_value());
                EXPECT_NEAR(
                    local_intercept,
                    annotation->gaussian.GetModel().GetIntercept(),
                    1e-12);
                EXPECT_TRUE(local_view.GetGaussianResult().fit_result.has_value());
            }
            const auto group_intercept{ rg::array_helper::ComputeMedian(member_intercepts) };
            EXPECT_NEAR(
                group_intercept,
                analysis_view.GetAtomGroupMean(group_key, class_key).GetIntercept(),
                1e-12);
            EXPECT_NEAR(
                group_intercept,
                analysis_view.GetAtomGroupMDPDE(group_key, class_key).GetIntercept(),
                1e-12);
            EXPECT_NEAR(
                group_intercept,
                analysis_view.GetAtomGroupPrior(group_key, class_key).GetIntercept(),
                1e-12);
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
            options.distance_max),
        rg::rhbm_helper::BuildMemberDataset(
            sample_entries_list.at(1),
            options.distance_min,
            options.distance_max)
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

TEST(GaussianEstimatorTest, RunLocalPotentialFittingPreservesInitialIntercepts)
{
    auto model{ MakeLocalFittingModel() };
    const auto options{ MakeOptions() };
    const auto tolerance{ rg::RHBMExecutionOptions{}.tolerance };
    std::vector<double> expected_intercepts;
    expected_intercepts.reserve(model->GetSelectedAtoms().size());
    for (const auto * atom : model->GetSelectedAtoms())
    {
        const auto local_view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
        const auto expected_result{
            ge::EstimateLocalGaussianWithIntercept(
                local_view.GetSamplingEntries(),
                local_view.GetAlphaR(),
                options)
        };
        expected_intercepts.emplace_back(expected_result.mdpde.GetModel().GetIntercept());
    }

    ge::RunLocalPotentialFitting(*model, options);

    for (std::size_t i = 0; i < model->GetSelectedAtoms().size(); i++)
    {
        const auto local_view{
            rg::AtomLocalPotentialView::RequireFor(*model->GetSelectedAtoms().at(i))
        };
        const auto & result{ local_view.GetGaussianResult() };
        EXPECT_NEAR(
            expected_intercepts.at(i),
            result.ols.GetModel().GetIntercept(),
            tolerance);
        EXPECT_NEAR(
            expected_intercepts.at(i),
            result.mdpde.GetModel().GetIntercept(),
            tolerance);
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

    EXPECT_NE(std::string::npos, output.find("\rLocal fitting iteration "));
    EXPECT_EQ(
        std::string::npos,
        output.find("stable streak = "));
    EXPECT_NE(
        std::string::npos,
        output.find("\nConverged after "));
    EXPECT_NE(
        std::string::npos,
        output.find("max parameter change"));
    EXPECT_NE(
        std::string::npos,
        output.find("90th percentile parameter change"));
    EXPECT_EQ(std::string::npos, output.find("Reached maximum iteration size"));
    EXPECT_EQ(std::string::npos, output.find("(20/20)"));
    for (const auto * atom : model->GetSelectedAtoms())
    {
        const auto local_view{ rg::AtomLocalPotentialView::RequireFor(*atom) };
        EXPECT_TRUE(local_view.GetGaussianResult().fit_result.has_value());
        EXPECT_EQ(expected_sample_size, local_view.GetSamplingEntries(false).size());
    }
}

TEST(GaussianEstimatorTest, RunLocalPotentialFittingStopsAfterCoupledMaxConvergence)
{
    auto model{ MakeCoupledLocalFittingModel() };
    const auto options{ MakeOptions() };
    const auto previous_log_level{ Logger::GetLogLevel() };
    std::vector<LocalPotentialSampleList> original_sample_entries_list;
    original_sample_entries_list.reserve(model->GetSelectedAtoms().size());
    for (const auto * atom : model->GetSelectedAtoms())
    {
        original_sample_entries_list.emplace_back(
            rg::AtomLocalPotentialView::RequireFor(*atom).GetSamplingEntries(false));
    }

    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    ge::RunLocalPotentialFitting(*model, options);
    const auto error_output{ testing::internal::GetCapturedStderr() };
    const auto output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_NE(
        std::string::npos,
        output.find("max parameter change"));
    EXPECT_NE(
        std::string::npos,
        output.find("90th percentile parameter change"));
    EXPECT_EQ(std::string::npos, error_output.find("Reached maximum iteration size"));
    const auto next_max_step{
        EstimateNextLocalFittingMaxStep(
            *model,
            original_sample_entries_list,
            options)
    };
    EXPECT_LT(next_max_step * next_max_step, rg::RHBMExecutionOptions{}.tolerance);
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
            sample_entries, -std::numeric_limits<double>::min(), options, 0.0),
        std::invalid_argument);
}

TEST(GaussianEstimatorTest, EstimateLocalGaussianRejectsInvalidIntercept)
{
    const auto options{ MakeOptions() };
    const auto sample_entries{ MakeSampleEntries() };

    EXPECT_THROW(
        ge::EstimateLocalGaussian(
            sample_entries,
            0.2,
            options,
            std::numeric_limits<double>::quiet_NaN()),
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

TEST(GaussianEstimatorTest, EstimateLocalGaussianMatchesDirectHelperPath)
{
    const auto options{ MakeOptions() };
    const auto sample_entries{ MakeSampleEntries() };
    constexpr double alpha_r{ 0.2 };
    const auto actual{
        ge::EstimateLocalGaussian(sample_entries, alpha_r, options, 0.0)
    };
    const auto dataset{
        rg::rhbm_helper::BuildMemberDataset(
            sample_entries,
            options.distance_min,
            options.distance_max)
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

TEST(GaussianEstimatorTest, EstimateLocalGaussianAppliesProvidedIntercept)
{
    const auto options{ MakeOptions() };
    const auto sample_entries{ MakeSampleEntries() };
    constexpr double alpha_r{ 0.2 };
    constexpr double intercept{ 0.25 };
    const auto actual{
        ge::EstimateLocalGaussian(sample_entries, alpha_r, options, intercept)
    };
    const auto shifted_sample_entries{
        sf::BuildResponseShiftedSampleEntries(sample_entries, intercept)
    };
    const auto dataset{
        rg::rhbm_helper::BuildMemberDataset(
            shifted_sample_entries,
            options.distance_min,
            options.distance_max)
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

TEST(GaussianEstimatorTest, EstimateLocalGaussianWithInterceptMatchesHelperPath)
{
    const auto options{ MakeOptions() };
    const auto sample_entries{ MakeSampleEntries() };
    constexpr double alpha_r{ 0.2 };
    const auto actual{
        ge::EstimateLocalGaussianWithIntercept(sample_entries, alpha_r, options)
    };
    const auto intercept{ actual.mdpde.GetModel().GetIntercept() };
    const auto shifted_sample_entries{
        sf::BuildResponseShiftedSampleEntries(sample_entries, intercept)
    };
    const auto dataset{
        rg::rhbm_helper::BuildMemberDataset(
            shifted_sample_entries,
            options.distance_min,
            options.distance_max)
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

TEST(GaussianEstimatorTest, EstimateLocalGaussianWithInterceptUsesRadiusMedianResiduals)
{
    const auto options{ MakeOptions() };
    constexpr double alpha_r{ 0.0 };
    const rg::GaussianModel3D signal_model{ 1.0, 0.5, 0.0 };
    const LocalPotentialSampleList sample_entries{
        {
            static_cast<float>(signal_model.SignalAtDistance(1.0) + 0.20),
            SamplingPoint{ 1.0f }
        },
        {
            static_cast<float>(signal_model.SignalAtDistance(1.0) + 0.30),
            SamplingPoint{ 1.0f }
        },
        {
            static_cast<float>(signal_model.SignalAtDistance(1.0) + 50.0),
            SamplingPoint{ 1.0f }
        },
        {
            static_cast<float>(signal_model.SignalAtDistance(1.2) + 0.50),
            SamplingPoint{ 1.2f }
        },
        {
            static_cast<float>(signal_model.SignalAtDistance(1.2) + 0.60),
            SamplingPoint{ 1.2f }
        },
        {
            static_cast<float>(signal_model.SignalAtDistance(1.2) - 40.0),
            SamplingPoint{ 1.2f }
        },
        {
            static_cast<float>(signal_model.SignalAtDistance(1.4) + 0.80),
            SamplingPoint{ 1.4f }
        },
        {
            static_cast<float>(signal_model.SignalAtDistance(1.4) + 0.90),
            SamplingPoint{ 1.4f }
        },
        {
            static_cast<float>(signal_model.SignalAtDistance(1.4) + 60.0),
            SamplingPoint{ 1.4f }
        }
    };

    const auto actual{
        ge::EstimateLocalGaussianWithIntercept(sample_entries, alpha_r, options)
    };

    EXPECT_NEAR(actual.mdpde.GetModel().GetIntercept(), 0.60, 0.10);
}

TEST(GaussianEstimatorTest, EstimateLocalGaussianWithInterceptDampingConverges)
{
    const auto options{ MakeOptions() };
    constexpr double alpha_r{ 0.2 };
    const auto sample_entries{ MakeInterceptIterationSampleEntries(3) };
    const auto intercept_initial{ EstimateInitialInterceptForTest(sample_entries) };
    const auto undamped{
        SimulateInterceptIteration(
            sample_entries, alpha_r, options, intercept_initial, 1.0)
    };
    const auto damped{
        SimulateInterceptIteration(
            sample_entries, alpha_r, options, intercept_initial, 0.5)
    };
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    const auto actual{
        ge::EstimateLocalGaussianWithIntercept(
            sample_entries, alpha_r, options, intercept_initial)
    };
    const auto error_output{ testing::internal::GetCapturedStderr() };
    const auto output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_FALSE(undamped.converged);
    EXPECT_TRUE(undamped.cycle_detected);
    ASSERT_TRUE(damped.converged);
    EXPECT_LT(damped.iteration_count, rg::RHBMExecutionOptions{}.max_iterations);
    EXPECT_EQ(std::string::npos, output.find("Cycle detected"));
    EXPECT_EQ(std::string::npos, error_output.find("Maximum iterations reached"));
    EXPECT_NEAR(
        damped.final_intercept,
        actual.mdpde.GetModel().GetIntercept(),
        1e-12);

    const auto expected{
        ge::EstimateLocalGaussian(
            sample_entries,
            alpha_r,
            options,
            actual.mdpde.GetModel().GetIntercept())
    };
    EXPECT_NEAR(
        expected.mdpde.GetModel().GetAmplitude(),
        actual.mdpde.GetModel().GetAmplitude(),
        1e-12);
    EXPECT_NEAR(
        expected.mdpde.GetModel().GetWidth(),
        actual.mdpde.GetModel().GetWidth(),
        1e-12);
    ASSERT_TRUE(expected.fit_result.has_value());
    ASSERT_TRUE(actual.fit_result.has_value());
    EXPECT_TRUE(actual.fit_result->beta_mdpde.isApprox(
        expected.fit_result->beta_mdpde, 1e-12));
}

TEST(GaussianEstimatorTest, EstimateLocalGaussianWithInterceptRefitsAtTwoCycleMean)
{
    VerifyInterceptCycleRefit(MakeOptions(), 17, 2);
}

TEST(GaussianEstimatorTest, EstimateLocalGaussianWithInterceptRefitsAtPeriodThreeMean)
{
    VerifyInterceptCycleRefit(MakeOptions(), 289, 3);
}

TEST(GaussianEstimatorTest, EstimateLocalGaussianWithInterceptDetectsLongCycle)
{
    VerifyInterceptCycleRefit(MakeOptions(), 57, 45);
}

TEST(GaussianEstimatorTest, EstimateLocalGaussianWithInterceptRefitsAtBestCandidateOnLimit)
{
    const auto options{ MakeOptions() };
    constexpr double alpha_r{ 0.2 };
    const auto sample_entries{ MakeInterceptIterationSampleEntries(56) };
    const auto intercept_initial{ EstimateInitialInterceptForTest(sample_entries) };
    const auto simulation{
        SimulateInterceptIteration(
            sample_entries, alpha_r, options, intercept_initial, 0.5)
    };
    ASSERT_TRUE(simulation.max_iterations_reached);
    ASSERT_FALSE(simulation.cycle_detected);
    const auto previous_log_level{ Logger::GetLogLevel() };

    Logger::SetLogLevel(LogLevel::Info);
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    const auto actual{
        ge::EstimateLocalGaussianWithIntercept(
            sample_entries, alpha_r, options, intercept_initial)
    };
    const auto error_output{ testing::internal::GetCapturedStderr() };
    const auto output{ testing::internal::GetCapturedStdout() };
    Logger::SetLogLevel(previous_log_level);

    EXPECT_EQ(std::string::npos, output.find("Cycle detected"));
    EXPECT_NE(std::string::npos, error_output.find("Maximum iterations reached"));
    EXPECT_NE(std::string::npos, error_output.find("best fixed-point candidate"));
    EXPECT_NEAR(
        simulation.final_intercept,
        actual.mdpde.GetModel().GetIntercept(),
        1e-12);

    const auto expected{
        ge::EstimateLocalGaussian(
            sample_entries, alpha_r, options, simulation.final_intercept)
    };
    EXPECT_NEAR(
        expected.mdpde.GetModel().GetAmplitude(),
        actual.mdpde.GetModel().GetAmplitude(),
        1e-12);
    EXPECT_NEAR(
        expected.mdpde.GetModel().GetWidth(),
        actual.mdpde.GetModel().GetWidth(),
        1e-12);
    ASSERT_TRUE(expected.fit_result.has_value());
    ASSERT_TRUE(actual.fit_result.has_value());
    EXPECT_TRUE(actual.fit_result->beta_mdpde.isApprox(
        expected.fit_result->beta_mdpde, 1e-12));
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
            sf::BuildResponseShiftedSampleEntries(sample_entries_list.at(i), intercept)
        };
        auto dataset{
            rg::rhbm_helper::BuildMemberDataset(
                shifted_sample_entries,
                options.distance_min,
                options.distance_max)
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
    std::vector<double> member_intercepts;
    member_intercepts.reserve(member_result_list.size());
    for (const auto & member_result : member_result_list)
    {
        member_intercepts.emplace_back(member_result.mdpde.GetModel().GetIntercept());
    }
    const auto expected_group_intercept{ rg::array_helper::ComputeMedian(member_intercepts) };
    EXPECT_NEAR(expected_mean.GetAmplitude(), actual.mean.GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_mean.GetWidth(), actual.mean.GetWidth(), 1e-12);
    EXPECT_NEAR(expected_group_intercept, actual.mean.GetIntercept(), 1e-12);
    EXPECT_NEAR(expected_mdpde.GetAmplitude(), actual.mdpde.GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_mdpde.GetWidth(), actual.mdpde.GetWidth(), 1e-12);
    EXPECT_NEAR(expected_group_intercept, actual.mdpde.GetIntercept(), 1e-12);
    EXPECT_NEAR(expected_prior.GetModel().GetAmplitude(), actual.prior.GetModel().GetAmplitude(), 1e-12);
    EXPECT_NEAR(expected_prior.GetModel().GetWidth(), actual.prior.GetModel().GetWidth(), 1e-12);
    EXPECT_NEAR(
        expected_group_intercept,
        actual.prior.GetModel().GetIntercept(),
        1e-12);
    ASSERT_EQ(sample_entries_list.size(), actual.member_results.size());
    for (std::size_t i = 0; i < actual.member_results.size(); i++)
    {
        EXPECT_NEAR(
            member_result_list.at(i).mdpde.GetModel().GetIntercept(),
            actual.member_results.at(i).mdpde.GetModel().GetIntercept(),
            1e-12);
    }
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
