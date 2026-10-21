#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/CommandTypes.hpp>
#include <rhbm_gem/core/JointComponentEstimator.hpp>
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include "core/command/detail/MapSimulation.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include "core/detail/JointCliProgressReporter.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "support/JointSyntheticWorkload.hpp"

namespace {
namespace core = rhbm_gem::core;
namespace joint = core::joint_component;

core::JointProblemInput MakeInput()
{
    core::JointProblemInput input;
    input.atom_ids = {"a", "b"};
    input.support.resize(2);
    for (int atom = 0; atom < 2; ++atom)
    {
        for (int k = 0; k < 40; ++k)
        {
            const double square = .003 * k * k;
            const auto basis = joint::EvaluateKernel(square, .5, 2.5);
            input.support[static_cast<std::size_t>(atom)].push_back({input.observations.size(), square});
            input.row_ids.push_back(std::to_string(input.observations.size()));
            input.observations.push_back(2 * basis.gaussian + .2 * basis.charge);
        }
    }
    return input;
}

core::JointProblemInput MakeObservableInput()
{
    core::JointProblemInput input;
    input.atom_ids = {"1"};
    input.support.resize(1);
    input.selection_domain.emplace();
    input.selection_domain->target_indices = {0};
    for (int k = 0; k < 40; ++k)
    {
        const double square = .003 * k * k;
        const auto basis = joint::EvaluateKernel(square, .5, 2.5);
        input.support[0].push_back({input.observations.size(), square});
        input.row_ids.push_back(std::to_string(input.observations.size()));
        input.observations.push_back(2 * basis.gaussian + .2 * basis.charge);
    }
    return input;
}

void ExpectSameChecks(const std::vector<rhbm_gem::JointCheck> & left,
    const std::vector<rhbm_gem::JointCheck> & right)
{
    ASSERT_EQ(left.size(), right.size());
    for (std::size_t k = 0; k < left.size(); ++k)
    {
        EXPECT_EQ(left[k].name, right[k].name);
        EXPECT_EQ(left[k].status, right[k].status);
        EXPECT_EQ(left[k].scope, right[k].scope);
        EXPECT_EQ(left[k].value, right[k].value);
        EXPECT_EQ(left[k].threshold, right[k].threshold);
        EXPECT_EQ(left[k].reason, right[k].reason);
    }
}

void ExpectSameState(const rhbm_gem::JointState & left, const rhbm_gem::JointState & right)
{
    EXPECT_EQ(left.ac, right.ac);
    EXPECT_EQ(left.b, right.b);
    EXPECT_EQ(left.log_b, right.log_b);
    EXPECT_EQ(left.width_gradient, right.width_gradient);
    EXPECT_DOUBLE_EQ(left.objective, right.objective);
    EXPECT_EQ(left.nuisance_amplitudes, right.nuisance_amplitudes);
}

void ExpectSameFitNumerics(const core::JointFitResult & left, const core::JointFitResult & right)
{
    EXPECT_EQ(left.initialization.valid, right.initialization.valid);
    EXPECT_EQ(left.initialization.reason, right.initialization.reason);
    EXPECT_EQ(left.initialization.b, right.initialization.b);
    EXPECT_EQ(left.search_completed, right.search_completed);
    EXPECT_EQ(left.available_row_mask, right.available_row_mask);
    EXPECT_EQ(left.prediction, right.prediction);
    EXPECT_EQ(left.objective, right.objective);
    EXPECT_EQ(left.RuntimeConvergence(), right.RuntimeConvergence());
    EXPECT_EQ(left.TargetRuntimeConvergence(), right.TargetRuntimeConvergence());
    ASSERT_EQ(left.assembled_state.has_value(), right.assembled_state.has_value());
    if (left.assembled_state) ExpectSameState(*left.assembled_state, *right.assembled_state);
    ExpectSameChecks(left.evidence, right.evidence);

    ASSERT_EQ(left.components.size(), right.components.size());
    for (std::size_t k = 0; k < left.components.size(); ++k)
    {
        const auto & a = left.components[k];
        const auto & b = right.components[k];
        EXPECT_EQ(a.id, b.id);
        EXPECT_EQ(a.stop_reason, b.stop_reason);
        EXPECT_EQ(a.atoms, b.atoms);
        EXPECT_EQ(a.rows, b.rows);
        EXPECT_EQ(a.search_completed, b.search_completed);
        EXPECT_EQ(a.profile_evaluations, b.profile_evaluations);
        EXPECT_EQ(a.reference_evaluations, b.reference_evaluations);
        EXPECT_EQ(a.accepted_updates, b.accepted_updates);
        EXPECT_EQ(a.native_status, b.native_status);
        ASSERT_EQ(a.state.has_value(), b.state.has_value());
        if (a.state) ExpectSameState(*a.state, *b.state);
        ExpectSameChecks(a.evidence, b.evidence);
    }
}

void ExpectLifecycle(const core::JointProblem & problem, const core::JointFitResult & fit,
    const std::vector<joint::JointProgressEvent> & events,
    const joint::FixedNeighborSearchPolicy & policy)
{
    const auto & data = core::JointProblemAccess::Get(problem);
    ASSERT_EQ(fit.components.size(), data.partition.components.size());
    ASSERT_FALSE(events.empty());
    EXPECT_EQ(events.front().phase, joint::JointProgressPhase::SolverConfigured);
    ASSERT_TRUE(events.front().solver_configuration);
    const auto expected_configuration = joint::ResolveJointSolverConfiguration(policy);
    EXPECT_EQ(events.front().solver_configuration->fixed_neighbor_core_atoms,
        expected_configuration.fixed_neighbor_core_atoms);
    EXPECT_EQ(std::count_if(events.begin(), events.end(), [](const auto & event) {
        return event.phase == joint::JointProgressPhase::SolverConfigured;
    }), 1);
    std::size_t cursor = 1;
    for (std::size_t k = 0; k < fit.components.size(); ++k)
    {
        const auto & view = data.partition.components[k];
        const auto & component = fit.components[k];
        ASSERT_LT(cursor, events.size());
        const auto & started = events[cursor++];
        EXPECT_EQ(started.phase, joint::JointProgressPhase::ComponentStarted);
        EXPECT_EQ(started.component_index, k + 1);
        EXPECT_EQ(started.component_count, fit.components.size());
        EXPECT_EQ(started.component_id, view.id);
        EXPECT_EQ(started.atom_count, view.atoms.size());
        EXPECT_EQ(started.row_count, view.rows.size());

        int last_evaluation = 0;
        int last_accepted = 0;
        std::optional<double> last_search_objective;
        std::optional<double> last_search_gradient;
        bool fixed_neighbor_progress{};
        while (cursor < events.size() && events[cursor].phase == joint::JointProgressPhase::SearchProgress)
        {
            const auto & progress = events[cursor++];
            EXPECT_EQ(progress.component_index, k + 1);
            EXPECT_EQ(progress.component_id, view.id);
            EXPECT_EQ(progress.atom_count, view.atoms.size());
            EXPECT_EQ(progress.row_count, view.rows.size());
            EXPECT_EQ(progress.profile_budget, data.context.profile_budget);
            if (progress.fixed_neighbor)
            {
                fixed_neighbor_progress = true;
                EXPECT_EQ(progress.update_budget, static_cast<int>(progress.fixed_neighbor->block_solves));
                EXPECT_EQ(progress.profile_evaluations,
                    static_cast<int>(progress.fixed_neighbor->profile_evaluations));
                EXPECT_EQ(progress.fixed_neighbor->accepted_blocks + progress.fixed_neighbor->unchanged_blocks +
                    progress.fixed_neighbor->failed_blocks, progress.fixed_neighbor->block_solves);
                EXPECT_LE(progress.fixed_neighbor->recovery_successes,
                    progress.fixed_neighbor->recovery_attempts);
                EXPECT_GE(progress.fixed_neighbor->maximum_accepted_log_width_step, 0.0);
                EXPECT_GT(progress.fixed_neighbor->maximum_target_width, 0.0);
            }
            else EXPECT_EQ(progress.update_budget, data.context.update_budget);
            EXPECT_GE(progress.profile_evaluations, last_evaluation);
            EXPECT_GE(progress.accepted_updates, last_accepted);
            EXPECT_GE(progress.elapsed_seconds, 0);
            last_evaluation = progress.profile_evaluations;
            last_accepted = progress.accepted_updates;
            if (progress.accepted_objective) last_search_objective = progress.accepted_objective;
            if (progress.accepted_gradient_inf_norm) last_search_gradient = progress.accepted_gradient_inf_norm;
        }
        EXPECT_GT(last_evaluation, 0);
        if (!fixed_neighbor_progress)
        {
            EXPECT_EQ(last_evaluation, component.profile_evaluations);
            EXPECT_EQ(last_accepted, component.accepted_updates);
        }

        ASSERT_LT(cursor, events.size());
        const auto & certifying = events[cursor++];
        EXPECT_EQ(certifying.phase, joint::JointProgressPhase::CertificationStarted);
        EXPECT_EQ(certifying.component_index, k + 1);
        EXPECT_EQ(certifying.profile_evaluations, component.profile_evaluations);
        EXPECT_EQ(certifying.accepted_updates, component.accepted_updates);
        EXPECT_EQ(certifying.profile_budget, data.context.profile_budget);
        EXPECT_EQ(certifying.update_budget, data.context.update_budget);
        EXPECT_EQ(certifying.accepted_objective, last_search_objective);
        EXPECT_EQ(certifying.accepted_gradient_inf_norm, last_search_gradient);

        ASSERT_LT(cursor, events.size());
        const auto & completed = events[cursor++];
        EXPECT_EQ(completed.phase, joint::JointProgressPhase::ComponentCompleted);
        EXPECT_EQ(completed.component_index, k + 1);
        EXPECT_EQ(completed.component_id, view.id);
        EXPECT_EQ(completed.atom_count, view.atoms.size());
        EXPECT_EQ(completed.row_count, view.rows.size());
        EXPECT_EQ(completed.profile_evaluations, component.profile_evaluations);
        EXPECT_EQ(completed.profile_budget, data.context.profile_budget);
        EXPECT_EQ(completed.accepted_updates, component.accepted_updates);
        EXPECT_EQ(completed.update_budget, data.context.update_budget);
        EXPECT_EQ(completed.stop_reason, component.stop_reason);
        EXPECT_EQ(completed.trusted_state, component.state.has_value());
        EXPECT_EQ(completed.accepted_objective, last_search_objective);
        EXPECT_EQ(completed.accepted_gradient_inf_norm, last_search_gradient);
    }
    ASSERT_EQ(events.size(), cursor + 2);
    EXPECT_EQ(events[cursor].phase, joint::JointProgressPhase::AssemblyStarted);
    EXPECT_EQ(events[cursor].component_count, fit.components.size());
    EXPECT_EQ(events[cursor + 1].phase, joint::JointProgressPhase::AssemblyCompleted);
    EXPECT_EQ(events[cursor + 1].component_count, fit.components.size());
    EXPECT_FALSE(events[cursor].accepted_objective);
    EXPECT_FALSE(events[cursor + 1].accepted_objective);
    EXPECT_FALSE(events[cursor].accepted_gradient_inf_norm);
    EXPECT_FALSE(events[cursor + 1].accepted_gradient_inf_norm);
    EXPECT_GE(events[cursor + 1].elapsed_seconds, 0);
}
}

TEST(JointProgressTest, PublicFitWithoutObserverIsSilent)
{
    const auto previous_level = Logger::GetLogLevel();
    Logger::SetLogLevel(LogLevel::Info);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    const auto fit = core::FitJointComponents(core::JointProblem(MakeInput()), {.55, .55});
    const auto stdout_text = testing::internal::GetCapturedStdout();
    const auto stderr_text = testing::internal::GetCapturedStderr();
    Logger::SetLogLevel(previous_level);
    ASSERT_TRUE(fit.assembled_state);
    EXPECT_TRUE(stdout_text.empty());
    EXPECT_TRUE(stderr_text.empty());
}

TEST(JointProgressTest, PublicFitUsesFixedNeighborProductionContract)
{
    const core::JointProblem problem(MakeInput());
    const std::vector<double> initial{.55, .55};
    const auto expected = joint::FitFixedNeighborComponents(problem, initial, {});
    const auto actual = core::FitJointComponents(problem, initial);

    EXPECT_EQ(actual.solver_provenance.search_method, "fixed-neighbor");
    EXPECT_EQ(actual.solver_provenance.fixed_neighbor_core_atoms, 12u);
    EXPECT_EQ(actual.solver_provenance.fixed_neighbor_policy_version,
        rhbm_gem::FixedNeighborPolicyContractVersion);
    EXPECT_EQ(actual.solver_provenance.fixed_neighbor_maximum_sweeps, 200u);
    EXPECT_EQ(actual.solver_provenance.fixed_neighbor_order, "forward");
    EXPECT_EQ(actual.solver_provenance.fixed_neighbor_local_search, "profile-lm");
    EXPECT_EQ(actual.RuntimeConvergence(), expected.RuntimeConvergence());
    ASSERT_TRUE(actual.assembled_state);
    ASSERT_TRUE(expected.assembled_state);
    ExpectSameFitNumerics(actual, expected);
}

TEST(JointProgressTest, FixedNeighborUsesConfiguredQualifiedDefaults)
{
    const core::JointProblem problem(MakeInput());
    const std::vector<double> initial{.55, .55};
    const joint::FixedNeighborSearchPolicy policy{};
    EXPECT_EQ(policy.core_atoms,12u);
    std::vector<joint::JointProgressEvent> events;
    const auto fit=joint::FitFixedNeighborComponents(problem,initial,policy,
        [&](const auto & event) { events.push_back(event); });

    ASSERT_EQ(fit.components.size(),2u);
    EXPECT_EQ(fit.solver_provenance.search_method,"fixed-neighbor");
    EXPECT_EQ(fit.solver_provenance.fixed_neighbor_core_atoms,12u);
    EXPECT_FALSE(fit.solver_provenance.fixed_neighbor_local_work);
    EXPECT_EQ(fit.solver_provenance.fixed_neighbor_policy_version,
        rhbm_gem::FixedNeighborPolicyContractVersion);
    EXPECT_EQ(fit.solver_provenance.fixed_neighbor_maximum_sweeps,200u);
    EXPECT_EQ(fit.solver_provenance.fixed_neighbor_order,"forward");
    EXPECT_EQ(fit.solver_provenance.fixed_neighbor_local_search,"profile-lm");
    ASSERT_TRUE(fit.assembled_state);
    EXPECT_EQ(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(events.front().solver_configuration);
    EXPECT_EQ(events.front().solver_configuration->fixed_neighbor_core_atoms,12u);
    EXPECT_FALSE(events.front().solver_configuration->bounded_width_profile);
    ASSERT_TRUE(std::any_of(events.begin(),events.end(),[](const auto & event) {
        return event.fixed_neighbor && event.fixed_neighbor->sweep>0;
    }));
    ExpectLifecycle(problem,fit,events,policy);
}

TEST(JointProgressTest, ObservableComponentsUseTheSameLifecycle)
{
    const core::JointProblem problem(MakeObservableInput());
    const std::vector<double> initial{.55};
    std::vector<joint::JointProgressEvent> events;
    const joint::JointProgressObserver observer = [&](const auto & event) { events.push_back(event); };
    const auto fit = joint::FitFixedNeighborComponents(problem, initial, {}, observer);
    EXPECT_TRUE(problem.ParameterLayout().groups.empty());
    ASSERT_TRUE(fit.assembled_state);
    ExpectLifecycle(problem, fit, events, {});
    EXPECT_TRUE(std::none_of(fit.components.begin(),fit.components.end(),[](const auto & component) {
        return component.stop_reason=="analytic-nuisance-only";
    }));
}

TEST(JointProgressTest, QuietWorkflowSuppressesProgressOutput)
{
    auto atom = std::make_unique<rhbm_gem::AtomObject>();
    atom->SetSerialID(1);
    atom->SetElement(Element::CARBON);
    atom->SetPosition(0, 0, 0);
    atom->SetChainID("A");
    atom->SetComponentID("ALA");
    atom->SetAtomID("C");
    std::vector<std::unique_ptr<rhbm_gem::AtomObject>> atoms;
    atoms.push_back(std::move(atom));
    rhbm_gem::ModelObject model(std::move(atoms));
    model.SelectAllAtoms();

    rhbm_gem::MapObject map({25, 25, 25}, {.3, .3, .3}, {-3.6, -3.6, -3.6});
    core::simulation::SimulationAtomPreparationResult generator;
    generator.atom_list.push_back(core::simulation::SimulationAtom{
        .serial_id = 1, .element = Element::CARBON, .position = {0, 0, 0}, .charge_used = .2});
    core::MapSimulationRequest request;
    request.job_count = 1;
    request.cutoff_distance = 2.5;
    request.potential_model_choice = core::PotentialModel::SINGLE_GAUS;
    core::simulation::PopulateMapValueArray(map, generator, request, .5);

    core::FitOptions options;
    options.estimator = core::PotentialEstimator::JOINT_COMPONENTS;
    options.thread_size = 1;
    options.quiet_mode = true;
    const auto previous_level = Logger::GetLogLevel();
    Logger::SetLogLevel(LogLevel::Info);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    core::RunPotentialFittingWorkflow(map, model, options);
    const auto stdout_text = testing::internal::GetCapturedStdout();
    const auto stderr_text = testing::internal::GetCapturedStderr();
    Logger::SetLogLevel(previous_level);

    EXPECT_TRUE(stdout_text.empty());
    EXPECT_TRUE(stderr_text.empty());
    EXPECT_TRUE(model.GetAnalysisView().GetJointResult());
}
