#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>
#include <rhbm_gem/core/CommandTypes.hpp>
#include "core/detail/joint_component/JointComponentEstimator.hpp"
#include <rhbm_gem/data/object/AtomObject.hpp>
#include <rhbm_gem/data/object/MapObject.hpp>
#include <rhbm_gem/data/object/ModelAnalysisView.hpp>
#include <rhbm_gem/data/object/ModelObject.hpp>
#include <rhbm_gem/utils/domain/Logger.hpp>

#include "core/command/detail/MapSimulation.hpp"
#include "core/command/detail/SimulationGeometry.hpp"
#include "core/detail/joint_component/JointCliProgressReporter.hpp"
#include "core/detail/joint_component/Problem.hpp"
#include "core/detail/joint_component/ResourceWork.hpp"
#include "support/JointSyntheticWorkload.hpp"
#include "support/JointComponentChecks.hpp"

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
        EXPECT_EQ(a.accepted_updates, b.accepted_updates);
        EXPECT_EQ(a.native_status, b.native_status);
        ASSERT_EQ(a.state.has_value(), b.state.has_value());
        if (a.state) ExpectSameState(*a.state, *b.state);
        ExpectSameChecks(a.evidence, b.evidence);
    }
}

void ExpectLifecycle(const core::JointProblem & problem, const core::JointFitResult & fit,
    const std::vector<joint::JointProgressEvent> & events)
{
    const auto & data = core::JointProblemAccess::Get(problem);
    ASSERT_EQ(fit.components.size(), data.partition.components.size());
    ASSERT_FALSE(events.empty());
    EXPECT_EQ(events.front().phase, joint::JointProgressPhase::SolverConfigured);
    ASSERT_TRUE(events.front().solver_configuration);
    EXPECT_EQ(events.front().solver_configuration->fixed_neighbor_core_atoms,
        joint::ProductionFixedNeighborCoreAtoms);
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
        std::optional<joint::JointFixedNeighborProgress> last_fixed_neighbor;
        std::size_t total_evaluations{}, total_local_accepted_updates{};
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
                last_fixed_neighbor = progress.fixed_neighbor;
                total_evaluations += progress.fixed_neighbor->profile_evaluations;
                total_local_accepted_updates += progress.fixed_neighbor->local_accepted_updates;
                EXPECT_EQ(progress.update_budget, static_cast<int>(progress.fixed_neighbor->block_solves));
                EXPECT_EQ(progress.profile_evaluations,
                    static_cast<int>(progress.fixed_neighbor->profile_evaluations));
                EXPECT_EQ(progress.fixed_neighbor->accepted_blocks + progress.fixed_neighbor->unchanged_blocks +
                    progress.fixed_neighbor->failed_blocks, progress.fixed_neighbor->block_solves);
                EXPECT_LE(progress.fixed_neighbor->recovery_successes,
                    progress.fixed_neighbor->recovery_attempts);
                EXPECT_GT(progress.fixed_neighbor->maximum_target_width, 0.0);
            }
            else
            {
                EXPECT_EQ(progress.update_budget, data.context.update_budget);
                EXPECT_GE(progress.profile_evaluations, last_evaluation);
                EXPECT_GE(progress.accepted_updates, last_accepted);
            }
            EXPECT_GE(progress.elapsed_seconds, 0);
            last_evaluation = progress.profile_evaluations;
            last_accepted = progress.accepted_updates;
            if (progress.accepted_objective) last_search_objective = progress.accepted_objective;
            if (progress.accepted_gradient_inf_norm) last_search_gradient = progress.accepted_gradient_inf_norm;
        }
        EXPECT_GT(last_evaluation, 0);
        if (last_fixed_neighbor)
        {
            EXPECT_EQ(total_evaluations, component.profile_evaluations);
            EXPECT_EQ(total_local_accepted_updates, component.accepted_updates);
        }
        else
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
        if (last_fixed_neighbor)
        {
            ASSERT_TRUE(certifying.fixed_neighbor);
            EXPECT_EQ(certifying.fixed_neighbor->sweep, last_fixed_neighbor->sweep);
            EXPECT_EQ(certifying.fixed_neighbor->maximum_sweeps, last_fixed_neighbor->maximum_sweeps);
            EXPECT_EQ(certifying.fixed_neighbor->profile_evaluations, last_fixed_neighbor->profile_evaluations);
            EXPECT_EQ(certifying.fixed_neighbor->accepted_blocks, last_fixed_neighbor->accepted_blocks);
            EXPECT_EQ(certifying.fixed_neighbor->local_accepted_updates, last_fixed_neighbor->local_accepted_updates);
            EXPECT_EQ(certifying.fixed_neighbor->failed_blocks, last_fixed_neighbor->failed_blocks);
            EXPECT_DOUBLE_EQ(certifying.fixed_neighbor->objective, last_fixed_neighbor->objective);
        }

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

std::size_t CountOccurrences(std::string_view text, std::string_view value)
{
    std::size_t count{};
    std::size_t position{};
    while ((position = text.find(value, position)) != std::string_view::npos)
    {
        ++count;
        position += value.size();
    }
    return count;
}

joint::JointProgressEvent MakeFixedNeighborEvent(joint::JointProgressPhase phase)
{
    joint::JointProgressEvent event;
    event.phase = phase;
    event.component_index = 1;
    event.component_count = 1;
    event.component_id = "component-A";
    event.atom_count = 12;
    event.row_count = 80;
    event.elapsed_seconds = 12.5;
    event.fixed_neighbor = joint::JointFixedNeighborProgress{
        .sweep = 18,
        .maximum_sweeps = 200,
        .block_solves = 8,
        .accepted_blocks = 3,
        .local_accepted_updates = 5,
        .profile_evaluations = 22,
        .unchanged_blocks = 3,
        .failed_blocks = 2,
        .recovery_attempts = 1,
        .recovery_successes = 0,
        .objective = 1.25,
        .global_ac_kkt = .003,
        .width_gradient_inf_norm = .004,
        .projected_width_kkt = 4e-7,
        .maximum_target_width = .7};
    return event;
}
}

TEST(JointProgressTest, BoundedSyntheticFitWithoutObserverIsSilent)
{
    const auto previous_level = Logger::GetLogLevel();
    Logger::SetLogLevel(LogLevel::Info);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    const auto fit=second_stage_test::matched::joint_abc::FitBoundedJointComponentsForTesting(
        core::JointProblem(MakeInput()),{.55,.55},2.2);
    const auto stdout_text = testing::internal::GetCapturedStdout();
    const auto stderr_text = testing::internal::GetCapturedStderr();
    Logger::SetLogLevel(previous_level);
    ASSERT_TRUE(fit.assembled_state);
    EXPECT_TRUE(stdout_text.empty());
    EXPECT_TRUE(stderr_text.empty());
}

TEST(JointProgressTest, BoundedSyntheticFitUsesFixedNeighborContract)
{
    const core::JointProblem problem(MakeInput());
    const std::vector<double> initial{.55, .55};
    const auto actual=second_stage_test::matched::joint_abc::FitBoundedJointComponentsForTesting(
        problem,initial,2.2);
    auto policy=joint::FixedNeighborSearchPolicy{};
    policy.eta_upper_bound=std::log(2.2);
    const auto expected=joint::FitFixedNeighborComponents(problem,initial,policy);

    EXPECT_EQ(actual.solver_provenance.search_method, "fixed-neighbor");
    EXPECT_EQ(actual.solver_provenance.fixed_neighbor_core_atoms, 12u);
    EXPECT_EQ(actual.solver_provenance.fixed_neighbor_policy_version,
        rhbm_gem::FixedNeighborPolicyContractVersion);
    EXPECT_EQ(actual.solver_provenance.fixed_neighbor_maximum_sweeps, 200u);
    EXPECT_EQ(actual.solver_provenance.fixed_neighbor_order, "forward");
    EXPECT_EQ(actual.solver_provenance.fixed_neighbor_local_search, "active-set-lm");
    EXPECT_EQ(actual.RuntimeConvergence(), expected.RuntimeConvergence());
    ASSERT_TRUE(actual.assembled_state);
    ASSERT_TRUE(expected.assembled_state);
    ExpectSameFitNumerics(actual, expected);
}

TEST(JointProgressTest, FixedNeighborUsesConfiguredQualifiedDefaults)
{
    const core::JointProblem problem(MakeInput());
    const std::vector<double> initial{.55, .55};
    auto policy=joint::FixedNeighborSearchPolicy{};
    policy.eta_upper_bound=std::log(2.2);
    std::vector<joint::JointProgressEvent> events;
    const auto fit=joint::FitFixedNeighborComponents(problem,initial,policy,
        [&](const auto & event) { events.push_back(event); });

    ASSERT_EQ(fit.components.size(),2u);
    EXPECT_EQ(fit.solver_provenance.search_method,"fixed-neighbor");
    EXPECT_EQ(fit.solver_provenance.fixed_neighbor_core_atoms,12u);
    EXPECT_EQ(fit.solver_provenance.fixed_neighbor_policy_version,
        rhbm_gem::FixedNeighborPolicyContractVersion);
    EXPECT_EQ(fit.solver_provenance.fixed_neighbor_maximum_sweeps,200u);
    EXPECT_EQ(fit.solver_provenance.fixed_neighbor_order,"forward");
    EXPECT_EQ(fit.solver_provenance.fixed_neighbor_local_search,"active-set-lm");
    ASSERT_TRUE(fit.assembled_state);
    EXPECT_EQ(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(events.front().solver_configuration);
    EXPECT_EQ(events.front().solver_configuration->fixed_neighbor_core_atoms,12u);
    ASSERT_TRUE(std::any_of(events.begin(),events.end(),[](const auto & event) {
        return event.fixed_neighbor && event.fixed_neighbor->sweep>0;
    }));
    ExpectLifecycle(problem,fit,events);
}

TEST(JointProgressTest, ObservableComponentsUseTheSameLifecycle)
{
    const core::JointProblem problem(MakeObservableInput());
    const std::vector<double> initial{.55};
    std::vector<joint::JointProgressEvent> events;
    const joint::JointProgressObserver observer = [&](const auto & event) { events.push_back(event); };
    auto policy=joint::FixedNeighborSearchPolicy{};
    policy.eta_upper_bound=std::log(2.2);
    const auto fit=joint::FitFixedNeighborComponents(problem,initial,policy,observer);
    ASSERT_TRUE(fit.assembled_state);
    ExpectLifecycle(problem,fit,events);
}

TEST(JointProgressTest, QuietWorkflowSuppressesProgressButKeepsFailureWarning)
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
    EXPECT_EQ(stderr_text,
        "[Warning] [Joint] Qualification failed; preserving diagnostics and skipping downstream processing.\n");
    EXPECT_TRUE(model.GetAnalysisView().GetJointResult());
}

TEST(JointProgressTest, InfoProgressIsCompactAndThrottlesAcceptedUpdateRefreshes)
{
    const auto previous_level = Logger::GetLogLevel();
    Logger::SetLogLevel(LogLevel::Info);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    core::detail::JointCliProgressReporter reporter;
    auto started = MakeFixedNeighborEvent(joint::JointProgressPhase::ComponentStarted);
    reporter.OnProgress(started);

    auto search = MakeFixedNeighborEvent(joint::JointProgressPhase::SearchProgress);
    search.fixed_neighbor->accepted_blocks = 4;
    reporter.OnProgress(search); // The search phase transition is immediate.
    search.fixed_neighbor->accepted_blocks = 5;
    reporter.OnProgress(search); // A changed accepted count alone is throttled.
    std::this_thread::sleep_for(std::chrono::milliseconds(520));
    search.fixed_neighbor->accepted_blocks = 6;
    reporter.OnProgress(search); // Long-running search refreshes after the interval.

    auto certifying = MakeFixedNeighborEvent(joint::JointProgressPhase::CertificationStarted);
    reporter.OnProgress(certifying); // Phase transitions are immediate.

    auto completed = MakeFixedNeighborEvent(joint::JointProgressPhase::ComponentCompleted);
    completed.stop_reason = "converged";
    completed.trusted_state = true;
    completed.accepted_objective = 1.25;
    reporter.OnProgress(completed);

    auto failed_assembly = MakeFixedNeighborEvent(joint::JointProgressPhase::AssemblyCompleted);
    failed_assembly.trusted_state = false;
    reporter.OnProgress(failed_assembly);

    const auto stdout_text = testing::internal::GetCapturedStdout();
    const auto stderr_text = testing::internal::GetCapturedStderr();
    Logger::SetLogLevel(previous_level);

    EXPECT_EQ(CountOccurrences(stdout_text, "\r"), 4u);
    EXPECT_EQ(CountOccurrences(stdout_text, "obj="), 5u); // Four progress rows and the summary.
    EXPECT_EQ(CountOccurrences(stdout_text, "KKT="), 4u);
    EXPECT_EQ(CountOccurrences(stdout_text, "width-grad="), 4u);
    EXPECT_EQ(CountOccurrences(stdout_text, "width-kkt=4.000e-07"), 4u);
    EXPECT_EQ(CountOccurrences(stdout_text, "accepted="), 4u);
    EXPECT_EQ(CountOccurrences(stdout_text, "failed="), 4u);
    EXPECT_NE(stdout_text.find("Component 1/1 completed | state_trusted=yes | stop=converged"),
        std::string::npos);
    EXPECT_EQ(stdout_text.back(), '\n');
    EXPECT_NE(stderr_text.find(
        "[Warning] [Joint] Global assembly did not produce a trusted state | elapsed=12.5s\n"),
        std::string::npos);
    EXPECT_EQ(stdout_text.find("unchanged_blocks"), std::string::npos);
    EXPECT_EQ(stdout_text.find("local_accepted_updates"), std::string::npos);
    EXPECT_EQ(stdout_text.find("recovery="), std::string::npos);
}

TEST(JointProgressTest, DebugProgressRetainsDetailedDiagnostics)
{
    const auto previous_level = Logger::GetLogLevel();
    Logger::SetLogLevel(LogLevel::Debug);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    core::detail::JointCliProgressReporter reporter;
    reporter.OnProgress(MakeFixedNeighborEvent(joint::JointProgressPhase::ComponentStarted));
    reporter.OnProgress(MakeFixedNeighborEvent(joint::JointProgressPhase::SearchProgress));
    auto completed = MakeFixedNeighborEvent(joint::JointProgressPhase::ComponentCompleted);
    completed.stop_reason = "converged";
    completed.trusted_state = true;
    completed.accepted_objective = 1.25;
    completed.profile_evaluations = 220;
    completed.profile_budget = 80;
    completed.accepted_updates = 50;
    completed.update_budget = 40;
    reporter.OnProgress(completed);

    const auto stdout_text = testing::internal::GetCapturedStdout();
    const auto stderr_text = testing::internal::GetCapturedStderr();
    Logger::SetLogLevel(previous_level);

    EXPECT_TRUE(stderr_text.empty());
    EXPECT_NE(stdout_text.find("[Debug] [Joint] Progress detail"), std::string::npos);
    EXPECT_NE(stdout_text.find("unchanged_blocks=3"), std::string::npos);
    EXPECT_NE(stdout_text.find("local_accepted_updates=5"), std::string::npos);
    EXPECT_NE(stdout_text.find("profile_evaluations=22"), std::string::npos);
    EXPECT_NE(stdout_text.find("recovery=1/0"), std::string::npos);
    EXPECT_NE(stdout_text.find("max_B=7.000e-01"), std::string::npos);
    EXPECT_EQ(stdout_text.find("max_log_width_step="), std::string::npos);
    EXPECT_EQ(stdout_text.find("eta_change="), std::string::npos);
    EXPECT_NE(stdout_text.find("Component detail"), std::string::npos);
    EXPECT_NE(stdout_text.find("evaluations=220 | accepted_updates=50 | stop=converged"), std::string::npos);
    EXPECT_EQ(CountOccurrences(stdout_text, "obj="), 3u);
    EXPECT_EQ(CountOccurrences(stdout_text, "KKT="), 2u);
    EXPECT_EQ(CountOccurrences(stdout_text, "width-grad="), 2u);
    EXPECT_EQ(stdout_text.find("width_grad="), std::string::npos);
}

TEST(JointProgressTest, MemoryPreflightLoggingPreservesWarningsAndDebugBreakdown)
{
    core::joint_component::MemoryEstimate estimate;
    estimate.stage = core::joint_component::MemoryStage::ProblemConstruction;
    estimate.status = core::joint_component::MemoryStatus::Pass;
    estimate.budget_bytes = 1000;
    estimate.estimated_peak_bytes = 900;
    estimate.shape.component_atoms = 12;

    const auto previous_level = Logger::GetLogLevel();
    Logger::SetLogLevel(LogLevel::Info);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    core::joint_component::LogMemoryEstimate(estimate);
    const auto info_text = testing::internal::GetCapturedStdout();
    const auto info_error = testing::internal::GetCapturedStderr();

    Logger::SetLogLevel(LogLevel::Debug);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    core::joint_component::LogMemoryEstimate(estimate);
    const auto debug_text = testing::internal::GetCapturedStdout();
    const auto debug_error = testing::internal::GetCapturedStderr();

    estimate.status = core::joint_component::MemoryStatus::Reject;
    estimate.reason = "resource-limit:test";
    Logger::SetLogLevel(LogLevel::Info);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    core::joint_component::LogMemoryEstimate(estimate);
    const auto warning_text = testing::internal::GetCapturedStdout();
    const auto warning_error = testing::internal::GetCapturedStderr();
    Logger::SetLogLevel(previous_level);

    EXPECT_TRUE(info_error.empty());
    EXPECT_NE(info_text.find("[Joint][Memory] problem-construction"), std::string::npos);
    EXPECT_EQ(info_text.find("component_atoms="), std::string::npos);
    EXPECT_NE(debug_text.find("component_atoms=12"), std::string::npos);
    EXPECT_NE(debug_text.find("estimated_peak="), std::string::npos);
    EXPECT_TRUE(debug_error.empty());
    EXPECT_TRUE(warning_text.empty());
    EXPECT_NE(warning_error.find("[Warning] [Joint][Memory]"), std::string::npos);
    EXPECT_NE(warning_error.find("resource-limit:test"), std::string::npos);
    EXPECT_NE(warning_error.find("component_atoms=12"), std::string::npos);
}

TEST(JointProgressTest, ReportingDoesNotChangeNumericsOrObserverCadence)
{
    const core::JointProblem problem(MakeInput());
    const std::vector<double> initial{.55, .55};
    const auto previous_level = Logger::GetLogLevel();
    auto policy=joint::FixedNeighborSearchPolicy{};
    policy.eta_upper_bound=std::log(2.2);
    const auto run = [&](LogLevel level, bool report) {
        Logger::SetLogLevel(level);
        Logger::FinishProgressLine();
        testing::internal::CaptureStdout();
        testing::internal::CaptureStderr();

        core::detail::JointCliProgressReporter reporter;
        std::vector<joint::JointProgressEvent> events;
        const auto observer = [&](const auto & event) {
            events.push_back(event);
            if (report) reporter.OnProgress(event);
        };
        const auto fit=joint::FitFixedNeighborComponents(problem,initial,policy,observer);
        const auto stdout_text = testing::internal::GetCapturedStdout();
        const auto stderr_text = testing::internal::GetCapturedStderr();
        return std::tuple{fit, std::move(events), stdout_text, stderr_text};
    };

    auto [baseline, baseline_events, baseline_stdout, baseline_stderr] = run(LogLevel::Info, false);
    auto [info_fit, info_events, info_stdout, info_stderr] = run(LogLevel::Info, true);
    auto [debug_fit, debug_events, debug_stdout, debug_stderr] = run(LogLevel::Debug, true);
    Logger::SetLogLevel(previous_level);

    EXPECT_TRUE(baseline_stdout.empty());
    EXPECT_TRUE(baseline_stderr.empty());
    ASSERT_TRUE(baseline.assembled_state);
    ASSERT_TRUE(info_fit.assembled_state);
    ASSERT_TRUE(debug_fit.assembled_state);
    ExpectSameFitNumerics(baseline, info_fit);
    ExpectSameFitNumerics(baseline, debug_fit);
    EXPECT_TRUE(info_stderr.empty());
    EXPECT_TRUE(debug_stderr.empty());
    EXPECT_EQ(CountOccurrences(info_stdout, " | certifying | sweep="), info_fit.components.size());
    EXPECT_EQ(CountOccurrences(debug_stdout, " | certifying | sweep="), debug_fit.components.size());
    EXPECT_EQ(info_stdout.find(" | certifying | eval="), std::string::npos);
    EXPECT_EQ(debug_stdout.find(" | certifying | eval="), std::string::npos);

    ASSERT_EQ(baseline_events.size(), info_events.size());
    ASSERT_EQ(baseline_events.size(), debug_events.size());
    for (std::size_t i = 0; i < baseline_events.size(); ++i)
    {
        const auto & baseline_event = baseline_events[i];
        const auto & info_event = info_events[i];
        const auto & debug_event = debug_events[i];
        EXPECT_EQ(baseline_event.phase, info_event.phase);
        EXPECT_EQ(baseline_event.phase, debug_event.phase);
        EXPECT_EQ(baseline_event.component_index, info_event.component_index);
        EXPECT_EQ(baseline_event.component_index, debug_event.component_index);
        EXPECT_EQ(baseline_event.profile_evaluations, info_event.profile_evaluations);
        EXPECT_EQ(baseline_event.profile_evaluations, debug_event.profile_evaluations);
        EXPECT_EQ(baseline_event.accepted_updates, info_event.accepted_updates);
        EXPECT_EQ(baseline_event.accepted_updates, debug_event.accepted_updates);
        ASSERT_EQ(baseline_event.fixed_neighbor.has_value(), info_event.fixed_neighbor.has_value());
        ASSERT_EQ(baseline_event.fixed_neighbor.has_value(), debug_event.fixed_neighbor.has_value());
        if (baseline_event.fixed_neighbor)
        {
            EXPECT_EQ(baseline_event.fixed_neighbor->sweep, info_event.fixed_neighbor->sweep);
            EXPECT_EQ(baseline_event.fixed_neighbor->sweep, debug_event.fixed_neighbor->sweep);
            EXPECT_DOUBLE_EQ(baseline_event.fixed_neighbor->objective, info_event.fixed_neighbor->objective);
            EXPECT_DOUBLE_EQ(baseline_event.fixed_neighbor->objective, debug_event.fixed_neighbor->objective);
        }
    }
}
