#include <gtest/gtest.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <rhbm_gem/core/GaussianEstimator.hpp>
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
#include "core/detail/joint_component/SparseFactor.hpp"

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
    input.atom_ids = {"1", "2", "3"};
    input.support.resize(3);
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
    for (const auto [atom, square, gaussian, charge] : {
        std::tuple{1, 4.4, .3, .1}, std::tuple{2, 5.1, .2, .15}})
    {
        input.support[static_cast<std::size_t>(atom)].push_back({39, square});
        const auto basis = joint::EvaluateKernel(square, .7, 2.5);
        input.observations[39] += gaussian * basis.gaussian + charge * basis.charge;
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

void ExpectSameRanks(const std::vector<rhbm_gem::JointRankEvidence> & left,
    const std::vector<rhbm_gem::JointRankEvidence> & right)
{
    ASSERT_EQ(left.size(), right.size());
    for (std::size_t k = 0; k < left.size(); ++k)
    {
        EXPECT_EQ(left[k].name, right[k].name);
        EXPECT_EQ(left[k].scope, right[k].scope);
        EXPECT_EQ(left[k].rank, right[k].rank);
        EXPECT_DOUBLE_EQ(left[k].threshold, right[k].threshold);
        EXPECT_EQ(left[k].singular_values, right[k].singular_values);
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
    ExpectSameRanks(left.ranks, right.ranks);

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
        ExpectSameRanks(a.ranks, b.ranks);
    }
}

void ExpectLifecycle(const core::JointProblem & problem, const core::JointFitResult & fit,
    const std::vector<joint::JointProgressEvent> & events, const joint::SearchPolicy & policy)
{
    const auto & data = core::JointProblemAccess::Get(problem);
    ASSERT_EQ(fit.components.size(), data.partition.components.size());
    ASSERT_FALSE(events.empty());
    EXPECT_EQ(events.front().phase, joint::JointProgressPhase::SolverConfigured);
    ASSERT_TRUE(events.front().solver_route);
    const auto expected_route = joint::ResolveJointSolverRoute(policy);
    EXPECT_EQ(events.front().solver_route->sparse_backend, expected_route.sparse_backend);
    EXPECT_EQ(events.front().solver_route->search_method, expected_route.search_method);
    EXPECT_EQ(events.front().solver_route->preconditioner, expected_route.preconditioner);
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
        while (cursor < events.size() && events[cursor].phase == joint::JointProgressPhase::SearchProgress)
        {
            const auto & progress = events[cursor++];
            EXPECT_EQ(progress.component_index, k + 1);
            EXPECT_EQ(progress.component_id, view.id);
            EXPECT_EQ(progress.atom_count, view.atoms.size());
            EXPECT_EQ(progress.row_count, view.rows.size());
            EXPECT_EQ(progress.profile_budget, data.context.profile_budget);
            EXPECT_EQ(progress.update_budget, data.context.update_budget);
            EXPECT_GE(progress.profile_evaluations, last_evaluation);
            EXPECT_GE(progress.accepted_updates, last_accepted);
            EXPECT_GE(progress.elapsed_seconds, 0);
            last_evaluation = progress.profile_evaluations;
            last_accepted = progress.accepted_updates;
        }
        EXPECT_GT(last_evaluation, 0);
        EXPECT_EQ(last_evaluation, component.profile_evaluations);
        EXPECT_EQ(last_accepted, component.accepted_updates);

        ASSERT_LT(cursor, events.size());
        const auto & certifying = events[cursor++];
        EXPECT_EQ(certifying.phase, joint::JointProgressPhase::CertificationStarted);
        EXPECT_EQ(certifying.component_index, k + 1);
        EXPECT_EQ(certifying.profile_evaluations, component.profile_evaluations);
        EXPECT_EQ(certifying.accepted_updates, component.accepted_updates);
        EXPECT_EQ(certifying.profile_budget, data.context.profile_budget);
        EXPECT_EQ(certifying.update_budget, data.context.update_budget);

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
    }
    ASSERT_EQ(events.size(), cursor + 2);
    EXPECT_EQ(events[cursor].phase, joint::JointProgressPhase::AssemblyStarted);
    EXPECT_EQ(events[cursor].component_count, fit.components.size());
    EXPECT_EQ(events[cursor + 1].phase, joint::JointProgressPhase::AssemblyCompleted);
    EXPECT_EQ(events[cursor + 1].component_count, fit.components.size());
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

TEST(JointProgressTest, ActiveSparseBackendMatchesCompiledFactorizationBackend)
{
    const auto expected = joint::SparseBackendEnabled() ? joint::SparseBackend::Spqr : joint::SparseBackend::Eigen;
    EXPECT_EQ(joint::ActiveSparseBackend(), expected);
}

TEST(JointProgressTest, LegacyCompactRouteOmitsUnusedPreconditioner)
{
    joint::SearchPolicy policy;
    policy.method = joint::SearchMethod::LegacyCompact;
    policy.preconditioner = joint::PreconditionerKind::Schwarz;
    const auto route = joint::ResolveJointSolverRoute(policy);
    EXPECT_EQ(route.sparse_backend, joint::ActiveSparseBackend());
    EXPECT_EQ(route.search_method, joint::SearchMethod::LegacyCompact);
    EXPECT_FALSE(route.preconditioner);
}

TEST(JointProgressTest, OperatorPcgRoutePreservesRequestedPreconditioner)
{
    joint::SearchPolicy policy;
    policy.method = joint::SearchMethod::OperatorPcg;
    for (const auto preconditioner : {joint::PreconditionerKind::Identity,
        joint::PreconditionerKind::Diagonal, joint::PreconditionerKind::Schwarz})
    {
        policy.preconditioner = preconditioner;
        const auto route = joint::ResolveJointSolverRoute(policy);
        EXPECT_EQ(route.sparse_backend, joint::ActiveSparseBackend());
        EXPECT_EQ(route.search_method, joint::SearchMethod::OperatorPcg);
        ASSERT_TRUE(route.preconditioner);
        EXPECT_EQ(*route.preconditioner, preconditioner);
    }
}

TEST(JointProgressTest, SolverRouteNamesAreCanonical)
{
    EXPECT_EQ(joint::SparseBackendName(joint::SparseBackend::Eigen), "EIGEN");
    EXPECT_EQ(joint::SparseBackendName(joint::SparseBackend::Spqr), "SPQR");
    EXPECT_EQ(joint::SearchMethodName(joint::SearchMethod::LegacyCompact), "LegacyCompact");
    EXPECT_EQ(joint::SearchMethodName(joint::SearchMethod::OperatorPcg), "OperatorPcg");
    EXPECT_EQ(joint::PreconditionerName(joint::PreconditionerKind::Identity), "Identity");
    EXPECT_EQ(joint::PreconditionerName(joint::PreconditionerKind::Diagonal), "Diagonal");
    EXPECT_EQ(joint::PreconditionerName(joint::PreconditionerKind::Schwarz), "Schwarz");
}

TEST(JointProgressTest, LifecycleAndNumericsMatchForBothSearchPolicies)
{
    const core::JointProblem problem(MakeInput());
    const std::vector<double> initial{.55, .55};
    for (const auto method : {joint::SearchMethod::LegacyCompact, joint::SearchMethod::OperatorPcg})
    {
        joint::SearchPolicy policy;
        policy.method = method;
        const auto without_observer = joint::FitWithSearchPolicy(problem, initial, policy);
        std::vector<joint::JointProgressEvent> events;
        const joint::JointProgressObserver observer = [&](const auto & event) { events.push_back(event); };
        const auto with_observer = joint::FitWithSearchPolicy(problem, initial, policy, observer);
        ExpectSameFitNumerics(without_observer, with_observer);
        ExpectLifecycle(problem, with_observer, events, policy);
    }
}

TEST(JointProgressTest, ObservableComponentsUseTheSameLifecycle)
{
    const core::JointProblem problem(MakeObservableInput());
    const std::vector<double> initial{.55, std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};
    std::vector<joint::JointProgressEvent> events;
    const joint::JointProgressObserver observer = [&](const auto & event) { events.push_back(event); };
    const auto fit = joint::FitWithSearchPolicy(problem, initial, {}, observer);
    ASSERT_EQ(problem.ParameterLayout().groups.size(), 1);
    ASSERT_TRUE(fit.assembled_state);
    ExpectLifecycle(problem, fit, events, {});
}

TEST(JointProgressTest, CliReporterFormatsResolvedSolverRoutes)
{
    const auto previous_level = Logger::GetLogLevel();
    Logger::SetLogLevel(LogLevel::Info);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();

    core::detail::JointCliProgressReporter reporter;
    joint::JointProgressEvent event;
    event.phase = joint::JointProgressPhase::SolverConfigured;
    joint::SearchPolicy policy;
    policy.method = joint::SearchMethod::LegacyCompact;
    event.solver_route = joint::ResolveJointSolverRoute(policy);
    reporter.OnProgress(event);

    policy.method = joint::SearchMethod::OperatorPcg;
    policy.preconditioner = joint::PreconditionerKind::Schwarz;
    event.solver_route = joint::ResolveJointSolverRoute(policy);
    reporter.OnProgress(event);

    const auto output = testing::internal::GetCapturedStdout();
    Logger::SetLogLevel(previous_level);
    const std::string sparse(joint::SparseBackendName(joint::ActiveSparseBackend()));
    const auto legacy = output.find("[Joint] Solver route: sparse=" + sparse + " | width-search=LegacyCompact");
    ASSERT_NE(legacy, std::string::npos);
    const auto legacy_line_end = output.find('\n', legacy);
    EXPECT_EQ(output.substr(legacy, legacy_line_end - legacy).find("preconditioner="), std::string::npos);
    EXPECT_NE(output.find("[Joint] Solver route: sparse=" + sparse
        + " | width-search=OperatorPcg | preconditioner=Schwarz"), std::string::npos);
}

TEST(JointProgressTest, CliReporterRefreshesPhasesAndFinishesLines)
{
    const auto previous_level = Logger::GetLogLevel();
    Logger::SetLogLevel(LogLevel::Info);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();

    core::detail::JointCliProgressReporter reporter;
    joint::JointProgressEvent event;
    event.phase = joint::JointProgressPhase::ComponentStarted;
    event.component_index = 1;
    event.component_count = 2;
    event.component_id = "a";
    event.atom_count = 1;
    event.row_count = 40;
    event.profile_budget = 200;
    event.update_budget = 100;
    reporter.OnProgress(event);

    event.phase = joint::JointProgressPhase::SearchProgress;
    event.profile_evaluations = 1;
    event.elapsed_seconds = .1;
    reporter.OnProgress(event);
    event.profile_evaluations = 2;
    event.accepted_updates = 1;
    reporter.OnProgress(event);

    event.phase = joint::JointProgressPhase::CertificationStarted;
    event.elapsed_seconds = .2;
    reporter.OnProgress(event);
    event.phase = joint::JointProgressPhase::ComponentCompleted;
    event.stop_reason = "native-lm-stop";
    event.trusted_state = true;
    reporter.OnProgress(event);

    event = {};
    event.phase = joint::JointProgressPhase::AssemblyStarted;
    event.component_count = 2;
    reporter.OnProgress(event);
    event.phase = joint::JointProgressPhase::AssemblyCompleted;
    event.trusted_state = true;
    event.elapsed_seconds = .3;
    reporter.OnProgress(event);

    const auto output = testing::internal::GetCapturedStdout();
    Logger::SetLogLevel(previous_level);
    EXPECT_EQ(std::count(output.begin(), output.end(), '\r'), 4);
    EXPECT_NE(output.find("searching | atoms=1 rows=40 | eval=2/200 | accepted=1/100"), std::string::npos);
    EXPECT_NE(output.find("certifying endpoint"), std::string::npos);
    EXPECT_NE(output.find("\n[Joint] Component 1/2 completed | atoms=1 rows=40"), std::string::npos);
    EXPECT_NE(output.find("stop=native-lm-stop | trusted=yes"), std::string::npos);
    EXPECT_NE(output.find("[Joint] Assembling global state..."), std::string::npos);
    EXPECT_NE(output.find("[Joint] Global assembly completed | trusted=yes"), std::string::npos);
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
