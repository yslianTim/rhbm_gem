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
#include "support/JointOperatorWorkload.hpp"

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

core::JointProblemInput MakeObservableInputWithAnalyticNuisance()
{
    auto input = MakeObservableInput();
    const auto row = input.observations.size();
    const double square = .3;
    const auto basis = joint::EvaluateKernel(square, .7, 2.5);
    input.atom_ids.push_back("4");
    input.support.emplace_back();
    input.support.back().push_back({row, square});
    input.row_ids.push_back(std::to_string(row));
    input.observations.push_back(.7 * basis.gaussian + .1 * basis.charge);
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

double NormalizedEndpointObjective(const joint::Endpoint & endpoint,double scale)
{
    return endpoint.certificate.objective/(scale*scale);
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

TEST(JointProgressTest, NormalizedProfileObjectiveUsesCertificateAndScale)
{
    joint::Evaluation evaluation;
    evaluation.valid = true;
    evaluation.certificate.available = true;
    evaluation.residual.resize(2);
    evaluation.residual << 3, 4;
    evaluation.certificate.objective = .5 * evaluation.residual.squaredNorm();

    const auto objective = joint::NormalizedProfileObjective(evaluation, 5);
    ASSERT_TRUE(objective);
    EXPECT_DOUBLE_EQ(*objective, evaluation.certificate.objective / 25);
    EXPECT_DOUBLE_EQ(*objective, .5 * (evaluation.residual / 5).squaredNorm());
    EXPECT_GE(*objective, 0);

    evaluation.valid = false;
    EXPECT_FALSE(joint::NormalizedProfileObjective(evaluation, 5));
    evaluation.valid = true;
    evaluation.certificate.available = false;
    EXPECT_FALSE(joint::NormalizedProfileObjective(evaluation, 5));
    evaluation.certificate.available = true;
    EXPECT_FALSE(joint::NormalizedProfileObjective(evaluation, 0));
    EXPECT_FALSE(joint::NormalizedProfileObjective(evaluation, std::numeric_limits<double>::infinity()));
    evaluation.certificate.objective = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(joint::NormalizedProfileObjective(evaluation, 5));
}

TEST(JointProgressTest, ProfileGradientInfinityNormValidatesAndUsesEvaluationGradient)
{
    joint::Evaluation evaluation;
    evaluation.valid = true;
    evaluation.gradient.resize(3);
    evaluation.gradient << .25, -4, 2;

    const auto gradient_inf_norm = joint::ProfileGradientInfinityNorm(evaluation);
    ASSERT_TRUE(gradient_inf_norm);
    EXPECT_DOUBLE_EQ(*gradient_inf_norm, evaluation.gradient.lpNorm<Eigen::Infinity>());
    EXPECT_DOUBLE_EQ(*gradient_inf_norm, 4);

    evaluation.valid = false;
    EXPECT_FALSE(joint::ProfileGradientInfinityNorm(evaluation));
    evaluation.valid = true;
    evaluation.gradient.resize(0);
    EXPECT_FALSE(joint::ProfileGradientInfinityNorm(evaluation));
    evaluation.gradient = joint::Vector::Constant(2, std::numeric_limits<double>::infinity());
    EXPECT_FALSE(joint::ProfileGradientInfinityNorm(evaluation));
}

TEST(JointProgressTest, RejectedCandidateKeepsAcceptedObjective)
{
    joint::Evaluation accepted;
    accepted.valid = true;
    accepted.certificate.available = true;
    accepted.certificate.objective = 10;
    std::optional<double> accepted_objective;
    joint::UpdateAcceptedProfileObjective(accepted_objective, accepted, 2, true);
    ASSERT_TRUE(accepted_objective);
    const double previous = *accepted_objective;

    auto rejected = accepted;
    rejected.certificate.objective = 40;
    joint::UpdateAcceptedProfileObjective(accepted_objective, rejected, 2, false);
    ASSERT_TRUE(accepted_objective);
    EXPECT_DOUBLE_EQ(*accepted_objective, previous);

    joint::UpdateAcceptedProfileObjective(accepted_objective, rejected, 2, true);
    ASSERT_TRUE(accepted_objective);
    EXPECT_DOUBLE_EQ(*accepted_objective, 10);
}

TEST(JointProgressTest, SearchEventsTrackOnlyAcceptedMetrics)
{
    const core::JointProblem problem(second_stage_test::OperatorWorkload("chain", 4));
    const auto & data = core::JointProblemAccess::Get(problem);
    const auto & view = data.partition.components.front();
    const auto y = joint::SelectValues(data.y, view.rows);
    auto context = joint::ChildContext(data.context, view, true);
    const joint::JointProgressComponent progress_component{
        1, data.partition.components.size(), view.id, view.atoms.size(), view.rows.size()};

    for (const double initial_width : {.05, 2.})
    {
        const auto initial = joint::Vector::Constant(static_cast<Eigen::Index>(view.atoms.size()), initial_width);
        const auto without_observer = joint::SearchProfile(view.domain, y, initial, context);
        std::vector<joint::JointProgressEvent> events;
        const joint::JointProgressObserver observer = [&](const auto & event) { events.push_back(event); };
        const auto observed = joint::SearchProfile(view.domain, y, initial, context, observer, &progress_component);

        EXPECT_EQ(observed.evaluations, without_observer.evaluations);
        EXPECT_EQ(observed.references, without_observer.references);
        EXPECT_EQ(observed.accepted, without_observer.accepted);
        EXPECT_EQ(observed.derivatives, without_observer.derivatives);
        EXPECT_EQ(observed.stop_reason, without_observer.stop_reason);
        EXPECT_EQ(observed.eta, without_observer.eta);
        EXPECT_EQ(observed.accepted_objective, without_observer.accepted_objective);
        EXPECT_EQ(observed.accepted_gradient_inf_norm, without_observer.accepted_gradient_inf_norm);
        EXPECT_EQ(observed.trials.size(), without_observer.trials.size());
        ASSERT_EQ(observed.trials.size(), without_observer.trials.size());
        for (std::size_t k = 0; k < observed.trials.size(); ++k)
        {
            EXPECT_EQ(observed.trials[k].endpoint.eta, without_observer.trials[k].endpoint.eta);
            EXPECT_EQ(observed.trials[k].endpoint.beta, without_observer.trials[k].endpoint.beta);
            EXPECT_EQ(observed.trials[k].accepted, without_observer.trials[k].accepted);
            EXPECT_EQ(observed.trials[k].accepted_update, without_observer.trials[k].accepted_update);
        }
        ASSERT_TRUE(observed.initial_accepted);
        ASSERT_TRUE(observed.accepted_objective);
        ASSERT_TRUE(observed.accepted_gradient_inf_norm);
        EXPECT_TRUE(std::isfinite(*observed.accepted_objective));
        EXPECT_GE(*observed.accepted_objective, 0);
        EXPECT_TRUE(std::isfinite(*observed.accepted_gradient_inf_norm));
        EXPECT_GE(*observed.accepted_gradient_inf_norm, 0);

        std::vector<std::optional<double>> accepted_objectives(
            static_cast<std::size_t>(observed.accepted + 1));
        std::vector<std::optional<double>> accepted_gradients(
            static_cast<std::size_t>(observed.accepted + 1));
        accepted_objectives[0] = NormalizedEndpointObjective(observed.initial, context.scale);
        accepted_gradients[0] = observed.initial.gradient.lpNorm<Eigen::Infinity>();
        int accepted_update{};
        bool observed_rejection{};
        for (const auto & trial : observed.trials)
        {
            if (trial.accepted && trial.accepted_update && *trial.accepted_update > 0)
            {
                accepted_update = *trial.accepted_update;
                ASSERT_LT(static_cast<std::size_t>(accepted_update), accepted_objectives.size());
                accepted_objectives[static_cast<std::size_t>(accepted_update)] =
                    NormalizedEndpointObjective(trial.endpoint, context.scale);
                accepted_gradients[static_cast<std::size_t>(accepted_update)] =
                    trial.endpoint.gradient.lpNorm<Eigen::Infinity>();
            }
            else if (!trial.accepted && trial.evaluation > 1)
            {
                observed_rejection = true;
                const auto event = std::find_if(events.begin(), events.end(), [&](const auto & value) {
                    return value.phase == joint::JointProgressPhase::SearchProgress
                        && value.profile_evaluations == trial.evaluation
                        && value.accepted_updates == accepted_update;
                });
                ASSERT_NE(event, events.end());
                ASSERT_TRUE(event->accepted_objective);
                EXPECT_DOUBLE_EQ(*event->accepted_objective,
                    *accepted_objectives[static_cast<std::size_t>(accepted_update)]);
                ASSERT_TRUE(event->accepted_gradient_inf_norm);
                EXPECT_DOUBLE_EQ(*event->accepted_gradient_inf_norm,
                    *accepted_gradients[static_cast<std::size_t>(accepted_update)]);
            }
        }
        EXPECT_TRUE(observed_rejection);
        for (const auto & objective : accepted_objectives) ASSERT_TRUE(objective);
        for (const auto & gradient : accepted_gradients) ASSERT_TRUE(gradient);
        EXPECT_DOUBLE_EQ(*observed.accepted_objective,
            *accepted_objectives[static_cast<std::size_t>(observed.accepted)]);
        EXPECT_DOUBLE_EQ(*observed.accepted_gradient_inf_norm,
            *accepted_gradients[static_cast<std::size_t>(observed.accepted)]);

        bool reported_initial_objective = false;
        for (const auto & event : events)
        {
            ASSERT_EQ(event.phase, joint::JointProgressPhase::SearchProgress);
            if (!event.accepted_objective)
            {
                EXPECT_FALSE(event.accepted_gradient_inf_norm);
                EXPECT_EQ(event.profile_evaluations, 1);
                EXPECT_EQ(event.accepted_updates, 0);
                continue;
            }
            ASSERT_TRUE(event.accepted_gradient_inf_norm);
            ASSERT_GE(event.accepted_updates, 0);
            ASSERT_LT(static_cast<std::size_t>(event.accepted_updates), accepted_objectives.size());
            EXPECT_DOUBLE_EQ(*event.accepted_objective,
                *accepted_objectives[static_cast<std::size_t>(event.accepted_updates)]);
            EXPECT_DOUBLE_EQ(*event.accepted_gradient_inf_norm,
                *accepted_gradients[static_cast<std::size_t>(event.accepted_updates)]);
            if (event.accepted_updates == 0) reported_initial_objective = true;
        }
        EXPECT_TRUE(reported_initial_objective);
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
    EXPECT_EQ(actual.solver_provenance.fixed_neighbor_maximum_sweeps, 30u);
    EXPECT_EQ(actual.solver_provenance.fixed_neighbor_order, "forward");
    EXPECT_EQ(actual.solver_provenance.fixed_neighbor_local_search, "legacy-compact");
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
    EXPECT_EQ(fit.solver_provenance.fixed_neighbor_maximum_sweeps,30u);
    EXPECT_EQ(fit.solver_provenance.fixed_neighbor_order,"forward");
    EXPECT_EQ(fit.solver_provenance.fixed_neighbor_local_search,"legacy-compact");
    ASSERT_TRUE(fit.assembled_state);
    EXPECT_EQ(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
    ASSERT_TRUE(events.front().solver_configuration);
    EXPECT_EQ(events.front().solver_configuration->fixed_neighbor_core_atoms,12u);
    ASSERT_TRUE(std::any_of(events.begin(),events.end(),[](const auto & event) {
        return event.fixed_neighbor && event.fixed_neighbor->sweep>0;
    }));
}

TEST(JointProgressTest, FixedNeighborObservableComponentsShareAssemblyContract)
{
    const core::JointProblem problem(MakeObservableInput());
    const std::vector<double> initial{.55,std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};
    const auto fit=joint::FitFixedNeighborComponents(problem,initial,{});
    ASSERT_EQ(problem.ParameterLayout().groups.size(),1u);
    ASSERT_EQ(fit.components.size(),1u);
    ASSERT_TRUE(fit.components.front().state);
    ASSERT_TRUE(fit.assembled_state);
    EXPECT_EQ(fit.components.front().layout->groups.size(),1u);
    EXPECT_EQ(fit.RuntimeConvergence(),rhbm_gem::JointCheckStatus::Passed);
}

TEST(JointProgressTest, ObservableComponentsUseTheSameLifecycle)
{
    const core::JointProblem problem(MakeObservableInput());
    const std::vector<double> initial{.55, std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};
    std::vector<joint::JointProgressEvent> events;
    const joint::JointProgressObserver observer = [&](const auto & event) { events.push_back(event); };
    const auto fit = joint::FitFixedNeighborComponents(problem, initial, {}, observer);
    ASSERT_EQ(problem.ParameterLayout().groups.size(), 1);
    ASSERT_TRUE(fit.assembled_state);
    ExpectLifecycle(problem, fit, events, {});
    for (std::size_t k = 0; k < fit.components.size(); ++k)
    {
        if (fit.components[k].stop_reason != "analytic-nuisance-only") continue;
        const auto completed = std::find_if(events.begin(), events.end(), [&](const auto & event) {
            return event.phase == joint::JointProgressPhase::ComponentCompleted
                && event.component_index == k + 1;
        });
        ASSERT_NE(completed, events.end());
        EXPECT_FALSE(completed->accepted_objective);
        EXPECT_FALSE(completed->accepted_gradient_inf_norm);
    }
}

TEST(JointProgressTest, CliReporterFormatsFixedNeighborConfiguration)
{
    const auto previous_level = Logger::GetLogLevel();
    Logger::SetLogLevel(LogLevel::Info);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();

    core::detail::JointCliProgressReporter reporter;
    joint::JointProgressEvent event;
    event.phase = joint::JointProgressPhase::SolverConfigured;
    event.solver_configuration = joint::ResolveJointSolverConfiguration({});
    reporter.OnProgress(event);

    const auto output = testing::internal::GetCapturedStdout();
    Logger::SetLogLevel(previous_level);
    const auto fixed = output.find(
        "[Joint] Solver: FixedNeighbor | sparse=EIGEN | core=12 | local-search=legacy-compact");
    ASSERT_NE(fixed, std::string::npos);
    const auto fixed_line_end = output.find('\n', fixed);
    EXPECT_EQ(output.substr(fixed, fixed_line_end - fixed).find("local-work="), std::string::npos);
    EXPECT_EQ(output.find("LegacyCompact (reference)"), std::string::npos);
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
    event.accepted_objective = 2.183e-4;
    event.accepted_gradient_inf_norm = 4.271e-6;
    reporter.OnProgress(event);

    event.phase = joint::JointProgressPhase::CertificationStarted;
    event.elapsed_seconds = .2;
    reporter.OnProgress(event);
    event.phase = joint::JointProgressPhase::ComponentCompleted;
    event.stop_reason = "native-lm-stop";
    event.trusted_state = true;
    event.accepted_objective = 2.183e-4;
    event.accepted_gradient_inf_norm = 4.271e-6;
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
    EXPECT_NE(output.find("searching | atoms=1 rows=40 | eval=2/200 | accepted=1/100 | obj=2.183e-04 | grad-inf=4.271e-06"), std::string::npos);
    EXPECT_NE(output.find("certifying endpoint | atoms=1 rows=40 | eval=2/200 | accepted=1/100 | accepted-obj=2.183e-04 | accepted-grad-inf=4.271e-06"), std::string::npos);
    EXPECT_NE(output.find("\n[Joint] Component 1/2 completed | atoms=1 rows=40"), std::string::npos);
    EXPECT_NE(output.find("accepted-obj=2.183e-04 | accepted-grad-inf=4.271e-06 | stop=native-lm-stop | trusted=yes"), std::string::npos);
    EXPECT_EQ(output.find("trusted-obj="), std::string::npos);
    EXPECT_EQ(output.find("trusted-grad="), std::string::npos);
    EXPECT_NE(output.find("[Joint] Assembling global state..."), std::string::npos);
    EXPECT_NE(output.find("[Joint] Global assembly completed | trusted=yes"), std::string::npos);
}

TEST(JointProgressTest, CliReporterOmitsUnavailableAcceptedMetrics)
{
    const auto previous_level = Logger::GetLogLevel();
    Logger::SetLogLevel(LogLevel::Info);
    Logger::FinishProgressLine();
    testing::internal::CaptureStdout();

    core::detail::JointCliProgressReporter reporter;
    joint::JointProgressEvent event;
    event.phase = joint::JointProgressPhase::ComponentStarted;
    event.component_index = 1;
    event.component_count = 1;
    event.profile_budget = 100;
    event.update_budget = 50;
    reporter.OnProgress(event);
    event.phase = joint::JointProgressPhase::SearchProgress;
    event.profile_evaluations = 1;
    reporter.OnProgress(event);
    event.phase = joint::JointProgressPhase::CertificationStarted;
    reporter.OnProgress(event);
    event.phase = joint::JointProgressPhase::ComponentCompleted;
    reporter.OnProgress(event);

    const auto output = testing::internal::GetCapturedStdout();
    Logger::SetLogLevel(previous_level);
    EXPECT_EQ(output.find(" | obj="), std::string::npos);
    EXPECT_EQ(output.find("grad-inf="), std::string::npos);
    EXPECT_EQ(output.find("accepted-obj="), std::string::npos);
    EXPECT_EQ(output.find("accepted-grad-inf="), std::string::npos);
    EXPECT_EQ(output.find("nan"), std::string::npos);
    EXPECT_EQ(output.find("N/A"), std::string::npos);
}

TEST(JointProgressTest, CliReporterRefreshesWhenAcceptedMetricsBecomeAvailable)
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
    event.profile_budget = 100;
    event.update_budget = 50;
    reporter.OnProgress(event);

    event.phase = joint::JointProgressPhase::SearchProgress;
    event.profile_evaluations = 1;
    event.elapsed_seconds = .1;
    reporter.OnProgress(event);
    event.accepted_objective = .1;
    event.accepted_gradient_inf_norm = .2;
    reporter.OnProgress(event);

    event = {};
    event.phase = joint::JointProgressPhase::ComponentStarted;
    event.component_index = 2;
    event.component_count = 2;
    event.profile_budget = 100;
    event.update_budget = 50;
    reporter.OnProgress(event);
    event.phase = joint::JointProgressPhase::SearchProgress;
    event.profile_evaluations = 1;
    event.elapsed_seconds = .1;
    reporter.OnProgress(event);
    event.accepted_objective = .3;
    event.accepted_gradient_inf_norm = .4;
    reporter.OnProgress(event);

    const auto output = testing::internal::GetCapturedStdout();
    Logger::SetLogLevel(previous_level);
    EXPECT_EQ(std::count(output.begin(), output.end(), '\r'), 6);
    EXPECT_NE(output.find("obj=1.000e-01 | grad-inf=2.000e-01"), std::string::npos);
    EXPECT_NE(output.find("obj=3.000e-01 | grad-inf=4.000e-01"), std::string::npos);
}

TEST(JointProgressTest, AnalyticNuisanceComponentHasNoAcceptedMetrics)
{
    const core::JointProblem problem(MakeObservableInputWithAnalyticNuisance());
    const std::vector<double> initial{.55, std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
    std::vector<joint::JointProgressEvent> events;
    const joint::JointProgressObserver observer = [&](const auto & event) { events.push_back(event); };
    const auto fit = joint::FitFixedNeighborComponents(problem, initial, {}, observer);

    const auto analytic = std::find_if(fit.components.begin(), fit.components.end(), [](const auto & component) {
        return component.stop_reason == "analytic-nuisance-only";
    });
    ASSERT_NE(analytic, fit.components.end());
    const auto component_index = static_cast<std::size_t>(analytic - fit.components.begin()) + 1;
    const auto completed = std::find_if(events.begin(), events.end(), [&](const auto & event) {
        return event.phase == joint::JointProgressPhase::ComponentCompleted
            && event.component_index == component_index;
    });
    ASSERT_NE(completed, events.end());
    EXPECT_FALSE(completed->accepted_objective);
    EXPECT_FALSE(completed->accepted_gradient_inf_norm);
    EXPECT_EQ(std::find_if(events.begin(), events.end(), [&](const auto & event) {
        return event.component_index == component_index
            && (event.phase == joint::JointProgressPhase::SearchProgress
                || event.phase == joint::JointProgressPhase::CertificationStarted);
    }), events.end());
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
