#pragma once

#include "SolverRoute.hpp"

#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <string>

namespace rhbm_gem::core::joint_component {

enum class JointProgressPhase
{
    SolverConfigured,
    ComponentStarted,
    SearchProgress,
    CertificationStarted,
    ComponentCompleted,
    AssemblyStarted,
    AssemblyCompleted
};

struct JointFixedNeighborProgress
{
    std::size_t sweep{},maximum_sweeps{},block_solves{},accepted_blocks{},local_accepted_updates{};
    double objective{std::numeric_limits<double>::quiet_NaN()};
    double global_ac_kkt{std::numeric_limits<double>::quiet_NaN()};
    double width_gradient_inf_norm{std::numeric_limits<double>::quiet_NaN()};
    double eta_change_inf{std::numeric_limits<double>::quiet_NaN()};
};

struct JointProgressEvent
{
    JointProgressPhase phase{};
    std::optional<JointSolverRoute> solver_route;
    std::size_t component_index{};
    std::size_t component_count{};
    std::string component_id;
    std::size_t atom_count{};
    std::size_t row_count{};
    int profile_evaluations{};
    int profile_budget{};
    int accepted_updates{};
    int update_budget{};
    std::optional<double> accepted_objective;
    std::optional<double> accepted_gradient_inf_norm;
    double elapsed_seconds{};
    std::string stop_reason;
    bool trusted_state{};
    std::optional<JointFixedNeighborProgress> fixed_neighbor;
};

using JointProgressObserver = std::function<void(const JointProgressEvent &)>;

inline void NotifyJointSolverConfigured(
    const JointProgressObserver & observer, const JointSolverRoute & route)
{
    if (!observer) return;
    JointProgressEvent event;
    event.phase = JointProgressPhase::SolverConfigured;
    event.solver_route = route;
    observer(event);
}

struct JointProgressComponent
{
    std::size_t index{};
    std::size_t count{};
    std::string id;
    std::size_t atom_count{};
    std::size_t row_count{};
};

inline JointProgressEvent MakeJointProgressEvent(
    JointProgressPhase phase, const JointProgressComponent & component)
{
    JointProgressEvent event;
    event.phase = phase;
    event.component_index = component.index;
    event.component_count = component.count;
    event.component_id = component.id;
    event.atom_count = component.atom_count;
    event.row_count = component.row_count;
    return event;
}

inline void NotifyJointProgress(
    const JointProgressObserver & observer, const JointProgressEvent & event)
{
    if (observer) observer(event);
}

inline void NotifyJointProgress(
    const JointProgressObserver & observer,
    JointProgressPhase phase,
    const JointProgressComponent & component,
    int profile_evaluations = 0,
    int profile_budget = 0,
    int accepted_updates = 0,
    int update_budget = 0,
    double elapsed_seconds = 0,
    const std::string & stop_reason = {},
    bool trusted_state = false,
    std::optional<double> accepted_objective = std::nullopt,
    std::optional<double> accepted_gradient_inf_norm = std::nullopt)
{
    if (!observer) return;
    auto event = MakeJointProgressEvent(phase, component);
    event.profile_evaluations = profile_evaluations;
    event.profile_budget = profile_budget;
    event.accepted_updates = accepted_updates;
    event.update_budget = update_budget;
    event.elapsed_seconds = elapsed_seconds;
    event.stop_reason = stop_reason;
    event.trusted_state = trusted_state;
    event.accepted_objective = accepted_objective;
    event.accepted_gradient_inf_norm = accepted_gradient_inf_norm;
    observer(event);
}

inline void NotifyJointAssemblyProgress(
    const JointProgressObserver & observer,
    JointProgressPhase phase,
    std::size_t component_count,
    double elapsed_seconds = 0,
    bool trusted_state = false)
{
    if (!observer) return;
    JointProgressEvent event;
    event.phase = phase;
    event.component_count = component_count;
    event.elapsed_seconds = elapsed_seconds;
    event.trusted_state = trusted_state;
    observer(event);
}

} // namespace rhbm_gem::core::joint_component
