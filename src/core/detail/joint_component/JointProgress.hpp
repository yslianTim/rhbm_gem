#pragma once

#include <cstddef>
#include <functional>
#include <string>

namespace rhbm_gem::core::joint_component {

enum class JointProgressPhase
{
    ComponentStarted,
    SearchProgress,
    CertificationStarted,
    ComponentCompleted,
    AssemblyStarted,
    AssemblyCompleted
};

struct JointProgressEvent
{
    JointProgressPhase phase{};
    std::size_t component_index{};
    std::size_t component_count{};
    std::string component_id;
    std::size_t atom_count{};
    std::size_t row_count{};
    int profile_evaluations{};
    int profile_budget{};
    int accepted_updates{};
    int update_budget{};
    double elapsed_seconds{};
    std::string stop_reason;
    bool trusted_state{};
};

using JointProgressObserver = std::function<void(const JointProgressEvent &)>;

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

} // namespace rhbm_gem::core::joint_component
