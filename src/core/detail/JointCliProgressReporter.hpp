#pragma once

#include "joint_component/JointProgress.hpp"

#include <chrono>
#include <optional>

namespace rhbm_gem::core::detail {

class JointCliProgressReporter
{
public:
    joint_component::JointProgressObserver Observer();
    void OnProgress(const joint_component::JointProgressEvent & event);

private:
    using Clock = std::chrono::steady_clock;
    bool m_solver_started{};
    bool m_has_output{};
    bool m_has_accepted_objective{};
    bool m_has_accepted_gradient{};
    std::size_t m_component_index{};
    int m_accepted_updates{};
    Clock::time_point m_last_output{};
    std::optional<joint_component::JointProgressPhase> m_last_phase;

    void ShowProgress(const joint_component::JointProgressEvent & event, bool certifying);
};

} // namespace rhbm_gem::core::detail
