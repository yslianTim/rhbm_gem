#include "JointCliProgressReporter.hpp"

#include <rhbm_gem/utils/domain/Logger.hpp>

#include <iomanip>
#include <sstream>
#include <string>

namespace rhbm_gem::core::detail {
namespace {
using joint_component::JointProgressEvent;

void LogJointInfo(const std::string & message)
{
    Logger::FinishProgressLine();
    Logger::Log(LogLevel::Info, message);
}

std::string Elapsed(double seconds)
{
    std::ostringstream output;
    output << std::fixed << std::setprecision(1) << seconds << 's';
    return output.str();
}

std::string Scientific(double value)
{
    std::ostringstream output;
    output << std::scientific << std::setprecision(3) << value;
    return output.str();
}

std::string ComponentSummary(const JointProgressEvent & event)
{
    std::ostringstream output;
    output << "[Joint] Component " << event.component_index << '/' << event.component_count
        << " completed | atoms=" << event.atom_count << " rows=" << event.row_count
        << " | eval=" << event.profile_evaluations
        << " accepted=" << event.accepted_updates;
    if (event.accepted_objective)
        output << " | accepted-obj=" << Scientific(*event.accepted_objective);
    if (event.accepted_gradient_inf_norm)
        output << " | accepted-grad-inf=" << Scientific(*event.accepted_gradient_inf_norm);
    output << " | stop=" << event.stop_reason
        << " | trusted=" << (event.trusted_state ? "yes" : "no")
        << " | elapsed=" << Elapsed(event.elapsed_seconds);
    return output.str();
}
}

joint_component::JointProgressObserver JointCliProgressReporter::Observer()
{
    return [this](const auto & event) { OnProgress(event); };
}

void JointCliProgressReporter::ShowProgress(const JointProgressEvent & event, bool certifying)
{
    std::ostringstream output;
    output << "Joint component " << event.component_index << '/' << event.component_count << " | "
        << (certifying ? "certifying endpoint" : "searching")
        << " | atoms=" << event.atom_count << " rows=" << event.row_count;
    if (event.fixed_neighbor)
    {
        const auto & progress=*event.fixed_neighbor;
        output << " | sweep=" << progress.sweep << '/' << progress.maximum_sweeps
            << " | blocks=" << progress.block_solves << " | accepted-blocks=" << progress.accepted_blocks
            << " | local-accepted=" << progress.local_accepted_updates
            << " | kkt=" << Scientific(progress.global_ac_kkt)
            << " | width-grad=" << Scientific(progress.width_gradient_inf_norm)
            << " | eta-change=" << Scientific(progress.eta_change_inf)
            << " | objective=" << Scientific(progress.objective);
    }
    else output << " | eval=" << event.profile_evaluations << '/' << event.profile_budget
        << " | accepted=" << event.accepted_updates << '/' << event.update_budget;
    if (event.accepted_objective)
        output << " | " << (certifying ? "accepted-obj=" : "obj=") << Scientific(*event.accepted_objective);
    if (event.accepted_gradient_inf_norm)
        output << " | " << (certifying ? "accepted-grad-inf=" : "grad-inf=")
            << Scientific(*event.accepted_gradient_inf_norm);
    output << " | " << Elapsed(event.elapsed_seconds);
    Logger::ProgressLine(output.str());
    m_component_index = event.component_index;
    m_accepted_updates = event.accepted_updates;
    m_has_accepted_objective = event.accepted_objective.has_value();
    m_has_accepted_gradient = event.accepted_gradient_inf_norm.has_value();
    m_last_phase = event.phase;
    m_last_output = Clock::now();
    m_has_output = true;
}

void JointCliProgressReporter::OnProgress(const JointProgressEvent & event)
{
    using joint_component::JointProgressPhase;
    switch (event.phase)
    {
    case JointProgressPhase::SolverConfigured:
    {
        if (!event.solver_configuration) return;
        const auto & configuration = *event.solver_configuration;
        std::ostringstream output;
        output << "[Joint] Solver: FixedNeighbor | sparse=EIGEN"
            << " | core=" << configuration.fixed_neighbor_core_atoms
            << " | local-search=legacy-compact";
        LogJointInfo(output.str());
        return;
    }
    case JointProgressPhase::ComponentStarted:
        if (!m_solver_started)
        {
            LogJointInfo("[Joint] Solving " + std::to_string(event.component_count) + " structural components");
            m_solver_started = true;
        }
        m_has_accepted_objective = false;
        m_has_accepted_gradient = false;
        ShowProgress(event, false);
        return;
    case JointProgressPhase::SearchProgress:
    {
        const auto now = Clock::now();
        const bool accepted_changed = event.accepted_updates != m_accepted_updates;
        const bool metrics_became_available = (!m_has_accepted_objective && event.accepted_objective.has_value())
            || (!m_has_accepted_gradient && event.accepted_gradient_inf_norm.has_value());
        const bool phase_changed = !m_last_phase || *m_last_phase != event.phase;
        const bool elapsed = !m_has_output || now - m_last_output >= std::chrono::milliseconds(500);
        if (event.component_index != m_component_index || accepted_changed || phase_changed
            || metrics_became_available || elapsed)
            ShowProgress(event, false);
        return;
    }
    case JointProgressPhase::CertificationStarted:
        ShowProgress(event, true);
        return;
    case JointProgressPhase::ComponentCompleted:
        Logger::FinishProgressLine();
        LogJointInfo(ComponentSummary(event));
        m_has_output = false;
        m_last_phase = event.phase;
        return;
    case JointProgressPhase::AssemblyStarted:
        LogJointInfo("[Joint] Assembling global state...");
        return;
    case JointProgressPhase::AssemblyCompleted:
    {
        std::ostringstream output;
        output << "[Joint] Global assembly completed | trusted=" << (event.trusted_state ? "yes" : "no")
            << " | elapsed=" << Elapsed(event.elapsed_seconds);
        LogJointInfo(output.str());
        return;
    }
    }
}

} // namespace rhbm_gem::core::detail
