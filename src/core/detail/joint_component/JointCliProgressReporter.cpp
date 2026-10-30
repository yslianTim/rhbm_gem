#include "JointCliProgressReporter.hpp"

#include <rhbm_gem/data/object/JointAnalysisResult.hpp>
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
        << " completed | state_trusted=" << (event.trusted_state ? "yes" : "no")
        << " | stop=" << event.stop_reason;
    if (event.accepted_objective)
        output << " | obj=" << Scientific(*event.accepted_objective);
    output << " | elapsed=" << Elapsed(event.elapsed_seconds);
    return output.str();
}

std::string ComponentDetail(const JointProgressEvent & event)
{
    // Component totals span blocks and sweeps; local budgets are not component limits.
    std::ostringstream output;
    output << "[Joint] Component detail | " << event.component_index << '/' << event.component_count
        << " | atoms=" << event.atom_count << " | rows=" << event.row_count
        << " | evaluations=" << event.profile_evaluations
        << " | accepted_updates=" << event.accepted_updates;
    output << " | stop=" << event.stop_reason
        << " | state_trusted=" << (event.trusted_state ? "yes" : "no");
    return output.str();
}
}

joint_component::JointProgressObserver JointCliProgressReporter::Observer()
{
    return [this](const auto & event) { OnProgress(event); };
}

void JointCliProgressReporter::ShowProgress(const JointProgressEvent & event, bool certifying)
{
    std::ostringstream output, detail;
    output << "[Joint] Component " << event.component_index << '/' << event.component_count << " | ";
    detail << "[Joint] Progress detail | component=" << event.component_index << '/' << event.component_count
        << " | phase=" << (certifying ? "certification" : "search")
        << " | atoms=" << event.atom_count << " | rows=" << event.row_count;
    if (event.fixed_neighbor)
    {
        const auto & progress=*event.fixed_neighbor;
        output << (certifying ? "certifying" : "searching")
            << " | sweep=" << progress.sweep << '/' << progress.maximum_sweeps
            << " | obj=" << Scientific(progress.objective)
            << " | KKT=" << Scientific(progress.global_ac_kkt)
            << " | width-grad=" << Scientific(progress.width_gradient_inf_norm)
            << " | width-kkt=" << Scientific(progress.projected_width_kkt)
            << " | accepted=" << progress.accepted_blocks
            << " | failed=" << progress.failed_blocks;
        detail << " | sweep=" << progress.sweep << '/' << progress.maximum_sweeps
            << " | blocks=" << progress.block_solves
            << " | accepted_blocks=" << progress.accepted_blocks
            << " | unchanged_blocks=" << progress.unchanged_blocks
            << " | failed_blocks=" << progress.failed_blocks
            << " | local_accepted_updates=" << progress.local_accepted_updates
            << " | profile_evaluations=" << progress.profile_evaluations
            << " | recovery=" << progress.recovery_attempts << '/' << progress.recovery_successes
            << " | max_B=" << Scientific(progress.maximum_target_width);
    }
    else
    {
        output << (certifying ? "certifying" : "searching")
            << " | eval=" << event.profile_evaluations << '/' << event.profile_budget
            << " | accepted=" << event.accepted_updates << '/' << event.update_budget;
        if (event.accepted_objective) output << " | obj=" << Scientific(*event.accepted_objective);
        if (event.accepted_gradient_inf_norm)
            output << " | gradient-inf=" << Scientific(*event.accepted_gradient_inf_norm);
        detail << " | eval=" << event.profile_evaluations << '/' << event.profile_budget
            << " | accepted=" << event.accepted_updates << '/' << event.update_budget;
        if (event.accepted_objective) detail << " | obj=" << Scientific(*event.accepted_objective);
        if (event.accepted_gradient_inf_norm)
            detail << " | gradient_inf=" << Scientific(*event.accepted_gradient_inf_norm);
    }
    output << " | elapsed=" << Elapsed(event.elapsed_seconds);
    if (Logger::GetLogLevel() == LogLevel::Debug)
    {
        Logger::FinishProgressLine();
        Logger::Log(LogLevel::Debug, detail.str());
    }
    Logger::ProgressLine(output.str());
    m_component_index = event.component_index;
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
        LogJointInfo("[Joint] Solver | FixedNeighbor / Active-set LM");
        std::ostringstream detail;
        detail << "[Joint] Solver detail | sparse=EIGEN | core=" << configuration.fixed_neighbor_core_atoms
            << " | local-search=" << FixedNeighborBoundedLocalSearchName;
        Logger::Log(LogLevel::Debug, detail.str());
        return;
    }
    case JointProgressPhase::ComponentStarted:
        if (!m_solver_started)
        {
            LogJointInfo("[Joint] Solving " + std::to_string(event.component_count) + " structural components");
            m_solver_started = true;
        }
        ShowProgress(event, false);
        return;
    case JointProgressPhase::SearchProgress:
    {
        const auto now = Clock::now();
        const bool phase_changed = !m_last_phase || *m_last_phase != event.phase;
        const bool elapsed = !m_has_output || now - m_last_output >= std::chrono::milliseconds(500);
        if (event.component_index != m_component_index || phase_changed || elapsed)
            ShowProgress(event, false);
        return;
    }
    case JointProgressPhase::CertificationStarted:
        ShowProgress(event, true);
        return;
    case JointProgressPhase::ComponentCompleted:
        Logger::FinishProgressLine();
        Logger::Log(LogLevel::Debug, ComponentDetail(event));
        Logger::Log(LogLevel::Info, ComponentSummary(event));
        m_has_output = false;
        m_last_phase = event.phase;
        return;
    case JointProgressPhase::AssemblyStarted:
        LogJointInfo("[Joint] Assembling global state...");
        return;
    case JointProgressPhase::AssemblyCompleted:
    {
        if (!event.trusted_state)
        {
            Logger::FinishProgressLine();
            Logger::Log(LogLevel::Warning, "[Joint] Global assembly did not produce a trusted state | elapsed=" +
                Elapsed(event.elapsed_seconds));
            return;
        }
        std::ostringstream output;
        output << "[Joint] Global assembly | trusted=yes | elapsed=" << Elapsed(event.elapsed_seconds);
        LogJointInfo(output.str());
        return;
    }
    }
}

} // namespace rhbm_gem::core::detail
