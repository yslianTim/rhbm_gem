#include "support/SecondStageNumericalProbe.hpp"
#include "support/EndpointRefinementExperiment.hpp"
#include "core/detail/second_stage/CandidateTransaction.hpp"
#include "core/detail/second_stage/IterationResult.hpp"
#include <atomic>
#include <iomanip>
#include <mutex>
#include <new>
#include <sstream>

namespace second_stage_test {
namespace {
std::atomic<bool> enabled{false};
std::array<std::atomic<std::size_t>,4> counts{};
std::mutex mutex;
std::vector<std::string> commits;
std::string terminal;
std::vector<std::vector<std::array<double,3>>> backgrounds;
std::atomic<AuditFault> fault{AuditFault::None};
}
void BeginNumericalCapture()
{
    terminal.clear(); backgrounds.clear(); commits.clear(); for(auto & count:counts) count=0; enabled=true;
}
NumericalCapture EndNumericalCapture()
{
    enabled=false; NumericalCapture result;
    for(std::size_t i=0;i<counts.size();++i) result.work[i]=counts[i].load();
    result.terminal=terminal; result.backgrounds=backgrounds; result.commits=commits; return result;
}
void CountWork(Work work) noexcept
{
    if(enabled.load(std::memory_order_relaxed) && !IsEndpointOperatorProbe()) counts[static_cast<std::size_t>(work)].fetch_add(1,std::memory_order_relaxed);
}
void CaptureCommit(const rhbm_gem::core::detail::CandidateCommitResult & result,
    const rhbm_gem::core::detail::QuarantineState & quarantine, const rhbm_gem::core::detail::TrustRegionStateSet & radii)
{
    if(!enabled) return;
    std::ostringstream out; out << std::setprecision(17);
    const auto keys=[&](const auto & list) {
        for(const auto & key:list) { out << '['; for(auto atom:key) out << atom << ','; out << "]r=" << radii.GetRadius(key) << ';'; }
    };
    out << "accepted:"; keys(result.accepted_key_list); out << "rejected:"; keys(result.rejected_key_list);
    out << "changed:"; keys(result.trust_region_update.changed_key_list);
    out << "q=" << quarantine.entered_target_count << '/' << quarantine.released_target_count << '/' << quarantine.failed_retry_count;
    for(const auto & [target,state]:quarantine.state_by_target)
    {
        out << '|' << static_cast<int>(target.kind) << ':'; for(auto atom:target.atom_index_list) out << atom << ',';
        out << ':' << static_cast<int>(state.lifecycle) << ':' << state.stable_iteration_count << ':' << state.last_recovery_revision;
        std::visit([&](const auto & reason) {
            using T=std::decay_t<decltype(reason)>;
            if constexpr(std::is_same_v<T,rhbm_gem::core::detail::JointOffsetSolveStatus>) out << ':' << static_cast<int>(reason);
            else out << ':' << static_cast<int>(reason.category) << ':' << (reason.guard_reason ? static_cast<int>(*reason.guard_reason) : -1);
        },state.reason);
    }
    std::lock_guard lock(mutex); commits.push_back(out.str());
}
void CaptureBackground(const rhbm_gem::core::detail::FrozenBackground & background)
{
    if (!enabled) return;
    std::vector<std::array<double,3>> models;
    for (const auto & model : background.model_by_atom)
        models.push_back({model.GetAmplitude(), model.GetWidth(), model.GetOffset()});
    backgrounds.push_back(std::move(models));
}
void CaptureTerminal(const rhbm_gem::core::detail::IterationResult & result)
{
    if (!enabled) return;
    terminal=std::to_string(static_cast<int>(result.stop_reason))+":"+
        std::to_string(result.attempt_number)+":"+std::to_string(result.accepted_iteration_count);
}
void SetAuditFault(AuditFault point) noexcept { fault=point; }
void CheckAuditFault(AuditFault point) { if(fault==point) throw std::bad_alloc{}; }
}
