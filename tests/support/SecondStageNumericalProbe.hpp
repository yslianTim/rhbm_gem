#pragma once

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#include <array>
#include <cstddef>
#include <string>
#include <vector>
namespace rhbm_gem::core::detail {
struct CandidateCommitResult;
struct IterationResult;
struct FrozenBackground;
struct QuarantineState;
class TrustRegionStateSet;
}
namespace second_stage_test {
enum class Work { Solver, Operator, Objective, Snapshot };
struct NumericalCapture { std::string terminal{}; std::array<std::size_t,4> work{}; std::vector<std::string> commits{}; std::vector<std::vector<std::array<double,3>>> backgrounds{}; };
void BeginNumericalCapture();
NumericalCapture EndNumericalCapture();
void CountWork(Work) noexcept;
void CaptureTerminal(const rhbm_gem::core::detail::IterationResult &);
void CaptureBackground(const rhbm_gem::core::detail::FrozenBackground &);
void CaptureCommit(const rhbm_gem::core::detail::CandidateCommitResult &,
    const rhbm_gem::core::detail::QuarantineState &, const rhbm_gem::core::detail::TrustRegionStateSet &);
// Failure injection is linked only into testing builds, never a runtime mode.
enum class AuditFault { None, Allocation, Collection, Writer };
void SetAuditFault(AuditFault) noexcept;
void CheckAuditFault(AuditFault);
}
#define RHBM_TEST_WORK(kind) ::second_stage_test::CountWork(::second_stage_test::Work::kind)
#define RHBM_TEST_COMMIT(result, quarantine, radii) ::second_stage_test::CaptureCommit(result, quarantine, radii)
#define RHBM_TEST_AUDIT_FAULT(point) ::second_stage_test::CheckAuditFault(::second_stage_test::AuditFault::point)
#else
#define RHBM_TEST_WORK(kind) ((void)0)
#define RHBM_TEST_COMMIT(result, quarantine, radii) ((void)0)
#define RHBM_TEST_AUDIT_FAULT(point) ((void)0)
#endif

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#define RHBM_TEST_BACKGROUND(value) ::second_stage_test::CaptureBackground(value)
#else
#define RHBM_TEST_BACKGROUND(value) ((void)0)
#endif

#ifdef RHBM_GEM_TEST_INSTRUMENTATION
#define RHBM_TEST_TERMINAL(value) ::second_stage_test::CaptureTerminal(value)
#else
#define RHBM_TEST_TERMINAL(value) ((void)0)
#endif
